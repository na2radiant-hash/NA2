%% ============================================================
%  STEP 2: RAY TRACING THROUGH THE CONSTRUCTED LENS
%           (solid-body fill + heatmap + emission profile)
%
%  RAY EMPHASIS MODE fixes:
%   1) grid equal (square grid cells)
%   2) extend rays to the right boundary even if they go "backwards"
%   3) lens MUCH more faded (fill + border lines)
%   4) remove heatmap (and its colorbar) in this mode
% ============================================================
clear; clc; close all;

%% -------------------------
%  User Parameters
% -------------------------
input_filename   = 'lens_data_gradual_N2003_test.mat';  % <- your file
%input_filename  = 'lens_data.mat';

N_rays           = 10;
seed             = 42;
draw_mirror      = true;
emission_profile = 'uniform';                 % 'uniform' or 'sin2'
rng(seed);

% Uniform emission params
uniform_min_angle = 0;                        % [deg]
uniform_max_angle = 178;                      % [deg] closer to 180 (avoid exact 180)
uniform_eps       = 1e-6;                     % [deg]

% =========================
%  NEW: Ray emphasis switch
% =========================
RAY_EMPHASIS_MODE = false;                     % <-- ON/OFF

%% -------------------------
%  Styling
% -------------------------
rayColor = [0.85 0.65 0.13];

if RAY_EMPHASIS_MODE
    % Rays: stronger
    rayLW = 1.6;

    % Lens: VERY faded (fill + edges)
    lensFace  = [0.68 0.92 0.68];
    lensAlpha = 0.015;                        % almost invisible fill
    edgeColor = [0.80 0.80 0.80];             % very light borders
    edgeLW    = 0.35;                         % very thin borders

    % Heatmap disabled in this mode
    ENABLE_HEATMAP = false;

    % Force extension to the right boundary (fix #2)
    EXTEND_TO_RIGHT_ALWAYS = true;
else
    rayLW = 0.8;

    lensFace  = [0.68 0.92 0.68];
    lensAlpha = 0.45;
    edgeColor = [0 0 0];
    edgeLW    = 1.25;

    ENABLE_HEATMAP = true;
    EXTEND_TO_RIGHT_ALWAYS = false;
end

%% -------------------------
%  Load lens geometry
% -------------------------
S = load(input_filename);
required = {'final_lens','n1','n2','layers_count','source_pos'};
for f = required
    if ~isfield(S, f{1}), error('Missing "%s" in %s', f{1}, input_filename); end
end
final_lens    = S.final_lens;
n1_global     = S.n1;
n2_global     = S.n2;
L             = S.layers_count;
origin_global = S.source_pos(:).';

% Ensure usable layers
L = min(L, numel(final_lens));
has_points = cellfun(@(c) ~isempty(c) && size(c,1) >= 2, final_lens);
if ~any(has_points), error('No usable layers found in final_lens.'); end

% Plot extents
allX   = cell2mat(cellfun(@(c) c(:,1), final_lens(has_points), 'uni', 0));
xmin   = min(allX);
xmax   = max(allX);
x_span = max(1e-9, xmax - xmin);
x_end  = xmax + 0.30 * x_span;

%% -------------------------
%  Build a solid-body lens polyshape (no leaking)
% -------------------------
layerPoly = cell(L,1);
for j = 1:L
    if isempty(final_lens{j}), continue; end
    up = final_lens{j}; lo = up; lo(:,2) = -up(:,2);
    ring = [up; flipud(lo)];
    layerPoly{j} = polyshape(ring(:,1), ring(:,2), 'Simplify', true);
end

lensBody = polyshape;
for j = 2:2:L
    if isempty(layerPoly{j-1}) || isempty(layerPoly{j}), continue; end
    shell = subtract(layerPoly{j}, layerPoly{j-1});
    lensBody = union(lensBody, shell);
end

%% -------------------------
%  Plot setup (MAIN)
% -------------------------
figure('Name','Ray Tracing Through Constructed Lens'); clf; hold on; axis equal;
set(gcf, 'Color', 'w');
xlabel('X (mm)'); ylabel('Y (mm)');
if RAY_EMPHASIS_MODE
    title(sprintf('Ray Tracing (%s emission) — Ray emphasis', emission_profile));
else
    title(sprintf('Ray Tracing (%s emission)', emission_profile));
end

% Fix #1: square grid cells
grid on;
set(gca,'DataAspectRatio',[1 1 1]);   % keep equal geometry
set(gca,'PlotBoxAspectRatio',[1 1 1]);% keep square plot box (square grid appearance)

% Lens fill (very faded in ray emphasis)
if ~isempty(lensBody.Vertices)
    plot(lensBody, 'FaceColor', lensFace, 'FaceAlpha', lensAlpha, 'EdgeColor','none');
end

% Lens edges (very faded in ray emphasis)
for j = 1:L
    if isempty(final_lens{j}), continue; end
    up = final_lens{j}; lo = up; lo(:,2) = -up(:,2);
    plot(up(:,1), up(:,2), 'Color', edgeColor, 'LineWidth', edgeLW);
    if draw_mirror
        plot(lo(:,1), lo(:,2), 'Color', edgeColor, 'LineWidth', edgeLW);
    end
end

% Source
plot(origin_global(1), origin_global(2), 'ro', 'MarkerFaceColor','r');

%% -------------------------
%  Heatmap logger (only when enabled)
% -------------------------
y_end_hits = [];

%% -------------------------
%  Sample launch angles, then trace
% -------------------------
alphas0 = sample_ray_angles(N_rays, emission_profile, seed, ...
                            uniform_min_angle, uniform_max_angle, uniform_eps);

for r = 1:N_rays
    alpha0    = alphas0(r);
    P0        = origin_global;
    alpha_cur = alpha0;
    n_in      = n1_global;
    n_out     = n2_global;

    path = P0;
    hit_all_layers = true;

    for j = 1:L
        [hit, P_hit, alpha_next] = first_hit_with_layer(P0, alpha_cur, j, n_in, n_out, final_lens);
        if ~hit, hit_all_layers = false; break; end

        path = [path; P_hit]; %#ok<AGROW>
        P0 = P_hit; alpha_cur = alpha_next;

        tmp = n_in; n_in = n_out; n_out = tmp;
    end

    % Fix #2: extend to x_end even if Df(1) is negative (or we are left of x_end)
    if hit_all_layers
        Df = [cosd(alpha_cur), sind(alpha_cur)];
        if abs(Df(1)) > 1e-12
            t_ext = (x_end - P0(1)) / Df(1);

            if EXTEND_TO_RIGHT_ALWAYS
                % draw extension regardless of sign (so "top ray" reaches the right)
                Pend = P0 + t_ext * Df;
                path = [path; Pend]; %#ok<AGROW>

                if ENABLE_HEATMAP
                    y_end_hits(end+1,1) = Pend(2); %#ok<AGROW>
                end
            else
                % original: only forward and only if x_end is ahead
                if (x_end > P0(1)) && (t_ext > 0)
                    Pend = P0 + t_ext * Df;
                    path = [path; Pend]; %#ok<AGROW>
                    if ENABLE_HEATMAP
                        y_end_hits(end+1,1) = Pend(2); %#ok<AGROW>
                    end
                end
            end
        end
    end

    plot(path(:,1), path(:,2), '-', 'Color', rayColor, 'LineWidth', rayLW);

    if draw_mirror
        path_m = path; path_m(:,2) = -path_m(:,2);
        plot(path_m(:,1), path_m(:,2), '-', 'Color', rayColor, 'LineWidth', rayLW);
    end
end

% View limits
if ~isempty(lensBody.Vertices)
    y_min_plot = min(lensBody.Vertices(:,2));
    y_max_plot = max(lensBody.Vertices(:,2));
else
    y_all_layers = cell2mat(cellfun(@(c) [c(:,2); -c(:,2)], final_lens(~cellfun('isempty',final_lens)), 'uni', 0));
    y_min_plot = min(y_all_layers); y_max_plot = max(y_all_layers);
end
y_lim = max(abs([y_min_plot, y_max_plot]));
ylim([-y_lim, y_lim]);
xlim([xmin - 0.1*x_span, x_end + 0.15*x_span]);

fprintf('Launched %d rays (%s profile)\n', N_rays, emission_profile);

%% -------------------------
%  Heatmap (DISABLED in ray-emphasis mode)
% -------------------------
if ENABLE_HEATMAP && ~isempty(y_end_hits)
    % Include mirrored rays
    y_det = [y_end_hits(:); -y_end_hits(:)];
    total_hits = numel(y_det);

    yl = ylim;
    y_min_hm = yl(1);
    y_max_hm = yl(2);

    nbins   = 200;
    edges   = linspace(y_min_hm, y_max_hm, nbins+1);
    centers = (edges(1:end-1) + edges(2:end)) / 2;
    dy      = edges(2) - edges(1);

    counts  = histcounts(y_det, edges).';
    frac    = counts / max(total_hits, 1);
    percent = 100 * frac;

    maxP = max(percent);
    if maxP == 0, maxP = 1; end

    C      = [percent percent];
    bar_w  = 0.03 * (xmax - xmin + eps);
    x0     = x_end;

    imagesc('XData',[x0, x0+bar_w], ...
            'YData',[centers(1), centers(end)], ...
            'CData', C);
    set(gca,'YDir','normal');

    ncm   = 256;
    yellow = [1 1 0.6];
    red    = [0.55 0 0];
    cmap   = [linspace(yellow(1),red(1),ncm)', ...
              linspace(yellow(2),red(2),ncm)', ...
              linspace(yellow(3),red(3),ncm)'];
    colormap(gca, cmap);
    caxis([0 maxP]);
    alpha(0.95);

    plot([x_end x_end],[y_min_hm y_max_hm],'k-','LineWidth',0.8);
    plot([x0 x0+bar_w x0+bar_w x0 x0], ...
         [y_min_hm y_min_hm y_max_hm y_max_hm y_min_hm], ...
         'k-','LineWidth',0.75);

    cb = colorbar('eastoutside');
    cb.Label.String = sprintf('Ray density (%% of rays per Δy bin)\\nΔy = %.3g mm', dy);
    cb.TickDirection = 'out';
end

%% ============================================================
%  EXTRA FIGURE: CLOSE-UP VIEW AROUND THE SOURCE (±5 mm window)
% ============================================================
figure('Name','Ray Tracing – Close Up'); clf; hold on; axis equal;
set(gcf, 'Color', 'w');
xlabel('X (mm)'); ylabel('Y (mm)');
title('Ray Tracing – Close-Up View');
grid on;
set(gca,'DataAspectRatio',[1 1 1]);
set(gca,'PlotBoxAspectRatio',[1 1 1]);

if ~isempty(lensBody.Vertices)
    plot(lensBody, 'FaceColor', lensFace, 'FaceAlpha', lensAlpha, 'EdgeColor','none');
end

for j = 1:L
    if isempty(final_lens{j}), continue; end
    up = final_lens{j}; lo = up; lo(:,2) = -up(:,2);
    plot(up(:,1), up(:,2), 'Color', edgeColor, 'LineWidth', edgeLW);
    if draw_mirror
        plot(lo(:,1), lo(:,2), 'Color', edgeColor, 'LineWidth', edgeLW);
    end
end

for r = 1:N_rays
    alpha0    = alphas0(r);
    P0        = origin_global;
    alpha_cur = alpha0;
    n_in      = n1_global;
    n_out     = n2_global;

    path = P0;
    hit_all = true;

    for j = 1:L
        [hit, P_hit, alpha_next] = first_hit_with_layer(P0, alpha_cur, j, n_in, n_out, final_lens);
        if ~hit, hit_all = false; break; end
        path = [path; P_hit]; %#ok<AGROW>
        P0 = P_hit; alpha_cur = alpha_next;
        tmp = n_in; n_in = n_out; n_out = tmp;
    end

    plot(path(:,1), path(:,2), '-', 'Color', rayColor, 'LineWidth', rayLW);

    if draw_mirror
        path_m = path; path_m(:,2) = -path_m(:,2);
        plot(path_m(:,1), path_m(:,2), '-', 'Color', rayColor, 'LineWidth', rayLW);
    end
end

plot(origin_global(1), origin_global(2), 'ro', 'MarkerFaceColor','r');
xlim([-5 5]); ylim([-5 5]);
axis square;

%% =========================
%  Local helpers
% =========================

function [alpha_new, intersection_point, theta_in, theta_out] = refract_ray_2d(origin_a, alpha, surface_a, surface_b, n1, n2)
D = [cosd(alpha), sind(alpha)]; O = origin_a(:).';
A = surface_a(:).'; B = surface_b(:).'; S = B - A;
M = [D(:), -S(:)]; rhs = (A - O).'; detM = det(M);
if abs(detM) < 1e-12, error('Ray and segment are parallel or nearly parallel; no unique intersection.'); end
sol = M \ rhs; t = sol(1); u = sol(2);
if t < -1e-12 || u < -1e-12 || u > 1+1e-12, error('No valid intersection.'); end
t = max(t,0); u = min(max(u,0),1);
intersection_point = O + t*D;

if norm(S) < 1e-12, error('Surface endpoints are identical.'); end
t_hat = S / norm(S);
n_hat_ccw = [-t_hat(2),  t_hat(1)];
n_hat_cw  = [ t_hat(2), -t_hat(1)];
if dot(D, n_hat_ccw) >= dot(D, n_hat_cw), n_hat = n_hat_ccw; else, n_hat = n_hat_cw; end

angN   = atan2d(n_hat(2), n_hat(1));
angRay = atan2d(D(2), D(1));
theta_in = wrapTo180(angRay - angN);
if abs(theta_in) > 90
    angN = wrapTo180(angN + 180);
    theta_in = wrapTo180(angRay - angN);
end

sin_out = (n1/n2) * sind(theta_in);
if abs(sin_out) > 1 + 1e-12
    alpha_new = NaN; theta_out = NaN; warning('Total internal reflection.'); return;
end
sin_out = max(min(sin_out,1),-1);
theta_out = asind(sin_out);
alpha_new = wrapTo360(angN + theta_out);
end

function a = wrapTo180(a), a = mod(a + 180, 360) - 180; end
function a = wrapTo360(a), a = mod(a, 360); if a < 0, a = a + 360; end, end

function [hit, P_hit, alpha_out] = first_hit_with_layer(P0, alpha_in, j, n_in, n_out, final_lens)
hit=false; P_hit=[NaN,NaN]; alpha_out=NaN;
if j>numel(final_lens)||isempty(final_lens{j})||size(final_lens{j},1)<2, return; end
pts=final_lens{j}; best_d=inf;
for k=1:size(pts,1)-1
    A=pts(k,:); B=pts(k+1,:);
    try
        [alpha_new,P_inter] = refract_ray_2d(P0,alpha_in,A,B,n_in,n_out);
    catch
        continue;
    end
    if any(isnan(P_inter))||isnan(alpha_new), continue; end
    d=hypot(P_inter(1)-P0(1),P_inter(2)-P0(2));
    if d>1e-9 && d<best_d
        best_d=d; P_hit=P_inter; alpha_out=alpha_new; hit=true;
    end
end
end

function alphas = sample_ray_angles(N, profile, seed, a_min, a_max, eps_ang)
if nargin >= 3 && ~isempty(seed), rng(seed); end
if nargin < 4 || isempty(a_min),   a_min   = 0;   end
if nargin < 5 || isempty(a_max),   a_max   = 180; end
if nargin < 6 || isempty(eps_ang), eps_ang = 0;   end

profile = lower(strrep(profile,'^',''));
switch profile
    case 'uniform'
        if a_max <= a_min
            error('uniform_max_angle must be > uniform_min_angle.');
        end
        lo = a_min + eps_ang;
        hi = a_max - eps_ang;
        if hi <= lo
            error('uniform_eps is too large for the [uniform_min_angle, uniform_max_angle] interval.');
        end
        alphas = linspace(lo, hi, N).';

    case 'sin2'
        alphas=zeros(N,1); k=0;
        while k<N
            M=max(1000,2*(N-k));
            c=180*rand(M,1); u=rand(M,1);
            acc=u<=sind(c).^2;
            nacc=sum(acc);
            if nacc>0
                take=min(nacc,N-k);
                idx=find(acc,take,'first');
                alphas(k+1:k+take)=c(idx);
                k=k+take;
            end
        end
    otherwise
        error('Unknown profile "%s". Use "uniform" or "sin2".',profile);
end
end
