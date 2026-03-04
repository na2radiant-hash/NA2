%% ============================================================
%  STEP 2: RAY TRACING THROUGH THE CONSTRUCTED LENS
%           + OPTIONAL MULTI-GLASS GAP REGIONS (n3) BETWEEN ELEMENTS
%
%  Updated to match STEP 1 (multi-glass + correct heights):
%    - Reads S.glass_info.gaps if present (preferred).
%    - Each gap uses x_left/x_right and h_glass already computed in STEP 1
%      as 1.3 * height(next element).
%    - Ray path: inserts entry/exit points for *all* gaps crossed by each
%      straight segment (ordered along propagation).
% ============================================================
clear; clc; close all;

%% -------------------------
%  User Parameters
%  -------------------------
input_filename   = 'lens_data_4_elements_with_glass.mat';

N_rays           = 10;                   % number of rays
seed             = 42;
draw_mirror      = true;
emission_profile = 'uniform';            % 'uniform' or 'sin2'
rng(seed);

% Uniform emission params
uniform_min_angle = 0;                   % [deg]
uniform_max_angle = 148;                 % [deg]
uniform_eps       = 1e-6;                % [deg]

% =========================
%  Multi-glass parameters (fallback if .mat lacks glass_info.gaps)
% =========================
ENABLE_GLASS = true;          % set false to ignore glass even if present in .mat
glass_pos    = [2 3 4];       % scalar or vector: insert glass BEFORE these elements
n3           = 1.5;           % refractive index in glass region(s)
glass_width  = 5;             % [mm] width of each glass region along x

% Colors
rayColor   = [0.85 0.65 0.13];          % dark yellow for rays
lensFace   = [0.68 0.92 0.68];          % light green
edgeColor  = [0 0 0];
glassFace  = [0.65 0.85 1.00];          % light blue fill
glassEdge  = [0.20 0.45 0.75];

%% -------------------------
%  Load lens geometry (+ optional glass_info)
%  -------------------------
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

% Common x_end = xmax + 30% span
allX   = cell2mat(cellfun(@(c) c(:,1), final_lens(has_points), 'uni', 0));
xmin   = min(allX);
xmax   = max(allX);
x_span = max(1e-9, xmax - xmin);
x_end  = xmax + 0.30 * x_span;

%% -------------------------
%  Resolve MULTI-glass region geometry list: gaps(k) has fields xL,xR,h,n
% -------------------------
[gaps, glass_enabled] = resolve_glass_gaps( ...
    ENABLE_GLASS, S, final_lens, L, has_points, glass_pos, n3, glass_width);

%% -------------------------
%  Build a *solid-body* polyshape of the lens (no leaking)
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
%  Plot setup (main)
% -------------------------
figure('Name','Ray Tracing Through Constructed Lens'); clf; hold on; grid on; axis equal;
set(gcf, 'Color', 'w');
xlabel('X (mm)');
ylabel('Y (mm)');
title(sprintf('Ray Tracing (%s emission)', emission_profile));

% Lens fill
if ~isempty(lensBody.Vertices)
    plot(lensBody, 'FaceColor', lensFace, 'FaceAlpha', 0.45, 'EdgeColor','none');
end

% Glass rectangles (light blue) BEFORE edges/rays
if glass_enabled
    for k = 1:numel(gaps)
        patch([gaps(k).xL gaps(k).xR gaps(k).xR gaps(k).xL], ...
              [gaps(k).h  gaps(k).h -gaps(k).h -gaps(k).h], ...
              glassFace, 'FaceAlpha', 0.25, 'EdgeColor', glassEdge, 'LineWidth', 1.2);
    end
end

% Lens edges
for j = 1:L
    if isempty(final_lens{j}), continue; end
    up = final_lens{j}; lo = up; lo(:,2) = -up(:,2);
    plot(up(:,1), up(:,2), 'Color', edgeColor, 'LineWidth', 1.25);
    if draw_mirror, plot(lo(:,1), lo(:,2), 'Color', edgeColor, 'LineWidth', 1.25); end
end
plot(origin_global(1), origin_global(2), 'ro', 'MarkerFaceColor','r');

%% -------------------------
%  Angle logger + detector-plane y logger (for heatmap)
% -------------------------
alpha_log  = NaN(N_rays, L);
y_end_hits = [];

%% -------------------------
%  Sample launch angles then trace
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
        alpha_log(r, j) = alpha_next;

        % Insert all glass crossings for the straight segment P0 -> P_hit
        if glass_enabled
            path = insert_all_glass_crossings(path, P0, alpha_cur, P_hit, gaps);
        end

        % Append hit point, then refract
        path = [path; P_hit]; %#ok<AGROW>
        P0 = P_hit; alpha_cur = alpha_next;

        % Toggle medium for next surface (same as original)
        tmp = n_in; n_in = n_out; n_out = tmp;
    end

    % Extend to detector plane (only if ray exits all layers & goes forward)
    if hit_all_layers
        Df = [cosd(alpha_cur), sind(alpha_cur)];
        if (x_end > P0(1)) && (abs(Df(1)) > 1e-12)
            t_ext = (x_end - P0(1)) / Df(1);
            if t_ext > 0
                Pend = P0 + t_ext * Df;

                % Insert glass crossings for final free-space segment
                if glass_enabled
                    path = insert_all_glass_crossings(path, P0, alpha_cur, Pend, gaps);
                end

                path = [path; Pend]; %#ok<AGROW>
                y_end_hits(end+1,1) = Pend(2); %#ok<AGROW>
            end
        end
    end

    plot(path(:,1), path(:,2), '-', 'Color', rayColor, 'LineWidth', 0.8);

    if draw_mirror
        path_m = path; path_m(:,2) = -path_m(:,2);
        plot(path_m(:,1), path_m(:,2), '-', 'Color', rayColor, 'LineWidth', 0.8);
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

%% ============================================================
%  EXTRA FIGURE: CLOSE-UP VIEW AROUND THE SOURCE (±5 mm window)
% ============================================================
figure('Name','Ray Tracing – Close Up');
clf; hold on; grid on; axis equal;
set(gcf, 'Color', 'w');
xlabel('X (mm)');
ylabel('Y (mm)');
title('Ray Tracing – Close-Up View');

if ~isempty(lensBody.Vertices)
    plot(lensBody, 'FaceColor', lensFace, 'FaceAlpha', 0.45, 'EdgeColor','none');
end

% Glass rectangles again
if glass_enabled
    for k = 1:numel(gaps)
        patch([gaps(k).xL gaps(k).xR gaps(k).xR gaps(k).xL], ...
              [gaps(k).h  gaps(k).h -gaps(k).h -gaps(k).h], ...
              glassFace, 'FaceAlpha', 0.25, 'EdgeColor', glassEdge, 'LineWidth', 1.2);
    end
end

for j = 1:L
    if isempty(final_lens{j}), continue; end
    up = final_lens{j}; lo = up; lo(:,2) = -up(:,2);
    plot(up(:,1), up(:,2), 'Color', edgeColor, 'LineWidth', 1);
    if draw_mirror, plot(lo(:,1), lo(:,2), 'Color', edgeColor, 'LineWidth', 1); end
end

% Replay rays (same logic)
for r = 1:N_rays
    alpha0    = alphas0(r);
    P0        = origin_global;
    alpha_cur = alpha0;
    n_in      = n1_global;
    n_out     = n2_global;

    path = P0; hit_all = true;

    for j = 1:L
        [hit, P_hit, alpha_next] = first_hit_with_layer(P0, alpha_cur, j, n_in, n_out, final_lens);
        if ~hit, hit_all = false; break; end

        if glass_enabled
            path = insert_all_glass_crossings(path, P0, alpha_cur, P_hit, gaps);
        end

        path = [path; P_hit]; %#ok<AGROW>
        P0 = P_hit; alpha_cur = alpha_next;
        tmp = n_in; n_in = n_out; n_out = tmp;
    end

    plot(path(:,1), path(:,2), '-', 'Color', rayColor, 'LineWidth', 0.8);

    if draw_mirror
        path_m = path; path_m(:,2) = -path_m(:,2);
        plot(path_m(:,1), path_m(:,2), '-', 'Color', rayColor, 'LineWidth', 0.8);
    end
end

plot(origin_global(1), origin_global(2), 'ro', 'MarkerFaceColor','r');
xlim([-5 5]);
ylim([-5 5]);
grid on; axis square;

%% =========================
%  Local helpers
% =========================

% Resolve gaps list.
% Prefer STEP1's saved glass_info.gaps (already has correct height).
function [gaps, enabled] = resolve_glass_gaps(ENABLE_GLASS, S, final_lens, L, has_points, glass_pos, n3, glass_width)
    enabled = false;
    gaps = struct('xL', {}, 'xR', {}, 'h', {}, 'n', {}, 'element', {});

    if ~ENABLE_GLASS
        return;
    end

    % Preferred: read from STEP 1 output
    if isfield(S,'glass_info') && isstruct(S.glass_info) && isfield(S.glass_info,'enabled') && S.glass_info.enabled ...
       && isfield(S.glass_info,'gaps') && ~isempty(S.glass_info.gaps)

        gi = S.glass_info;
        enabled = true;

        for k = 1:numel(gi.gaps)
            gaps(k).xL      = gi.gaps(k).x_left;
            gaps(k).xR      = gi.gaps(k).x_right;
            gaps(k).h       = gi.gaps(k).h_glass;  % already = 1.3 * height(next element)
            gaps(k).n       = gi.n3;
            gaps(k).element = gi.gaps(k).element;
        end

        % Optional sanity warning
        bad = find(isnan([gaps.h]), 1);
        if ~isempty(bad)
            warning('Some glass heights are NaN in file glass_info.gaps (did STEP 1 complete fully?).');
        end
        return;
    end

    % Fallback inference if file lacks gaps
    glass_pos = unique(glass_pos(:).','sorted');
    max_elem = floor(L/2);
    if any(glass_pos < 2) || any(glass_pos > max_elem)
        warning('glass_pos out of valid range for layers_count; disabling glass fallback.');
        return;
    end

    for k = 1:numel(glass_pos)
        elem = glass_pos(k);
        j_insert = 2*(elem-1);
        if j_insert < 1 || j_insert > L || isempty(final_lens{j_insert})
            continue;
        end

        xL = max(final_lens{j_insert}(:,1));
        xR = xL + glass_width;

        % Approx height (same rule): 1.3 * height(next element) => use layer 2*elem
        j_elem2 = 2*elem;
        if j_elem2 <= numel(final_lens) && ~isempty(final_lens{j_elem2})
            h_elem = max(abs(final_lens{j_elem2}(:,2)));
        else
            h_elem = max(cellfun(@(c) max(abs(c(:,2))), final_lens(has_points)));
        end
        hG = 1.3 * h_elem;

        gaps(end+1).xL = xL; %#ok<AGROW>
        gaps(end).xR   = xR;
        gaps(end).h    = hG;
        gaps(end).n    = n3;
        gaps(end).element = elem;
    end

    enabled = ~isempty(gaps);
end

% Insert crossing points for ALL gaps a segment crosses (in correct order).
% This just adds points on the straight line; no refraction occurs inside gaps.
function path_out = insert_all_glass_crossings(path_in, P_start, alpha, P_end, gaps)
    path_out = path_in;

    D = [cosd(alpha), sind(alpha)];
    if abs(D(1)) < 1e-12
        return; % vertical direction: skip
    end

    x0 = P_start(1);
    y0 = P_start(2);

    % Segment parameter along x
    tEnd = (P_end(1) - x0) / D(1);
    tmin = min(0, tEnd);
    tmax = max(0, tEnd);

    t_list = [];
    P_list = [];

    for k = 1:numel(gaps)
        xL = gaps(k).xL;
        xR = gaps(k).xR;

        % quick reject by x-range overlap
        if (max(P_start(1), P_end(1)) < xL) || (min(P_start(1), P_end(1)) > xR)
            continue;
        end

        tL = (xL - x0) / D(1);
        tR = (xR - x0) / D(1);
        t_in  = min(tL, tR);
        t_out = max(tL, tR);

        if t_out < tmin || t_in > tmax
            continue;
        end

        t_in_cl  = max(t_in,  tmin);
        t_out_cl = min(t_out, tmax);

        if (t_out_cl - t_in_cl) < 1e-9
            continue; % touch-only
        end

        P_entry = [x0, y0] + t_in_cl  * D;
        P_exit  = [x0, y0] + t_out_cl * D;

        t_list = [t_list; t_in_cl; t_out_cl]; %#ok<AGROW>
        P_list = [P_list; P_entry; P_exit]; %#ok<AGROW>
    end

    if isempty(t_list)
        return;
    end

    [t_sorted, idx] = sort(t_list, 'ascend');
    P_sorted = P_list(idx, :);

    % remove near duplicates
    keep = true(size(t_sorted));
    for i = 2:numel(t_sorted)
        if abs(t_sorted(i) - t_sorted(i-1)) < 1e-8 && norm(P_sorted(i,:) - P_sorted(i-1,:)) < 1e-6
            keep(i) = false;
        end
    end
    P_sorted = P_sorted(keep, :);

    path_out = [path_out; P_sorted]; %#ok<AGROW>
end

% --- Refraction helper (unchanged) ---
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

% --- Find first intersection with layer j ---
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

% --- Emission sampler ---
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
