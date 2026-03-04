%% ============================================================
%  STEP 2: RAY TRACING THROUGH THE CONSTRUCTED LENS (solid-body fill + quantitative heatmap)
%  ============================================================
%clear; clc; close all;

%% -------------------------
%  User Parameters
%  -------------------------
input_filename = 'lens_data_close_up.mat';
N_rays         = 10;               % adjust as needed
seed           = 42;
draw_mirror    = true;
rng(seed);

% Colors
rayColor  = [0.85 0.65 0.13];      % dark yellow for rays
lensFace  = [0.68 0.92 0.68];      % slightly darker green
edgeColor = [0 0 0];

%% -------------------------
%  Refraction helper (unchanged)
%  -------------------------
function [alpha_new, intersection_point, theta_in, theta_out] = refract_ray_2d(origin_a, alpha, surface_a, surface_b, n1, n2)
D = [cosd(alpha), sind(alpha)]; O = origin_a(:).';
A = surface_a(:).'; B = surface_b(:).'; S = B - A;
M = [D(:), -S(:)]; rhs = (A - O).'; detM = det(M);
if abs(detM) < 1e-12, error('Ray and segment are parallel or nearly parallel; no unique intersection.'); end
sol = M \ rhs; t = sol(1); u = sol(2);
if t < -1e-12 || u < -1e-12 || u > 1+1e-12, error('No valid intersection between the ray and the segment.'); end
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

%% -------------------------
%  Load lens geometry
%  -------------------------
S = load(input_filename);
required = {'final_lens','n1','n2','layers_count','source_pos'};
for f = required
    if ~isfield(S, f{1})
        error('Missing "%s" in %s. Re-run the construction step.', f{1}, input_filename);
    end
end
final_lens    = S.final_lens;
n1_global     = S.n1;
n2_global     = S.n2;
L             = S.layers_count;
origin_global = S.source_pos(:).';

% Ensure we have usable layers
L = min(L, numel(final_lens));
has_points = cellfun(@(c) ~isempty(c) && size(c,1) >= 2, final_lens);
if ~any(has_points), error('No usable layers found in final_lens.'); end

% Common x_end = xmax + 30% span
allX = [];
for j = 1:L
    if isempty(final_lens{j}), continue; end
    allX = [allX; final_lens{j}(:,1)];
end
xmin = min(allX); xmax = max(allX);
x_span = max(1e-9, xmax - xmin);
x_end  = xmax + 0.30 * x_span;

%% -------------------------
%  Build a *solid-body* polyshape of the lens (no leaking)
%  -------------------------
layerPoly = cell(L,1);
for j = 1:L
    if isempty(final_lens{j}), continue; end
    up = final_lens{j};
    lo = up; lo(:,2) = -up(:,2);
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
%  Plot setup
%  -------------------------
figure('Name','Ray Tracing Through Constructed Lens'); clf; hold on; grid on; axis equal;
set(gcf, 'Color', 'w'); xlabel('X'); ylabel('Y');
title(sprintf('Ray Tracing with sin^2(\\alpha) acceptance (N = %d)', N_rays));

if ~isempty(lensBody.Vertices)
    plot(lensBody, 'FaceColor', lensFace, 'FaceAlpha', 0.45, 'EdgeColor','none');
end

for j = 1:L
    if isempty(final_lens{j}), continue; end
    up = final_lens{j};
    lo = up; lo(:,2) = -up(:,2);
    plot(up(:,1), up(:,2), 'Color', edgeColor, 'LineWidth', 1.25);
    if draw_mirror
        plot(lo(:,1), lo(:,2), 'Color', edgeColor, 'LineWidth', 1.25);
    end
end
plot(origin_global(1), origin_global(2), 'ro', 'MarkerFaceColor','r');

%% -------------------------
%  Angle logger + detector-plane y logger (for heatmap)
%  -------------------------
alpha_log  = NaN(N_rays, L);
y_end_hits = [];

%% -------------------------
%  Sample rays with sin^2(alpha) acceptance and trace
%  -------------------------
accepted = 0; attempts = 0; max_attempts = 1e7;

while accepted < N_rays && attempts < max_attempts
    attempts = attempts + 1;

    alpha0 = 180 * rand();
    p      = min(max(sind(alpha0)^2,0),1);
    if rand() >= p, continue; end

    P0        = origin_global;
    alpha_cur = alpha0;
    n_in      = n1_global;
    n_out     = n2_global;

    path = P0;
    hit_all_layers = true;
    row = accepted + 1;

    for j = 1:L
        [hit, P_hit, alpha_next] = first_hit_with_layer( ...
            P0, alpha_cur, j, n_in, n_out, final_lens);
        if ~hit, hit_all_layers = false; break; end

        alpha_log(row, j) = alpha_next;

        path      = [path; P_hit];
        P0        = P_hit;
        alpha_cur = alpha_next;
        tmp = n_in; n_in = n_out; n_out = tmp;
    end

    if hit_all_layers
        Df = [cosd(alpha_cur), sind(alpha_cur)];
        if (x_end > P0(1)) && (abs(Df(1)) > 1e-12)
            t_ext = (x_end - P0(1)) / Df(1);
            if t_ext > 0
                Pend = P0 + t_ext * Df;
                path = [path; Pend];
                y_end_hits(end+1,1) = Pend(2);
            end
        end
        plot(path(:,1), path(:,2), '-', 'Color', rayColor, 'LineWidth', 0.8);
        if draw_mirror
            path_m = path; path_m(:,2) = -path_m(:,2);
            plot(path_m(:,1), path_m(:,2), '-', 'Color', rayColor, 'LineWidth', 0.8);
        end
    end

    accepted = accepted + 1;
end

if ~isempty(lensBody.Vertices)
    y_min_plot = min(lensBody.Vertices(:,2));
    y_max_plot = max(lensBody.Vertices(:,2));
else
    y_all_layers = cell2mat(cellfun(@(c) [c(:,2); -c(:,2)], final_lens(~cellfun('isempty',final_lens)), 'uni', 0));
    y_min_plot = min(y_all_layers);
    y_max_plot = max(y_all_layers);
end
y_lim = max(abs([y_min_plot, y_max_plot]));
ylim([-y_lim, y_lim]);
xlim([xmin - 0.1*x_span, x_end + 0.15*x_span]);

fprintf('Accepted rays: %d (attempts: %d), x_end = %.3f\n', accepted, attempts, x_end);
legend({'Lens fill (solid body)','Lens edges','Source'}, 'Location','best');

%% -------------------------
%  Output vertical heatmap with quantitative colorbar
%  -------------------------
if ~isempty(y_end_hits)
    % mirrored data
    y_det = [y_end_hits(:); -y_end_hits(:)];

    yl = ylim; y_min_hm = yl(1); y_max_hm = yl(2);

    nbins   = 200;
    edges   = linspace(y_min_hm, y_max_hm, nbins+1);
    centers = (edges(1:end-1) + edges(2:end))/2;
    dy      = edges(2) - edges(1);

    % density = rays per unit Y
    counts   = histcounts(y_det, edges).';
    density  = counts / max(dy, eps);
    maxD     = max(density);
    if maxD == 0, maxD = 1; end

    C = [density density];
    bar_w = 0.03 * (xmax - xmin + eps);
    x0    = x_end;

    imagesc('XData',[x0, x0+bar_w], 'YData',[centers(1), centers(end)], 'CData', C);
    set(gca,'YDir','normal');

    ncm   = 256;
    yellow= [1.00, 1.00, 0.60];
    red   = [0.55, 0.00, 0.00];
    cmap  = [linspace(yellow(1), red(1), ncm)', ...
             linspace(yellow(2), red(2), ncm)', ...
             linspace(yellow(3), red(3), ncm)'];
    colormap(gca, cmap);
    caxis([0 maxD]);
    alpha(0.95);

    % border
    plot([x_end x_end],[y_min_hm y_max_hm],'k-','LineWidth',0.8);
    plot([x0 x0+bar_w x0+bar_w x0 x0],[y_min_hm y_min_hm y_max_hm y_max_hm y_min_hm], ...
         'k-','LineWidth',0.75);

    % quantitative colorbar
    cb = colorbar('eastoutside');
    cb.Label.String = sprintf('Ray density (rays / unit Y)\\nΔy = %.3g', dy);
    cb.TickDirection = 'out';
end

%% =========================
%  Local helper
%  =========================
function [hit, P_hit, alpha_out] = first_hit_with_layer(P0, alpha_in, j, n_in, n_out, final_lens)
    hit = false; P_hit = [NaN, NaN]; alpha_out = NaN;
    if j > numel(final_lens) || isempty(final_lens{j}) || size(final_lens{j},1) < 2, return; end
    pts = final_lens{j};
    best_d = inf; best_P = []; best_alpha = NaN;
    for k = 1:size(pts,1)-1
        A = pts(k,   :); B = pts(k+1, :);
        try
            [alpha_new, P_inter] = refract_ray_2d(P0, alpha_in, A, B, n_in, n_out);
        catch
            continue;
        end
        if any(isnan(P_inter)) || any(isinf(P_inter)) || isnan(alpha_new), continue; end
        d = hypot(P_inter(1)-P0(1), P_inter(2)-P0(2));
        if d > 1e-9 && d < best_d
            best_d = d; best_P = P_inter; best_alpha = alpha_new;
        end
    end
    if ~isempty(best_P), hit = true; P_hit = best_P; alpha_out = best_alpha; end
end
