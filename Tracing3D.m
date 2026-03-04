%% ============================================================
%  STEP 4: 3-D RAY TRACING THROUGH REVOLVED LENS
%  Publication-quality figure for research paper
%  ============================================================
clear; clc; close all;

%% ── User parameters ─────────────────────────────────────────
mat_filename     = 'lens_data_gradual_N2003_test.mat';
N_rays           = 10;
N_phi            = 8;
emission_profile = 'uniform';
seed             = 42;
rng(seed);

%% ── Color palette ───────────────────────────────────────────
RAY_COLOR  = [1.00  0.65  0.00];   % deep amber

%% ── 1) Load data ────────────────────────────────────────────
S = load(mat_filename);
req = {'final_lens','n1','n2','layers_count','source_pos'};
for k = 1:numel(req)
    if ~isfield(S, req{k})
        error('Missing field "%s" in %s', req{k}, mat_filename);
    end
end

final_lens    = S.final_lens;
n1_global     = S.n1;
n2_global     = S.n2;
L             = S.layers_count;
origin_global = S.source_pos(:).';

% X-axis extent
has_pts = cellfun(@(c) ~isempty(c) && size(c,1)>=2, final_lens);
allX    = cell2mat(cellfun(@(c) c(:,1), final_lens(has_pts), 'uni', 0));
xmin    = min(allX);   xmax = max(allX);
xspan   = max(1e-9, xmax - xmin);
x_end   = xmax + 0.55 * xspan;

%% ── 2) Figure & axes ────────────────────────────────────────
fig = figure('Name','3D Ray Tracing — GRIN Lens', ...
             'Units','centimeters', ...
            'Position',[2 2 24 20]);
set(fig, 'Color', 'w');

ax = axes('Parent', fig);
hold(ax, 'on');
grid(ax, 'on');
axis(ax, 'equal');
set(ax, ...
    'Color',         'w', ...
    'XColor',        'k', ...
    'YColor',        'k', ...
    'ZColor',        'k', ...
    'GridColor',     [0.75 0.75 0.75], ...
    'GridAlpha',     0.5, ...
    'GridLineStyle', ':', ...
    'FontSize',      10, ...
    'FontName',      'Times New Roman', ...
    'LineWidth',     0.8, ...
    'TickDir',       'out', ...
    'Box',           'off');

xlabel(ax, '\itX\rm (mm)', 'FontSize',11, 'FontName','Times New Roman');
ylabel(ax, '\itY\rm (mm)', 'FontSize',11, 'FontName','Times New Roman');
zlabel(ax, '\itZ\rm (mm)', 'FontSize',11, 'FontName','Times New Roman');
title(ax, '3-D Ray Tracing Through Gradient-Index Lens', ...
          'FontSize',12, 'FontName','Times New Roman', 'FontWeight','bold');

view(ax, 18, 14);
ax.CameraViewAngle = 9;

%% ── 3) Discrete layer colormap (one distinct color per layer) ──
% Each layer gets its own unambiguous color so reader can count them
layer_colors = [0.15 0.35 0.80;   % Layer 1 — deep blue
                0.05 0.70 0.75;   % Layer 2 — cyan
                0.10 0.70 0.30;   % Layer 3 — green
                0.85 0.50 0.05;   % Layer 4 — amber
                0.75 0.15 0.15;   % Layer 5 — red
                0.55 0.10 0.70;   % Layer 6 — purple
                0.10 0.50 0.50;   % Layer 7 — teal
                0.40 0.40 0.40];  % Layer 8 — grey
% Truncate or extend to actual number of layers
if L <= size(layer_colors,1)
    cmap_lens = layer_colors(1:L, :);
else
    cmap_lens = [layer_colors; parula(L - size(layer_colors,1))];
end

%% ── 4) Revolve lens surfaces ────────────────────────────────
n_theta   = 150;
theta_rev = linspace(0, 2*pi, n_theta);

for j = 1:L
    curve = final_lens{j};
    if isempty(curve), continue; end
    curve = curve(curve(:,2) >= 0, :);
    if isempty(curve), continue; end

    x = curve(:,1);
    r = curve(:,2);
    [T, Xg] = meshgrid(theta_rev, x);
    R = repmat(r, 1, n_theta);

    % Inner layers more opaque so core is visible; outer layers ghost-like
    alpha_val = 0.05 + 0.12 * (j / L);

    surf(ax, Xg, R.*cos(T), R.*sin(T), ...
        'FaceColor', cmap_lens(j,:), ...
        'FaceAlpha', alpha_val, ...
        'EdgeColor', cmap_lens(j,:) * 0.6, ...
        'EdgeAlpha', 0.12, ...
        'LineWidth',  0.3, ...
        'AmbientStrength',  0.40, ...
        'DiffuseStrength',  0.70, ...
        'SpecularStrength', 0.30);
end

%% ── 5) 2-D ray tracing (meridional plane) ───────────────────
alphas0 = sample_ray_angles(N_rays, emission_profile, seed);
paths2D = cell(N_rays, 1);

for ri = 1:N_rays
    alpha_cur = alphas0(ri);
    P0        = origin_global;
    n_in      = n1_global;
    n_out     = n2_global;
    path      = P0;
    survived  = true;

    for j = 1:L
        [hit, P_hit, alpha_next] = first_hit_with_layer( ...
            P0, alpha_cur, j, n_in, n_out, final_lens);
        if ~hit, survived = false; break; end
        path      = [path; P_hit]; %#ok<AGROW>
        P0        = P_hit;
        alpha_cur = alpha_next;
        tmp = n_in; n_in = n_out; n_out = tmp;
    end

    if survived
        D = [cosd(alpha_cur), sind(alpha_cur)];
        if (x_end > P0(1)) && (abs(D(1)) > 1e-12)
            t_ext = (x_end - P0(1)) / D(1);
            if t_ext > 0
                path = [path; P0 + t_ext*D]; %#ok<AGROW>
            end
        end
    end
    paths2D{ri} = path;
end

%% ── 6) Revolve rays and plot ────────────────────────────────
phi_vals = linspace(0, 360, N_phi+1);
phi_vals(end) = [];

for ri = 1:N_rays
    path = paths2D{ri};
    if isempty(path), continue; end
    x2d = path(:,1);
    r2d = abs(path(:,2));

    for p = 1:numel(phi_vals)
        ph = phi_vals(p);
        if ph == 0
            ray_col = RAY_COLOR;           % fully opaque meridional ray
            lw      = 1.8;
        else
            ray_col = [RAY_COLOR, 0.25];   % transparent azimuthal copies
            lw      = 0.35;
        end
        Y3 = r2d * cosd(ph);
        Z3 = r2d * sind(ph);
        plot3(ax, x2d, Y3, Z3, '-', 'Color', ray_col, 'LineWidth', lw);
    end
end

%% ── 7) Source marker — single occurrence at true origin ─────
plot3(ax, origin_global(1), origin_global(2), 0, ...
      'p', ...
      'MarkerSize',      14, ...
      'MarkerFaceColor', [1.0 0.85 0.0], ...
      'MarkerEdgeColor', [0.6 0.35 0.0], ...
      'LineWidth',       1.2);
text(ax, origin_global(1) + 2, origin_global(2), 8, ...
     'Source', ...
     'FontName',   'Times New Roman', ...
     'FontSize',   10, ...
     'FontWeight', 'bold', ...
     'Color',      [0.5 0.25 0.0]);

%% ── 8) Colorbar ─────────────────────────────────────────────
colormap(ax, cmap_lens);
cb = colorbar(ax, 'eastoutside');
cb.Label.String   = 'Layer index (inner \rightarrow outer)';
cb.Label.FontName = 'Times New Roman';
cb.Label.FontSize = 9;
cb.FontName       = 'Times New Roman';
cb.FontSize       = 8;
cb.TickDirection  = 'out';
clim(ax, [1 L]);

%% ── 9) Axis limits — manual control to show full lens ───────

% First let MATLAB render everything, then expand limits
axis(ax, 'tight');          % auto-fit to data first
pause(0.01);                % allow render

% Then manually pad all six directions
xl = xlim(ax); yl = ylim(ax); zl = zlim(ax);
xlim(ax, [xl(1) - 0.15*(xl(2)-xl(1)),  xl(2) + 0.15*(xl(2)-xl(1))]);
ylim(ax, [yl(1) - 0.20*(yl(2)-yl(1)),  yl(2) + 0.20*(yl(2)-yl(1))]);
zlim(ax, [zl(1) - 0.25*(zl(2)-zl(1)),  zl(2) + 0.25*(zl(2)-zl(1))]);  % <-- extra Z room

% Remove axis equal — it fights with manual limits in 3D
% axis(ax, 'equal');   % <-- comment this out or delete it

%% ── 10) Export ──────────────────────────────────────────────
exportgraphics(fig, 'ray_tracing_3D.tiff', ...
               'Resolution', 300, 'BackgroundColor', 'white');
exportgraphics(fig, 'ray_tracing_3D.pdf', ...
               'ContentType', 'vector', 'BackgroundColor', 'white');

fprintf('Exported: ray_tracing_3D.tiff + .pdf\n');
fprintf('Rays: N_rays=%d, N_phi=%d, profile=%s\n', N_rays, N_phi, emission_profile);

%% ============================================================
%  LOCAL FUNCTIONS
% ============================================================
function [alpha_new, P_inter, theta_in, theta_out] = ...
        refract_ray_2d(O, alpha, A, B, n1, n2)
    D  = [cosd(alpha), sind(alpha)];
    O  = O(:).'; A = A(:).'; B = B(:).';
    S  = B - A;
    M  = [D(:), -S(:)];
    if abs(det(M)) < 1e-12, error('Ray parallel to segment.'); end
    sol = M \ (A - O).';
    t = sol(1); u = sol(2);
    if t < -1e-12 || u < -1e-12 || u > 1+1e-12
        error('No valid intersection.');
    end
    t = max(t,0); u = min(max(u,0),1);
    P_inter = O + t*D;

    t_hat     = S / norm(S);
    n_hat_ccw = [-t_hat(2),  t_hat(1)];
    n_hat_cw  = [ t_hat(2), -t_hat(1)];
    if dot(D, n_hat_ccw) >= dot(D, n_hat_cw)
        n_hat = n_hat_ccw;
    else
        n_hat = n_hat_cw;
    end

    angN     = atan2d(n_hat(2), n_hat(1));
    angRay   = atan2d(D(2),     D(1));
    theta_in = wrapTo180(angRay - angN);
    if abs(theta_in) > 90
        angN     = wrapTo180(angN + 180);
        theta_in = wrapTo180(angRay - angN);
    end

    sin_out = (n1/n2) * sind(theta_in);
    if abs(sin_out) > 1 + 1e-12
        alpha_new = NaN; theta_out = NaN;
        warning('Total internal reflection.');
        return;
    end
    sin_out   = max(min(sin_out,1),-1);
    theta_out = asind(sin_out);
    alpha_new = wrapTo360(angN + theta_out);
end

function a = wrapTo180(a)
    a = mod(a + 180, 360) - 180;
end

function a = wrapTo360(a)
    a = mod(a, 360);
    if a < 0, a = a + 360; end
end

function [hit, P_hit, alpha_out] = ...
        first_hit_with_layer(P0, alpha_in, j, n_in, n_out, final_lens)
    hit = false; P_hit = [NaN, NaN]; alpha_out = NaN;
    if j > numel(final_lens) || isempty(final_lens{j}) || ...
            size(final_lens{j},1) < 2
        return;
    end
    pts    = final_lens{j};
    best_d = inf;
    for k = 1:size(pts,1)-1
        try
            [a_new, P_i] = refract_ray_2d( ...
                P0, alpha_in, pts(k,:), pts(k+1,:), n_in, n_out);
        catch
            continue;
        end
        if any(isnan(P_i)) || isnan(a_new), continue; end
        d = hypot(P_i(1)-P0(1), P_i(2)-P0(2));
        if d > 1e-9 && d < best_d
            best_d    = d;
            P_hit     = P_i;
            alpha_out = a_new;
            hit       = true;
        end
    end
end

function alphas = sample_ray_angles(N, profile, seed)
    if nargin >= 3 && ~isempty(seed), rng(seed); end
    switch lower(strrep(profile,'^',''))
        case 'uniform'
            alphas = linspace(0, 180, N).';
        case 'sin2'
            alphas = zeros(N,1); k = 0;
            while k < N
                M   = max(1000, 2*(N-k));
                c   = 180*rand(M,1); u = rand(M,1);
                acc = u <= sind(c).^2;
                n_acc = sum(acc);
                if n_acc > 0
                    take = min(n_acc, N-k);
                    idx  = find(acc, take, 'first');
                    alphas(k+1:k+take) = c(idx);
                    k = k + take;
                end
            end
        otherwise
            error('Unknown profile: "%s". Use "uniform" or "sin2".', profile);
    end
end