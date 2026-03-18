%% ============================================================
%  LIVE PLOT GENERATOR
%  (10/20/30 rays) x (Blue/DarkYellow/ScientificBlue) x (Regular/Zoom) x (Grid/No Grid)
% ============================================================
clear; clc; close all;

%% -------------------------
%  User Parameters
% -------------------------
input_filename   = 'lens_data_gradual_N2003.mat';  % <- your file
seed             = 42;
draw_mirror      = true;
emission_profile = 'uniform';
rng(seed);

uniform_min_angle = 0;
uniform_max_angle = 178;
uniform_eps       = 1e-6;

RAY_EMPHASIS_MODE = false;

%% -------------------------
%  Variation Switches  (set to false to skip)
% -------------------------
% SVG export
EXPORT_SVG     = true;

% Ray counts
USE_N10        = true;
USE_N20        = false;
USE_N30        = false;

% Colors
USE_BLUE           = false;
USE_DARKYELLOW     = false;
USE_SCIENTIFICBLUE = true;   % light-blue lens + deep-navy rays, scientific paper style

% View modes
USE_REGULAR     = true;
USE_ZOOM        = true;

% Grid modes
USE_NOGRID      = true;
USE_WITHGRID    = false;

%% -------------------------
%  Batch Styling Configurations
% -------------------------
% 1. Ray counts
n_rays_list = [10, 20, 30];
n_rays_list = n_rays_list([USE_N10, USE_N20, USE_N30]);

% 2. Colors
colors.blue            = [0.15 0.35 0.75];
colors.darkyellow      = [0.75 0.55 0.00];
colors.scientificblue  = [0.08 0.28 0.62];   % deep blue for rays in scientific style
all_color_names        = fieldnames(colors);
color_names            = all_color_names([USE_BLUE, USE_DARKYELLOW, USE_SCIENTIFICBLUE]);

% 3. Views and Grids
view_modes = {'regular', 'zoom'};
grid_modes = {'nogrid', 'withgrid'};
view_modes = view_modes([USE_REGULAR, USE_ZOOM]);
grid_modes = grid_modes([USE_NOGRID, USE_WITHGRID]);

% Lens / Heatmap styles based on emphasis mode
if RAY_EMPHASIS_MODE
    rayLW = 0.9;
    lensFace  = [0.68 0.92 0.68]; lensAlpha = 0.015;
    edgeColor = [0.80 0.80 0.80]; edgeLW    = 0.35;
    ENABLE_HEATMAP = false;
    EXTEND_TO_RIGHT_ALWAYS = true;
else
    rayLW = 0.65;
    lensFace  = [0.68 0.92 0.68]; lensAlpha = 0.45;
    edgeColor = [0 0 0];          edgeLW    = 1.5;
    ENABLE_HEATMAP = false;
    EXTEND_TO_RIGHT_ALWAYS = true;
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

L = min(L, numel(final_lens));
has_points = cellfun(@(c) ~isempty(c) && size(c,1) >= 2, final_lens);
if ~any(has_points), error('No usable layers found in final_lens.'); end

allX   = cell2mat(cellfun(@(c) c(:,1), final_lens(has_points), 'uni', 0));
xmin   = min(allX);
xmax   = max(allX);
x_span = max(1e-9, xmax - xmin);
x_end  = xmax + 0.30 * x_span;

%% -------------------------
%  Build solid-body lens polyshape
% -------------------------

% Precompute tip intersection points for each outer element (e >= 2)
% so they can be baked into the layerPoly rings (fill matches outline).
n_elements = floor(L/2);
tip_pts = cell(n_elements, 1);   % tip_pts{e} = [x, y] of crescent tip (upper side)
for e = 2:n_elements
    li  = 2*e - 1;  loe = 2*e;
    if li > L || loe > L || isempty(final_lens{li}) || isempty(final_lens{loe}), continue; end
    if size(final_lens{li},1) < 2 || size(final_lens{loe},1) < 2, continue; end
    p1 = final_lens{li}(end,:);
    a1 = atan2d(p1(2)-final_lens{li}(end-1,2),  p1(1)-final_lens{li}(end-1,1));
    p2 = final_lens{loe}(end,:);
    a2 = atan2d(p2(2)-final_lens{loe}(end-1,2), p2(1)-final_lens{loe}(end-1,1));
    tpt = tip_intersection(p1, a1, p2, a2);
    if ~any(isnan(tpt)), tip_pts{e} = tpt; end
end

layerPoly = cell(L,1);
for j = 1:L
    if isempty(final_lens{j}), continue; end
    up = final_lens{j}; lo = up; lo(:,2) = -up(:,2);
    e  = ceil(j/2);
    if e >= 2 && e <= n_elements && ~isempty(tip_pts{e})
        tpt  = tip_pts{e};
        tptm = [tpt(1), -tpt(2)];           % mirrored tip
        ring = [up; tpt; tptm; flipud(lo)]; % V-shaped tip baked in
    else
        ring = [up; flipud(lo)];
    end
    layerPoly{j} = polyshape(ring(:,1), ring(:,2), 'Simplify', true);
end

% Precompute each glass shell (annular region between adjacent layer pairs)
glassShells = cell(L,1);
lensBody = polyshape;
for j = 2:2:L
    if isempty(layerPoly{j-1}) || isempty(layerPoly{j}), continue; end
    shell = subtract(layerPoly{j}, layerPoly{j-1});
    glassShells{j} = shell;
    lensBody = union(lensBody, shell);
end

if ~isempty(lensBody.Vertices)
    y_min_plot = min(lensBody.Vertices(:,2));
    y_max_plot = max(lensBody.Vertices(:,2));
else
    y_all_layers = cell2mat(cellfun(@(c) [c(:,2); -c(:,2)], final_lens(has_points), 'uni', 0));
    y_min_plot = min(y_all_layers); y_max_plot = max(y_all_layers);
end
y_lim_base = max(abs([y_min_plot, y_max_plot]));

%% -------------------------
%%  Export lens geometry to SVG
%% -------------------------
if EXPORT_SVG
    fig_svg = figure('Color','w','Units','centimeters','Position',[2 2 13 9]);
    hold on; axis equal;
    set(gca, 'DataAspectRatio', [1 1 1]);
    if ~isempty(lensBody.Vertices)
        plot(lensBody, 'FaceColor', lensFace, 'FaceAlpha', lensAlpha, 'EdgeColor', 'none');
    end
    for j = 1:L
        if isempty(final_lens{j}), continue; end
        up = final_lens{j}; lo = up; lo(:,2) = -up(:,2);
        plot(up(:,1), up(:,2), 'Color', edgeColor, 'LineWidth', edgeLW);
        plot(lo(:,1), lo(:,2), 'Color', edgeColor, 'LineWidth', edgeLW);
    end
    plot(origin_global(1), origin_global(2), 'ro', 'MarkerFaceColor', 'r', 'MarkerSize', 5);
    xlim([xmin - 0.1*x_span, x_end + 0.15*x_span]);
    ylim([-y_lim_base, y_lim_base]);
    xlabel('$x$ (mm)', 'Interpreter', 'latex');
    ylabel('$y$ (mm)', 'Interpreter', 'latex');
    set(gca, 'FontSize', 9, 'FontName', 'Helvetica', 'TickDir', 'out', 'Box', 'on', 'LineWidth', 0.75);
    set(fig_svg, 'Renderer', 'painters');
    print(fig_svg, 'lens_geometry.svg', '-dsvg');
    close(fig_svg);
    fprintf('Lens geometry exported to lens_geometry.svg\n');
end



%% ============================================================
%  START LIVE PLOTTING
% ============================================================

total_figs = length(n_rays_list) * length(color_names) * length(view_modes) * length(grid_modes);
current_fig = 1;

fprintf('Generating %d live figures...\n', total_figs);

for n_idx = 1:length(n_rays_list)
    N_rays = n_rays_list(n_idx);
    
    % --- 1. Pre-calculate Rays for this N_rays ---
    alphas0 = sample_ray_angles(N_rays, emission_profile, seed, ...
                                uniform_min_angle, uniform_max_angle, uniform_eps);
    all_paths = cell(N_rays, 1);
    
    for r = 1:N_rays
        alpha0    = alphas0(r); P0 = origin_global;
        alpha_cur = alpha0;     n_in = n1_global; n_out = n2_global;
        
        path = P0;
        hit_all_layers = true;

        for j = 1:L
            [hit, P_hit, alpha_next] = first_hit_with_layer(P0, alpha_cur, j, n_in, n_out, final_lens);
            if ~hit, hit_all_layers = false; break; end
            path = [path; P_hit]; %#ok<AGROW>
            P0 = P_hit; alpha_cur = alpha_next;
            tmp = n_in; n_in = n_out; n_out = tmp;
        end

        % Always extend in final direction so every ray reaches the plot edge
        Df = [cosd(alpha_cur), sind(alpha_cur)];
        if abs(Df(1)) > 1e-12
            t_ext = (x_end - P0(1)) / Df(1);
            Pend = P0 + t_ext * Df; path = [path; Pend]; %#ok<AGROW>
        end
        all_paths{r} = path;
    end
    
    % --- 2. Iterate through Visual Styles ---
    for c_idx = 1:length(color_names)
        cName = color_names{c_idx};
        current_ray_color = colors.(cName);
        
        for v_idx = 1:length(view_modes)
            vName = view_modes{v_idx};
            
            for g_idx = 1:length(grid_modes)
                gName = grid_modes{g_idx};
                
                % Create visible figure
                fig_name = sprintf('N=%d | %s | %s | %s', N_rays, cName, vName, gName);
                fig = figure('Name', fig_name, 'NumberTitle', 'off'); 
                hold on; axis equal;
                set(fig, 'Units', 'centimeters', 'Position', [2 2 13 9], 'Color', 'w');
                set(gca,'DataAspectRatio',[1 1 1]);
                set(gca,'PlotBoxAspectRatio',[1 1 1]);
                
                % Per-figure lens and ray style
                is_sci = strcmp(cName, 'scientificblue');
                if is_sci
                    cur_facecol   = [0.72 0.88 0.94];   % light cyan-blue
                    cur_facealpha = 0.72;
                    cur_edgelw    = 1.0;
                    cur_glowcol   = [0.86 0.94 0.97];   % pale cyan glow
                    cur_raylw     = 0.70;
                else
                    cur_facecol   = lensFace;
                    cur_facealpha = lensAlpha;
                    cur_edgelw    = edgeLW;
                    cur_raylw     = rayLW;
                end

                % Draw each glass shell — fill only (EdgeColor='none' avoids bridge lines)
                for j = 2:2:L
                    if j == 2
                        % Element 1 (innermost): explicit patch polygon so the flat-tip
                        % diagonal closure matches the outline exactly.
                        % Polygon traces: outer curve → flat tip (top) → inner curve reversed
                        %                 → flat tip (bottom) → outer curve reversed (bottom)
                        if isempty(final_lens{1}) || isempty(final_lens{2}), continue; end
                        up1 = final_lens{1}; lo1 = [up1(:,1), -up1(:,2)];
                        up2 = final_lens{2}; lo2 = [up2(:,1), -up2(:,2)];
                        verts = [up2; flipud(up1); lo1; flipud(lo2)];
                        patch(verts(:,1), verts(:,2), cur_facecol, ...
                              'FaceAlpha', cur_facealpha, 'EdgeColor', 'none');
                        if is_sci
                            patch(verts(:,1), verts(:,2), [0.52 0.78 0.91], ...
                                  'FaceAlpha', 0.25, 'EdgeColor', 'none');
                        end
                    else
                        % Elements 2+: polyshape with tip extensions baked in
                        if isempty(glassShells{j}) || isempty(glassShells{j}.Vertices), continue; end
                        sh = glassShells{j};
                        plot(sh, 'FaceColor', cur_facecol, 'FaceAlpha', cur_facealpha, 'EdgeColor', 'none');
                        if is_sci
                            plot(sh, 'FaceColor', [0.52 0.78 0.91], 'FaceAlpha', 0.25, 'EdgeColor', 'none');
                        end
                    end
                end

                % Edge color (used in both edge loop and tip closure below)
                if is_sci, ec = [0 0 0]; else, ec = edgeColor; end

                % White gap width under each black edge — creates visible separation
                % between adjacent filled shells regardless of fill opacity.
                gapLW = cur_edgelw + 3.0;

                % Draw lens edges (white gap first, then black on top)
                for j = 1:L
                    if isempty(final_lens{j}), continue; end
                    up = final_lens{j}; lo = up; lo(:,2) = -up(:,2);
                    % White gap
                    plot(up(:,1), up(:,2), 'Color', 'w', 'LineWidth', gapLW);
                    if draw_mirror
                        plot(lo(:,1), lo(:,2), 'Color', 'w', 'LineWidth', gapLW);
                    end
                end
                for j = 1:L
                    if isempty(final_lens{j}), continue; end
                    up = final_lens{j}; lo = up; lo(:,2) = -up(:,2);
                    if is_sci
                        plot(up(:,1), up(:,2), 'Color', cur_glowcol, 'LineWidth', cur_edgelw + 1.8);
                        plot(up(:,1), up(:,2), 'Color', ec, 'LineWidth', cur_edgelw);
                        if draw_mirror
                            plot(lo(:,1), lo(:,2), 'Color', cur_glowcol, 'LineWidth', cur_edgelw + 1.8);
                            plot(lo(:,1), lo(:,2), 'Color', ec, 'LineWidth', cur_edgelw);
                        end
                    else
                        plot(up(:,1), up(:,2), 'Color', ec, 'LineWidth', cur_edgelw);
                        if draw_mirror
                            plot(lo(:,1), lo(:,2), 'Color', ec, 'LineWidth', cur_edgelw);
                        end
                    end
                end

                % Close crescent tips (same method as LensConstructionGradual)
                % Element 1 (layers 1&2): flat/truncated tip
                if L >= 2 && ~isempty(final_lens{1}) && ~isempty(final_lens{2})
                    p1t = final_lens{1}(end,:);
                    p2t = final_lens{2}(end,:);
                    plot([p1t(1),p2t(1)],[p1t(2),p2t(2)],'Color','w','LineWidth',gapLW);
                    if is_sci
                        plot([p1t(1),p2t(1)],[p1t(2),p2t(2)],'Color',cur_glowcol,'LineWidth',cur_edgelw+1.8);
                    end
                    plot([p1t(1),p2t(1)],[p1t(2),p2t(2)],'Color',ec,'LineWidth',cur_edgelw);
                    if draw_mirror
                        plot([p1t(1),p2t(1)],[-p1t(2),-p2t(2)],'Color','w','LineWidth',gapLW);
                        if is_sci
                            plot([p1t(1),p2t(1)],[-p1t(2),-p2t(2)],'Color',cur_glowcol,'LineWidth',cur_edgelw+1.8);
                        end
                        plot([p1t(1),p2t(1)],[-p1t(2),-p2t(2)],'Color',ec,'LineWidth',cur_edgelw);
                    end
                end

                % Elements 2+ (layers 3&4, 5&6, ...): pointy tip via tangent intersection
                for e = 2:n_elements
                    li  = 2*e - 1;
                    loe = 2*e;
                    if isempty(tip_pts{e}), continue; end
                    tpt = tip_pts{e};
                    p1  = final_lens{li}(end,:);
                    p2  = final_lens{loe}(end,:);
                    if ~any(isnan(tpt))
                        % White gap
                        plot([p1(1),tpt(1)],[p1(2),tpt(2)],'Color','w','LineWidth',gapLW);
                        plot([p2(1),tpt(1)],[p2(2),tpt(2)],'Color','w','LineWidth',gapLW);
                        if is_sci
                            plot([p1(1),tpt(1)],[p1(2),tpt(2)],'Color',cur_glowcol,'LineWidth',cur_edgelw+1.8);
                            plot([p2(1),tpt(1)],[p2(2),tpt(2)],'Color',cur_glowcol,'LineWidth',cur_edgelw+1.8);
                        end
                        plot([p1(1),tpt(1)],[p1(2),tpt(2)],'Color',ec,'LineWidth',cur_edgelw);
                        plot([p2(1),tpt(1)],[p2(2),tpt(2)],'Color',ec,'LineWidth',cur_edgelw);
                    end

                    if draw_mirror
                        p1m = [p1(1), -p1(2)];
                        p2m = [p2(1), -p2(2)];
                        tpb = [tpt(1), -tpt(2)];
                        if ~any(isnan(tpb))
                            plot([p1m(1),tpb(1)],[p1m(2),tpb(2)],'Color','w','LineWidth',gapLW);
                            plot([p2m(1),tpb(1)],[p2m(2),tpb(2)],'Color','w','LineWidth',gapLW);
                            if is_sci
                                plot([p1m(1),tpb(1)],[p1m(2),tpb(2)],'Color',cur_glowcol,'LineWidth',cur_edgelw+1.8);
                                plot([p2m(1),tpb(1)],[p2m(2),tpb(2)],'Color',cur_glowcol,'LineWidth',cur_edgelw+1.8);
                            end
                            plot([p1m(1),tpb(1)],[p1m(2),tpb(2)],'Color',ec,'LineWidth',cur_edgelw);
                            plot([p2m(1),tpb(1)],[p2m(2),tpb(2)],'Color',ec,'LineWidth',cur_edgelw);
                        end
                    end
                end

                % Draw Source
                plot(origin_global(1), origin_global(2), 'ro', 'MarkerFaceColor','r');

                % Arrowhead size (in data units = mm), scaled per view
                if strcmp(vName, 'zoom')
                    ah_len = 0.5;   ah_wid = 0.07;
                else
                    ah_len = 4.0;   ah_wid = 0.55;
                end

                % ---- Draw Stored Rays (per-segment; arrows on air segments) ----
                for r = 1:N_rays
                    path = all_paths{r};
                    if size(path,1) < 2, continue; end
                    do_mirror = draw_mirror && max(abs(path(:,2))) > 1e-3;
                    n_segs = size(path,1) - 1;
                    for k = 1:n_segs
                        P1 = path(k,:); P2 = path(k+1,:);
                        is_air_seg = (mod(k-1,2) == 0); % seg 1,3,5,... = air (0-indexed even)
                        % --- upper half ---
                        if is_sci
                            plot([P1(1),P2(1)],[P1(2),P2(2)],'-', ...
                                'Color',[current_ray_color,0.15],'LineWidth',cur_raylw+1.2);
                        end
                        plot([P1(1),P2(1)],[P1(2),P2(2)],'-', ...
                            'Color',current_ray_color,'LineWidth',cur_raylw);
                        if is_air_seg
                            ray_arrowhead(P1,P2,current_ray_color,ah_len,ah_wid);
                        end
                        % --- mirror half ---
                        if do_mirror
                            P1m=[P1(1),-P1(2)]; P2m=[P2(1),-P2(2)];
                            if is_sci
                                plot([P1m(1),P2m(1)],[P1m(2),P2m(2)],'-', ...
                                    'Color',[current_ray_color,0.15],'LineWidth',cur_raylw+1.2);
                            end
                            plot([P1m(1),P2m(1)],[P1m(2),P2m(2)],'-', ...
                                'Color',current_ray_color,'LineWidth',cur_raylw);
                            if is_air_seg
                                ray_arrowhead(P1m,P2m,current_ray_color,ah_len,ah_wid);
                            end
                        end
                    end
                end

                % Apply View Limits
                if strcmp(vName, 'zoom')
                    xlim([-10, 10]);
                    ylim([-10, 10]);
                else
                    xlim([xmin - 0.1*x_span, x_end + 0.15*x_span]);
                    ylim([-y_lim_base, y_lim_base]);
                end

                % Apply Grid
                if strcmp(gName, 'withgrid')
                    grid on;
                else
                    grid off;
                end

                % No axis decorations — clean scientific style
                set(gca,'DataAspectRatio',[1 1 1]);
                axis off;

                % Scale bar (bottom-left corner, in data coordinates)
                % Regular view: 5 mm;  Zoom view: 1 mm
                if strcmp(vName, 'zoom')
                    sb_len = 1;  sb_label = '1 mm';
                else
                    sb_len = 10;  sb_label = '10 mm';
                end
                xl = xlim; yl = ylim;
                sb_x = xl(1) + 0.04*(xl(2)-xl(1));
                sb_y = yl(1) + 0.04*(yl(2)-yl(1));
                sb_h = 0.010*(yl(2)-yl(1));   % end-cap height
                line([sb_x, sb_x+sb_len],         [sb_y, sb_y],           'Color','k','LineWidth',1.5);
                line([sb_x,       sb_x],           [sb_y-sb_h, sb_y+sb_h],'Color','k','LineWidth',1.5);
                line([sb_x+sb_len, sb_x+sb_len],   [sb_y-sb_h, sb_y+sb_h],'Color','k','LineWidth',1.5);
                text(sb_x+sb_len/2, sb_y+2.5*sb_h, sb_label, ...
                    'HorizontalAlignment','center','VerticalAlignment','bottom', ...
                    'FontSize',8,'FontName','Helvetica','Color','k');

                set(fig, 'Renderer', 'painters');
                
                % Force MATLAB to draw the figure immediately
                drawnow;
                
                fprintf('Plotted %d/%d: %s\n', current_fig, total_figs, fig_name);
                current_fig = current_fig + 1;
            end
        end
    end
end

fprintf('All %d figures have been generated!\n', total_figs);


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
    alpha_new = NaN; theta_out = NaN; return; 
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

function ray_arrowhead(P1, P2, color, ah_len, ah_wid)
% Draw a thin filled triangle at the midpoint of segment P1→P2.
% ah_len = arrowhead length (mm), ah_wid = half-width at base (mm).
dx = P2(1)-P1(1); dy = P2(2)-P1(2);
seg_len = hypot(dx,dy);
if seg_len < 1.5*ah_len, return; end
ux = dx/seg_len; uy = dy/seg_len;   % unit direction
px = -uy;        py =  ux;          % unit perpendicular
tx = P2(1);            ty = P2(2);             % tip at ray endpoint
bx = P2(1) - ah_len*ux; by = P2(2) - ah_len*uy; % base behind tip
fill([tx,  bx+ah_wid*px,  bx-ah_wid*px], ...
     [ty,  by+ah_wid*py,  by-ah_wid*py], ...
     color, 'EdgeColor', color, 'LineWidth', 0.4);
end

function pt = tip_intersection(p1, a1, p2, a2)
% Find intersection of two lines: each defined by a point and an angle (degrees).
if mod(a1,180) == 90
    x = p1(1); m2 = tand(a2); y = m2*(x-p2(1)) + p2(2);
elseif mod(a2,180) == 90
    x = p2(1); m1 = tand(a1); y = m1*(x-p1(1)) + p1(2);
else
    m1 = tand(a1); m2 = tand(a2);
    if abs(m1-m2) < 1e-6, pt = [NaN,NaN]; return; end
    x = (p2(2)-p1(2) + m1*p1(1) - m2*p2(1)) / (m1-m2);
    y = m1*(x-p1(1)) + p1(2);
end
if abs(x) > 1e6 || abs(y) > 1e6, pt = [NaN,NaN]; return; end
pt = [x, y];
end

function alphas = sample_ray_angles(N, profile, seed, a_min, a_max, eps_ang)
if nargin >= 3 && ~isempty(seed), rng(seed); end
if nargin < 4 || isempty(a_min),   a_min   = 0;   end
if nargin < 5 || isempty(a_max),   a_max   = 180; end
if nargin < 6 || isempty(eps_ang), eps_ang = 0;   end

profile = lower(strrep(profile,'^',''));
switch profile
    case 'uniform'
        if a_max <= a_min, error('uniform_max_angle must be > uniform_min_angle.'); end
        lo = a_min + eps_ang; hi = a_max - eps_ang;
        if hi <= lo, error('uniform_eps is too large.'); end
        alphas = linspace(lo, hi, N).';
    case 'sin2'
        alphas=zeros(N,1); k=0;
        while k<N
            M=max(1000,2*(N-k));
            c=180*rand(M,1); u=rand(M,1);
            acc=u<=sind(c).^2; nacc=sum(acc);
            if nacc>0
                take=min(nacc,N-k); idx=find(acc,take,'first');
                alphas(k+1:k+take)=c(idx); k=k+take;
            end
        end
    otherwise
        error('Unknown profile "%s".',profile);
end
end