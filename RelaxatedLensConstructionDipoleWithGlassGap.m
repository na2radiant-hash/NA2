%% ============================================================
%  STEP 1: LENS CONSTRUCTION (ANGLE-CUTOFF + MULTI-GLASS GAPS)
%  ------------------------------------------------------------
%  Feature:
%    glass_pos can be a vector of element indices, e.g. [2 3 4]
%    -> inserts a rectangular glass gap BEFORE those elements
%    -> i.e. before layers (2*glass_pos - 1): 3,5,7,...
%
%  Construction rule (forced spacer, same algorithm):
%    - after boundary layer j_insert = 2*(elem-1), we enforce an x-gap
%      of width glass_width BEFORE constructing layer (j_insert+1)
%    - the first layer after each gap must be entirely to the right of
%      x_right (no intrusion). If not valid, we shift start-x right by
%      jump_size until valid (same “move-right” algorithm).
%
%  Glass height rule (FIXED):
%    - Each glass rectangle height = 1.3 * height of the NEXT element
%      (the element immediately after the gap).
%    - This height is computed AFTER that element is constructed
%      (because before that, its geometry doesn't exist yet).
%
%  Output:
%    Saves final_lens + glass_info (with per-gap x_left/x_right/h_glass)
% ============================================================
clear; clc; close all;

%% -------------------------
%  Construction Parameters
%  -------------------------
n1            = 1;        % Outside medium refractive index
n2            = 2.003;    % Lens medium refractive index
layers_count  = 8;        % Number of layers/surfaces (must be even = 2*#elements)
rays_count    = 10;      % Number of rays to *design* the lens with
source_pos    = [0, 0];   % Source position
output_filename = 'lens_data_4_elements_no_glass.mat';

% NEW: cutoff for the design fan (degrees)
MAX_FAN_ANGLE_DEG = 175;

% =========================
% NEW: Multi-glass insertion params
% =========================
ENABLE_GLASS = false;

% Can be scalar or vector (elements are 1..layers_count/2)
glass_pos   = [];    % insert glass BEFORE these elements (e.g. [2 3 4])
n3          = 1.5;    % refractive index in the gap(s)
glass_width = 5;      % [mm] gap width (same for all gaps)

fprintf('--- Lens Construction Phase ---\n');
fprintf('Constructing %d-layer lens for %d uniform rays up to %.1f deg...\n', ...
        layers_count, rays_count, MAX_FAN_ANGLE_DEG);

if ENABLE_GLASS
    fprintf('Glass insertion BEFORE elements: [%s] (width %.3g mm, n3=%.4g)\n', ...
        num2str(glass_pos), glass_width, n3);
end

try
    [final_lens, gamma_matrix, glass_info] = construct_lens_with_multi_glass_forced( ...
        layers_count, rays_count, n1, n2, source_pos, MAX_FAN_ANGLE_DEG, ...
        ENABLE_GLASS, glass_pos, n3, glass_width);

    fprintf('Lens construction complete.\n');

    save(output_filename, 'final_lens', 'n1', 'n2', 'layers_count', 'source_pos', ...
         'MAX_FAN_ANGLE_DEG', 'glass_info');

    fprintf('Lens data successfully saved to %s\n', output_filename);

    figure('Name', 'Lens Construction Output');
    clf; hold on;
    set(gcf, 'Color', 'w');
    xlabel('X-axis'); ylabel('Y-axis');
    axis equal; grid on;
    title(sprintf('Constructed Lens Geometry (fan up to %.1f°)', MAX_FAN_ANGLE_DEG));

    % ---- Plot all glass rectangles first (light blue) ----
    if ENABLE_GLASS && glass_info.enabled && ~isempty(glass_info.gaps)
        glassFace = [0.65 0.85 1.00];
        glassEdge = [0.20 0.45 0.75];

        for k = 1:numel(glass_info.gaps)
            g = glass_info.gaps(k);
            if isnan(g.h_glass)
                warning('Gap before element %d has NaN height (should not happen).', g.element);
                continue;
            end
            xL = g.x_left;  xR = g.x_right;  hG = g.h_glass;

            patch([xL xR xR xL], [ hG  hG -hG -hG], ...
                  glassFace, 'FaceAlpha', 0.25, 'EdgeColor', glassEdge, 'LineWidth', 1.2);
            plot([xL xL xR xR xL], [ hG -hG -hG hG hG], ...
                 'Color', glassEdge, 'LineWidth', 1.2);
        end
    end

    mirrored_lens = mirror_lens(final_lens);
    for j = 1:layers_count
        plot(final_lens{j}(:,1), final_lens{j}(:,2), 'k', 'LineWidth', 1.5);
        plot(mirrored_lens{j}(:,1), mirrored_lens{j}(:,2), 'k', 'LineWidth', 1.5);
    end
    plot(source_pos(1), source_pos(2), 'ro', 'MarkerSize', 8, 'MarkerFaceColor', 'r');

catch ME
    fprintf('ERROR during lens construction: %s\n', ME.message);
end

%% ============================================================
%  ===============  LENS CONSTRUCTION FUNCTIONS  ==============
% ============================================================

function [alpha, theta_in, theta_out, gamma] = calculate_angles(L, R, n1, n2, max_fan_angle_deg)
    epsilon   = 0.1;
    alpha     = zeros(R, L+1);
    theta_in  = zeros(R, L);
    theta_out = zeros(R, L);
    gamma     = zeros(R, L);

    % Uniform fan in [0, max_fan_angle_deg - epsilon]
    for i = 1:R
        alpha(i,1) = (i - 1) * (max_fan_angle_deg - epsilon) / (R - 1);
    end

    for j = 1:L
        for i = 1:R
            if mod(j,2) == 1
                n_in  = n1; n_out = n2;
            else
                n_in  = n2; n_out = n1;
            end

            delta_alpha = alpha(i,1) / L;
            alpha_j     = alpha(i,j);
            alpha_next  = alpha_j - delta_alpha;

            % n_in*sin(alpha_j - phi) == n_out*sin(alpha_next - phi)
            syms phi_sym
            f = n_in * sind(alpha_j - phi_sym) - n_out * sind(alpha_next - phi_sym);

            phi_guess = 0;
            try
                phi_val = double(vpasolve(f == 0, phi_sym, [-180, 180]));
                if isempty(phi_val); phi_val = phi_guess; end
            catch
                phi_val = phi_guess;
            end

            th_in  = wrapTo180(alpha_j    - phi_val);
            th_out = wrapTo180(alpha_next - phi_val);

            theta_in(i,j)  = th_in;
            theta_out(i,j) = th_out;
            gamma(i,j)     = phi_val + 90;   % surface tangent angle

            alpha(i,j+1) = alpha_next;
        end
    end
end

function a = wrapTo180(a), a = mod(a + 180, 360) - 180; end

function intersection = find_intersection(p1, a1, p2, a2)
    if mod(a1, 180) == 90
        x  = p1(1); m2 = tand(a2); y  = m2 * (x - p2(1)) + p2(2);
    elseif mod(a2, 180) == 90
        x  = p2(1); m1 = tand(a1); y  = m1 * (x - p1(1)) + p1(2);
    else
        m1 = tand(a1); m2 = tand(a2);
        if abs(m1 - m2) < 1e-6, intersection = [NaN, NaN]; return; end
        x = (p2(2) - p1(2) + m1 * p1(1) - m2 * p2(1)) / (m1 - m2);
        y = m1 * (x - p1(1)) + p1(2);
    end
    if abs(x) > 1e6 || abs(y) > 1e6, intersection = [NaN, NaN]; return; end
    intersection = [x, y];
end

function valid = is_valid_layer(new_layer, final_lens, j)
    valid = true;
    if j == 1; return; end

    for i = 1:size(new_layer, 1)-1
        for k = 1:j
            if k > length(final_lens) || isempty(final_lens{k}), continue; end
            prev_layer = final_lens{k};
            for m = 1:size(prev_layer, 1)-1
                if check_intersection(new_layer(i, :), new_layer(i+1, :), ...
                                      prev_layer(m, :), prev_layer(m+1, :))
                    valid = false; return;
                end
            end
        end
    end
end

function intersects = check_intersection(p1, p2, q1, q2)
    intersects = false;
    if any(isinf(p1)) || any(isnan(p1)) || any(isinf(p2)) || any(isnan(p2)) || ...
       any(isinf(q1)) || any(isnan(q1)) || any(isinf(q2)) || any(isnan(q2))
        return;
    end
    denom = (p1(1) - p2(1)) * (q1(2) - q2(2)) - (p1(2) - p2(2)) * (q1(1) - q2(1));
    if abs(denom) < 1e-9, return; end
    t = ((p1(1) - q1(1)) * (q1(2) - q2(2)) - (p1(2) - q1(2)) * (q1(1) - q2(1))) / denom;
    u = -((p1(1) - p2(1)) * (p1(2) - q1(2)) - (p1(2) - p2(2)) * (p1(1) - q1(1))) / denom;
    intersects = (t >= 0 && t <= 1 && u >= 0 && u <= 1);
end

% ============================================================
%  MULTI-GLASS: forced gaps + same move-right algorithm
%  + FIXED GLASS HEIGHT: computed AFTER next element exists
% ============================================================
function [final_lens, gamma_matrix, glass_info] = construct_lens_with_multi_glass_forced( ...
    L, R, n1, n2, source_pos, max_fan_angle_deg, ...
    enable_glass, glass_pos, n3, glass_width)

    if mod(L,2) ~= 0
        error('layers_count must be even (2 layers per element). Got L=%d.', L);
    end

    [alpha, ~, ~, gamma] = calculate_angles(L, R, n1, n2, max_fan_angle_deg);
    gamma_matrix = gamma;

    d          = 1;
    jump_size  = 0.1;
    final_lens = cell(L+1, 1);

    % Normalize glass_pos to sorted unique row vector
    if enable_glass
        if isempty(glass_pos)
            enable_glass = false;
            glass_pos = [];
        else
            glass_pos = unique(glass_pos(:).','sorted');
            max_elem = L/2;
            if any(glass_pos < 2) || any(glass_pos > max_elem)
                error('glass_pos values must be in [2 .. %d] for layers_count=%d.', max_elem, L);
            end
        end
    else
        glass_pos = [];
    end

    % For each element e: boundary after element (e-1) is j_insert=2*(e-1),
    % and the first layer after the gap is j_enforce=j_insert+1
    if enable_glass
        j_inserts = 2*(glass_pos - 1);
        j_enforce = j_inserts + 1;
    else
        j_inserts = [];
        j_enforce = [];
    end

    % glass_info struct
    glass_info = struct();
    glass_info.enabled     = enable_glass;
    glass_info.glass_pos   = glass_pos;
    glass_info.n3          = n3;
    glass_info.glass_width = glass_width;
    glass_info.gaps        = struct('element', {}, 'j_insert', {}, 'j_enforce', {}, ...
                                    'x_left', {}, 'x_right', {}, 'h_glass', {});

    % ---- Layer 1 ----
    lens_temp = [d, 0];
    for i = 1:R-1
        point1 = lens_temp(i, :);
        angle1 = gamma(i, 1);
        point2 = source_pos;
        angle2 = alpha(i+1, 1);

        intersection_point = find_intersection(point1, angle1, point2, angle2);
        if any(isnan(intersection_point))
            intersection_point = [d, tand(angle2)*d];
        end
        lens_temp = [lens_temp; intersection_point]; %#ok<AGROW>
    end
    final_lens{1} = lens_temp;

    tol = 1e-6;

    % ---- Layers 2..L ----
    for j = 2:L
        layer_valid = false;
        attempts    = 0;
        max_attempts= 100000;

        base_start = final_lens{j-1}(1,1) + jump_size;

        % Is this layer right after a glass gap?
        enforce_idx = [];
        if enable_glass
            enforce_idx = find(j_enforce == j, 1, 'first');
        end

        enforce_no_intrusion = false;
        x_right = NaN;

        if ~isempty(enforce_idx)
            % x_left from last layer before gap (j-1 == j_insert)
            x_left  = max(final_lens{j-1}(:,1));
            x_right = x_left + glass_width;

            % initial forced jump over gap
            base_start = base_start + glass_width;

            enforce_no_intrusion = true;

            % Record this gap now; height will be filled later (after element exists)
            elem = glass_pos(enforce_idx); % "next element" after the gap

            g.element   = elem;
            g.j_insert  = j_inserts(enforce_idx);
            g.j_enforce = j_enforce(enforce_idx);
            g.x_left    = x_left;
            g.x_right   = x_right;
            g.h_glass   = NaN;            % <-- fill after element elem is built
            glass_info.gaps(end+1) = g; %#ok<AGROW>
        end

        current_layer_x_start_guess = base_start;

        while ~layer_valid && attempts < max_attempts
            attempts = attempts + 1;
            lens_temp = [current_layer_x_start_guess, 0];

            for i = 1:R-1
                point1 = lens_temp(i, :);
                angle1 = gamma(i, j);

                point2 = final_lens{j-1}(i+1, :);
                angle2 = alpha(i+1, j);

                intersection_point = find_intersection(point1, angle1, point2, angle2);
                if any(isnan(intersection_point))
                    break;
                end
                lens_temp = [lens_temp; intersection_point]; %#ok<AGROW>
            end

            if size(lens_temp, 1) == R
                % Hard constraint: do not intrude into the gap region
                if enforce_no_intrusion
                    if min(lens_temp(:,1)) < (x_right + tol)
                        current_layer_x_start_guess = current_layer_x_start_guess + jump_size;
                        continue;
                    end
                end

                % Original: avoid intersections with previous layers
                if is_valid_layer(lens_temp, final_lens(1:j-1), j)
                    final_lens{j} = lens_temp;
                    layer_valid   = true;

                    % ===== FIX: fill glass heights AFTER the next element exists =====
                    % When j is even, it completes element elem_built = j/2
                    if enable_glass && mod(j,2) == 0 && ~isempty(glass_info.gaps)
                        elem_built = j/2;

                        % For any gap whose next element is elem_built and height not filled:
                        for kk = 1:numel(glass_info.gaps)
                            if glass_info.gaps(kk).element == elem_built && isnan(glass_info.gaps(kk).h_glass)
                                h_elem = max(abs(final_lens{j}(:,2)));   % outer surface of that element
                                glass_info.gaps(kk).h_glass = 1.3 * h_elem;
                            end
                        end
                    end
                    % ===============================================================

                else
                    current_layer_x_start_guess = current_layer_x_start_guess + jump_size;
                end
            else
                current_layer_x_start_guess = current_layer_x_start_guess + jump_size;
            end
        end

        if ~layer_valid
            error('Failed to construct layer %d after %d attempts.', j, max_attempts);
        end
    end

    % Safety: ensure all glass heights filled
    if enable_glass && ~isempty(glass_info.gaps)
        for kk = 1:numel(glass_info.gaps)
            if isnan(glass_info.gaps(kk).h_glass)
                warning('Glass height for gap before element %d was not filled.', glass_info.gaps(kk).element);
            end
        end
    end
end

function mirrored_lens = mirror_lens(lens_cell_array)
    mirrored_lens = cell(size(lens_cell_array));
    for j = 1:length(lens_cell_array)
        if isempty(lens_cell_array{j}), continue; end
        mirrored_lens{j}      = lens_cell_array{j};
        mirrored_lens{j}(:,2) = -lens_cell_array{j}(:,2);
    end
end
