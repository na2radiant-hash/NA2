%% ============================================================
%  STEP 1: LENS CONSTRUCTION
%  180° fan + weighted bending, but STRICT collimation
%  ============================================================
clear; clc; close all;

%% -------------------------
%  Construction Parameters
%  -------------------------
n1            = 1;        % Outside medium refractive index
n2            = 2.003;     % Lens medium refractive index
layers_count  = 8;       % Number of layers/surfaces to construct
rays_count    = 10;       % Number of rays to *design* the lens with
source_pos    = [0, 0];   % Source position
output_filename = 'lens_data_gradual_N2003_test.mat';
enable_filled_overlay = false;   % toggle: true/false


fprintf('--- Lens Construction Phase ---\n');

try
    fprintf('Constructing %d-layer lens for %d rays (full 0–180° fan, collimated)...\n', ...
        layers_count, rays_count);

    % ---- MOD: get theta_in/out matrices too (like gamma_matrix) ----
    [final_lens, gamma_matrix, theta_in_matrix, theta_out_matrix] = construct_lens( ...
        layers_count, rays_count, n1, n2, source_pos);

    fprintf('Lens construction complete.\n');

    % --- Save the lens data ---
    % (Optional) Save theta matrices too, since they’re useful for debugging/tracing
    save(output_filename, 'final_lens', 'n1', 'n2', ...
         'layers_count', 'source_pos', 'gamma_matrix', ...
         'theta_in_matrix', 'theta_out_matrix');
    fprintf('Lens data successfully saved to %s\n', output_filename);

    % --- Plot the constructed lens ---
    figure('Name', 'Lens Construction Output');
    clf; hold on;
    set(gcf, 'Color', 'w');
    xlabel('X-axis'); ylabel('Y-axis');
    axis equal; grid on;
    title('Constructed Lens Geometry (weighted bending, collimated output)');

    mirrored_lens = mirror_lens(final_lens);
    for j = 1:layers_count
        if isempty(final_lens{j}), continue; end
        plot(final_lens{j}(:,1),    final_lens{j}(:,2),    'k', 'LineWidth', 1.5);
        plot(mirrored_lens{j}(:,1), mirrored_lens{j}(:,2), 'k', 'LineWidth', 1.5);
    end
    plot(source_pos(1), source_pos(2), 'ro', 'MarkerSize', 8, 'MarkerFaceColor', 'r');
    fprintf('Plotted constructed lens.\n');
        % --------------------------------------------------------
    % OPTIONAL FILLED OVERLAY (easy to disable)
    % Creates closed polygons for each layer (top + mirrored bottom)
    % --------------------------------------------------------
    if enable_filled_overlay
        overlay_handles = gobjects(0);

        for j = 1:layers_count
            if isempty(final_lens{j}), continue; end

            top = final_lens{j};

            % Make a closed polygon by adding mirrored bottom and closing
            poly_xy = [top; flipud([top(:,1), -top(:,2)])];
            if any(isnan(poly_xy(:))) || any(isinf(poly_xy(:)))
                continue;
            end
            if ~isequal(poly_xy(1,:), poly_xy(end,:))
                poly_xy = [poly_xy; poly_xy(1,:)];
            end

            h = patch(poly_xy(:,1), poly_xy(:,2), [0.2 0.7 0.35], ...
                      'FaceAlpha', 0.25, 'EdgeColor', 'none');
            overlay_handles(end+1) = h; %#ok<AGROW>
        end

        % Put overlay behind black outlines (so outlines stay crisp)
        uistack(overlay_handles, 'bottom');

        % Store handles so you can toggle later via command line if you want
        setappdata(gcf, 'lens_overlay_handles', overlay_handles);

        fprintf('Filled overlay enabled. (You can disable by setting enable_filled_overlay=false)\n');
    else
        fprintf('Filled overlay disabled.\n');
    end

        % --------------------------------------------------------
    % EXPORT FIGURE TO SVG (vector) for Inkscape
    % --------------------------------------------------------
    set(gcf, 'Renderer', 'painters');              % ensure vector output (no rasterization)
    svg_filename = 'lens_geometry.svg';

    % Optional: make it cleaner for publication
    grid off; box off;
    set(gca, 'TickDir','out');

    set(gcf,'Renderer','painters');   % חשוב לוקטור אמיתי

    svg_filename = 'lens_geometry.svg';
    print(gcf, svg_filename, '-dsvg');


    fprintf('Figure exported to SVG: %s\n', svg_filename);

    % --------------------------------------------------------
    % EXPORT OPTIONS
    % --------------------------------------------------------
    thickness_mm  = 10;            % info only for 2D exports
    export_format = 'step';        % 'xt' or 'step'

    switch lower(export_format)
        case 'xt'
            % 2D top-half closed polygon from layers 1&2
            export_lens_to_xt(final_lens, thickness_mm, 'lens_export.xt');
        case 'step'
            % 2D top-half closed polygon from layers 1&2
            export_lens_to_step(final_lens, thickness_mm, 'lens_export.step');
        otherwise
            warning('Unknown export format "%s". No export done.', export_format);
    end

catch ME
    fprintf('ERROR during lens construction: %s\n', ME.message);
    fprintf('Please check the construction functions.\n');
end


%% ============================================================
%  ===============  LENS CONSTRUCTION FUNCTIONS  ==============
%% ============================================================

function [alpha, theta_in, theta_out, gamma] = calculate_angles(L, R, n1, n2)
    % Angles in degrees. alpha(:,1) is the initial global angle of each ray.
    epsilon   = 0.1;
    alpha     = zeros(R, L+1);
    theta_in  = zeros(R, L);
    theta_out = zeros(R, L);
    gamma     = zeros(R, L);

    % ---- Initial fan of design rays in [0, 180 - eps] (full span) ----
    for i = 1:R
        alpha(i,1) = (i - 1) * (180 - epsilon) / (R - 1);
    end

    % ------------ enforce common final angle (COLLIMATION) -------------
    alpha0        = alpha(:,1);               % initial angle per ray
    theta_target  = 0;                        % all rays must end at this angle
    alpha_final   = theta_target * ones(R,1);
    total_bend    = alpha0 - alpha_final;     % total bend per ray

    % ------------- per-layer/per-ray weights (distribution) ------------
    J_front      = min(10, L);
    k_shape      = 4.0;
    base_w       = 0.8;
    hardness_exp = 2.0;
    layer_exp    = 1.2;

    alpha_max = max(alpha0);
    w    = zeros(R, L);
    Wsum = zeros(R, 1);

    for i = 1:R
        hardness = (alpha0(i) / alpha_max)^hardness_exp;

        for j = 1:L
            if j <= J_front
                frac_layer = (J_front - (j-1)) / J_front;
                frac_layer = frac_layer^layer_exp;

                extra   = k_shape * hardness * frac_layer;
                w(i,j) = base_w + extra;
            else
                w(i,j) = base_w;
            end
        end
        Wsum(i) = sum(w(i,:));
    end

    % ---------------- compute angles layer by layer ----------------
    for j = 1:L
        for i = 1:R
            % Alternating media across successive surfaces
            if mod(j,2) == 1
                n_in  = n1;
                n_out = n2;
            else
                n_in  = n2;
                n_out = n1;
            end

            % distribute the per-ray total bend with the weights
            delta_alpha = total_bend(i) * ( w(i,j) / Wsum(i) );

            alpha_j    = alpha(i,j);
            alpha_next = alpha_j - delta_alpha;

            % Snell in the correct frame:
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

            % Local incidence & refraction relative to the chosen normal
            th_in  = wrapTo180(alpha_j    - phi_val);
            th_out = wrapTo180(alpha_next - phi_val);

            % Store results
            theta_in(i,j)  = th_in;
            theta_out(i,j) = th_out;
            gamma(i,j)     = phi_val + 90;      % surface tangent angle

            % Propagate global angle
            alpha(i,j+1) = alpha_next;
        end
    end
end

% ---- helpers ----
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

function [final_lens, gamma_matrix, theta_in_matrix, theta_out_matrix] = construct_lens(L, R, n1, n2, source_pos)
    % ---- MOD: keep theta_in/theta_out as matrices (like gamma_matrix) ----
    [alpha, theta_in, theta_out, gamma] = calculate_angles(L, R, n1, n2);

    gamma_matrix     = gamma;
    theta_in_matrix  = theta_in;
    theta_out_matrix = theta_out;

    d          = 1;
    jump_size  = 0.1;
    final_lens = cell(L+1, 1);

    % ---- Construct the first layer (j=1) ----
    lens_temp = [d, 0]; % Start at (d, 0)
    for i = 1:R-1
        point1 = lens_temp(i, :);
        angle1 = gamma(i, 1);
        point2 = source_pos;
        angle2 = alpha(i+1, 1);
        intersection_point = find_intersection(point1, angle1, point2, angle2);
        if any(isnan(intersection_point))
            intersection_point = [d, tand(angle2)*d]; % Failsafe: simple projection
        end
        lens_temp = [lens_temp; intersection_point];
    end
    final_lens{1} = lens_temp;

    % ---- Construct subsequent layers (j=2..L) ----
    for j = 2:L
        layer_valid = false;
        attempts    = 0;
        max_attempts= 100000;
        current_layer_x_start_guess = final_lens{j-1}(1,1) + jump_size;

        while ~layer_valid && attempts < max_attempts
            attempts = attempts + 1;
            lens_temp = [current_layer_x_start_guess, 0];

            for i = 1:R-1
                point1 = lens_temp(i, :);
                angle1 = gamma(i, j);

                % Intersect with ray from previous layer point
                point2 = final_lens{j-1}(i+1, :);
                angle2 = alpha(i+1, j);

                intersection_point = find_intersection(point1, angle1, point2, angle2);
                if any(isnan(intersection_point))
                    break;
                end
                lens_temp = [lens_temp; intersection_point];
            end

            if size(lens_temp, 1) == R
                if is_valid_layer(lens_temp, final_lens(1:j-1), j)
                    final_lens{j} = lens_temp;
                    layer_valid   = true;
                else
                    current_layer_x_start_guess = current_layer_x_start_guess + jump_size;
                end
            else
                current_layer_x_start_guess = current_layer_x_start_guess + jump_size;
            end
        end

        if ~layer_valid
            error(['Failed to construct layer ', num2str(j), ...
                   ' after ', num2str(max_attempts), ' attempts.']);
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


%% ============================================================
%  XT EXPORT (2D, single element, closed polygon from layers 1&2)
%% ============================================================
function export_lens_to_xt(final_lens, thickness_mm, filename)
    layers = final_lens(~cellfun(@isempty, final_lens));
    if numel(layers) < 2
        error('export_lens_to_xt:NotEnoughLayers', ...
              'Need at least 2 non-empty layers (front/back) to build a single element.');
    end

    front = layers{1};
    back  = layers{2};

    front_top = front(front(:,2) >= 0, :);
    back_top  = back(back(:,2)  >= 0, :);

    if size(front_top,1) < 2 || size(back_top,1) < 2
        error('export_lens_to_xt:TooFewPoints', ...
              'Not enough points in top half of layers 1 and 2.');
    end

    [~, idxF] = sort(front_top(:,1), 'ascend');
    front_top = front_top(idxF,:);

    [~, idxB] = sort(back_top(:,1), 'ascend');
    back_top  = back_top(idxB,:);

    poly_xy = [front_top; flipud(back_top)];

    if ~isequal(poly_xy(1,:), poly_xy(end,:))
        poly_xy = [poly_xy; poly_xy(1,:)];
    end

    N = size(poly_xy,1);

    fid = fopen(filename, 'w');
    if fid == -1
        error('export_lens_to_xt:FileOpenError', 'Could not open file "%s" for writing.', filename);
    end

    fprintf(fid, '# 2D LENS PROFILE XT EXPORT\n');
    fprintf(fid, '# thickness_mm (for info only): %.6f\n', thickness_mm);
    fprintf(fid, '# One closed polygon built from layers 1 & 2 (top half only)\n\n');

    fprintf(fid, 'PROFILE N %d\n', N);
    for i = 1:N
        x = poly_xy(i,1);
        y = poly_xy(i,2);
        z = 0.0;
        fprintf(fid, '%.9g %.9g %.9g\n', x, y, z);
    end

    fclose(fid);
    fprintf('2D XT export complete: %s (points: %d, thickness info: %.3f mm)\n', ...
            filename, N, thickness_mm);
end

%% ============================================================
%  STEP EXPORT (placeholder wrapper if you already have it)
%% ============================================================
function export_lens_to_step(final_lens, thickness_mm, filename) %#ok<INUSD>
    % If you already have a global STEP exporter, keep using it.
    % This stub avoids breaking your script if export_lens_to_step exists elsewhere.
    warning('export_lens_to_step is not defined here. Using your existing implementation (if on path).');
    if exist('export_lens_to_step', 'file') == 2
        % If you have a separate file with this function, MATLAB will call it.
        feval('export_lens_to_step', final_lens, thickness_mm, filename);
    end
end