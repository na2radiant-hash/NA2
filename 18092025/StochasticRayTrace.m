%% ============================================================
%  STEP 2: RAY TRACING THROUGH THE CONSTRUCTED LENS
%  ============================================================
clear; clc; close all;

%% -------------------------
%  User Parameters
%  -------------------------
input_filename = 'lens_data.mat';  % must exist from step 1
N_rays         = 200;              % number of accepted rays to trace
seed           = 42;               % RNG seed for reproducibility
draw_mirror    = true;             % also plot mirrored lens (y -> -y)
max_ray_extent = 2.5;              % how far to extend the final outgoing segment (in x)

rng(seed);


%% -------------------------

function [alpha_new, intersection_point, theta_in, theta_out] = refract_ray_2d(origin_a, alpha, surface_a, surface_b, n1, n2)
%REFRACT_RAY_2D  Refract a 2D ray on a line segment using Snell's law.
%
% Inputs
%   origin_a   : [x,y] start point of ray (in medium n1)
%   alpha      : ray angle (deg) relative to +x axis (0..180 allowed; wraps OK)
%   surface_a  : [x,y] first endpoint of the segment (interface)
%   surface_b  : [x,y] second endpoint of the segment (interface)
%   n1         : refractive index of incident medium
%   n2         : refractive index of transmitted medium
%
% Outputs
%   alpha_new        : refracted ray angle (deg) relative to +x axis in medium n2
%   intersection_point : [x,y] intersection point on the segment
%   theta_in         : signed incidence angle (deg) relative to surface normal
%   theta_out        : signed refraction angle (deg) relative to same normal
%
% Notes
% - Sign convention: positive angles are counter-clockwise. theta_in/theta_out are
%   measured from the chosen normal to the ray; signs are preserved by Snell.
% - If total internal reflection occurs, alpha_new = NaN and theta_out = NaN.

% ---- Ray and segment parametric forms ----
% Ray: R(t) = O + t*D, t >= 0
D = [cosd(alpha), sind(alpha)];
O = origin_a(:).';  % row vec

% Segment: S(u) = A + u*(B-A), u in [0,1]
A = surface_a(:).';
B = surface_b(:).';
S = B - A;

% Solve for intersection: O + t*D = A + u*S
M = [D(:), -S(:)]; % 2x2
rhs = (A - O).';

detM = det(M);
if abs(detM) < 1e-12
    error('Ray and segment are parallel or nearly parallel; no unique intersection.');
end

sol = M \ rhs;  % [t; u]
t = sol(1); u = sol(2);

% Verify intersection lies forward on the ray and within the segment
if t < -1e-12 || u < -1e-12 || u > 1+1e-12
    error('No valid intersection between the ray and the segment.');
end

% Clamp small numerical drift
t = max(t, 0);
u = min(max(u, 0), 1);

intersection_point = O + t*D;

% ---- Surface local frame (tangent & normals) ----
% Unit tangent along the segment
if norm(S) < 1e-12
    error('Surface endpoints are identical.');
end
t_hat = S / norm(S);

% Two unit normals (±90° rotation of tangent)
% n_hat_ccw rotates t_hat by +90°, n_hat_cw by -90°.
n_hat_ccw = [-t_hat(2),  t_hat(1)];
n_hat_cw  = [ t_hat(2), -t_hat(1)];

% Choose the normal that "faces" the incoming ray (dot(D, n_hat) > 0)
% This ensures the incidence angle is in [-90°, +90°].
if dot(D, n_hat_ccw) >= dot(D, n_hat_cw)
    n_hat = n_hat_ccw;
else
    n_hat = n_hat_cw;
end

% ---- Compute signed incidence angle theta_in (normal -> ray) ----
angN   = atan2d(n_hat(2), n_hat(1));      % normal angle
angRay = atan2d(D(2), D(1));              % ray angle
theta_in = wrapTo180(angRay - angN);      % signed, should be within [-90,90]
if abs(theta_in) > 90
    % If numerical drift pushes us outside, flip normal once.
    angN = wrapTo180(angN + 180);
    theta_in = wrapTo180(angRay - angN);
end

% ---- Snell's law ----
% n1 * sin(theta_in) = n2 * sin(theta_out), preserve sign of theta_in.
sin_out = (n1/n2) * sind(theta_in);

% Total internal reflection check
if abs(sin_out) > 1 + 1e-12
    alpha_new = NaN;
    theta_out = NaN;
    warning('Total internal reflection: no refracted ray. Returning NaN for alpha_new and theta_out.');
    return;
end

% Clamp small numerical error
sin_out = max(min(sin_out, 1), -1);

theta_out = asind(sin_out);                 % signed
alpha_new = wrapTo360(angN + theta_out);    % global refracted direction

end

% ---- helpers ----
function a = wrapTo180(a)
% Wrap angle in degrees to (-180, 180]
a = mod(a + 180, 360) - 180;
end

function a = wrapTo360(a)
% Wrap angle in degrees to [0, 360)
a = mod(a, 360);
if a < 0, a = a + 360; end
end



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
if ~any(has_points)
    error('No usable layers found in final_lens.');
end

%% -------------------------
%  Plot setup
%  -------------------------
figure('Name','Ray Tracing Through Constructed Lens'); clf; hold on; grid on; axis equal;
set(gcf, 'Color', 'w');
xlabel('X'), ylabel('Y');
title(sprintf('Ray Tracing with sin^2(\\alpha) acceptance (N = %d)', N_rays));

% Plot constructed lens (upper)
for j = 1:L
    if isempty(final_lens{j}), continue; end
    plot(final_lens{j}(:,1), final_lens{j}(:,2), 'k', 'LineWidth', 1.25);
end
% Optional mirror
if draw_mirror
    for j = 1:L
        if isempty(final_lens{j}), continue; end
        mirror = final_lens{j}; mirror(:,2) = -mirror(:,2);
        plot(mirror(:,1), mirror(:,2), 'k', 'LineWidth', 1.25);
    end
end
plot(origin_global(1), origin_global(2), 'ro', 'MarkerFaceColor','r', 'DisplayName', 'Source');

%% -------------------------
%  Sample rays with sin^2(alpha) acceptance and trace
%  -------------------------
accepted = 0;
attempts = 0;
max_attempts = 1e7;  % safety

while accepted < N_rays && attempts < max_attempts
    attempts = attempts + 1;

    % 3) Draw alpha ~ Uniform(0,180); compute sin^2(alpha)
    alpha = 180 * rand();
    p     = sind(alpha)^2;

    % 4) Accept/reject
    if rand() >= p
        continue; % reject
    end

    % 5) If accepted, trace through the L layers
    P0        = origin_global;
    alpha_cur = alpha;
    n_in      = n1_global;
    n_out     = n2_global;

    hit_all_layers = true;

    for j = 1:L
        % 6) Find first intersection with layer j
        [hit, P_hit, alpha_next] = first_hit_with_layer( ...
            P0, alpha_cur, j, n_in, n_out, final_lens);

        if ~hit
            hit_all_layers = false;
            break;
        end

        % 7–8) Plot segment and advance
        plot([P0(1), P_hit(1)], [P0(2), P_hit(2)], 'b-');
        P0        = P_hit;
        alpha_cur = alpha_next;

        % Alternate indices for next surface
        tmp = n_in; n_in = n_out; n_out = tmp;
    end

    % Draw the final outgoing leg (just for visualization)
    if hit_all_layers
        Df = [cosd(alpha_cur), sind(alpha_cur)];
        if abs(Df(1)) > 1e-6
            t_ext = (max_ray_extent) / abs(Df(1));
        else
            t_ext = max_ray_extent;
        end
        Pend = P0 + t_ext * Df;
        plot([P0(1), Pend(1)], [P0(2), Pend(2)], 'b-');
    end

    accepted = accepted + 1;
end

fprintf('Accepted rays: %d (attempts: %d)\n', accepted, attempts);
legend({'Lens (upper)','Lens (mirror)','Source'}, 'Location','best');

%% =========================
%  Local helper (explicitly passed data)
%  =========================
function [hit, P_hit, alpha_out] = first_hit_with_layer(P0, alpha_in, j, n_in, n_out, final_lens)
%FIRST_HIT_WITH_LAYER Find the nearest forward intersection with layer j and refract.
    hit = false; P_hit = [NaN, NaN]; alpha_out = NaN;

    if j > numel(final_lens) || isempty(final_lens{j}) || size(final_lens{j},1) < 2
        return;
    end

    best_d = inf; best_P = []; best_alpha = NaN;
    pts = final_lens{j};

    for k = 1:size(pts,1)-1
        A = pts(k,   :);
        B = pts(k+1, :);
        try
            [alpha_new, P_inter, ~, ~] = refract_ray_2d(P0, alpha_in, A, B, n_in, n_out);
        catch
            % Parallel / no unique intersection / numerical issue
            continue;
        end
        if any(isnan(P_inter)) || any(isinf(P_inter)) || isnan(alpha_new)
            % Reject invalid intersections or TIR
            continue;
        end
        d = hypot(P_inter(1)-P0(1), P_inter(2)-P0(2));
        if d > 1e-9 && d < best_d
            best_d     = d;
            best_P     = P_inter;
            best_alpha = alpha_new;
        end
    end

    if ~isempty(best_P)
        hit = true; P_hit = best_P; alpha_out = best_alpha;
    end
end
