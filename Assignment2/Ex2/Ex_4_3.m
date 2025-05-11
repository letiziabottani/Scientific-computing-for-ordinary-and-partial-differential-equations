% Parameters
e = 0.01 / pi;                    % Viscosity
N = 2^9;                         % Number of intervals (increase if needed)
h = 2 / N;                        % Spatial step size
x = -1 + h * (0:N);              % Spatial grid

% Time stepping parameters
k_scale = 1 / 100;                % CFL-satisfying constant
k = h^2 / (2 * e) * k_scale;      % Time step
t_target = 1.6037 / pi;          % Target time
max_j = round(t_target / k);     % Number of time steps

% Initial condition: eta(x) = -sin(pi x)
U = zeros(max_j+1, N+1);
U(1, :) = -sin(pi * x);

% Time stepping
for j = 1:max_j
    U(j+1, :) = Burger_iteration(U(j, :), 0, 0, k, e);
end

% Estimate ∂u/∂x at x = 0 using centered difference
[~, ix0] = min(abs(x)); % Index closest to x = 0
ux_approx = (U(end, ix0 + 1) - U(end, ix0 - 1)) / (2 * h);

% Output result
fprintf('Estimated ∂u/∂x at x = 0, t = %.5f: %.5f\n', t_target, ux_approx);
