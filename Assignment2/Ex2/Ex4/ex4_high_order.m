clc; clear; close all

% Parameters
e = 0.01 / pi;
t_target = 1.6037 / pi;

% Uniform grid
N = 128;
x = linspace(-1, 1, N+1);
h = x(2) - x(1);

% Time step
k_scale = 1/100;
k = h^2 / (2*e) * k_scale;
max_j = ceil(t_target / k);

fprintf("N = %d, k = %.2e, steps = %d\n", N, k, max_j);

% Initial condition and allocation
U = zeros(2, N+1);
U(1,:) = -sin(pi * x);

% Time loop
t = 0; j = 0;
while t < t_target
    j = j + 1;
    t = j * k;
    U(2,:) = Burger_iteration_highorder(U(1,:), 0, 0, k, e, x);
    U(1,:) = U(2,:);
end

% Derivative estimation at x = 0 using 5-point stencil
[~, idx0] = min(abs(x));
if idx0 < 3 || idx0 > N-1
    error("Not enough space for 5-point stencil at x = 0");
end
x_stencil = x(idx0-2:idx0+2);
u_stencil = U(2, idx0-2:idx0+2);
coeffs_d1 = fdcoeffV(1, x(idx0), x_stencil);
ux0_est = coeffs_d1 * u_stencil';

% Display
fprintf("\n--- Estimation Result ---\n");
fprintf("Estimated ∂u/∂x at x = 0, t = %.5f: %.6f\n", t, ux0_est);
fprintf("Reference value: -152.00516\n");
fprintf("Absolute error: %.5e\n", abs(ux0_est + 152.00516));

% Plot
figure;
plot(x, U(2,:), 'b', 'LineWidth', 2);
title('Stable High-Order Solution at t = 1.6037/π');
xlabel('x'); ylabel('u(x,t)'); grid on;
