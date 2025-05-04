clear all; clc; close all;

% Given parameters
a = 0.5;
lambda = 1;                         % wavelength
dx = lambda / 100;                  % 100 points per wavelength
x = 0:dx:1-dx;                      % periodic domain [0,1)
N = length(x);
cr = 0.8;                           % CFL number
dt = cr * dx / a;
T_period = 1 / a;                   % one wave period
t_final = 40 * T_period;           % 40 wave periods
Nt = round(t_final / dt);
dt = t_final / Nt;                 % recompute dt for integer steps

% Initial condition
u = sin(2*pi*x);
u_initial = u;

% Time-stepping using FTBS
for n = 1:Nt
    u_prev = [u(end), u(1:end-1)];
    u = u - cr * (u - u_prev);
end

% Exact solution after 40 periods (should be identical to initial)
u_exact = u_initial;

% Compute errors
amp_error = norm(u) / norm(u_initial);         % total amplitude decay
phase_error = finddelay(u_exact, u) * dx;      % crude phase shift estimate

fprintf('Total amplitude ratio after 40 periods: %.6f\n', amp_error);
fprintf('Approximate phase shift (in space units): %.6f\n', phase_error);

% Plot final vs exact
figure;
plot(x, u, 'b', 'LineWidth', 2); hold on;
plot(x, u_exact, 'r--', 'LineWidth', 2);
xlabel('x'); ylabel('u');
legend('Numerical after 40T', 'Exact');
title('FTBS after 40 wave periods (c_r = 0.8)');
grid on;
