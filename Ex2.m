clc; clear; close all;

% Parameters
epsilon = 0.1;
alpha = -1; 
beta = 1.5;
N = 100; % Number of grid points
h = 1/N;
tol = 1e-6;
max_iter = 100;

% Compute w_0 and x_bar from equations (2.103) and (2.104)
a = 0; b = 1;
w0 = 0.5 * (a - b + beta - alpha);
x_bar = 0.5 * (a + b - alpha - beta);

% Define grid points
x = linspace(a, b, N+1)';

% Use the corrected initial guess
%u = alpha + (x - a) * (beta - alpha) / (b - a) + w0 * tanh(w0 * (x - x_bar) / (2 * epsilon));
u = x-x_bar + w0 * tanh(w0 * (x - x_bar) / (2 * epsilon));

% Apply Dirichlet boundary conditions
u(1) = alpha;  % Enforce u(0) = alpha
u(end) = beta; % Enforce u(1) = beta

for iter = 1:max_iter
    % Compute Residual F
    F = zeros(N-1,1);
    for j = 2:N
        F(j-1) = epsilon * (u(j+1) - 2*u(j) + u(j-1)) / h^2 ...
                 + (u(j) * (u(j+1) - u(j-1)) / (2*h)) - u(j);
    end

    % Compute Jacobian J
    J = zeros(N-1, N-1);
    for j = 2:N
        row_idx = j-1; % Since MATLAB indices start from 1

        % Handle left diagonal (J(row_idx, row_idx-1)) safely
        if row_idx > 1
            J(row_idx, row_idx-1) = epsilon / h^2 - u(j) / (2*h);
        end

        % Main diagonal
        J(row_idx, row_idx) = -2*epsilon / h^2 - (u(j+1) - u(j-1)) / (2*h) - 1;

        % Handle right diagonal (J(row_idx, row_idx+1)) safely
        if row_idx < N-1
            J(row_idx, row_idx+1) = epsilon / h^2 + u(j) / (2*h);
        end
    end

    % Solve linear system J * delta_u = -F
    delta_u = J \ (-F);
    
    % Update solution (keeping boundary conditions fixed)
    u(2:N) = u(2:N) + delta_u;

    % Convergence check
    if norm(delta_u, inf) < tol
        fprintf('Converged in %d iterations\n', iter);
        break;
    end
end

% Plot the solution
figure;
plot(x, u, 'b', 'LineWidth', 2);
xlabel('x'); ylabel('u(x)');
title('Solution of the Nonlinear BVP using Newton''s Method');
grid on;
