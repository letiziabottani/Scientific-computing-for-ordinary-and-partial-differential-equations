clc; clear; close all;

% Parameters
epsilon = 0.1;
alpha = -1; 
beta = 1.5;
N_exact = 2^8; % Number of grid points
h = 1/N_exact;
tol = 1e-6;
max_iter = 100;

% Compute w_0 and x_bar from equations (2.103) and (2.104)
a = 0; b = 1;
w0 = 0.5 * (a - b + beta - alpha);
x_bar = 0.5 * (a + b - alpha - beta);

% Define grid points
x = linspace(a, b, N_exact+1)';

% Use the corrected initial guess
%u = alpha + (x - a) * (beta - alpha) / (b - a) + w0 * tanh(w0 * (x - x_bar) / (2 * epsilon));
u = x-x_bar + w0 * tanh(w0 * (x - x_bar) / (2 * epsilon));

% Apply Dirichlet boundary conditions
u(1) = alpha;  % Enforce u(0) = alpha
u(end) = beta; % Enforce u(1) = beta

for iter = 1:max_iter
    % Compute Residual F
    F = zeros(N_exact-1,1);
    N=N_exact
    for j = 2:N
        F(j-1) = epsilon * (u(j+1) - 2*u(j) + u(j-1)) / h^2 ...
                 + (u(j) * (u(j+1) - u(j-1)) / (2*h)) - u(j);
    end

    % Compute Jacobian J
   
    main_diag  = -2 * epsilon / h^2 - (u(3:N+1) - u(1:N-1)) / (2*h) - 1;
    upper_diag = epsilon / h^2 + u(2:N) / (2*h);
    lower_diag = epsilon / h^2 - u(2:N) / (2*h);
    
    % Construct sparse Jacobian matrix J
    J = spdiags([lower_diag, main_diag, upper_diag], [-1, 0, 1], N-1, N-1);

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
u_exact=u;

%%

N = 2.^[2:6]; % Number of grid points
h = 1./N;
err=[]
% Compute w_0 and x_bar from equations (2.103) and (2.104)

w0 = 0.5 * (a - b + beta - alpha);
x_bar = 0.5 * (a + b - alpha - beta);

for i=1:length(N)
    % Define grid points
    x = linspace(a, b, N(i)+1)';
    
    % Use the corrected initial guess
    %u = alpha + (x - a) * (beta - alpha) / (b - a) + w0 * tanh(w0 * (x - x_bar) / (2 * epsilon));
    u = x-x_bar + w0 * tanh(w0 * (x - x_bar) / (2 * epsilon));
    
    % Apply Dirichlet boundary conditions
    u(1) = alpha;  % Enforce u(0) = alpha
    u(end) = beta; % Enforce u(1) = beta
    
    for iter = 1:max_iter
        % Compute Residual F
        F = zeros(N(i)-1,1);
        for j = 2:N(i)
            F(j-1) = epsilon * (u(j+1) - 2*u(j) + u(j-1)) / h(i)^2 ...
                     + (u(j) * (u(j+1) - u(j-1)) / (2*h(i))) - u(j);
        end
    
        % Compute Jacobian J
       
        main_diag  = -2 * epsilon / h(i)^2 - (u(3:N(i)+1) - u(1:N(i)-1)) / (2*h(i)) - 1;
        upper_diag = epsilon / h(i)^2 + u(2:N(i)) / (2*h(i));
        lower_diag = epsilon / h(i)^2 - u(2:N(i)) / (2*h(i));
        
        % Construct sparse Jacobian matrix J
        J = spdiags([lower_diag, main_diag, upper_diag], [-1, 0, 1], N(i)-1, N(i)-1);
    
        % Solve linear system J * delta_u = -F
        delta_u = J \ (-F);
        
        % Update solution (keeping boundary conditions fixed)
        u(2:N(i)) = u(2:N(i)) + delta_u;
    
        % Convergence check
        if norm(delta_u, inf) < tol
            fprintf('Converged in %d iterations\n', iter);
            break;
        end
    end
    steps=[1:(N_exact/N(i)):N_exact+1]
    err(i)=max(abs(u-u_exact(steps)))

end

figure()
loglog(h,err,h,h.^2)
legend("error","h^2")
xlabel("h")
ylabel("error")