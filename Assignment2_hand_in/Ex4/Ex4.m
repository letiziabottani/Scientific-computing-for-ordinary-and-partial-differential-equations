clc; clear; close all

% exact solution and parameters
e = 0.1;
u = @(x,t) -tanh((x + 0.5 - t)/(2*e)) + 1;

% simulation parameters
N = 2^8;
h = 2/N;
k_scale = 1/100;
k = h^2 / (2*e) * k_scale;

max_time = 2;
max_j = round(max_time / k);

% check for memory issues
if max_j > 1e5
    warning("Time step too small, reducing simulation steps.");
    max_j = 1e5;
end

x = -1 + h*(0:N);
U = zeros(2, N+1); 
U(1,:) = u(x, 0); 

for j = 1:max_j
    t = j * k;
    U(2,:) = Burger_iteration(U(1,:), u(-1,t), u(1,t), k, e);
    U(1,:) = U(2,:); 
end
%%
% Plot final result
figure;
plot(x, U(1,:), 'b', 'DisplayName', 'Numerical'); hold on;
plot(x, u(x, max_j * k), 'r--', 'DisplayName', 'Exact');
xlabel("x")
legend; title('Final solution');

% Full 2D comparison
[X, T] = meshgrid(x, (0:max_j)*k);
U_full = zeros(max_j+1, N+1);
U_full(1,:) = u(x, 0);

for j = 1:max_j
    t = j * k;
    U_full(j+1,:) = Burger_iteration(U_full(j,:), u(-1,t), u(1,t), k, e);
end

figure;
surf(x, (0:max_j)*k, U_full, 'EdgeColor', 'none');
title('Numerical Solution U(x,t)'); xlabel('x'); ylabel('t');

% Compare with exact
u_exact = u(X, T);
err = norm(U_full - u_exact, 2);
fprintf("Global L2 error: %.2e\n", err);

%% convergence test in space
clc; clear; close all
% Parameters
e = 0.1;
u = @(x,t) -tanh((x + 0.5 - t)/(2*e)) + 1;
max_time = 2;
k_scale = [0.1 0.01 0.001 0.0001 ];

% Grids to test
N = 128;               
h = 2 ./ N;
L2_errors = zeros(size(k_scale));
Linf_errors = zeros(size(k_scale));

for idx = 1:length(k_scale)
    
    k = h^2 / (2*e) * k_scale(idx);
    max_j = round(max_time / k);
    
    x = -1 + h*(0:N);
    U_old = u(x, 0);  
    U_new = zeros(1, N+1);
    
    % Time integration (only two time levels)
    for j = 1:max_j
        t = j * k;
        U_new = Burger_iteration(U_old, u(-1,t), u(1,t), k, e);
        U_old = U_new;
    end
    
    % Compute error at final time
    u_exact = u(x, max_time);
    residual = U_new - u_exact;
    L2_errors(idx) = sqrt(h) * norm(residual, 2);         
    Linf_errors(idx) = norm(residual, Inf);               
end



%% convergence test in time
clc; clear; close all

% Parameters
e = 0.1;
u = @(x,t) -tanh((x + 0.5 - t)/(2*e)) + 1;
max_time = 2;
k_scale = [0.1 0.01 0.001 0.0001];

% Fixed spatial grid
N = 128;
h = 2 / N;
x = -1 + h*(0:N);

% Allocate errors
L2_errors = zeros(size(k_scale));
Linf_errors = zeros(size(k_scale));

% Compute Reference Solution 
k_ref = h^2 / (2*e) * k_scale(end)/2;  
max_j_ref = round(max_time / k_ref);
U_ref = u(x, 0);  

for j = 1:max_j_ref
    t = j * k_ref;
    U_ref = Burger_iteration(U_ref, u(-1,t), u(1,t), k_ref, e);
end

% Loop over coarser time steps 
for idx = 1:length(k_scale)
    k = h^2 / (2*e) * k_scale(idx);
    max_j = round(max_time / k);
    
    U_num = u(x, 0);  
    for j = 1:max_j
        t = j * k;
        U_num = Burger_iteration(U_num, u(-1,t), u(1,t), k, e);
    end
    
    % Compute error vs reference solution
    residual = U_num - U_ref;
    L2_errors(idx) = sqrt(h) * norm(residual, 2);
    Linf_errors(idx) = norm(residual, Inf);
end


%%
% Plot convergence results
figure;
loglog(k_scale, L2_errors, 'o-', k_scale, k_scale, '--', k_scale,k_scale.^2, ':', 'LineWidth', 1.5);
legend('L^2 error','O(h)','O(h^2)', 'Location','northwest');
xlabel('h'); ylabel('Error'); title('Convergence in L^2 norm');

figure;
loglog(k_scale, Linf_errors, 'o-', k_scale, k_scale, '--', k_scale,k_scale.^2, ':', 'LineWidth', 1.5);
legend('L^\infty error','O(h)','O(h^2)', 'Location','northwest');
xlabel('h'); ylabel('Error'); title('Convergence in L^\infty norm');

% Optional: Print convergence rates
rates_L2 = log2(L2_errors(1:end-1) ./ L2_errors(2:end));
rates_Inf = log2(Linf_errors(1:end-1) ./ Linf_errors(2:end));
disp('L2 convergence rates:'); disp(rates_L2);
disp('Linf convergence rates:'); disp(rates_Inf);



