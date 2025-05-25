clc; clear; close all

% exact solution and parameters
e = 0.1;
u = @(x,t) -tanh((x + 0.5 - t)/(2*e)) + 1;

% simulation parameters
N = 2^8;
h = 2/N;
k_scale = 1/10;
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



%% Spatial Convergence

clc; clear; close all;

% Parameters
e = 0.1;
T = 2;
u_exact = @(x,t) -tanh((x + 0.5 - t) / (2*e)) + 1;

% Space resolutions
N_values = [25, 50, 100, 150, 200, 400];
L2_errors = zeros(size(N_values));
Linf_errors = zeros(size(N_values));


k_fixed = 0.000005;  

for idx = 1:length(N_values)
    N = N_values(idx);
    h = 2 / N;
    x = -1 + h*(0:N);
    
    k = k_fixed;                     
    max_j = round(T / k);           
    k = T / max_j;                  

    % Initial condition
    U = zeros(2, N+1);
    U(1,:) = u_exact(x, 0);

    % Time stepping
    for j = 1:max_j
        t = j * k;
        U(2,:) = Burger_iteration(U(1,:), u_exact(-1, t), u_exact(1, t), k, e);
        U(1,:) = U(2,:);
    end

    % Error 
    u_ex = u_exact(x, T);
    res = U(1,:) - u_ex;
    L2_errors(idx) = sqrt(h) * norm(res, 2);
    Linf_errors(idx) = norm(res, Inf);
end

% Plot
figure;
loglog(N_values, L2_errors, 'o-', ...
       N_values, 1./N_values, '--', ...
       N_values, 1./(N_values.^2), ':', 'LineWidth', 1.5);
legend('L^2 error','O(h)','O(h^2)', 'Location','northeast');
xlabel('N (grid points)'); ylabel('Error');
title('Burgers - L^2 Spatial Convergence'); grid on;

figure;
loglog(N_values, Linf_errors, 'o-', ...
       N_values, 1./N_values, '--', ...
       N_values, 1./(N_values.^2), ':', 'LineWidth', 1.5);
legend('L^\infty error','O(h)','O(h^2)', 'Location','northeast');
xlabel('N (grid points)'); ylabel('Error');
title('Burgers - L^\infty Spatial Convergence'); grid on;

% Convergence rates
rates_L2 = log2(L2_errors(1:end-1) ./ L2_errors(2:end));
rates_Inf = log2(Linf_errors(1:end-1) ./ Linf_errors(2:end));
disp('L2 convergence rates:'); disp(rates_L2);
disp('Linf convergence rates:'); disp(rates_Inf);



%% convergence test in time
clc; clear; close all

% Parameters
e = 0.1;
u = @(x,t) -tanh((x + 0.5 - t)/(2*e)) + 1;
max_time = 2;
k_scale = [0.1 0.01 0.001 0.0001];


N = 128;
h = 2 / N;
x = -1 + h*(0:N);

L2_errors = zeros(size(k_scale));
Linf_errors = zeros(size(k_scale));

% Reference solution
k_ref = h^2 / (2*e) * k_scale(end)/2; 
max_j_ref = round(max_time / k_ref);
U_ref = u(x, 0);  

for j = 1:max_j_ref
    t = j * k_ref;
    U_ref = Burger_iteration(U_ref, u(-1,t), u(1,t), k_ref, e);
end

% Coarser time grids
for idx = 1:length(k_scale)
    k = h^2 / (2*e) * k_scale(idx);
    max_j = round(max_time / k);
    
    U_num = u(x, 0);  
    for j = 1:max_j
        t = j * k;
        U_num = Burger_iteration(U_num, u(-1,t), u(1,t), k, e);
    end
    
    % Compute error 
    residual = U_num - U_ref;
    L2_errors(idx) = sqrt(h) * norm(residual, 2);
    Linf_errors(idx) = norm(residual, Inf);
end


%%
% Plot 
figure;
loglog(k_scale, L2_errors, 'o-', k_scale, k_scale, '--', k_scale,k_scale.^2, ':', 'LineWidth', 1.5);
legend('L^2 error','O(h)','O(h^2)', 'Location','northwest');
xlabel('h'); ylabel('Error'); title('Convergence in L^2 norm');

figure;
loglog(k_scale, Linf_errors, 'o-', k_scale, k_scale, '--', k_scale,k_scale.^2, ':', 'LineWidth', 1.5);
legend('L^\infty error','O(h)','O(h^2)', 'Location','northwest');
xlabel('h'); ylabel('Error'); title('Convergence in L^\infty norm');

% Convergence rates
rates_L2 = log2(L2_errors(1:end-1) ./ L2_errors(2:end));
rates_Inf = log2(Linf_errors(1:end-1) ./ Linf_errors(2:end));
disp('L2 convergence rates:'); disp(rates_L2);
disp('Linf convergence rates:'); disp(rates_Inf);