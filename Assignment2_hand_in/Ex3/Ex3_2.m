clear all
clc
close all
a = 0.5;                 
L = 2;                   
T = 1;                  

grid_points = [50, 100, 200, 300, 400]; 
errors = zeros(size(grid_points));


for idx = 1:length(grid_points)
    N = grid_points(idx);
    dx = L / N;
    x = linspace(-1, 1, N+1); 
    x(end) = [];              

    
    CFL = 0.75;                  
    dt = CFL * dx / a;
    Nt = round(T/dt);         
    dt = T / Nt;              

    
    u = sin(2*pi*x);


    for n = 1:Nt
        u_previous = [u(end), u(1:end-1)]; 
        u = u - (a*dt/dx) * (u - u_previous);
    end

    u_exact = sin(2*pi*(x - a*T));

    if idx == length(grid_points)
        figure;
        plot(x, u, 'b-', 'LineWidth', 2); hold on;
        plot(x, u_exact, 'r--', 'LineWidth', 2);
        legend('FTBS', 'Exact solution');
        xlabel('x'); ylabel('u(x, T)');
        title(sprintf('FTBS vs exact solution at T = %.2f (N = %d)', T, N));
        grid on;
    end

    % u_exact is continuous function
    errors(idx) = sqrt(sum((u - u_exact).^2) * dx);
end



%%

% Temporal Convergence

clc; clear; close all;


a = 0.5;


L = 2;               
T = 1;               

% Fixed spatial grid
N = 400;             
dx = L / N;
x = linspace(-1, 1, N+1); 
x(end) = [];         

cfl_scales = [0.8, 0.4, 0.2, 0.1];  
L2_errors = zeros(size(cfl_scales));
Linf_errors = zeros(size(cfl_scales));

% Reference solution 
cfl_ref = 0.01;                               
dt_ref = cfl_ref * dx / a;
Nt_ref = round(T / dt_ref);
dt_ref = T / Nt_ref;

u_ref = sin(2*pi*x);                           

for n = 1:Nt_ref
    u_prev = [u_ref(end), u_ref(1:end-1)];     
    u_ref = u_ref - (a*dt_ref/dx) * (u_ref - u_prev);
end

% Loop over time step sizes 
for idx = 1:length(cfl_scales)
    cfl = cfl_scales(idx);
    dt = cfl * dx / a;
    Nt = round(T / dt);
    dt = T / Nt;

    u = sin(2*pi*x);

    for n = 1:Nt
        u_prev = [u(end), u(1:end-1)];
        u = u - (a*dt/dx) * (u - u_prev);
    end

    residual = u - u_ref;
    L2_errors(idx) = sqrt(dx) * norm(residual, 2);
    Linf_errors(idx) = norm(residual, Inf);
end

% Plot convergence
figure;
loglog(cfl_scales, L2_errors, 'o-', ...
       cfl_scales, cfl_scales, '--', ...
       cfl_scales, cfl_scales.^2, ':', 'LineWidth', 1.5);
legend('L^2 error','O(k)','O(k^2)', 'Location','northwest');
xlabel('Time Step'); ylabel('Error'); title('FTBS - L^2 Time Convergence');

figure;
loglog(cfl_scales, Linf_errors, 'o-', ...
       cfl_scales, cfl_scales, '--', ...
       cfl_scales, cfl_scales.^2, ':', 'LineWidth', 1.5);
legend('L^\infty error','O(k)','O(k^2)', 'Location','northwest');
xlabel('Time Step'); ylabel('Error'); title('FTBS - L^\infty Time Convergence');

% Convergence rates
rates_L2 = log2(L2_errors(1:end-1) ./ L2_errors(2:end));
rates_Inf = log2(Linf_errors(1:end-1) ./ Linf_errors(2:end));
disp('L2 convergence rates:'); disp(rates_L2);
disp('Linf convergence rates:'); disp(rates_Inf);




%%

% Spatial Convergence

clc; clear; close all;

a = 0.5;
T = 1;                   
dt_fixed = 0.00125;      
Nt = round(T / dt_fixed);
dt = T / Nt;             

% Spatial resolutions
grid_points = [50, 100, 200, 400, 800];
L2_errors = zeros(size(grid_points));
Linf_errors = zeros(size(grid_points));

% Reference solution 
N_ref = 2500;
dx_ref = 2 / N_ref;
x_ref = linspace(-1, 1, N_ref+1); x_ref(end) = [];
u_ref = sin(2*pi*x_ref);

for n = 1:Nt
    u_prev = [u_ref(end), u_ref(1:end-1)];
    u_ref = u_ref - (a * dt / dx_ref) * (u_ref - u_prev);
end

% Loop over spatial resolutions
for idx = 1:length(grid_points)
    N = grid_points(idx);
    dx = 2 / N;
    x = linspace(-1, 1, N+1); x(end) = [];
    u = sin(2*pi*x);

    for n = 1:Nt
        u_prev = [u(end), u(1:end-1)];
        u = u - (a * dt / dx) * (u - u_prev);
    end

    % Interpolate reference onto coarse grid
    u_ref_interp = interp1(x_ref, u_ref, x, 'spline');

    residual = u - u_ref_interp;
    L2_errors(idx) = sqrt(dx) * norm(residual, 2);
    Linf_errors(idx) = norm(residual, Inf);
end

% Plot convergence 
figure;
loglog(grid_points, L2_errors, 'o-', ...
       grid_points, 1./grid_points, '--', ...
       grid_points, 1./(grid_points.^2), ':', 'LineWidth', 1.5);
legend('L^2 error','O(h)','O(h^2)', 'Location','northeast');
xlabel('N (spatial grid points)'); ylabel('Error');
title('FTBS - L^2 Spatial Convergence'); grid on;

figure;
loglog(grid_points, Linf_errors, 'o-', ...
       grid_points, 1./grid_points, '--', ...
       grid_points, 1./(grid_points.^2), ':', 'LineWidth', 1.5);
legend('L^\infty error','O(h)','O(h^2)', 'Location','northeast');
xlabel('N (spatial grid points)'); ylabel('Error');
title('FTBS - L^\infty Spatial Convergence'); grid on;

% Print observed convergence rates ===
rates_L2 = log2(L2_errors(1:end-1) ./ L2_errors(2:end));
rates_Inf = log2(Linf_errors(1:end-1) ./ Linf_errors(2:end));
disp('L2 convergence rates:'); disp(rates_L2);
disp('Linf convergence rates:'); disp(rates_Inf);

