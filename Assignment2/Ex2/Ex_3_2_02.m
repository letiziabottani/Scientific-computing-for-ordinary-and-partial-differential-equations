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

    
    ctx = 0.7;                  
    dt = ctx * dx / a;
    Nt = round(T/dt);         
    dt = T / Nt;              

    
    u = sin(2*pi*x);


    for n = 1:Nt
        u_prev = [u(end), u(1:end-1)]; 
        u = u - (a*dt/dx) * (u - u_prev);
    end

    u_exact = sin(2*pi*(x - a*T));

    if idx == length(N_values)
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



figure;
loglog(N_values, errors, '-o', 'LineWidth', 2); hold on;


N_ref = N_values(1);            
error_ref = errors(1);          
theoretical = error_ref * (N_ref ./ N_values);  

loglog(N_values, theoretical, 'k--', 'LineWidth', 2);  

legend('Numerical Error', 'Theoretical 1. Order');
xlabel('Nr of grid points');
ylabel('Error');
title('Convergence of FTBS vs. Theoretical 1. Order');
grid on;