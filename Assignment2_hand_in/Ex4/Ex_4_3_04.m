% Parameters
e = 0.01 / pi;                    
N = 2^9;                         
h = 2 / N;                       
x = -1 + h * (0:N);              

k_scale = 1 / 100;                
k = h^2 / (2 * e) * k_scale;      
t_target = 1.6037 / pi;          
max_j = round(t_target / k);     

% Initial condition
U = zeros(max_j+1, N+1);
U(1, :) = -sin(pi * x);

% Burger Iteration
for j = 1:max_j
    U(j+1, :) = Burger_iteration(U(j, :), 0, 0, k, e);
end

% Centered difference
[~, ix0] = min(abs(x)); 
ux_approx = (U(end, ix0 + 1) - U(end, ix0 - 1)) / (2 * h);

% Result
fprintf('Estimated ∂u/∂x at x = 0, t = %.5f: %.5f\n', t_target, ux_approx);




figure;
plot(x, U(end, :), 'b-', 'LineWidth', 2);
title(sprintf('Solution u(x,t) at t = %.5f', t_target));
xlabel('x');
ylabel('u(x,t)');
grid on;
