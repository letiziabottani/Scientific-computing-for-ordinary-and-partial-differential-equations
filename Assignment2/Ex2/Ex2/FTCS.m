clc
clear all
close all

% exact solution and RHS

e = 0.1; %epsilon
Nodes = 3;

al = [1, 4, 16];
a = ones(1,3);
b = zeros(1,3);

u=@(x,t) exp(-e*al(1).^2.*t).*(a(1).*cos(al(1)*x)+b(1).*sin(al(1).*x)) + ...
         exp(-e*al(2).^2.*t).*(a(2).*cos(al(2)*x)+b(2).*sin(al(2).*x)) + ...
         exp(-e*al(3).^2.*t).*(a(3).*cos(al(3)*x)+b(3).*sin(al(3).*x));

%approximation perameters
N = 100;
h = 2/N;
k = h^2/(2*e);
max_time = 2;

%% run
x = -1+h*(0:N);

max_j = round(max_time/k);

U = zeros(max_j+1,N+1);

U(1,:) = u(x,0); %initial values


for j = 1:max_j
    
    t = j*k;
    
    Unew = FTCS_iteration(U(j,:),u(-1,t),u(1,t),k,e);
    
    U(j+1,:)=Unew;

    

end
%% analysis
surf(U,EdgeColor="none")
ylabel("t")
xlabel("x")


%% test comp

Time = (0:max_j)*k;

[X,Y]=meshgrid(x,Time);

surf(X,Y,u(X,Y),EdgeColor="none")




%% Compare solutions at a fixed time point

% Final time
t_final = max_time;


idx_final = size(U, 1);  

u_num_final = U(idx_final, :);
u_exact_final = u(x, t_final);

% Compare
figure
plot(x, u_exact_final, 'k--', 'LineWidth', 2)
hold on
plot(x, u_num_final, 'b-', 'LineWidth', 2)
hold off

% Labels and formatting
xlabel('x')
ylabel('u(x, t})')
title(['Comparison at Final Time t = ', num2str(t_final)])
legend('Exact solution', 'Numerical solution', 'Location', 'best')
grid on





%% Spatial convergence

N_s = 2.^(3:6);

h_vals = 2./N_s;



errors_L2 = zeros(1, length(N_s));
errors_Inf = zeros(1, length(N_s));

max_time=2;

k_scale = 1/10;
k = (h_vals(end))^2/(2*e)*k_scale;

max_j = round(max_time/k);
    
for i = 1:length(N_s)
    N= N_s(i);
    h = h_vals(i);
    
    

    x = -1+h*(0:N);
    
    U = zeros(max_j+1,N+1);
    
    U(1,:) = u(x,0); %initial values
    
    for j = 1:max_j
        t = j*k;
        Unew = FTCS_iteration(U(j,:),u(-1,t),u(1,t),k,e);
        U(j+1,:)=Unew;
    end

    Time = (0:max_j)*k;

    [X,Y] = meshgrid(x,Time);

    true_u = u(X,Y);
    residuals = true_u - U;
    errors_L2(i) = norm(residuals);               
    errors_Inf(i) = norm(residuals, Inf); 

end

%% 
figure;

loglog(h_vals, errors_L2, 'o-', 'LineWidth', 1.5, 'DisplayName', 'L^2 error'); hold on;
loglog(h_vals, h_vals, '--', 'LineWidth', 1.5, 'DisplayName', 'O(h)');
loglog(h_vals, h_vals.^2, ':', 'LineWidth', 1.5, 'DisplayName', 'O(h^2)');

xlabel('h (spatial step size)');
ylabel('Error');
title('FTCS - L^2 Spatial Convergence');
legend('Location', 'northwest');
grid on;

%% 
figure;

loglog(h_vals, errors_Inf, 's-', 'LineWidth', 1.5, 'DisplayName', 'L^\infty error'); hold on;
loglog(h_vals, h_vals, '--', 'LineWidth', 1.5, 'DisplayName', 'O(h)');
loglog(h_vals, h_vals.^2, ':', 'LineWidth', 1.5, 'DisplayName', 'O(h^2)');

xlabel('h (spatial step size)');
ylabel('Error');
title('FTCS - L^\infty Spatial Convergence');
legend('Location', 'northwest');
grid on;




%% Temporal convergence 



% Spatial setup (fixed grid)
N = 100;
h = 2/N;
x = -1 + h*(0:N);



% Time convergence setup
k_scales = [1/4, 1/8, 1/16, 1/32];   
L2_errors = zeros(size(k_scales));
Linf_errors = zeros(size(k_scales));

% Reference solution 
k_ref = h^2 / (2*e) * (1/128); 
max_j_ref = round(max_time / k_ref);
U_ref = zeros(max_j_ref+1, N+1);
U_ref(1,:) = u(x, 0);

for j = 1:max_j_ref
    t = j * k_ref;
    Unew = FTCS_iteration(U_ref(j,:), u(-1,t), u(1,t), k_ref, e);
    U_ref(j+1,:) = Unew;
end

u_ref_final = U_ref(end, :);

for idx = 1:length(k_scales)
    k = h^2 / (2*e) * k_scales(idx);
    max_j = round(max_time / k);
    
    U = zeros(max_j+1, N+1);
    U(1,:) = u(x, 0);

    for j = 1:max_j
        t = j * k;
        Unew = FTCS_iteration(U(j,:), u(-1,t), u(1,t), k, e);
        U(j+1,:) = Unew;
    end

    u_num_final = U(end, :);
    residual = u_num_final - u_ref_final;

    L2_errors(idx) = sqrt(h) * norm(residual, 2);
    Linf_errors(idx) = norm(residual, Inf);
end

%%
% Plot time convergence (L2 and Linf)
figure;
loglog(k_scales, L2_errors, 'o-', ...
       k_scales, k_scales, '--', ...
       k_scales, k_scales.^2, ':', 'LineWidth', 1.5);
legend('L^2 error','O(k)','O(k^2)', 'Location','northwest');
xlabel('Time step scale'); ylabel('Error'); title('FTCS - L^2 Time Convergence');

figure;
loglog(k_scales, Linf_errors, 'o-', ...
       k_scales, k_scales, '--', ...
       k_scales, k_scales.^2, ':', 'LineWidth', 1.5);
legend('L^\infty error','O(k)','O(k^2)', 'Location','northwest');
xlabel('Time step scale'); ylabel('Error'); title('FTCS - L^\infty Time Convergence');

% Convergence rates
rates_L2 = log2(L2_errors(1:end-1) ./ L2_errors(2:end));
rates_Inf = log2(Linf_errors(1:end-1) ./ Linf_errors(2:end));
disp('L2 convergence rates:'); disp(rates_L2);
disp('Linf convergence rates:'); disp(rates_Inf);


