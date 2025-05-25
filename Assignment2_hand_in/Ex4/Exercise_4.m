% Exercise 4

clc
clear all
close all

% exact solution and RHS

%e=0.01/pi;
e=0.1;

u=@(x,t) -tanh((x+0.5-t)/(2*e))+1;

%approximation perameters
N = 2^6;
h = 2/N;
k_scale = 1/100; %problems when e<<k_scale?
k = h^2/(2*e)*k_scale;

max_time = 2;

max_j = round(max_time/k);




%% run
x = -1+h*(0:N);

U = zeros(max_j+1,N+1); %define matrix for U

U(1,:) = u(x,0); %initial values


for j = 1:max_j
    
    t = j*k;
    

    Unew = Burger_iteration(U(j,:),u(-1,t),u(1,t),k,e);
    
    U(j+1,:)=Unew;

end
%% plot result
surf(x,(0:max_j).*k,U, EdgeColor="none")
ylabel("t")
xlabel("x")

%% compare to true solution

Time = (0:max_j)*k;

[X,Y]=meshgrid(x,Time);

surf(X,Y,u(X,Y),EdgeColor="none")


norm(U-u(X,Y),2)

%% demonstrate convergence (space)

N_s = 2.^(4:7);

h_vals = 2./N_s;

errors = zeros(1,length(N_s));

max_time=2;

k_scale = 1/100;

k = (h_vals(end))^2/(2*e)*k_scale;

max_j = round(max_time/k);

for i = 1:length(N_s)
    N= N_s(i);
    h = h_vals(i);
    
    %k = (h)^2/(2*e)*k_scale;

    %max_j = round(max_time/k);

    x = -1+h*(0:N);
    
    U = zeros(max_j+1,N+1);
    
    Unew = u(x,0); %initial values
    
    for j = 1:max_j
        t = j*k;
        Unew = Burger_iteration(Unew,u(-1,t),u(1,t),k,e);
        
    end

    %Time = (0:max_j)*k;

    %[X,Y] = meshgrid(x,Time);

    residuals = Unew - u(x,max_time);

    %errors(i) = norm(residuals(2:end,2:(end-1)));
    errors(i) = norm(residuals,inf);

end
%% plot convergance rate
%loglog(h_vals,errors,h_vals,h_vals,h_vals,h_vals.^2)

loglog(h_vals, errors, 'o-', h_vals, h_vals, '--', h_vals,h_vals.^2, ':', 'LineWidth', 1.5);
xlim([min(h_vals),max(h_vals)])
legend('L^\infty error','O(h)','O(h^2)', 'Location','northwest');
xlabel('h'); ylabel('Error'); title('Convergence in L^\infty norm');
%legend('L^2 error','O(h)','O(h^2)', 'Location','northwest');
%xlabel('h'); ylabel('Error'); title('Convergence in L^2 norm');

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