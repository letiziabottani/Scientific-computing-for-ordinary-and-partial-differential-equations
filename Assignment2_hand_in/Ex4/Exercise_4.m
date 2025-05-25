% Exercise 4

clc
clear all
close all

% exact solution and RHS

%e=0.01/pi;
e=0.1;

u=@(x,t) -tanh((x+0.5-t)/(2*e))+1;

%approximation perameters
N = 2^8;
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

%% demonstrate convergence

N_s = 2.^(3:6);

h_vals = 2./N_s;

errors = zeros(1,length(N_s));

max_time=2;

k_scale = 1/100;

k = (h_vals(end))^2/(2*e)*k_scale;

max_j = round(max_time/k);

for i = 1:length(N_s)
    N= N_s(i);
    h = h_vals(i);
    
    k = (h)^2/(2*e)*k_scale;

    max_j = round(max_time/k);

    x = -1+h*(0:N);
    
    U = zeros(max_j+1,N+1);
    
    U(1,:) = u(x,0); %initial values
    
    for j = 1:max_j
        t = j*k;
        Unew = Burger_iteration(U(j,:),u(-1,t),u(1,t),k,e);
        U(j+1,:)=Unew;
    end

    Time = (0:max_j)*k;

    [X,Y] = meshgrid(x,Time);

    residuals = U - u(X,Y);

    %errors(i) = norm(residuals(2:end,2:(end-1)));
    errors(i) = norm(residuals);

end
%% plot convergance rate
loglog(h_vals,errors,h_vals,h_vals,h_vals,h_vals.^2)

