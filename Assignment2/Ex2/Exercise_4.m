% Exercise 4

clc
clear all
close all

% exact solution and RHS

%e=0.01/pi;
e=0.01;

u=@(x,t) -tanh((x+0.5-t)/(2*e))+1;

%approximation perameters
N = 2^8;
h = 2/N;
k_scale = 1/1000; %problems when e<<k_scale
k = h^2/(2*e)*k_scale;

max_time = 2;

max_j = round((max_time*N^2*e/2)/k_scale);




%% run
x = -1+h*(0:N);

U = zeros(max_j+1,N+1); %define matrix for U

U(1,:) = u(x,0); %initial values


for j = 1:max_j
    
    t = j*k;
    

    Unew = Burger_iteration(U(j,:),u(-1,t),u(1,t),k,e);
    
    

    U(j+1,:)=Unew;

    

end
%% analysis
surf(x,(0:max_j).*k,U, EdgeColor="none")
ylabel("t")
xlabel("x")


%% test comp

Time = (0:max_j)*k;

[X,Y]=meshgrid(x,Time);

surf(X,Y,u(X,Y),EdgeColor="none")

%% demonstrate convergenc

N_s = 2.^(3:6);

h_vals = 2./N_s;

e=0.1;

errors = zeros(1,length(N_s));

max_time=2;

k_scale = 1/2;

for i = 1:length(N_s)
    N= N_s(i);
    h = h_vals(i);
    
    k = h^2/(2*e)*k_scale;

    max_j = (max_time*N^2*e/2)/k_scale;

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
    %errors(i) = norm(residuals(2:end,2:(end-1)));
    errors(i) = norm(residuals);

end
%% plot
loglog(h_vals,errors,h_vals,h_vals.^2)

