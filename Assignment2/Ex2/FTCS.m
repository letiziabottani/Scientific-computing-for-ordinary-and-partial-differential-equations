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



%% demonstrate convergence (space)

N_s = 2.^(3:6);

h_vals = 2./N_s;



errors = zeros(1,length(N_s));

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
    %errors(i) = norm(residuals(2:end,2:(end-1)));
    errors(i) = norm(residuals);

end
%% plot
loglog(h_vals,errors,h_vals,h_vals.^2)

xlabel("h")
ylabel("error (2-norm)")


%% demonstrate convergenc (of time)


N= 2.^(6);
h = 2./N;


k_scale = 1./((2).^(2:8));

k_s = (h)^2/(2*e).*k_scale;


errors = zeros(1,length(k_s));
inf_norm = zeros(1,length(k_s));
max_e = zeros(1,length(k_s));

max_time=2;

x = -1+h*(0:N);
    
for i = 1:length(errors)
    
    k=k_s(i);

    max_j = round(max_time/k);
    
    
    U = zeros(max_j+1,N+1);
    
    U(1,:) = u(x,0); %initial values
    
    for j = 1:max_j
        t = j*k;
         
        U(j+1,:) = FTCS_iteration(U(j,:),u(-1,t),u(1,t),k,e);
    end

    Time = (0:max_j)*k;

    [X,Y] = meshgrid(x,Time);

    true_u = u(X,Y);
    residuals = true_u - U;

    %errors(i) = norm(residuals(2:end,2:(end-1)));
    errors(i) = norm(residuals);
    %inf_norm(i) = norm(residuals,inf);
    max_e(i) = max(abs(residuals(end,:)),[],"all");
    

end

%% plot
loglog(k_s,errors,k_s,k_s)
%loglog(k_s,errors,k_s,k_s,k_s,inf_norm)

xlabel("k")
ylabel("error")

%% analysis
figure(1)
surf(x,Time,U,EdgeColor="none")
ylabel("t")
xlabel("x")

%% test comp

figure(2)
surf(X,Y,true_u,EdgeColor="none")

%% resiluals?

surf(x,Time,residuals,EdgeColor="none")