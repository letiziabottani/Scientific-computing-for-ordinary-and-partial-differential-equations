% Exercise 4

clc
clear all
close all
%%

% exact solution and RHS

e=0.01/pi;


u_initial=@(x) -sin(pi*x);

%approximation perameters
N = 2^10;
h = 2/N;
k_scale = 1/100; 
k = h^2/(2*e)*k_scale;

max_time = 1.6037/pi;

max_j = ceil(max_time/k);

k = max_time/max_j;


%% run
x = -1+h*(0:N);

U = zeros(max_j+1,N+1); %define matrix for U

U(1,:) = u_initial(x); %initial values

i_0 =N/2+1;


for j = 1:max_j
    
    t = j*k;
    
    U(j+1,:) = Burger_iteration(U(j,:),0,0,k,e);
end
dx_uni_0 = (U(:,i_0+1)-U(:,i_0-1))/(2*h); %uniform grid slope at x=0

%% delta x of u
i_0 =N/2+1;
dx0 = (U(end,i_0+1)-U(end,i_0-1))/(2*h); 
[dx0,dx0+152.00516]

%% plot result
plot = surf(x,(0:100:max_j).*k,U(1:100:end,:), EdgeColor="none");
set(gca,'xdir','reverse','ydir','reverse')
ylabel("t")
xlabel("x")
title('Numerical Solution U(x,t) for N=2^{10}'); xlabel('x'); ylabel('t');

%%
plot(x,U(end,:))
xlabel("x")
xlim([-1,1])

%% plot dx of u at t_end
dx = (U(end,(2:N)+1)-U(end,(2:N)-1))/(2*h);
plot(x(2:N),dx)
title('Numerical Solution dx(U(x,t_end))'); xlabel('x'); ylabel('dx');
min(dx)
%% non- uniform space grid
% keep k reduse N

N_non=2^7;

x_initialise = -1+2*(0:N_non)/N_non;

%alpha = 0.72; % N 8 (local j)

%alpha = 0.5; % N 8 (global j)

%alpha = 0.51; % N 7 (local j)

alpha = 0.51;

x_non = (1-alpha).*x_initialise.^3 + alpha.*x_initialise; %new x values

h_non = x_non(2:end)-x_non(1:end-1); %calculate step size

% new t
k_scale=1/100;

k_non = min(h_non)^2/(2*e)*k_scale; %use smallest h to calcualte k

max_j_non = ceil(max_time/k_non);

%max_j_non = max_j;

k_non = max_time/max_j_non; % round to match final value


%% non uniform solver

U_non = zeros(max_j_non+1,N_non+1); %define matrix for U

U_non(1,:) = u_initial(x_non); %initial values

i_0 =N_non/2+1;
dx_non_0 = zeros(max_j_non,1);

for j = 1:max_j_non
    
    t = j*k_non;
    
    [U_non(j+1,:),dx_non] = Burger_NU_grid(U_non(j,:),0,0,k_non,e,h_non);
    dx_non_0(j) = dx_non(N_non/2);
end
%%

[dx_non_0(end),dx_non_0(end)+152.00516]
plot((1:max_j_non)/max_j_non,dx_non_0,(1:max_j+1)/max_j,dx_uni_0)

%% plot result
plot(x_non,U_non(end,:),x,U(end,:))
xlim([-1,1])
xlabel("x")

%%
plot(x_non(2:N_non),dx_non)
min(dx_non)


%% 3d plot
surf(x_non,((max_j_non-100):max_j_non).*k_non,U_non((end-100):end,:), EdgeColor="none")
ylabel("t")
xlabel("x")
title('Numerical Solution U(x,t)'); xlabel('x'); ylabel('t');

%%
surf(x_non,(1:10000:max_j_non).*k_non,U_non(1:10000:end,:), EdgeColor="none")
ylabel("t")
xlabel("x")
title('Numerical Solution U(x,t)'); xlabel('x'); ylabel('t');

%% test different alpha

N_non=2^9;

x_initialise = -1+2*(0:N_non)/N_non;

alphas = (1:10)./10;

dx_alphas = zeros(1,length(alphas));

for i = 1:length(alphas)
    alpha = alphas(i);
    
    x_non = (1-alpha).*x_initialise.^3 + alpha.*x_initialise;
    
    h_non = x_non(2:end)-x_non(1:end-1);
    
    % new t
    k_scale=1/100;
    
    k_non = min(h_non)^2/(2*e)*k_scale;
    
    max_j_non = ceil(max_time/k_non);
    
    %max_j_non = max_j;
    
    k_non = max_time/max_j_non;
    
    U_non = u_initial(x_non); %initial values
    
    i_0 =N_non/2+1;
    
    
    for j = 1:max_j_non
        
        t = j*k_non;
        
        [U_non,dx_non] = Burger_NU_grid(U_non,0,0,k_non,e,h_non);
        
    end
    dx_alphas(i) = dx_non(N_non/2);
end
%% plot error

plot(alphas,abs(dx_alphas+152.00516),'-o',alphas,ones(1,length(alphas))*0.0028)
%ylim([0,0.04])
xlabel("\alpha");
ylabel("|error|");
legend('non-uniform N^{9}','target accuracy (N^{10})', 'Location','northeast');
title('Error of non-uniform solver for varying \alpha');