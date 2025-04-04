clear all
clc
close all
%% 1_a_a
syms h 
A=[1,1,1,1,1; 4, 3, 2, 1, 0; 8 9/2 2 1/2 0;-32/3 -9/2 -4/3 -1/6 0; 32/3 27/8 2/3 1/24 0];
b=[0 0 1/h.^2 0 0 ]';
x=A\b
%% 1_a_b
syms h 
A=[1,1,1,1,1; -2 -1 0 1 2; 2 1/2 0 1/2 2; -4/3 -1/6 0 1/6 4/3; 2/3 1/24 0 1/24 2/3];
b=[0 0 1/h.^2 0 0 ]';
x=A\b

%% 1_c
u=@(x) exp(cos(x))

h = 0.2
% Define parameters
k = 2; % Second derivative
xbar = 0; % Point where we approximate the derivative
%x = [-h 0 h]; % Equidistant stencil (should be ordered symmetrically)
x_1 = [-2*h -h 0 h +2*h]
n_1 = length(x_1);

% Compute finite difference coefficients
%c = fdcoeffV(k, xbar, x);
c_1 = fdcoeffF(k, xbar, x_1);

% Compute the approximation using the stencil
u_values = u(x_1); % Evaluate function at stencil points
%u_approx = c * u_values'; % Apply finite difference formula
u_approx_1 = c_1 * u_values';

% Display result
%fprintf('Approximation of u''''(0): %f\n', u_approx);
fprintf('Approximation of u''''(0): %f\n', u_approx_1);

%% 1_d
u=@(x)exp(cos(x))
u_prime_prime=@(x)(sin(x).^2).*exp(cos(x))-cos(x).*exp(cos(x))
p=[2:15];
h = 1./(2.^p);
% Define parameters
k = 2; % Second derivative
xbar = 0; % Point where we approximate the derivative
err=[]
for ii= 1:length(p)
    %x = [-2*h(ii) -h(ii) 0 h(ii) +2*h(ii)]
    x=[-4*h(ii) -3*h(ii) -2*h(ii) -h(ii) 0]
    n = length(x);

    % Compute finite difference coefficients
    c= fdcoeffF(k, xbar, x);

    % Compute the approximation using the stencil
    u_values = u(x); % Evaluate function at stencil point
    u_approx = c * u_values';
    err(ii)=max(abs(u_prime_prime(xbar)-u_approx));

end 

figure()
loglog(h,err,h,h.^3,h,h.^4)
legend("error","h^3","h^4")
grid on


