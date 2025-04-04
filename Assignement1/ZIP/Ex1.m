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
u=@(x)exp(cos(x));
u_prime_prime=@(x)(sin(x).^2).*exp(cos(x))-cos(x).*exp(cos(x));
h = 0.2;

k = 2; % Second derivative
xbar = 0; % Point where we approximate the derivative
%x = [-2*h -h 0 h +2*h];
x=[-4*h -3*h -2*h -h 0];
n = length(x);

% Compute finite difference coefficients
c = fdcoeffF(k, xbar, x);

% Compute the approximation using the stencil
u_values = u(x); 
u_approx = c * u_values';

% Display result
fprintf('Approximation of u''''(0): %f\n', u_approx);
err=u_prime_prime(xbar)-u_approx

%% 1_d
u=@(x)exp(cos(x))
u_prime_prime=@(x)(sin(x).^2).*exp(cos(x))-cos(x).*exp(cos(x))
p=[2:12];
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
loglog(h,err,'Linewidth',1.5)
hold on
loglog(h,h.^3,h,h.^4)
legend("Error","h^3","h^4")
xlabel("h")
ylabel("Error")
grid on

%% 1_e
%% calculate for p=2
u=@(x) exp(cos(x));

k=0;
h=0.1;
p=4;
x=linspace(-h/2,h/2,p);

c = [fdcoeffF(k,0,x)];

c*[u(x)]'

%% for loop p=2
it=15;
tau=zeros(it,1);
k=0;
p=2;
h=ones(1,it)./(2.^(1:it));
for i = 1:it
    x=linspace(-h(i)/2,h(i)/2,p);

    c = [fdcoeffF(k,0,x)];

    tau(i)=abs(c*[u(x)]'-exp(1));
end
figure()
loglog(h,tau,h,h.^2,h,h.^3)
legend("tau","h^2","h^3")
grid on

%% for loop p=3
it=26;
tau=zeros(it,1);
k=0;
p=3;
h=ones(1,it)./(2.^((1:it)/2));
for i = 1:it
    x=linspace(-h(i)/2,h(i)*3/2,p);

    c = [fdcoeffF(k,0,x)];

    tau(i)=abs(c*[u(x)]'-exp(1));
end
figure()
loglog(h,tau,h,h.^2,h,h.^3,h,h.^4)
legend("\tau","h^2","h^3","h^4")
xlabel('h') 
ylabel('\tau') 
grid on

%% for p=4

it=15;
tau=zeros(it,1);
k=0;
p=4;
h=ones(1,it)./(2.^(1:it));
for i = 1:it
    x=linspace(-h(i)*3/2,h(i)*3/2,p);

    c = [fdcoeffF(k,0,x)];

    tau(i)=abs(c*[u(x)]'-exp(1));
end
figure;
loglog(h,tau,h,h.^2,h,h.^3,h,h.^4)
legend("tau","h^2","h^3","h^4")
grid on

