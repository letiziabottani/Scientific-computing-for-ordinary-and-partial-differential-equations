clear all
clc
%close all

u = @(x,y) sin(4*pi*(x+y))+cos(4*pi*x.*y);
f = @(x,y) -16*pi^2*(2*sin(4*pi*(x+y))+cos(4*pi*x.*y).*(x.^2+y.^2));

m = [10 15 20 40 50  ];
err=zeros(length(m),1);
h=1./(m+1);
rho=zeros(length(m),1);
for i=1: length(m)
	A = poisson5(m(i));
    

	F = form_rhs(m(i),f,u);
    
    rho(i)=2*(sqrt(condest(A))-1)/(sqrt(condest(A))+1);
    
	u_values = pcg(-A,-F,1e-6,1000); 

	u_matrix=reshape(u_values,m(i),m(i));

	x = linspace(h(i),1-h(i),m(i));
	y = linspace(h(i),1-h(i),m(i));
	[X,Y]=meshgrid(x,y);
	err(i)=max(max((abs(u_matrix-u(X,Y)))));
end 

figure;surf(X,Y,u_matrix)

figure; surf(X,Y,u(X,Y))
rho

figure; 
loglog(h,err,h,h,h,h.^2,h,h.^3)
legend("Err","h","h^2","h^3")
title("Convergence test for 5-point scheme")
grid on 
xlabel("h")
ylabel("Error")

% we estimate to have order 2, we also see that we have a very high modulus
% of the eigenvalues => A is bad conditioned 
%% convergence history

figure;
hold on
colors = lines(length(m));
for i = 1:length(m)
    A = poisson5(m(i));
    F = form_rhs(m(i),f,u);
    
    % Solve using PCG and store residuals
    [u_values,~,~,~,resvec] = pcg(-A,-F,1e-6,100);
    
    % Plot convergence history
    semilogy(0:length(resvec)-1, resvec, 'o-', 'DisplayName', ['m=' num2str(m(i))]);
    
    % Error computation
    u_matrix = reshape(u_values, m(i), m(i));
    x = linspace(h(i),1-h(i),m(i));
    y = linspace(h(i),1-h(i),m(i));
    [X,Y] = meshgrid(x,y);
    err(i) = max(max(abs(u_matrix - u(X,Y))));
end

xlabel('Iteration')
ylabel('Residual Norm')
title('PCG Convergence History')
legend show
grid on

