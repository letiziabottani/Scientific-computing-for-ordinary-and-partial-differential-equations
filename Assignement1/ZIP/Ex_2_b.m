clear all
clc
close all

u = @(x,y) sin(4*pi*(x+y))+cos(4*pi*x.*y);
f = @(x,y) -16*pi^2*(2*sin(4*pi*(x+y))+cos(4*pi*x.*y).*(x.^2+y.^2));

%m=80;
m = [10 20 50 80 100 150];
err=[];
h=1./(m+1);

for i=1: length(m)
	A = poisson5(m(i));
	F = form_rhs(m(i),f,u);
	u_values = A\F;

	u_matrix=reshape(u_values,m(i),m(i));

	x = linspace(h(i),1-h(i),m(i));
	y = linspace(h(i),1-h(i),m(i));
	[X,Y]=meshgrid(x,y);
	err(i)=max(max((abs(u_matrix-u(X,Y)))));
end 

figure;surf(X,Y,u_matrix)

figure; surf(X,Y,u(X,Y))

figure; 
 surf(X,Y,u_matrix -u(X,Y))

figure; 
loglog(h,err,h,h,h,h.^2,h,h.^3)
legend("Err","h","h^2","h^3")
title("Convergence test for 5-point scheme")
grid on 
xlabel("h")
ylabel("Error")

