u = @(x,y) sin(4*pi*(x+y))+cos(4*pi*x*y);
f = @(x,y) -16*pi^2*(2*sin(4*pi*(x+y))+cos(4*pi*x*y)*(x^2+y^2));


m = 10;

A = poisson5(m);

F = form_rhs(m,f,u);


u_values = A\F;

u_matrix=reshape(u_values,m,m);

h=1/(m+1);
x = linspace(h,1-h,m);
y = linspace(h,1-h,m);
[X,Y]=meshgrid(x,y);

figure;surf(X,Y,u_matrix)

figure; surf(X,Y,u(X,Y))