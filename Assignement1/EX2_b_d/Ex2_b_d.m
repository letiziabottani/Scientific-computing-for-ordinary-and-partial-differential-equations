%Ex2b

%% c)

u=@(x,y) sin(4*pi*(x+y)) + cos(4*pi*x.*y);
f=@(x,y) -16*pi^2*(2*sin(4*pi*(x + y)) + cos(4*pi*x.*y).*(x.^2 + y.^2));

m=50;

A = poisson5(m);

F=form_rhs(m,f,u);

u_est=reshape(A\F,[m,m]);

h=1/(m+1);
interval=linspace(h,1-h,m);
    
[X,Y]=meshgrid(interval,interval);

figure; surf(X,Y,u_est)


%% d)

u=@(x,y) sin(4*pi*(x+y)) + cos(4*pi*x.*y);
f=@(x,y) -16*pi^2*(2*sin(4*pi*(x + y)) + cos(4*pi*x.*y).*(x.^2 + y.^2));

m=50;

A = poisson9(m);


F=form_rhs_9(m,f,u);

u_est=reshape(A\F,[m,m]);

h=1/(m+1);
interval=linspace(h,1-h,m);
    
[X,Y]=meshgrid(interval,interval);

figure; surf(X,Y,u_est);


figure; surf(X,Y,u(X,Y)-u_est)

%% convergence test 
clear all 
clc 
close all

u=@(x,y) sin(4*pi*(x+y)) + cos(4*pi*x.*y);
f=@(x,y) -16*pi^2*(2*sin(4*pi*(x + y)) + cos(4*pi*x.*y).*(x.^2 + y.^2));
laplacian_f=@(x,y) 256*(4*pi^2*sin(4*pi*(x+y))+(-1/4+(x.^2+y.^2).^2*pi^2).*cos(4*pi*x.*y)+2*pi.*y.*sin(4*pi*x.*y).*x)*pi^2;

m=[10 50 80 100 150 200];

for i=1:length(m)
    A = poisson9(m(i));
    F=form_rhs_9(m(i),f,u,laplacian_f);

    u_est=reshape(A\F,[m(i),m(i)]);

    h(i)=1/(m(i)+1);
    interval=linspace(h(i),1-h(i),m(i));
    
    [X,Y]=meshgrid(interval,interval);
    err(i)=max(max(abs(u_est-u(X,Y))));
end
figure; surf(X,Y,u_est);


figure; 
loglog(h,err,h,h.^2,h,h.^3,h,h.^4)
legend("Err",",h^2","h^3","h^4")
title("Convergence test for 9-point scheme")
grid on 
xlabel("h")
ylabel("Error")

