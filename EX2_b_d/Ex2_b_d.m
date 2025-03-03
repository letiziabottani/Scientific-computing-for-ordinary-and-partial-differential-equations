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

%% 