function [y_high,y_low,err] = Burger_RK_iteration(U,t,gl,gr,k,e)

% gl and gr given as functions of t

%initialize
N = length(U)-1;
h = 2/N;
i=(1:(N-1))+1;

% a coefficients for intermediate stages
a21 = 1/2;    a31 = -1;    a32 = 2;

% b coefficients (3rd order solution)
b1 = 1/6;    b2 = 2/3;    b3 = 1/6;

% b_hat coefficients (2nd order solution)
b1_hat = 1/4;    b2_hat = 1/2;    b3_hat = 1/4;


% d_i = b_i - b_hat_i (used to estimate the local error)
% d1=b1-b1_hat;    d2=b2-b2_hat;    d3=b3-b3_hat;

% c coefficients
c1=0;     c2=1/2;     c3=1; 

%spatial methods
    function f = f(t,U)
        %dxx part
        dxx = U(i-1) - 2*U(i) + U(i+1);
        
        % dx(u)*u
        dx_u = (U(i+1)-U(i-1)).*U(i);

        f = zeros(1,N+1);

        f(i) = (e*(N^2)/4)*dxx - (N/4)*dx_u;
        %bounduarry conditions
        f([1,end]) = [gl(t),gr(t)];
    end

% calculate Y[i] = f(t[i],xi[i]) (insted of just xi[i])
Y1 = f(t,U);
Y2 = f(t + c2*k, U + a21*k.*Y1);
Y3 = f(t + c3*k, U + a31*k*Y1 + a32*k*Y2);

% 3rd-order solution
y_high = U + k * (b1*Y1 + b2*Y2 + b3*Y3);

% 2nd-order solution
y_low = U + k * (b1_hat*Y1 + b2_hat*Y2 + b3_hat*Y3);

% fix 
y_high([1,end]) = [gl(t+k),gr(t+k)];
y_low([1,end]) = [gl(t+k),gr(t+k)];


err=norm(y_high-y_low);

end