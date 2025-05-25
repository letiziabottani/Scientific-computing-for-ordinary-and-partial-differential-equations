function [Unew,dx] = Burger_NU_grid(U,gl,gr,k,e,h)

N = length(U)-1;

%initialize values
i=(1:(N-1))+1;

h_minus = h(1:end-1);
h_plus = h(2:end);
%repeating values
h_add = h_plus + h_minus;

h_mult = h_minus.*h_plus;

h_am_plus = h_add.*h_plus;

h_am_minus = h_add.*h_minus;

%dxx(u) part
dxx = U(i-1)./h_am_minus - U(i)./h_mult + U(i+1)./h_am_plus;

% dx(u)
dx = -U(i-1).*h_plus./h_am_minus + ...
    U(i).*(h_plus-h_minus)./h_mult + ...
    U(i+1).*h_minus./h_am_plus;

Unew = zeros(1,N+1);

Unew(i) = U(i) + (k*e*2)*dxx - (k)*dx.*U(i); %combine

Unew([1,end]) = [gl,gr]; %add boundary conditions 

end