function Unew = Burger_iteration(U,gl,gr,k,e)

N = length(U)-1;

h = 2/N;

i=(1:(N-1))+1;

%dxx part
AC = U(i-1) - 2*U(i) + U(i+1);

% dx(u)*u
AX = (U(i+1)-U(i-1)).*U(i);

Unew = zeros(1,N+1);

Unew(i) = U(i) + (k*e*(N^2)/4)*AC - (k*N/4)*AX;

%AX = (U(i)-U(i-1)).*U(i);
%Unew(i) = U(i) + (k*e*(N^2)/4)*AC - (k*N/2)*AX;

Unew([1,end]) = [gl,gr];

end