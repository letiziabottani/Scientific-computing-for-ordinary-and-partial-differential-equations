function Unew = FTCS_iteration(U,gl,gr,k,e)

N = length(U)-1;

h = 2/N;

i=(1:(N-1))+1;

AU = U(i-1)-2.*U(i)+U(i+1);

Unew = zeros(1,N+1);

Unew(i) = U(i) + (k*e/(h^2)).*AU;

Unew([1,end]) = [gl,gr];

end
