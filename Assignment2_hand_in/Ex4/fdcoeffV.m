function c = fdcoeffV(k, xbar, x)
n = length(x);
A = ones(n,n);
xrow = (x(:)' - xbar);     % displacements
for i = 2:n
    A(i,:) = xrow.^(i-1) / factorial(i-1);
end
b = zeros(n,1);
b(k+1) = 1;                 % derivative order
c = A\b;
c = c';
end
