% Define the function u(x)
%function y = u(x)
    %y = exp(cos(x));
%end

u_1=@(x) exp(cos(x))

h = 0.2
% Define parameters
k = 2; % Second derivative
xbar = 0; % Point where we approximate the derivative
%x = [-h 0 h]; % Equidistant stencil (should be ordered symmetrically)
x_1 = [-2*h -h 0 h +2*h]
n_1 = length(x_1);

% Compute finite difference coefficients
%c = fdcoeffV(k, xbar, x);
c_1 = fdcoeffV(k, xbar, x_1);

% Compute the approximation using the stencil
u_values = u_1(x_1); % Evaluate function at stencil points
%u_approx = c * u_values'; % Apply finite difference formula
u_approx_1 = c_1 * u_values';

% Display result
%fprintf('Approximation of u''''(0): %f\n', u_approx);
fprintf('Approximation of u''''(0): %f\n', u_approx_1);

% Function to compute finite difference coefficients
function c = fdcoeffV(k, xbar, x)
    n = length(x); % Ensure n is defined here
    A = ones(n, n);
    xrow = (x(:) - xbar)'; % Displacements as a row vector.
    
    for i = 2:n
        A(i, :) = (xrow .^ (i-1)) ./ factorial(i-1);
    end
    
    b = zeros(n, 1); % Right-hand side vector
    b(k+1) = 1; % Enforce the k-th derivative condition
    c = A \ b; % Solve system for coefficients
    c = c'; % Return as row vector
end