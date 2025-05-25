function Unew = Burger_iteration(U, gl, gr, k, e)
% Performs one explicit update step for the viscous Burgers' equation
% using central differences and Forward Euler time integration.

N = length(U) - 1;          % number of intervals
i = 2:N;                    % interior indices

% Precompute constants
diffusion_coeff = (k * e * N^2) / 4;   % corresponds to k*e/h^2
advection_coeff = (k * N) / 4;         % corresponds to k/(2h)

% Compute second spatial derivative (diffusion term)
d2u = U(i-1) - 2 * U(i) + U(i+1);

% Compute advection term: u * (du/dx)
du = (U(i+1) - U(i-1)) .* U(i);

% Update interior points
Unew = zeros(1, N+1);
Unew(i) = U(i) + diffusion_coeff * d2u - advection_coeff * du;

% Apply Dirichlet boundary conditions
Unew([1, end]) = [gl, gr];
end