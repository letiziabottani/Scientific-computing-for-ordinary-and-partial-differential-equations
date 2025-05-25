function Unew = Burger_iteration_highorder(U, gl, gr, k, e, x)
% High-order finite difference time step for viscous Burgers' equation
% using a 5-point (4th-order accurate) central difference stencil
% on a uniform grid

N = length(U) - 1;
Unew = zeros(1, N+1);
Unew(1) = gl;
Unew(end) = gr;

stencil_size = 7;
half = floor(stencil_size / 2);

for i = (1 + half):(N + 1 - half)
    x_stencil = x(i - half : i + half);
    u_stencil = U(i - half : i + half);

    % 1st derivative (advection term): u * du/dx
    coeffs_d1 = fdcoeffV(1, x(i), x_stencil);
    ux = coeffs_d1 * u_stencil';

    % 2nd derivative (diffusion term)
    coeffs_d2 = fdcoeffV(2, x(i), x_stencil);
    uxx = coeffs_d2 * u_stencil';

    % Burgers update
    Unew(i) = U(i) + k * (e * uxx - U(i) * ux);
end
end
