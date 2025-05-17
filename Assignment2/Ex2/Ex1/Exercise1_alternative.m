function rk23_flame_solver()
    % Initial condition
    delta = 0.02;
    y0 = delta;
    t0 = 0;
    tf = 2 / delta;

    % Time stepping setup
    N = 10000; % increase this if needed
    h = (tf - t0) / N;
    t = linspace(t0, tf, N+1);
    y = zeros(1, N+1);
    y_hat = zeros(1, N+1); % 2nd order approx
    error_est = zeros(1, N+1);

    y(1) = y0;
    y_hat(1) = y0;

    % RHS
    f = @(t, y) y.^2 - y.^3;

    for n = 1:N
        tn = t(n);
        yn = y(n);

        % RK stages
        k1 = f(tn, yn);
        k2 = f(tn + h/2, yn + h * (1/2) * k1);
        k3 = f(tn + h, yn - h * k1 + 2 * h * k2);

        % Third-order solution
        y(n+1) = yn + h * (1/6 * k1 + 2/3 * k2 + 1/6 * k3);

        % Second-order embedded solution
        y_hat(n+1) = yn + h * (1/4 * k1 + 1/2 * k2 + 1/4 * k3);

        % Error estimate (norm if vector valued)
        error_est(n+1) = abs(y(n+1) - y_hat(n+1));
    end

    % Plot solution
    figure;
    plot(t, y, 'b', 'DisplayName', 'RK23 - 3rd Order'); hold on;
    plot(t, y_hat, 'g--', 'DisplayName', 'RK23 - 2nd Order');
    xlabel('t'); ylabel('y(t)');
    title(['RK23 Solution, \delta = ', num2str(delta)]);
    legend(); grid on;

    % Plot error
    figure;
    semilogy(t, error_est, 'r');
    xlabel('t'); ylabel('Estimated Local Error');
    title('Local Error Estimate (|y^{(3)} - y^{(2)}|)');
    grid on;
end





