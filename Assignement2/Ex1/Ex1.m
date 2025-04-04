clear all
clc
close all

%% point1
%parameters and initialization 
%delta=0.02;
delta_val=[0.02, 0.01, 0.005, 0.001];
t0=0;
f = @(t, y) y^2 - y^3;

for i=1:length(delta_val)
    delta=delta_val(i);
    tf=2/delta; 
    y0=delta; 
      
    h=0.001; %initial step size (guess)
    t=t0; 
    y=y0; 
    
    T = t;
    Y = y;
    
    while t < tf
        if t + h > tf
            h = tf - t;
        end
        
        % RK23 step and error estimate
        [y_high,~,~] = rk23_step(f, t, y, h);
    
        t = t + h;
        y = y_high;

        T(end+1) = t;
        Y(end+1) = y;
    end
    
    [t_builtin, y_builtin] = ode23(f, [0 tf], delta);

    % Plot
    figure;
    plot(T, Y, 'b.-', 'DisplayName', 'Custom RK23');
    hold on;
    plot(t_builtin, y_builtin, 'r--', 'DisplayName', 'MATLAB ode23');
    xlabel('t');
    ylabel('y(t)');
    legend;
    title(['Custom RK23 vs MATLAB ode23, \delta= ', num2str(delta) ]);
    grid on;

end



% f is coninuous everywhere + df/dy=2y-3y^2 id coninuous, then f is locally
% Lipschitz in y (SC for LC) => on any compact interval around y0=delta, the derivative is bounded =>
% Lipschitz
%=> the ODE admits a unique local solution around t=0, but for this ODE,
%the existance is guaranteed only locally for small delta since for large
%delta the solution may blow up in finite time. For small delta we get
%global existance up to t=2/delta.


%% point 2
% question about b_hat

clear all
clc
close all
%parameters and initialization 
%delta=0.02;
delta_val=[0.02, 0.01, 0.005, 0.001];
t0=0;
f = @(t, y) y^2 - y^3;


for i=1:length(delta_val)
    delta=delta_val(i);
    tf=2/delta; 
    y0=delta; 
    
    aeps = 1e-6;   % Absolute error tolerance
    reps = 1e-4;   % Relative error tolerance
    
    h=0.01; %initial step size (guess)
    t=t0; 
    y=y0; 
    
    T = t;
    Y = y;
    
    % Efficiency counters
    func_evals = 0;
    accepted_steps = 0;
    rejected_steps = 0;

    while t < tf
        if t + h > tf
            h = tf - t;
        end
        
        % RK23 step and error estimate
        [y_high,y_low,err] = rk23_step(f, t, y, h);

        tol = reps * norm(y_high) + aeps;
    
        func_evals = func_evals + 3; % 3 evaluations per step (k1, k2, k3)

        if err <= tol
            t = t + h;
            y = y_high;
            T(end+1) = t;
            Y(end+1) = y;
            accepted_steps = accepted_steps + 1;
        else 
            rejected_steps = rejected_steps + 1;
        end
        
        % Update step size
        h = h * min(5, max(0.1, 0.9 * (tol / norm(err))^(1/3)));
    end
    
    [t_builtin, y_builtin] = ode23(f, [0 tf], delta);

    % Plot
    figure;
    plot(T, Y, 'b.-', 'DisplayName', 'Custom RK23');
    hold on;
    plot(t_builtin, y_builtin, 'r--', 'DisplayName', 'MATLAB ode23');
    xlabel('t');
    ylabel('y(t)');
    legend;
    title(['Comparison: Custom RK23 vs MATLAB ode23, \delta= ' , num2str(delta)]);
    grid on;

    % Efficiency output
    fprintf('\n--- Efficiency Report for δ = %.4f ---\n', delta);
    fprintf('Total accepted steps    : %d\n', accepted_steps);
    fprintf('Total rejected steps    : %d\n', rejected_steps);
    fprintf('Total function evaluations: %d\n', func_evals);
    fprintf('Average f evals per unit t: %.2f\n', func_evals / tf);
end
