% Exercise 4 3 RK

clc
clear all
close all

%% General initialization
e=0.01/pi;

u_initial=@(x) -sin(pi*x);

% left and right boundaries
gl = @(t) 0;
gr = @(t) 0;

max_time = 1.6037/pi;

%% run numerical scheme

aeps = 1e-8;   % Absolute error tolerance
reps = 1e-6;   % Relative error tolerance

%aeps = 1e-6;   % Absolute error tolerance
%reps = 1e-4;   % Relative error tolerance


%precition

N = 2^10;

h = 2/N;

x = -1+h*(0:N);

k_scale = 1/10; 
k = (4/N^2)/(2*e)*k_scale; %initial step size (guess)

t=0; 
U=u_initial(x); %initial condition

%save values for plotting
t_save = zeros(1,N*100);
t_save(1) = t;
U_save = zeros(N*100,N+1);
U_save(1,:) = U;

% Efficiency counters
func_evals = 0;
accepted_steps = 0;
rejected_steps = 0;

while t < max_time 
    if t + k > max_time
        k = max_time - t;
    end
    
    % RK23 step and error estimate
    
    [y_high,y_low,err] = Burger_RK_iteration(U,t,gl,gr,k,e);

    tol = reps * norm(y_high) + aeps;

    func_evals = func_evals + 3; % 3 evaluations per step (Y1, Y2, Y3)

    if err <= tol 
        
        t = t + k;
        U = y_high;

        accepted_steps = accepted_steps + 1;

        if accepted_steps <= length(t_save)  
            t_save(accepted_steps) = t;
            U_save(accepted_steps,:) = U;
        else
            t_save(end+1) = t;
            U_save(end+1,:) = U;
        end
    else 
        rejected_steps = rejected_steps + 1;
    end
    
    % Update step size
    k = k * min(5, max(0.1, 0.9 * (tol / norm(err))^(1/3)));

    if abs(k) < 1e-8
        fprintf("Using a stepsize less than 1E-8*(tfinal-t0) at t = %g\n",t);
    end

end
%clear tol err y_low y_low 

if accepted_steps < N*100 
    t_save((accepted_steps+1):end) = [];
    U_save((accepted_steps+1):end,:) = [];
end

%% DX at time 1.6037/π

i_0=N/2+1;

dx_uni_0 = (U_save(:,i_0+1)-U_save(:,i_0-1))/(2*h);

%plot(t_save,dx_uni_0)

[dx_uni_0(end),dx_uni_0(end)+152.00516]

%% plot solution
surf(x,t_save,U_save, EdgeColor="none")
ylabel("t")
xlabel("x")
title('Numerical Solution U(x,t)'); xlabel('x'); ylabel('t');


%% how does the t_steps vary
t_steps = t_save(2:end)-t_save(1:end-1);

plot(t_steps)

%% plot dx of u

dx = (U((2:N)+1)-U((2:N)-1))/(2*h);
plot(x(2:N),dx)
title('Numerical Solution dx(U(x,t_end))'); xlabel('x'); ylabel('dx');
min(dx)

%% convergence test (varying N)

N_s = 2.^(3:12);

h_vals = 2./N_s;

aeps = 1e-8;   % Absolute error tolerance
reps = 1e-6;   % Relative error tolerance

k_scale = 1/10; 

dx_log = zeros(1,length(N_s));

for i = 1:length(N_s)

    N = N_s(i);

    h = h_vals(i);

    x = -1+h*(0:N);
    
   
    k = (4/N^2)/(2*e)*k_scale; %initial step size (guess)

    t=0; 
    U_non=u_initial(x); 


    % Efficiency counters
    func_evals = 0;
    accepted_steps = 0;
    rejected_steps = 0;
    
    while t < max_time 
        if t + k > max_time
            k = max_time - t;
        end
        
        % RK23 step and error estimate
        
        [y_high,y_low,err] = Burger_RK_iteration(U_non,t,gl,gr,k,e);
    
        tol = reps * norm(y_high) + aeps;
    
        func_evals = func_evals + 3; % 3 evaluations per step (Y1, Y2, Y3)
    
        if err <= tol
            
            t = t + k;
            U_non = y_high;
    
            accepted_steps = accepted_steps + 1;
        else 
            rejected_steps = rejected_steps + 1;
        end
        
        % Update step size
        k = k * min(5, max(0.1, 0.9 * (tol / norm(err))^(1/3)));
    
        if abs(k) < 1e-8
            fprintf("Using a stepsize less than 1E-8*(tfinal-t0) at t = %g\n",t);
        end
    
    end
    
    i_0=N/2+1;

    dx_uni_0 = (U_non(i_0+1)-U_non(i_0-1))/(2*h);
    
    dx_log(i) = dx_uni_0(end);


end

%% plot errors
loglog(h_vals,dx_log+152.00516,'-o',h_vals,h_vals.^2)

dx_log

%% non- uniform space grid

N_non=2^9;

x_initialise = -1+2*(0:N_non)/N_non;

alpha = 0.8;

x_non = (1-alpha).*x_initialise.^3 + alpha.*x_initialise; %non uniform x 

h_non = x_non(2:end)-x_non(1:end-1); % calculate new step size

%% Run 
k_non_scale=1/10;
k_non = min(h_non)^2/(2*e)*k_non_scale; %initial step size (guess)

aeps = 1e-8;   % Absolute error tolerance
reps = 1e-6;   % Relative error tolerance

t=0; 
U_non=u_initial(x_non); 

t_non_save = zeros(1,N_non*100);
t_non_save(1) = t;
U_non_save = zeros(N_non*100,N_non+1);
U_non_save(1,:) = U_non;

% Efficiency counters
func_evals = 0;
accepted_steps = 0;
rejected_steps = 0;

while t < max_time 
    if t + k_non > max_time
        k_non = max_time - t;
    end
    
    % RK23 step and error estimate
    
    [y_high,y_low,err] = Burger_RK_NU(U_non,t,gl,gr,k_non,e,h_non); %now includes h_non

    tol = reps * norm(y_high) + aeps;

    func_evals = func_evals + 3; % 3 evaluations per step (Y1, Y2, Y3)

    if err <= tol
        
        t = t + k_non;
        U_non = y_high;        

        accepted_steps = accepted_steps + 1;

        if accepted_steps <= length(t_non_save)  
            t_non_save(accepted_steps) = t;
            U_non_save(accepted_steps,:) = U_non;
        else
            t_non_save(end+1) = t;
            U_non_save(end+1,:) = U_non;
        end
    else 
        rejected_steps = rejected_steps + 1;
    end
    
    % Update step size
    k_non = k_non * min(5, max(0.1, 0.9 * (tol / norm(err))^(1/3)));

    if abs(k_non) < 1e-8
        fprintf("Using a stepsize less than 1E-8*(tfinal-t0) at t = %g\n",t);
    end

end
clear y_low y_low 

if accepted_steps < length(t_non_save) 
    t_non_save((accepted_steps+1):end) = [];
    U_non_save((accepted_steps+1):end,:) = [];
end


%% analysis
i_0=N_non/2+1;

dx_non_0 = (U_non_save(:,i_0+1)-U_non_save(:,i_0-1))./(2*h_non(i_0));

%plot(t_non_save,dx_non_0)

%value and error comparison
[dx_uni_0(end),dx_uni_0(end)+152.00516]

[dx_non_0(end),dx_non_0(end)+152.00516]

%% plot result

surf(x_non,t_non_save,U_non_save,EdgeColor="none")

%% compare varying alpha

alpha_s = (1:10)/10;

dx_alphas = zeros(1,length(alpha_s));

for i = 1:length(alpha_s)
    N_non=2^9;

    x_initialise = -1+2*(0:N_non)/N_non;
    
    alpha = alpha_s(i);
    
    x_non = (1-alpha).*x_initialise.^5 + alpha.*x_initialise;
    
    h_non = x_non(2:end)-x_non(1:end-1);

    k_non_scale=1/10;
    k_non = min(h_non)^2/(2*e)*k_non_scale; %initial step size (guess)
    
    aeps = 1e-8;   % Absolute error tolerance
    reps = 1e-6;   % Relative error tolerance
    
    t=0; 
    U_non=u_initial(x_non); 
    
    % Efficiency counters
    func_evals = 0;
    accepted_steps = 0;
    rejected_steps = 0;
    
    while t < max_time 
        if t + k_non > max_time
            k_non = max_time - t;
        end
        
        % RK23 step and error estimate
        
        [y_high,y_low,err] = Burger_RK_NU(U_non,t,gl,gr,k_non,e,h_non);
    
        tol = reps * norm(y_high) + aeps;
    
        func_evals = func_evals + 3; % 3 evaluations per step (Y1, Y2, Y3)
    
        if err <= tol
            
            t = t + k_non;
            U_non = y_high;        
    
            accepted_steps = accepted_steps + 1;
        else 
            rejected_steps = rejected_steps + 1;
        end
        
        % Update step size
        k_non = k_non * min(5, max(0.1, 0.9 * (tol / norm(err))^(1/3)));
    
        if abs(k_non) < 1e-8
            fprintf("Using a stepsize less than 1E-8*(tfinal-t0) at t = %g\n",t);
        end
    
    end
    clear y_low y_low 
    
    
    i_0=N_non/2+1;

    dx_alphas(i) = (U_non(i_0+1)-U_non(i_0-1))./(2*h_non(i_0));
end

%% plot errors
plot(alpha_s,abs(dx_alphas+152.00516),'-o',alpha_s,ones(1,length(alpha_s))*0.0131,'--',alpha_s,ones(1,length(alpha_s))*0.0028)
%ylim([0,0.04])
xlabel("\alpha");
ylabel("|error|");
legend('RK23 non-uniform N^{9}','RK23 uniform N^{9}','target accuracy', 'Location','northeast');
title('Error of RK23 non-uniform for varying \alpha');

