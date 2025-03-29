clc
clear all
close all
% exact solution and RHS
u=@(x,y) exp(pi*x).*sin(pi*y)+0.5*(x.*y).^2;
f=@(x,y) x.^2+y.^2;
m=2^6-1;
U =zeros(m*m,1);
F =form_rhs(m,f,u); %% TODO: Form the right-hand side
epsilon = 1.0E-6;
omega=find_best_omega(m);
%omega=2/3;

%%
figure()
for i=1:100
    R =F+Amult(U,m);
    fprintf('*** Outer iteration: %3d, rel. resid.: %e\n', ...
        i, norm(R,2)/norm(F,2));
    % fprintf(' R:%e\n',norm(R,2))
    % fprintf(' F:%e\n',norm(F,2) )
    if(norm(R,2)/norm(F,2) < epsilon)
        break;
    end
    U=Vcycle(U,omega,3,m,F);
    
    plotU(m,U);
    pause(.5);
end

%%
function Unew=Vcycle(U,omega,nsmooth,m,F)
% Approximately solve: A*U = F
h=1/(m+1);
l2m=log2(m+1);
assert(l2m==round(l2m));
assert(length(U)==m*m);
if(m==1)
    % if we are at the coarsest level
    % TODO: solve the only remaining equation directly!
     A = poisson5(m); 
    Unew = A \ F;
    
    return;
else
    % 1. TODO: pre-smooth the error
    %    perform <nsmooth> Jacobi iterations
    for i=1:nsmooth
        U=smooth(U,omega,m,F);
    end 
    % 2. TODO: calculate the residual
    R=F+Amult(U,m);
    % 3. TODO: coarsen the residual
    Rcoarse=coarsen(R,m);
    % 4. recurse to Vcycle on a coarser grid
    mc=(m-1)/2;
    Ecoarse=Vcycle(zeros(mc*mc,1),omega,nsmooth,mc,-Rcoarse);
    % 5. TODO: interpolate the error
    E=interpolate(Ecoarse,m);
    % 6. TODO: update the solution given the interpolated error
    %disp(['||E|| = ', num2str(norm(E)), '   ||U|| = ', num2str(norm(U))]);
    %E = smooth(E, omega, m, zeros(size(E)));  
    U = U -  E;
    %disp(['||U updated|| = ', num2str(norm(U))]);

    % 7. TODO: post-smooth the error
    %    perform <nsmooth> Jacobi iterations
    for i=1:nsmooth
        U=smooth(U,omega,m,F);
    end 
    Unew=U;
end
end


function plotU(m,U)
h=1/(m+1);
x=linspace(h,1-h,m);
y=linspace(h,1-h,m);
[X,Y]=meshgrid(x,y);
surf(X, Y, reshape(U,[m,m])');
shading interp;
title('Computed solution');
xlabel('x');
ylabel('y');
zlabel('U');
end


function optimal_omega = find_best_omega(m_values)
    lambda = @(p,q,h,w) 1+w.*((cos(p.*pi.*h)-1)+(cos(q.*pi.*h)-1))/2;
    omega_values = linspace(0, 2, 1000);% Omega values from 0 to 2
    num_m = length(m_values);
    optimal_omegas = zeros(size(m_values));
    
    figure; hold on; % Create a single figure for all m values
    
    for idx = 1:num_m
        m = m_values(idx);
        
        h = 1 / (m + 1); % Grid spacing
        
        interval =  ceil(m/2):m;
        [p,q]=meshgrid(interval,interval);
        
        
        
        l=zeros(length(omega_values),1);
        %calculate max gamma for each omega:
        for i=1:length(omega_values)
            l(i)=max(abs(lambda(p,q,h,omega_values(i))),[],"all");
        end
        
        
        % Find the optimal omega that minimizes max |gamma_pq^omega|
        [max_gamma, min_idx] = min(l);

        optimal_omega = omega_values(min_idx);
        optimal_omegas(idx) = optimal_omega;
        
        % Plot for this m value
        plot(omega_values, l, 'LineWidth', 2, 'DisplayName', ['m = ', num2str(m)]);
        plot(optimal_omega, max_gamma, 'ro', 'MarkerSize', 10, 'MarkerFaceColor', 'r');
    end
    
    % Formatting the plot
    xlabel('\omega');
    ylabel('max |\gamma_{p,q}^\omega|');
    title('Effect of \omega on max |\gamma_{p,q}^\omega| for different m');
    legend show;
    grid on;
    hold off;

    % Print optimal omegas
    for i = 1:num_m
        fprintf('Optimal omega for m = %d: %.4f\n', m_values(i), optimal_omegas(i));
    end
end