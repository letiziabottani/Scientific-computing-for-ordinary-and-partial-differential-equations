function [y_high,y_low,err]=rk23_step(f,t,y,h)
    % RK23 coefficients (Bogacki–Shampine method), see table 5.2 in the
    % book
    
    % a coefficients for intermediate stages
    a21 = 1/2;
    a31 = -1;
    a32 = 2;
    
    % b coefficients (3rd order solution)
    b1 = 1/6;
    b2 = 2/3;
    b3 = 1/6;
    
    % b_hat coefficients (2nd order solution)
    b1_hat = 1/4;
    b2_hat = 1/2;
    b3_hat = 1/4;

    
    % d_i = b_i - b_hat_i (used to estimate the local error)
    d1=b1-b1_hat;
    d2=b2-b2_hat;
    d3=b3-b3_hat;

    c1=0; 
    c2=1/2; 
    c3=1; 
    
    k1 = f(t,y);
    k2 = f(t + c2*h, y + a21*h*k1);
    k3 = f(t + c3*h, y + a31*h*k1 + a32*h*k2);
    
    % 3rd-order solution
    y_high = y + h * (b1*k1 + b2*k2 + b3*k3);

    % 2nd-order solution
    y_low = y + h * (b1_hat*k1 + b2_hat*k2 + b3_hat*k3);
    err=norm(y_high-y_low);
end 


