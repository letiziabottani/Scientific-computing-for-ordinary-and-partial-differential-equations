function F = form_rhs_9(m,f,u)
    h=1/(m+1);
    interval=linspace(h,1-h,m);
    
    [X,Y]=meshgrid(interval,interval);

    f_val=f(X,Y);
    f_val(:,1)=f_val(:,1)-(1/6*(m+1)^2)*(u(interval'-h,0)+4*u(interval',0)+u(interval'+h,0));

    f_val(:,m)=f_val(:,m)-(1/6*(m+1)^2)*(u(interval'-h,1)+4*u(interval',1)+u(interval'+h,1));

    f_val(1,:)=f_val(1,:)-(1/6*(m+1)^2)*(u(0,interval-h)+4*u(0,interval)+u(0,interval+h));

    f_val(m,:)=f_val(m,:)-(1/6*(m+1)^2)*(u(1,interval-h)+4*u(1,interval)+u(1,interval+h));

    %correction for adding
    f_val([1,m],[1,m])=f_val([1,m],[1,m])+1/6*((m+1)^2)*u([0,0;1,1],[0,1;0,1]);

    F=reshape(f_val',[],1);
end
