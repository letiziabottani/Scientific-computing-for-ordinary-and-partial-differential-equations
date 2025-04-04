function F = form_rhs_9(m,f,u,laplacian_f)
    h=1/(m+1);
    interval=linspace(h,1-h,m);
    
    [X,Y]=meshgrid(interval,interval);

    F=f(X,Y);
    F(:,1)=F(:,1)-(1/6*(m+1)^2)*(u(interval'-h,0)+4*u(interval',0)+u(interval'+h,0));
    F(:,m)=F(:,m)-(1/6*(m+1)^2)*(u(interval'-h,1)+4*u(interval',1)+u(interval'+h,1));
    F(1,:)=F(1,:)-(1/6*(m+1)^2)*(u(0,interval-h)+4*u(0,interval)+u(0,interval+h));
    F(m,:)=F(m,:)-(1/6*(m+1)^2)*(u(1,interval-h)+4*u(1,interval)+u(1,interval+h));
    
    %correction of corners for dubble adding.
    F([1,m],[1,m])=F([1,m],[1,m])+1/6*((m+1)^2)*u([0,0;1,1],[0,1;0,1]);
    
    laplac_eval=laplacian_f(X,Y);
    F= F+1/12*h^2*laplac_eval;

    F=reshape(F',[],1);
end
