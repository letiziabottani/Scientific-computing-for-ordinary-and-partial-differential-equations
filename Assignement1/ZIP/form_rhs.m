function F=form_rhs(m,f,u)
    h=1/(m+1);
    interval=linspace(h,1-h,m);
    
    [X,Y]=meshgrid(interval,interval);

    F=f(X,Y);
    
    %incorporating the boundary conditions
    F(:,1)=F(:,1)-((m+1)^2)*u(interval',0);
    F(:,m)=F(:,m)-((m+1)^2)*u(interval',1);
    F(1,:)=F(1,:)-((m+1)^2)*u(0,interval);
    F(m,:)=F(m,:)-((m+1)^2)*u(1,interval);
    
    %reshape 
    F=reshape(F',[],1);
end