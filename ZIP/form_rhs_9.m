function F = form_rhs_9(m,f,u,laplacian_f)
    h=1/(m+1);
    interval=linspace(h,1-h,m);
    
    [X,Y]=meshgrid(interval,interval);

    F=f(X,Y);
   for i = 1:(m)
        F(i,1)=F(i,1)-(1/(6*h^2))*(4*u(0,i/(m+1))+u(0,(i-1)/(m+1))+u(0,(i+1)/(m+1)));
        F(i,m)=F(i,m)-(1/(6*h^2))*(4*u(1,i/(m+1))+u(1,(i-1)/(m+1))+u(1,(i+1)/(m+1)));
    end
    
    for i = 1:m
        F(1,i)=F(1,i)-(1/(6*h^2))*(4*u(i/(m+1),0)+u((i-1)/(m+1),0)+u((i+1)/(m+1),0));
        F(m,i)=F(m,i)-(1/(6*h^2))*(4*u(i/(m+1),1)+u((i-1)/(m+1),1)+u((i+1)/(m+1),1));
    end
    
    %correction for adding
    F(1,1)=F(1,1)+1/6*((m+1)^2)*u(0,0);
    F(1,m)=F(1,m)+1/6*((m+1)^2)*u(0,1);
    F(m,m)=F(m,m)+1/6*((m+1)^2)*u(1,1);
    F(m,1)=F(m,1)+1/6*((m+1)^2)*u(1,0);
    
    laplac_eval=laplacian_f(X,Y);
    F= F+1/12*h^2*laplac_eval;

    F=reshape(F',[],1);
end
