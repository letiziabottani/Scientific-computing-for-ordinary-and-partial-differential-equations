function F=form_rhs(m,f,u)
    h=1/(m+1);
    x = linspace(h,1-h,m);
    y = linspace(h,1-h,m)
    [X,Y]=meshgrid(x,y);
    F = f(X,Y);
    for i = 1:(m)
        F(i,1)=F(i,1)-(1/h^2)*u(0,i/(m+1));
        F(i,m)=F(i,m)-(1/h^2)*u(1,i/(m+1));
    end

    for i = 1:m
        F(1,i)=F(1,i)-(1/h^2)*u(i/(m+1),0);
        F(m,i)=F(m,i)-(1/h^2)*u(i/(m+1),1);
    end

    F = reshape(F',[m*m,1]);
    

end