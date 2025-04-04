function AU=Amult(U,m)
h=1/(m+1); %define h

if max(size(U))>m
    U=reshape(U,m,m);
end

%add zero padding to U
U0=zeros(m+2,m+2);
U0((1:m)+1,(1:m)+1)=U;

%sum of the neighboring U values
AU=-(U0((1:m)+1-1,(1:m)+1)+U0((1:m)+1+1,(1:m)+1)+U0((1:m)+1,(1:m)+1-1)+U0((1:m)+1,(1:m)+1+1)-4*U)/h^2;

AU=reshape(AU,[],1);
end 