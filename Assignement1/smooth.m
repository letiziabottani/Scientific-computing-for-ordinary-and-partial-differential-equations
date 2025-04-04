function Unew=smooth(U,omega,m,F)
h=1/(m+1); %define h

if max(size(F))>m % if F is a vector convert to matrix
    F=reshape(F,m,m);
end

%add zero padding to U
U0=zeros(m+2,m+2);
U0((1:m)+1,(1:m)+1)=U;

%sum of the neighboring U values
Unew=U0((1:m)+1-1,(1:m)+1)+U0((1:m)+1+1,(1:m)+1)+U0((1:m)+1,(1:m)+1-1)+U0((1:m)+1,(1:m)+1+1);

%Unew=omega*(Unew - h^2*F)/4;?

% Unew = Gw*U+c:
Unew=(1-omega)*U+(omega*Unew - h^2*F)/4; %complete calculation of next itteration



