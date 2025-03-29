function Unew = smooth(U, omega, m, F)
%     h = 1/(m+1);
%     U = reshape(U, m, m);  % reshape into 2D grid
%     F = reshape(F, m, m);
%     Unew = U;
% 
%     for i = 2:m-1
%         for j = 2:m-1
%             Unew(i,j) = (1 - omega) * U(i,j) + omega * 0.25 * ...
%                 (U(i+1,j) + U(i-1,j) + U(i,j+1) + U(i,j-1) - h^2 * F(i,j));
%         end
%     end
% 
%     Unew = reshape(Unew, m*m, 1);  % flatten back to vector
% end

h=1/(m+1); %define h

if max(size(F))>m % if F is a vector convert to matrix
    F=reshape(F,m,m);
end
if max(size(U))>m
    U=reshape(U,m,m);
end
%add zero padding to U
U0=zeros(m+2,m+2);
U0((1:m)+1,(1:m)+1)=U;

%sum of the neighboring U values
Unew=U0((1:m)+1-1,(1:m)+1)+U0((1:m)+1+1,(1:m)+1)+U0((1:m)+1,(1:m)+1-1)+U0((1:m)+1,(1:m)+1+1);


% Unew = Gw*U+c:
Unew=(1-omega)*U+omega*(Unew - h^2*F)/4; %complete calculation of next iteration

Unew=reshape(Unew,[],1);
