function [x,k]=GaussSeidelIteration(A,b,x0,c,N)

% k is the actual number of iterations
% N is the upper bound of the number of iterations
% x0 is the initial vector
% c is the precision

U=diag(diag(A))-triu(A);
G=tril(A)^(-1)*U;
f=tril(A)^(-1)*b;
x=G*x0+f;
k=1;

while norm(x-x0,2)>=c
    x0=x;
    x=G*x0+f;
    k=k+1;
    if k>=N
        break;
    end
end
