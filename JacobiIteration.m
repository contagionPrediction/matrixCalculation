function [x,k]=JacobiIteration(A,b,x0,c,N)

% k is the actual number of iterations
% N is the upper bound of the number of iterations
% x0 is the initial vector
% c is the precision

D=diag(diag(A));
J=eye(length(A))-D^(-1)*A;
f=D^(-1)*b;
x=J*x0+f;
k=1;

while norm(x-x0,2)>=c
    x0=x;
    x=J*x0+f;
    k=k+1;
    if k>=N
        break;
    end
end