function [x,k]=SORIteration(A,b,w,x0,c,N)

% k is the actual number of iterations
% N is the upper bound of the number of iterations
% x0 is the initial vector
% c is precision
% w is the relaxation factor (w>0)

D=diag(diag(A));
L=diag(diag(A))-tril(A);
R=eye(length(A))-w*((D-w*L)^(-1))*A;
f=w*((D-w*L)^(-1))*b;
x=R*x0+f;
k=1;

while norm(x-x0,2)>=c
    x0=x;
    x=R*x0+f;
    k=k+1;
    if k>=N
        break;
    end
end