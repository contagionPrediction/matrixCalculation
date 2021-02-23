function [x,k]=conjugateGradient(A,b,x0,c,N)

% k is the actual number of iterations
% N is the upper bound of the number of iterations
% x0 is the initial vector
% c is the precision

n=length(A);
r=zeros(n,N);
p=zeros(n,N);
x=x0;
r(:,1)=b-A*x;
alpha=zeros(N,1);
beta=zeros(N,1);
k=1;

while norm(r(:,k),2)>=c
    k=k+1;
    if k==2
        p(:,1)=r(:,1);
    else
        beta(k-2)=((r(:,k-1))'*r(:,k-1))/((r(:,k-2))'*r(:,k-2));
        p(:,k-1)=r(:,k-1)+beta(k-2)*p(:,k-2);
    end
    alpha(k-1)=((r(:,k-1))'*r(:,k-1))/((p(:,k-1))'*A*p(:,k-1));
    x=x+alpha(k-1)*p(:,k-1);
    r(:,k)=r(:,k-1)-alpha(k-1)*A*p(:,k-1);
    if k>=N
        break;
    end
end

k=k-1;
