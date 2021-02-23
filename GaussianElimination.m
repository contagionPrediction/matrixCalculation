function x=GaussianElimination(A,b)

% gaussian elimination with maximal column pivoting

% elimination process
n=length(A);
a=[A,b];
m=zeros(n);
for k=1:n-1
    maxa=max(abs(a(k:n,k)));
    if maxa==0
        return
    end
    for i=k:n
        if abs(a(i,k))==maxa
            y=a(i,k:n+1);
            a(i,k:n+1)=a(k,k:n+1);
            a(k,k:n+1)=y;
            break
        end
    end
    for i=k+1:n
        m(i,k)=a(i,k)/a(k,k);
        a(i,k+1:n+1)=a(i,k+1:n+1)-m(i,k).*a(k,k+1:n+1);
    end
end

% back substitution process
if a(n,n)==0
    return
end
x(n)=a(n,n+1)/a(n,n);
for i=n-1:-1:1
    x(i)=(a(i,n+1)-sum(a(i,i+1:n).*x(i+1:n)))/a(i,i);
end
x=x.';
