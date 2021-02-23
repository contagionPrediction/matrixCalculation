M0=csvread('2020北京市汇总_cleaned.csv',1,1);

L=length(M0);
q=10;
c=1e-7;
N=1e6;
x0=zeros(5,1);

r=zeros(1,5,q);
M=zeros(5,5,q);
R=zeros(5,1,q);
x=zeros(5,1,q);
k=zeros(1,1,q);

for i=1:q
    r(:,:,i)=randi(L,[1,5]);
    M(:,:,i)=M0(r(:,:,i),1:5);
    R(:,:,i)=M0(r(:,:,i),8);
    [x(:,:,i),k(:,:,i)]=conjugateGradient(M(:,:,i),R(:,:,i),x0,c,N);
end
