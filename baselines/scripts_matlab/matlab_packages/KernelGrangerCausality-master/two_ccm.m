function [xx,yy]=two_ccm(n,m,a,eps1,eps2,sig)
n1=n+1000;
x=zeros(n1,1);
y=zeros(n1,1);
x(1:m)=rand(m,1);
y(1:m)=rand(m,1);
test=1;
while test>0
    for i=m+1:n1
        bx=abs(mean(x(i-m:i-1)));
        by=abs(mean(y(i-m:i-1)));
        x(i)=(1-eps1)*(1-a*bx^2)+eps1*(1-a*by^2)+sig*randn;
        y(i)=(1-eps2)*(1-a*by^2)+eps2*(1-a*bx^2)+sig*randn;
    end
    X=[x y];
    test=isinf(max(max(abs(X)))^2);
end
xx=x(n1-n+1:n1);
yy=y(n1-n+1:n1);