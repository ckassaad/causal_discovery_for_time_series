function [xx,yy]=two_ccm_vareps(n,m,a,eps1,eps2,sig,modulation1,modulation2)
n1=n+1000;
modulation1=[ones(1000,1);modulation1];
modulation2=[ones(1000,1);modulation2];
x=zeros(n1,1);
y=zeros(n1,1);
x(1:m)=rand(m,1);
y(1:m)=rand(m,1);
epsmod1=eps1*modulation1;
epsmod2=eps2*modulation2;
test=1;
while test>0
    for i=m+1:n1
        bx=abs(mean(x(i-m:i-1)));
        by=abs(mean(y(i-m:i-1)));
        x(i)=(1-epsmod1(i))*(1-a*bx^2)+epsmod1(i)*(1-a*by^2)+sig*randn;
        y(i)=(1-epsmod2(i))*(1-a*by^2)+epsmod2(i)*(1-a*bx^2)+sig*randn;
    end
    X=[x y];
    test=isinf(max(max(abs(X)))^2);
end
xx=x(n1-n+1:n1);
yy=y(n1-n+1:n1);