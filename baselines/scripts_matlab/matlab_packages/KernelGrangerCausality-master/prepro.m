function [X,Y,xt,yt,L]=prepro(m,x,y)
n=size(x,1);
x=x-mean(x);
x=x./std(x);
y=y-mean(y);
y=y./std(y);
L=n-m;
X=zeros(m,L);
Y=zeros(m,L);
xt=zeros(L,1);
yt=zeros(L,1);
for i=1:L
    X(1:m,i)=x(i:i+m-1);
    Y(1:m,i)=y(i:i+m-1);
    xt(i)=x(i+m);
    yt(i)=y(i+m);
end


