function [X x]=init_causality(XX,m)
if nargin<2
    m=1;
end
[nn nvar]=size(XX);
for i=1:nvar;
    XX(:,i)=(XX(:,i)-mean(XX(:,i)))/std(XX(:,i));
end
n=nn-m;
X=zeros(nvar,m,n);
x=zeros(n,nvar);
for i=1:n
    for k=1:nvar
        for j=1:m
            X(k,j,i)=XX(m-j+i,k);
        end
        x(i,k)=XX(i+m,k);
    end
end
% media nulla
Xm=reshape(repmat(mean(X,3),1,n),nvar,m,n);
X=X-Xm;
xm=repmat(mean(x,1),n,1);
x=x-xm;
x=x./repmat(sqrt(diag(x'*x)'),n,1);
