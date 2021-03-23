function [V D ifail]=filtro(X,type,p,f,polycall)
global term
if nargin<5
    polycall=true;
end
if nargin<4
    f=1.e-6;
end
[m n]=size(X);
if type=='g'
    fast=false;
    rmax=n;
else
    rmax=min(n,nchoosek(p+m,m));
    fast=nchoosek(p+m,m)<rmax+1;
end
if fast
    %polynomial kernel (p+m)!/(p!m!)-1< rmax'
    [L ifail]=Leval(X,m,p,polycall);
    if ifail>0
        v=0;
        return
    end
    [VN D]=eig(L'*L);
    V=L*VN;
    ifail=0;
else
    [L P ifa]=cholesky(X,type,p,rmax,f);
    if ifa
        ifail=1;
        V=0;
        return
    end
    ifail=0;
    [VN D]=eig(L'*L);
    V=L(P,:)*VN;
end
xnorm=repmat(sqrt(dot(V,V)),n,1);
V=V./xnorm;
[s ind]=sort(abs(diag(D)),'descend');
r=sum(s>f*s(1));
ind=ind(1:r);
V=V(:,ind);
D=D(ind,ind);