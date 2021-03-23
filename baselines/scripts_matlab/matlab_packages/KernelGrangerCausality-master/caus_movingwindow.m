function cw=caus_movingwindow(data,m,type,par,dt,dts,delta,tau,nt)
if nargin<8
    tau=1;nt=1;
end
if nargin<9
    delta=0;
end
[np, nv]=size(data);
ts=1:dts:np-dt;
nk=length(ts);
cw=zeros(nv,nv,nk);
for k=1:nk
    t=ts(k):ts(k)+dt-1;
    XX=data(t,:);
    cw(:,:,k)=causality(XX,type,par,m);
end