function cw=corr_movingwindow(data,dt,dts)
[np, nv]=size(data);
ts=1:dts:np-dt;
nk=length(ts);
cw=zeros(nk,1);
for k=1:nk
    t=ts(k):ts(k)+dt-1;
    XX=data(t,:);
    cc=corrcoef(XX(:,1),XX(:,2));
    cw(k)=cc(1,2);
end