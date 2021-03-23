%load data   %load your data or comment and run the script if data are already loaded (n.b. data dimensions [npoints nvariables])
mtot=1:20;  %possible model orders
lam=0.05;   %regularization parameter
k=0;
[n nvar]=size(data);
type='polynomial';
par=2;
%%% In the following lines we introduce the lenght of the time window and the step, so the
%%% model order is averaged across all the sliding windows.
%%% If you are not interested in moving window, just put dt=n;
dt=60;  %window length
dts=1;  %step
ts=1:dts:n-dt;
nwind=length(ts);
empytot=zeros(length(mtot),length(ts),4);looytot=empytot;
for iw=1:nwind
    disp([iw nwind]);
    t=ts(iw):ts(iw)+dt-1;
    XX=data(t,:);
    for i=1:nvar
        for j=1:nvar
            k=k+1;
            for im=1:length(mtot);
                m=mtot(im);
                x=XX(:,i);
                y=XX(:,j);
                [X,Y,xt,yt,L]=prepro(m,x,y);
                T=X;t=yt;
                [empytot(im,iw,k),looytot(im,iw,k)]=compute_errors(lam,type,par,m,dt-m,T,t);
            end
        end
    end
end
empy=squeeze(mean(empytot,2));looy=squeeze(mean(looytot,2));
figure;plot(mtot,mean(empy,2),'ob');hold on;plot(mtot,mean(looy,2),'or');
title('leave one out cross validation');legend('empirical risk','generalization error');
xlabel('model order');ylabel('mean error')