function [cb cbt]=causality_biv_trials(y,xt,type,par,m,ntrials,th)
% Input:    y       : matrix (n x 1) of driver data;
%           xt       : matrix (n x 1) of target data
%                       nvar = number of conditioning time series
%           n    = number of samples; ntrials = number of trials
%           type    : kernel function 'p' polynomial 'g' gaussian ;
%           par     : parameter of the kernel function;
%           m       : order of the model
% Output:   cb      :  y->xt
%           ifail   :  0 no error
%                      1 cholesky algoritm fail
%                      2 complex eigensystem
%                      3 error in polypower
%Last version now
n=length(y);
% for i=1:nvar;
%     Y(:,i)=(Y(:,i)-mean(Y(:,i)))/std(Y(:,i));
% end
% y=(y-mean(y))/std(y);
% xt=(xt-mean(xt))/std(xt);


Y=xt;
%%% N=trial length, M=number of trials
N=n/ntrials;
order=m;
M=ntrials;
x=[];
y_past=[];
Y_past=[];
past_ind = repmat([1:order],N-order,1) + repmat([0:N-order-1]',1,order);
for i=1:M
    %now
    xc = xt((i-1)*N+order+1:i*N,:);
    %past
    Y_past_c = reshape(Y((i-1)*N+past_ind,:),N-order,order*size(Y,2));
    y_past_c = reshape(y((i-1)*N+past_ind,:),N-order,order);
    %%%%%
    %%% concatenate
    x=[x;xc];
    Y_past=[Y_past;Y_past_c];
    y_past=[y_past;y_past_c];
end

%%%%%%% x is the target

%%%%%%%%%% Y_past is the past of the system without the driver

XY_past=[Y_past y_past]; %%% past including the driver 

%%%%%%%%%%ora normalizzazione
x=(x-mean(x))/std(x);
for j=1:size(Y_past,2)
    Y_past(:,j)=(Y_past(:,j)-mean(Y_past(:,j)))/std(Y_past(:,j));
end

for j=1:size(XY_past,2)
    XY_past(:,j)=(XY_past(:,j)-mean(XY_past(:,j)))/std(XY_past(:,j));
end


f=1.e-6;
Xr=XY_past';
cb=0;
rr=0;
pp=0;
[VV D ifail]=filtro(Xr,type,par,f,true);
if ifail>0
    return
end
VT=VV*D.^0.5;
polycall=true;
Xr=Y_past';
[V, D, ifail]=filtro(Xr,type,par,f,polycall);
if ifail>0
    return
end
polycall=false;
[VN, ifail]=vnorma(VT,V,VV);
if ifail>0
    return
end
xv=x-V*V'*x;
[rrt, ppt]=corr(xv,VN);

rn=rrt.^2;
cbt=sum(rn);
thb=th/length(rrt);
indpr=find(ppt>thb);
rn(indpr)=0;
cb=sum(rn);
cb=-log(1-cb);
