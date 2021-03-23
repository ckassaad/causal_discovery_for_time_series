clear;clc
path('C:\Program Files\MATLAB\R2011b\toolbox\Kernel Causality Last', path);
tau=1;              %separation of the m points used for prediction (generally we leave it =1)
nt=1;               %starting point, normally always 1
delta=0;            %gap between training and test
type='p';           %type of kernel ('p' for polynomial, 'g' for gaussian)
par=2;              %parameter of the model (order for polynomial kernel, or width for the gaussian kernel
N=2100;             %number of points of simulated time series
msim=1;             %order of the model in the simulations
m=1;
%%%%%%%%%% influences %%%%%%%%%
eps1=0;     %influence 2->1
eps2=0.1;   %influence 1->2
modulation1=ones(N,1);
modulation2=zeros(N,1);modulation2(ceil(N/3):2*ceil(N/3))=1;  %step modulation
%modulation2=sin((1:N)/600).^10;      %sine-like modulation
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
windowlength=20;
step=1;
ts=windowlength+1:step:N;   %time vector in points
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
sig=0.03;   %noise
ntrials=100;
Gtot=zeros(ntrials,2,2,length(ts));
corrtot=zeros(ntrials,length(ts));
for itrials=1:ntrials
    disp(itrials)
    [xx,yy]=two_ccm_vareps(N,msim,1.8,eps1,eps2,sig,modulation1,modulation2);
    data=[xx yy];            %data matrix must have the dimensions [n_points n_variables]
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    [np nv]=size(data); %number of observations (time points) and number of variables
    Gtot(itrials,:,:,:)=caus_movingwindow(data,m,type,par,windowlength,step);            %compute causality matrix
    corrtot(itrials,:)=corr_movingwindow(data,windowlength,step);
end
G=squeeze(mean(Gtot,1));
C=squeeze(mean(corrtot,1));
figure(3);
hold on
subplot(3,1,1);hold on;plot(ts,squeeze(G(1,2,:)),'go');ylabel('causality 1->2');%, \epsilon = ',num2str(eps),' \sigma = ',num2str(sig)));
subplot(3,1,2);hold on;plot(ts,squeeze(G(2,1,:)),'go');ylabel('causality 2->1');
subplot(3,1,3);hold on;plot(ts,C,'go');ylabel('correlation');
suptitle(['\epsilon = ' num2str(eps2) ' ' ' \sigma = ' num2str(sig)]);
