clear;clc
%%% this section simulates a dataset of 5 variables
%%% couplings: 1->2, 1->3, 1->4, 4->5, 5->4
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
path('C:\Users\admin\Desktop\KGC_export\Kernel Causality Last',path);
nonlinear=true;    %set to true for nonlinear model
m=6;                %order of the model
tau=1;              %separation of the m points used for prediction (generally we leave it =1)
nt=1;               %starting point, normally always 1
delta=0;            %gap between training and test
type='p';           %type of kernel ('p' for polynomial, 'g' for gaussian)
par=2;              %parameter of the model (order for polynomial kernel, or width for the gaussian kernel
N=2000;             %number of points of simulated time series
settleTime=1000;    %settling time
%create a null distribution with randomized phases. this is not always
%necessary, since there is already the test on the coefficients, but it's good to have
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
N = N + settleTime;
X = randn(5,N);
%X = cca_normrnd(0,1,5,N);
r = sqrt(2);
if nonlinear
    %%%%%% nonlinear case, so you can see that it works better with p=2
    for i=4:N,
        X(1,i) = X(1,i) + 0.95.*r.*X(1,i-1) - 0.9025.*X(1,i-2);
        X(2,i) = X(2,i) + 0.5.*X(1,i-2)^2;
        X(3,i) = X(3,i) - 0.4.*X(1,i-3);
        X(4,i) = X(4,i) - 0.5.*X(1,i-2)^2 + 0.5.*r.*X(4,i-1) + 0.25.*r.*X(5,i-1);
        X(5,i) = X(5,i) - 0.5.*r.*X(4,i-1) + 0.5.*r.*X(5,i-1);
    end
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
else
    %%%%%% linear influences only
    for i=4:N,
        X(1,i) = X(1,i) + 0.95.*r.*X(1,i-1) - 0.9025.*X(1,i-2);
        X(2,i) = X(2,i) + 0.5.*X(1,i-2);
        X(3,i) = X(3,i) - 0.4.*X(1,i-3);
        X(4,i) = X(4,i) - 0.5.*X(1,i-2) + 0.25.*r.*X(4,i-1) + 0.25.*r.*X(5,i-1);
        X(5,i) = X(5,i) - 0.25.*r.*X(4,i-1) + 0.25.*r.*X(5,i-1);
    end
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
end
X = X(:,settleTime+1:end);
data=X';            %data matrix must have the dimensions [n_points n_variables]
%NTimepoints=[NTimepoints;data];
%end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
[np nv]=size(data); %number of observations (time points) and number of variables
n=np-1;             %endpoint
[GC_matrix]=causality(data,type,par,m);            %compute causality matrix
figure;imagesc(GC_matrix);colorbar;
xlabel('to');ylabel('from')
set(gca,'XTick',[1:5])
set(gca,'YTick',[1:5])