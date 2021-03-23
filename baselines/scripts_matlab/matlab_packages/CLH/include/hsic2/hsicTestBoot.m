%This function implements the HSIC independence test using a bootstrap approximation
%to the test threshold

%Inputs: 
%        X contains dx columns, m rows. Each row is an i.i.d sample
%        Y contains dy columns, m rows. Each row is an i.i.d sample
%        alpha is the level of the test
%        params.bootForce=1: force bootstrap. If 0: bootstrap if no
%        previous file found (recommended = 1)
%        params.shuff is number of shuffles to approximate null distribution
%        params.sigx is kernel size for x (set to median distance if -1)
%        params.sigy is kernel size for y (set to median distance if -1)

%Outputs: 
%        thresh: test threshold for level alpha test
%        testStat: test statistic

%Set kernel size to median distance between points, if no kernel specified


%Copyright (c) Arthur Gretton, 2007
%31/05/07

function [thresh,testStat,pval] = hsicTestBoot(X,Y,alpha,params);

    

%If this is 1, save the test threshold for later use
doSaveThresh = 0;


m=size(X,1);
dx=size(X,2);



%Set kernel size to median distance between points, if no kernel specified
%Use at most 100 points to save time.
if params.sigx == -1
    size1=size(X,1);
    if size1>100
      Xmed = X(1:100,:);
      size1 = 100;
    else
      Xmed = X;
    end
    G = sum((Xmed.*Xmed),2);
    Q = repmat(G,1,size1);
    R = repmat(G',size1,1);
    dists = Q + R - 2*Xmed*Xmed';
    dists = dists-tril(dists);
    dists=reshape(dists,size1^2,1);
    params.sigx = sqrt(0.5*median(dists(dists>0)));  %rbf_dot has factor of two in kernel
end

if params.sigy == -1
    size1=size(Y,1);
    if size1>100
      Ymed = Y(1:100,:);
      size1 = 100;
    else
      Ymed = Y;
    end    
    G = sum((Ymed.*Ymed),2);
    Q = repmat(G,1,size1);
    R = repmat(G',size1,1);
    dists = Q + R - 2*Ymed*Ymed';
    dists = dists-tril(dists);
    dists=reshape(dists,size1^2,1);
    params.sigy = sqrt(0.5*median(dists(dists>0)));
end





%Compute the Gram matrices
K = rbf_dot(X,X,params.sigx);
L = rbf_dot(Y,Y,params.sigy);


bone = ones(m,1);
H = eye(m)-1/m*ones(m,m);

Kc = H*K*H;
testStat = 1/m^2*sum(sum(Kc'.*L));


threshFileName = strcat('hsicTestThresh',num2str(m),'_',num2str(dx));



if ~exist(strcat(threshFileName,'.mat'),'file') || params.bootForce==1
  
%  disp(strcat('Generating new threshold: ',threshFileName))
  
  HSICarr = zeros(params.shuff,1);
  for whichSh=1:params.shuff
    
    [notUsed,indL] = sort(rand(m,1));
    HSICarr(whichSh) = 1/m^2*sum(sum(Kc'.*L(indL,indL)));
    
  end 

  HSICarr = sort(HSICarr);
  thresh = HSICarr(round((1-alpha)*params.shuff));
  pval = sum(HSICarr>testStat)/params.shuff;

  if doSaveThresh
    save(threshFileName,'thresh','HSICarr');  
  end
    
else
  load(threshFileName);
end
