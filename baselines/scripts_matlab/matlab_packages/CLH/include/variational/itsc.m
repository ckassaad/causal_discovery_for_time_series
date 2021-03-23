function [A_Hat, sigmaHat, wHat, fVal_all] = itsc(X, opts)
% Identifying structural coefficients of vector
% autoregressive processes with hidden components
% Inputs:
%       X = data points
%       opts = parameters
% Outputs:
%       A_Hat = estimated mixing matrix
%       sigmaHat = estimated mixture of gaussian parameters
%       wHat  = estimated prior parameters

%% initilize parameters
m = opts.m;
T = size(X,2);
n = size(X,1) + m;
A = opts.parsEM.A;
mu = opts.parsEM.mu;
sigma = opts.parsEM.sigma;
w = opts.parsEM.w;
updatePrior = opts.parsEM.updatePrior;
thres = opts.parsEM.thres;
maxIter = opts.parsEM.maxIter;
maxIter0 = opts.parsEM.maxIter0;

E = A(1:m,1:m);
D = A(1:m,m+1:end);
C = A(m+1:end,1:m);
B = A(m+1:end,m+1:end);
A_Hat = A;
wHat = w;
muHat = mu;
sigmaHat = sigma;

%% EM main iterations
% initialize posterior of l^x and l^t
qlzi = repmat(w(1:m,:),[1, 1, T-1]);
qlxi = repmat(w(m+1:n,:),[1, 1, T-1]);

iter = 0;
fVal_all = [];
fVal_prev = inf;
fVal = inf;
while iter < 10 || iter < maxIter && abs(fVal_prev-fVal) >= abs(fVal_prev)*thres
    fVal_prev = fVal;
    fprintf('iter %d, fval = %.5f\n', iter, fVal);
    A_Hat
    if updatePrior == true
        wHat
        muHat
        sigmaHat
        fprintf('------------------------------------\n');
    end
    %% E step
    iter0 = 0;
    fVal0_prev = inf;
    fVal0 = fVal;
    fVal0_all = [];
    while iter0 < 1 || iter0 < maxIter0 && abs(fVal0_prev-fVal0) >= abs(fVal0_prev)*thres
        fVal0_prev = fVal0;
        fVal0_all = [fVal0_all fVal0];
        % alternate between z and l
        [posMeanZ, posCovZ, posCovZZ, loglik] = evaluateZ(X, E, D, C, B, qlzi, qlxi, mu, sigma);
        fVal0 = istcBound(w, qlzi, qlxi, loglik);
        [qlzi, qlxi] = evaluateL(X, E, D, C, B, posMeanZ, posCovZ, posCovZZ, w, mu, sigma);
        iter0 = iter0 + 1;
%         plot(fVal0_all);
    end
    iter = iter+1;
    fVal = fVal0;
    fVal_all = [fVal_all fVal];
    
    numGauss = size(mu,2);
    sigmaz = sigma(1:m,:);
    sigmax = sigma(m+1:n,:);
    muz = mu(1:m,:);
    mux = mu(m+1:n,:);
    wz = w(1:m,:);
    wx = w(m+1:n,:);
    
    %% M step
    DX = D*X(:,1:T-2);
    
    % update E
    for i =1:m
        E(i,:) = (sum(repmat(sum(qlzi(i,:,2:T-1)./repmat(sigmaz(i,:).^2,[1 1 T-2]),2),...
            [m m 1]).*posCovZ(:,:,1:T-2),3)\...
            sum(sum(repmat(qlzi(i,:,2:T-1)./repmat(sigmaz(i,:).^2,[1 1 T-2]),[m 1 1])...
            .*(repmat(reshape(posCovZZ(i,:,2:T-1),[m 1 T-2]),[1 numGauss 1])...
            -repmat(reshape(posMeanZ(:,1:T-2),[m 1 T-2]).*repmat(reshape(DX(i,:),[1 1 T-2]),[m 1 1]),[1 numGauss 1])...
            -repmat(reshape(posMeanZ(:,1:T-2),[m 1 T-2]),[1 numGauss]).*repmat(muz(i,:),[m 1 T-2])),2),3))';
    end
    
    covX = zeros(n-m,n-m,T-2);
    for t=1:T-1
        covX(:,:,t) = X(:,t)*X(:,t)';
    end
    EZ = E*posMeanZ(:,1:T-2);
    
    % update D
    for i =1:m
        D(i,:) = (sum(repmat(sum(qlzi(i,:,2:T-1)./repmat(sigmaz(i,:).^2,[1 1 T-2]),2),...
            [n-m, n-m 1]).*covX(:,:,1:T-2),3)\...
            sum(sum(repmat(qlzi(i,:,2:T-1)./repmat(sigmaz(i,:).^2,[1 1 T-2]),[n-m 1 1])...
            .*repmat(reshape(X(:,1:T-2),[n-m 1 T-2]),[1 numGauss 1])...
            .*(repmat(reshape(posMeanZ(i,2:T-1),[1 1 T-2]),[n-m numGauss 1])-repmat(reshape(EZ(i,1:T-2),[1 1 T-2]),[n-m numGauss 1])...
            -repmat(muz(i,:),[n-m 1 T-2])),2),3))';
    end
    
    BX = B*X(:,1:T-1);
    
    % update C
    for i =1:n-m
        C(i,:) = (sum(repmat(sum(qlxi(i,:,1:T-1)./repmat(sigmax(i,:).^2,[1 1 T-1]),2),...
            [m m 1]).*posCovZ(:,:,1:T-1),3)\...
            sum(sum(repmat(qlxi(i,:,1:T-1)./repmat(sigmax(i,:).^2,[1 1 T-1]),[m 1 1])...
            .*repmat(reshape(posMeanZ(:,1:T-1),[m 1 T-1]),[1 numGauss 1])...
            .*(repmat(reshape(X(i,2:T),[1 1 T-1]),[m numGauss 1])-repmat(reshape(BX(i,1:T-1),[1 1 T-1]),[m numGauss 1])...
            -repmat(mux(i,:),[m 1 T-1])),2),3))';
    end
    
    CZ = C*posMeanZ(:,1:T-1);
    % update B
    for i =1:n-m
        B(i,:) = (sum(repmat(sum(qlxi(i,:,1:T-1)./repmat(sigmax(i,:).^2,[1 1 T-1]),2),...
            [n-m n-m 1]).*covX(:,:,1:T-1),3)\...
            sum(sum(repmat(qlxi(i,:,1:T-1)./repmat(sigmax(i,:).^2,[1 1 T-1]),[n-m 1 1])...
            .*repmat(reshape(X(:,1:T-1),[n-m 1 T-1]),[1 numGauss 1])...
            .*(repmat(reshape(X(i,2:T),[1 1 T-1]),[n-m numGauss 1])-repmat(reshape(CZ(i,1:T-1),[1 1 T-1]),[n-m numGauss 1])...
            -repmat(mux(i,:),[n-m 1 T-1])),2),3))';
    end
    
    if updatePrior == true
        % upate w
        w = updateW(qlzi, qlxi, numGauss, m, n);
        
        % update mu
        mu = updateMu(qlzi, qlxi, posMeanZ, numGauss, m, n, E, D, C, B, X, T);
        
        % update sigma
        sigma = updateSigma(qlzi, qlxi, mu, posMeanZ, posCovZ, posCovZZ, numGauss, m, n, E, D, C, B, X, T);
    end
    
    A_Hat = A;
    A_Hat(1:m,1:m) = E;
    A_Hat(1:m,m+1:end) = D;
    A_Hat(m+1:end,1:m) = C;
    A_Hat(m+1:end,m+1:end) = B;
    muHat = mu;
    sigmaHat = sigma;
    wHat = w;
end

function w = updateW(qlzi, qlxi, numGauss, m, n)
%% Update w
% the posterior of qi
w = zeros(n, numGauss);
w(1:m,:) = sum(qlzi,3)/size(qlzi,3);
w(m+1:n,:) = sum(qlxi,3)/size(qlxi,3);

function mu = updateMu(qlzi, qlxi, posMeanZ, numGauss, m, n, E, D, C, B, X, T)
% update mu

muz = sum(qlzi.*repmat(reshape(posMeanZ - [zeros(m,1), E*posMeanZ(:,1:end-1)+D*X(:,1:end-2)],...
    [m, 1 T-1]),[1 numGauss]),3)./sum(qlzi,3);
mux = sum(qlxi.*repmat(reshape(X(:,2:end) - (C*posMeanZ...
    +B*X(:,1:end-1)),[n-m, 1 T-1]),[1 numGauss]),3)./sum(qlxi,3);
mu = [muz; mux];


function sigma = updateSigma(qlzi, qlxi, mu, posMeanZ, posCovZ, posCovZZ, numGauss, m, n, E, D, C, B, X, T)
% update sigma
muz = mu(1:m,:);
mux = mu(m+1:n,:);
DX = D*X(:,1:T-2);
EZ = E*posMeanZ(:,1:T-2);
BX = B*X(:,1:T-1);
CZ = C*posMeanZ(:,1:T-1);
diagPosCovZ = zeros(m, T-1);
for i = 1:T-1
    diagPosCovZ(:,i) = diag(posCovZ(:,:,i));
end
ztEztm1 = zeros(m, T-2);
for i = 1:T-2
    ztEztm1(:,i) = diag(posCovZZ(:,:,i+1)*E');
end
ztm1EEztm1 = zeros(m, T-2);
for i = 1:T-2
    ztm1EEztm1(:,i) = diag(E*posCovZ(:,:,i)*E');
end
ztm1CCztm1 = zeros(n-m, T-1);
for i = 1:T-1
    ztm1CCztm1(:,i) = diag(C*posCovZ(:,:,i)*C');
end

sigmaz2 = sum(qlzi.*repmat(reshape(diagPosCovZ + ...
    [zeros(m,1) -2*ztEztm1-2*posMeanZ(:,2:end).*DX + ztm1EEztm1 + 2*EZ.*DX + DX.^2],...
    [m 1 T-1]),[1 numGauss]),3)./sum(qlzi,3) - muz.^2;
sigmax2 = sum(qlxi.*repmat(reshape(X(:,2:end).^2 - 2*CZ.*X(:,2:end)-2*BX.*X(:,2:end) ...
    + ztm1CCztm1 + 2*CZ.*BX + BX.^2,[n-m 1 T-1]),[1 numGauss]),3)./sum(qlxi,3) - mux.^2;
sigma = [sqrt(sigmaz2); sqrt(sigmax2)];


function [qlzi, qlxi] = evaluateL(X, E, D, C, B, posMeanZ, posCovZ, posCovZZ, w, mu, sigma)
%% posterior of l^x, l^z

m = size(E,1);
n = size(X,1)+m;
sigmaz = sigma(1:m,:);
sigmax = sigma(m+1:n,:);
muz = mu(1:m,:);
mux = mu(m+1:n,:);
wz = w(1:m,:);
wx = w(m+1:n,:);
condProbZ = conditionalProbZ(X, E, D, muz, sigmaz, posMeanZ, posCovZ, posCovZZ);
qlzi = condProbZ.*repmat(wz,[1 1 size(condProbZ,3)])./repmat(sum(condProbZ.*repmat(wz,[1 1 size(condProbZ,3)]),2), [1 size(wz,2)]);
condProbX = conditionalProbX(X, C, B, mux, sigmax, posMeanZ, posCovZ);
qlxi = condProbX.*repmat(wx,[1 1 size(condProbX,3)])./repmat(sum(condProbX.*repmat(wx,[1 1 size(condProbX,3)]),2), [1 size(wx,2)]);

function condProb = conditionalProbZ(X, E, D, muz, sigmaz, posMeanZ, posCovZ, posCovZZ)
%% Conditional probability, conditioned on l^z

T = size(X,2);
m = size(E,1);
numGauss = size(muz,2);
logCondProb = zeros(m, numGauss, T-1);
logCondProb(:,:,1) = -0.5* ((repmat(diag(posCovZ(:,:,1)),[1 numGauss]) - ...
    2*repmat(posMeanZ(:,1),[1 numGauss]).*muz + muz.^2)./(sigmaz.^2)+2*log(sigmaz));
for t=2:T-1
    logCondProb(:,:,t) = -0.5*((repmat(diag(posCovZ(:,:,t)),[1 numGauss]) -2*repmat(diag(posCovZZ(:,:,t)*E'),[1 numGauss])...
        + repmat(diag(E*posCovZ(:,:,t-1)*E'), [1 numGauss]) - 2*repmat((posMeanZ(:,t)-...
        E*posMeanZ(:,t-1)).*(D*X(:,t-1)),[1 numGauss]) + repmat((D*X(:,t-1)).^2,[1 numGauss]) ...
        - 2*repmat(posMeanZ(:,t)-E*posMeanZ(:,t-1)-D*X(:,t-1),[1 numGauss]).*muz + muz.^2)./(sigmaz.^2)+2*log(sigmaz));
end
condProb = exp(logCondProb);

for t = 1:T-1
    [maxRow, maxRowIdx] = max(logCondProb(:,:,t),[],2);
    logCondProb(maxRow>700,:,t) = logCondProb(maxRow>700,:,t) - repmat(maxRow(maxRow>700,:),[1 size(sigmaz,2)]) + 700;
    condProb(maxRow>700,:,t) = exp(logCondProb(maxRow>700,:,t));
    %     [minRow, minRowIdx] = min(logCondProb(:,:,t),[],2);
    %     logCondProb(minRow<-700,:,t) = logCondProb(minRow<-700,:,t) - repmat(minRow(minRow<-700,:),[1 size(sigmaz,2)]) - 700;
    %     condProb(minRow<-700,:,t) = exp(logCondProb(minRow<-700,:,t));
    %     condProb(logCondProb>700) = exp(700);
end


function condProb = conditionalProbX(X, C, B, mux, sigmax, posMeanZ, posCovZ)
%% Conditional probability, conditioned on l^x

T = size(X,2);
numGauss = size(mux,2);
logCondProb = zeros(size(B,1), numGauss, T-1);
% logCondProb(:,:,1) = -0.5* ((repmat(X(:,1).^2,[1 numGauss]) - ...
%     2*repmat(X(:,1),[1 numGauss]).*mux + mux.^2)./(sigmax.^2)+2*log(sigmax));
for t=2:T
    logCondProb(:,:,t-1) = -0.5*(((repmat(X(:,t)-B*X(:,t-1),[1 numGauss])-mux).^2-...
        2*(repmat(X(:,t)-B*X(:,t-1),[1 numGauss])-mux).*repmat(C*posMeanZ(:,t-1),[1 numGauss])...
        +repmat(diag(C*posCovZ(:,:,t-1)*C'),[1 numGauss]))./(sigmax.^2)+2*log(sigmax));
end
condProb = exp(logCondProb);

for t = 1:T-1
    [maxRow, maxRowIdx] = max(logCondProb(:,:,t),[],2);
    logCondProb(maxRow>700,:,t) = logCondProb(maxRow>700,:,t) - repmat(maxRow(maxRow>700,:),[1 size(sigmax,2)]) + 700;
    condProb(maxRow>700,:,t) = exp(logCondProb(maxRow>700,:,t));
    %     [minRow, minRowIdx] = min(logCondProb(:,:,t),[],2);
    %     logCondProb(minRow<-700,:,t) = logCondProb(minRow<-700,:,t) - repmat(minRow(minRow<-700,:),[1 size(sigmax,2)]) - 700;
    %     condProb(minRow<-700,:,t) = exp(logCondProb(minRow<-700,:,t));
    %     condProb(logCondProb>700) = exp(700);
end


function [posMeanZ, posCovZ, posCovZZ, loglik] = evaluateZ(X, E, D, C, B, qlzi, qlxi, mu, sigma)
%% posterior of Z via kalman smoothing
T = size(X,2);
m = size(qlzi,1);
n = size(X,1) + m;
sigmaz = sigma(1:m,:);
sigmax = sigma(m+1:n,:);
muz = mu(1:m,:);
mux = mu(m+1:n,:);

% varational parameters
SigmaHatZ = zeros(m, m, T-1);
sigmaTmp = 1/sum(1/repmat(sigmaz.^2,[1 1 T-1]) .* qlzi, 2);
for t = 1:T-1
    SigmaHatZ(:,:,t) = diag(sigmaTmp(:,t));
end

meanHatZ = sigmaTmp .* sum((repmat(muz,[1 1 T-1])+reshape(repmat([zeros(m,1),D*X(:,1:T-2)],...
    [size(muz,2) 1]),[size(muz,1) size(muz,2), T-1]))./repmat(sigmaz.^2,[1 1 T-1]) .* qlzi, 2);

SigmaHatX = zeros(n-m, n-m, T-1);
sigmaTmp = 1/sum(1/repmat(sigmax.^2,[1 1 T-1]) .* qlxi, 2);
for t = 1:T-1
    SigmaHatX(:,:,t) = diag(sigmaTmp(:,t));
end
meanHatX = sigmaTmp .* sum(repmat(mux,[1 1 T-1])./repmat(sigmax.^2,[1 1 T-1]) .* qlxi, 2);

% kalman smoother
initX = zeros(m,1);
initV = SigmaHatZ(:,:,1);
E1 = repmat(E,[1 1 T-1]);
C1 = repmat(C,[1 1 T-1]);

Y = X(:,2:end) - B*X(:,1:end-1) - reshape(meanHatX,[size(meanHatX,1) size(meanHatX,3)]);

[posMeanZ, posCovZ, posCovZZ, loglik] = kalman_smoother(Y, E1, C1, SigmaHatZ, SigmaHatX, ...
    initX, initV, 'model', 1:T-1, 'u', meanHatZ, 'B', repmat(eye(size(E,1)),[1 1 T-1]));

for t = 1:T-1
    posCovZ(:,:,t) = posCovZ(:,:,t) + posMeanZ(:,t)*posMeanZ(:,t)';
end
for t = 1:T-2
    posCovZZ(:,:,t+1) = posCovZZ(:,:,t+1) + posMeanZ(:,t+1)*posMeanZ(:,t)';
end
