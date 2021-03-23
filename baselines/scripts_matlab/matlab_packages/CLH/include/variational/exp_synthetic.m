addpath('../KalmanAll/Kalman');
addpath('../KalmanAll/KPMtools');
addpath('../KalmanAll/KPMstats');
rng('default')
rng(1);
opts.expDir = 'data/synthetic';
n = 3; % number of channels, dimension
m = 1; % number of missing channels
T = 500; % length of time series
noiseType = 3;
% gaussian mixture parameters
switch noiseType
    case 1 % superGauss
        mu = repmat([0 0],[n,1]);
        w = repmat([0.8 0.2],[n 1]);
        sigma = repmat([0.1 1], [n 1]);
    case 2 % subGauss
        mu = repmat([-1 1],[n,1]);
        w = repmat([0.5 0.5],[n 1]);
        sigma = repmat([1 1], [n 1]);
    case 3 % gauss
        mu = repmat([0 0],[n,1]);
        w = repmat([0.5 0.5],[n 1]);
        sigma = repmat([1 1],[n 1]);
end
A = (rand(n)-0.5); % causal matrix
% A = [0.7 0 0 0; 0.4 0.7 0 0 ; 0.4 0 0.7 0; 0 0.4 0 0.7]; 
if ~exist(opts.expDir,'dir')
    mkdir(opts.expDir);
end
%% generate synthetic data
p = zeros(n, T);
E = zeros(n, T);
for i = 1:n
    p(i,:) = random('binomial', 1, w(i,2), 1, size(E,2)) + 1;
    for j = 1:2
        EE = mu(i,j) + sigma(i,j) * randn(1,sum(p(i,:)==j));
        E(i,p(i,:)==j) = EE;
    end
end

x1 = ones(n,T);
x1(:,1) = E(:,1);
for i=2:T
    x1(:,i) = A * x1(:,i-1) + E(:,i);
end
A1 = x1(:,2:T) * x1(:,1:T-1)' * inv(x1(:,1:T-1) * x1(:,1:T-1)');
B1 = x1(m+1:end,2:T) * x1(m+1:end,1:T-1)' * inv(x1(m+1:end,1:T-1) * x1(m+1:end,1:T-1)');

opts.x = x1;
opts.m = m;
% estimation
opts.parsEM.thres = 1e-6;
opts.parsEM.maxIter = 200;
opts.parsEM.maxIter0 = 5;
opts.parsEM.A = zeros(n);
opts.parsEM.A(m+1:n,m+1:n) = B1;
opts.parsEM.A(1:m,1:m) = rand(m);
opts.parsEM.A(1:m,m+1:n) = rand(m,n-m);
opts.parsEM.A(m+1:n,1:m) = rand(n-m,m);
opts.parsEM.mu = mu;
opts.parsEM.sigma = sigma;
opts.parsEM.w = w;
opts.parsEM.updatePrior = 1;

[A_Hat, sigmaHat, wHat, loglAll] = itsc(x1(m+1:end,:), opts);
figure, plot(loglAll);
B = A(m+1:n,m+1:n);
B2 = A_Hat(m+1:n,m+1:n);
mseB1 = sum((B1(:)-B(:)).^2)/(n-m)^2
mseB2 = sum((B2(:)-B(:)).^2)/(n-m)^2