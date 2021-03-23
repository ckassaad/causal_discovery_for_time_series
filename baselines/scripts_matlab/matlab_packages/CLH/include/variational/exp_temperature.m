% data temprature
load('../data/temperature.mat');

x1 = temperature(2001:1:3200,[2 3 4])'; 

x1 = x1 - repmat(mean(x1')', 1,size(x1,2));
for i = 1:size(x1,1)
    x1(i,:) = x1(i,:)./std(x1(i,:));
end

% k1 = 2;
m = 1;
n = size(x1,1);
T = size(x1,2);
w = repmat([0.7 0.3],[n 1]);
mu = repmat([-0.25 0.75],[n 1]);
sigma = repmat([1 1], [n 1]);

A1 = x1(:,2:T) * x1(:,1:T-1)' * inv(x1(:,1:T-1) * x1(:,1:T-1)');
B1 = x1(m+1:end,2:T) * x1(m+1:end,1:T-1)' * inv(x1(m+1:end,1:T-1) * x1(m+1:end,1:T-1)');

opts.x = x1;
opts.m = m;
% estimation
opts.parsEM.thres = 1e-6;
opts.parsEM.maxIter = 300;
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