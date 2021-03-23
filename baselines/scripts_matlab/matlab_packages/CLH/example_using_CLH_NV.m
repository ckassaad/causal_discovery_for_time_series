%% Example for using CLH_NV
%
% Summary:
% An example for how to apply CLH_NV on a simulated data set.


%% Preliminaries:

%addpath(genpath('include'));


%% Simulate data:

% From which model:

K_X = 2; % two observed components
K_Z = 1; % one hidden
K = K_X + K_Z;
L = 4000; % sample length;

% Fix some true A:

A = [0.4 0.2 0.3; 0 0.5 0.4; 0 0 0.9]

% Generate super-Gaussian noise (code snippet taken from Mingming):

mu = repmat([0 0],[K,1]);
w = repmat([0.8 0.2],[K 1]);
sigma = repmat([0.1 1], [K 1]);
p_ming = zeros(K, L);
E = zeros(K, L);
for i = 1:K
    p_ming(i,:) = random('binomial', 1, w(i,2), 1, size(E,2)) + 1;
    for j = 1:2
        EE = mu(i,j) + sigma(i,j) * randn(1,sum(p_ming(i,:)==j));
        E(i,p_ming(i,:)==j) = EE;
    end
end

% Generate the complete process W and take first K_X components as X:

N = E; % take noise from above
W = zeros(K,L);
for i = (1 + 1):L
    W(:,i) = A * W(:,i-1) + N(:,i);
end
X = W(1:K_X,:); % ok, we are done with generating X


%% Model check:

[p_val_indep, p_val_Gauss] = CLH_MC(X)


%% Set the options and run CLH_NV:

% Note: generally it makes sense to start itsc with many different (random)
% initializations of the parameters to try to find a global maximum of the
% log likelihood (i.e. at the end pick that A_est with the highest
% log_likelihood).
% Here we only do one initialization, partially random, partially fixed. It
% seems to make sense to use linear Granger on X as initial value for A.

options = struct;

% Initial value for A:

options.additional.A_init = rand(K, K) - 0.5; % initialize A randomly, except ...
options.additional.A_init(1:K_X, 1:K_X) = X(:,2:L) * X(:,1:L-1)' * inv(X(:,1:L-1) * X(:,1:L-1)'); % ... for initialization of that part of A that underlies X take linear Granger on X

% Convergence threshold and max. number of iterations:

options.parsEM.thres = 1e-6; % threshold for stopping optimization ...
options.parsEM.maxIter = 300; % ... otherwise: max. number of iterations
options.parsEM.maxIter0 = 5;

% Mixture of Gaussians parameters:

options.parsEM.mu = repmat([-0.25 0.75],[K 1]);
options.parsEM.sigma = repmat([1 1], [K 1]);
options.parsEM.w = repmat([0.7 0.3],[K 1]);

% Forgot what the following option is for, maybe to automatically estimate certain initial paramerts:

options.parsEM.updatePrior = 1;

% Run CLH_NV:

[A_est, log_likelihood] = CLH_NV(X, K_Z, options)
