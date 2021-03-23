% Run all
clc
clear all
addpath('glmnet_matlab/')

%% Generate a simple synthetic dataset
N = 20;     % # of time series
T = 1000;    % length of time series
sig = 0.2;
genSynth(N, T, sig);

%% Run Forward Lasso-Granger
lambda = 1e-2;
L = 1;      % Only one lag for analysis
load synth.mat
causeF = zeros(N, N);
for in = 1:N
    index = [in, 1:(in-1), (in+1):N];
    [~, temp] = lassoGranger(series(index, :), L, lambda, 'l');
    causeF(in, :) = temp([2:in, 1, (in+1):N])';
    fprintf('%c%c%c%c%c%c', 8,8,8,8,8,8);
    fprintf('%5d ', in);
end

%% Run Backward Lasso-Granger

series = fliplr(series);
causeB = zeros(N, N);
for in = 1:N
    index = [in, 1:(in-1), (in+1):N];
    [~, temp] = lassoGranger(series(index, :), L, lambda, 'l');
    causeB(in, :) = temp([2:in, 1, (in+1):N])';
    fprintf('%c%c%c%c%c%c', 8,8,8,8,8,8);
    fprintf('%5d ', in);
end

%% Visualized the Matrix
subplot(2, 1, 1)
showMatrix(A)
title('Ground Truth')

subplot(2, 1, 2)
showMatrix(0.5.*(causeF+causeB'))
title('Inferred Causality: Forward Backward')