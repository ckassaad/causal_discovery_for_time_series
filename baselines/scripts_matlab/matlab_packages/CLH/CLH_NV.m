%% CLH_NV
%
% Summary:
% A Causal inference algorithm for Linear processes with Hidden confounders
% and Non-gaussian noise based on Variational expectation maximization.
% This is an implementation of Algorithm 1 in the paper "Causal Inference 
% by Identification of Vector Autoregressive Processes with Hidden
% Components". The parameter names are chosen in compliance with the
% paper's notation. Basically, this function, CLH_NV, is a wrapper for 
% "itsc.m", which uses different notation and thus may be confusing. 
% Authors are Mingming Gong and Philipp Geiger.
%
% Inputs:
% - X: observed time series with one column for each time step
% - K_Z: dimension of the assumed hidden process Z which together with X
%   forms a VAR process
%
% Output:
% - A_est: the estimated transition matrix A of the complete process (X,Z).
%   Note that there are only strong theretical gurantees for the entries of
%   A that capture the direct interaction of the components of X, i.e. the
%   upper left dim(X) times dim(X) submatrix of A.
% - log_likelihood: the log_likelihood for A_est and the remaining
%   parameters of the model, which are not returned.


function [A_est, log_likelihood] = CLH_NV(X, K_Z, options)


% Preliminaries:

addpath(genpath('include'));


% Infer the remaining parameters:

K_X = size(X,1);


% Add some parameters and switch from paper notation to Mingming's 
% notation:

opts = options;
opts.n = size(X) + K_Z;
opts.m = K_Z;
opts.parsEM.A = XZ_to_ZX_representation(options.additional.A_init, K_X);


% Run the variational EM estimator (which is based on Kevin Murphy's library):

[A_Hat_ZX_repr, sigmaHat, wHat, loglAll] = itsc(X, opts); % core estimation
A_Hat = ZX_to_XZ_representation(A_Hat_ZX_repr, K_X); % switch to the usual ordering (X,Z)


% Output:

A_est = A_Hat;
log_likelihood = loglAll(end);
