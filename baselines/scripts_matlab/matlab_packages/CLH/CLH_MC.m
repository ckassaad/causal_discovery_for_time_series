%% CLH_MC
%
% Summary:
% A Model Checking algorithm for the framework of Causal inference on 
% Linear processes with Hidden confounders.
% This is an implementation of the model check suggested in the paper 
% "Causal Inference by Identification of Vector Autoregressive Processes 
% with Hidden Components" in Section 7.3. This model check should be
% applied prior to using methods such as CLH_NV to potentially falsify the
% model assumptions. Two tests are performed: (1) a test for a certain
% independence and (2) a test for Gaussianity.
% Author: Philipp Geiger.
%
% Inputs:
% - X: observed time series with one column for each time step
%
% Output:
% - non_singular: true iff certain non-singularity condition is met
% - p_val_indep: p-value for the HSCI test statistic for the mentioned
%   independence hypothesis
% - p_val_Gauss: p-value for the test for Gaussianity


function [non_singular, p_val_indep, p_val_Gauss] = CLH_MC(X)

%% Preliminaries:

addpath(genpath('include'));


%% Calculate the generalized residual R (for details see the paper)

K_X = size(X,1);
sample_length = size(X,2);

R = @(U_1,U_2)( X(:, 4:end) - U_1 * X(:, 3:end-1) - U_2 * X(:, 2:end-2) );
R_single_arg = @(U)( R( U(:, 1:K_X), U(:, (K_X + 1):(2*K_X)) ) );

autocov = @(k,l)( 1/(sample_length-3-1) * X(:, (4-k):(end-k) ) * X(:, (4-l):(end-l) )' );

M = [ autocov(1,2), autocov(1,3) ; autocov(2,2) , autocov(2,3) ];
N = [  autocov(0,2), autocov(0,3) ];

% Check whether the autocovariance submatrix that we use is non-singular:
non_singular = eq(2*K_X, rank(M));

U = linsolve(M', N')';


%% Testing 

params.bootForce=1;
params.sigx = -1;
params.sigy = -1;

params.shuff = 1000;



% Independence test using HSIC:
R_single_arg_U = R_single_arg(U);
X_truncated = X(:, 4:end); % such that X_truncated(:,1) corresponds to R_single_arg_U(:1,)
[thresh,testStat,p_val_indep] = hsicTestBoot([X_truncated(:, 1:end-3)', X_truncated(:, 2:end-2)'], R_single_arg_U(:, 4:end)' , 0.05, params );

% Test for Gaussianity:
[h,p_val_Gauss] = kstest(X');

