function [vals2, cause, aic] = lassoGranger(series, P, lambda, type)
% This code finds the Granger causality relationship among the input time series
% For more information on Lasso-Granger, see the following paper:
% 	A. Arnold, Y. Liu, and N. Abe. Temporal causal modeling with graphical granger methods. In KDD, 2007.

% INPUTS:
% 'series': A NxT matrix 
% 'P':		Length of the lag
% 'lambda': 	Value of the penalization parameter in Lasso
% 'type': 	'l' Lasso with the given Lambda
%       	'r' Ridge regression with very small Lambda
% OUTPUTS:
% 'cause':	Sum of the Granger causality coefficient from one timeseries to anohter one
% 'vals2':	Unconverted cause
% 'aic':		To help select the value of Lambda
% Dependency: This code requires the GLMnet package to perform Lasso.

[N T] = size(series);
Am = zeros(T-P, P*N);
bm = zeros(T-P, 1);
for i = (P+1):T
    bm(i-P) = series(1, i);
    Am(i-P, :) = reshape(fliplr(series(:, (i-P):(i-1))), 1, N*P);
end

if type == 'r'
    vals2 = (Am'*Am+(1e-6)*eye(P*N))\(Am'*bm); 
else
    % Lasso using GLMnet
    opt = glmnetSet;
    opt.lambda = lambda;
    opt.alpha = 1;
    fit = glmnet(Am, bm, 'gaussian', opt);
    vals2 = fit.beta;
end

% Outputting AIC metric for variable tuning purpose
th = 0;
aic = norm(Am*vals2 - bm)^2/(T-P) + (sum(abs(vals2) > th))*2/(T-P);

% Reformatting the results into NxP matrix
n1Coeff = zeros(N, P);
for i = 1:N
    n1Coeff(i, :) = vals2( ((i-1)*P+1): (i*P));
end

sumCause = sum(abs(n1Coeff), 2);
if type == 'r'
    cause = sumCause;
else
    cause = (sumCause > th).*sumCause;
end