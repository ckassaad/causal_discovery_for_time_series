%% rsquared
%
% [R^2] and adjusted [R^2] statistics
%
% <matlab:open('rsquared.m') code>
%
%% Syntax
%
%     [RSQ,RSQADJ] = rsquared(X,E)
%
%% Arguments
%
% See also <mvgchelp.html#4 Common variable names and data structures>.
%
% _input_
%
%     X          multi-trial time series data
%     E          residuals time series
%
% _output_
%
%     RSQ        R^2 statistic
%     RSQADJ     adjusted R^2 statistic
%
%% Description
%
% Return [R^2] and adjusted [R^2] statistics. The time series data X and
% residual errors E may be single- or multi-trial time series.
%
%% See also
%
% <mvgc_demo_stats.html |mvgc_demo_stats|>
%
%% Copyright notice
%
% [(C)] _Lionel Barnett and Anil K. Seth, 2012. See file
% <matlab:open('license.txt') license.txt> in root directory for licensing
% terms._
%
%%

function [RSQ,RSQADJ] = rsquared(X,E)

[n,m,N] = size(X);
p = m-size(E,2);     % number of lags in model
assert(m >= n,'too few observations');
assert(size(E,1) == n && size(E,3) == N,'residuals don''t match data');
assert(p > 0,'bad number of lags');

X = demean(X);

if N > 1 % multi-trial
    X = X(:,p+1:m,:);
    X = X(:,:);              % stack data
    E = E(:,:);              % stack residuals
    s = N*(m-p);             % sample size (number of observations)
else
    X = X(:,p+1:m);
    s = m-p;                 % sample size (number of observations)
end
r = n*n*p;                   % number of regressors in model

SSR = sum(E.^2,2)';          % residuals sum of squares
SST = sum(X.^2,2)';          % total sum of squares

RSQ = 1 - SSR./SST;

RSQADJ = 1 - ((s-1)/(s-r-1))*(1-RSQ);
