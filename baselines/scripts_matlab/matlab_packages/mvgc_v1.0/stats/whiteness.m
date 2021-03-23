%% whiteness
%
% Durbin-Watson test for whiteness (no serial correlation) of VAR residuals
%
% <matlab:open('whiteness.m') code>
%
%% Syntax
%
%     [dw,pval] = whiteness(X,E)
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
%     dw         vector of Durbin-Watson statistics
%     pval       vector of p-values
%
%% Description
%
% Returns the Durbin-Watson test statistic |dw| along with p-values |pval|
% for a multivariate regression of time series data |X| with residuals |E|
% (may be single- or multi-trial). This routine computes the statistic
% separately for each individual variable in |X|.
%
% A standard rule of thumb is that |dw < 1| or |dw > 3| indicates a high
% chance of residuals serial correlation; this implies poor VAR model fit.
%
% NOTE: To test for significance you should correct for multiple null
% hypotheses, or test for false discovery rate; see <significance.html
% |significance|>.
%
%% References
%
% [1] J. Durbin and G. S. Watson, "Testing for Serial Correlation in Least
% Squares Regression I", _Biometrika_, 37, 1950.
%
% [2] A. Bhargava, L. Franzini and W. Narendranathan, "Serial Correlation and
% the Fixed Effects Model", _Review of Economic Studies_, 49, 1982.
%
%% See also
%
% <significance.html |significance|> |
% <mvgc_demo_stats.html |mvgc_demo_stats|>
%
%% Copyright notice
%
% [(C)] _Lionel Barnett and Anil K. Seth, 2012. See file
% <matlab:open('license.txt') license.txt> in root directory for licensing
% terms._
%
%%

function [dw,pval] = whiteness(X,E)

[n,m,N] = size(X);
assert(m >= n,'too few observations');
assert(size(E,1) == n && size(E,3) == N,'residuals don''t match time series data');

X = demean(X);

dw = zeros(1,n);
pval = zeros(1,n);
for i = 1:n
    Ei = squeeze(E(i,:,:));
    [dw(i),pval(i)] = durbinwatson(X(:,:),Ei(:)); % concatenate trials (ref. [2])
end

function [dw,pval] = durbinwatson(X,E)

[n,m] = size(X);

% calculate Durbin Watson (DW) statistic: rule of thumb: if < 1 or > 3 then
% high chance of residual correlation

dw = sum(diff(E).^2)/sum(E.^2);

% calculate critical values for the DW statistic using approx method (ref. [1])

A = X*X'; 
B = filter([-1,2,-1],1,X');
B([1,m],:) = (X(:,[1,m])-X(:,[2,m-1]))';
D = B/A;
C = X*D;
nu1 = 2*(m-1)-trace(C);
nu2 = 2*(3*m-4)-2*trace(B'*D)+trace(C*C);
mu = nu1/(m-n);
sigma = sqrt(2/((m-n)*(m-n+2))*(nu2-nu1*mu));

% evaluate p-value

pval = normcdf(dw,mu,sigma);
pval = 2*min(pval,1-pval); % two tailed test
