%% mvgc_kpss (EXPERIMENTAL)
%
% KPSS unit root stationarity test
%
% <matlab:open('mvgc_kpss.m') code>
%
%% Syntax
%
%     [ksstat,cval] = mvgc_kpss(X,alpha,q)
%
%% Arguments
%
% See also <mvgchelp.html#4 Common variable names and data structures>.
%
% _input_
%
%     X          multi-trial time series data
%     alpha      significance level
%     q          number of lags; default: sqrt(number of observations)
%
% _output_
%
%     ksstat     matrix of KPSS test statistics
%     cval       KPSS test critical value
%
%% Description
%
% Calculates the KPSS unit root test statistics |ksstat| and critical value
% |cval| for a (possibly multi-trial) multivariate time series |X| at
% significance level |alpha| [1]. The returned test statistics matrix |ksstat|
% is |N x n| for |N| trials of an |n|-variate time series.
%
% *_Note:_* Adapted from code provided by
% <http://morgana.unimore.it/forni_mario/ Mario Forni>.
%
%% References
%
% [1] D. Kwiatkowski, P. C. B. Phillips, P. Schmidt, and Y. Shin, "Testing the
% Null Hypothesis of Stationarity against the Alternative of a Unit Root",
% _J. Econometrics_ 54, 1992.
%
%% See also
%
% <mvgc_adf.html |mvgc_adf|> |
% <mvgc_demo_stats.html |mvgc_demo_stats|>
%
% (C) Lionel Barnett and Anil K. Seth, 2012. See file license.txt in
% installation directory for licensing terms.
%
%%

function [ksstat,cval] = mvgc_kpss(X,alpha,q)

ctable = [0.119 0.146 0.216]; % critical values for 0.1,0.05,0.01 respectively
pvals  = [0.100 0.050 0.010];
 
[n,m,N] = size(X);
assert(m >= n,'fewer observations than variables, exiting');

if nargin < 3 || isempty(q)
    q = floor(sqrt(m)); % default
end

[tf,pidx] = ismember(alpha,pvals);
if ~tf
    adiff = abs(alpha-pvals);
    pidx = find(min(adiff) == adiff);
    alpha = pvals(pidx);
    fprintf(2,'WARNING(mvgc_kpss): using nearest p-value = %f\n',alpha);
end
cval = ctable(pidx);

% calculate KPSS statistics

ksstat = zeros(N,n);
z = [ones(1,m); 1:m]; % deterministic trend
for r = 1:N
    for i = 1:n
        e = X(i,:,r)-(X(i,:,r)/z)*z; % residuals against trend
        prod = zeros(q,1);
        for j = 1:q
            prod(j) = e(j+1:m)*e(1:m-j)';
        end
        s2 = e*e'+2*(1-(1:q)/(q+1))*prod;
        S = cumsum(e);
        ksstat(r,i) = (S*S')/(s2*m);
    end
end
