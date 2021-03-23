%% var_to_tsdata_nonstat
%
% Generate random multi-trial non-stationary Gaussian VAR time series
%
% <matlab:open('var_to_tsdata_nonstat.m') code>
%
%% Syntax
%
%     X = var_to_tsdata_nonstat(A,SIG,N)
%
%% Arguments
%
% See also <mvgchelp.html#4 Common variable names and data structures>.
%
% _input_
%
%     A          time-varying VAR coefficients matrix
%     SIG        time-varying residuals covariance matrix
%     N          number of trials (default: 1)
%
% _output_
%
%     X          multi-trial Gaussian VAR time series
%     E          residuals time series
%     mtrunc     number of initial time steps truncated
%
%% Description
%
% Return |N| independent random Gaussian non-stationary VAR time series, with
% time-varying coefficients |A| and residuals covariances |SIG|. The last index
% in |A| and |SIG| is the time index (the number of observations is taken as the
% size of the last dimension). If |SIG| is 2 dimensional, then it is replicated
% at each time.
%
% Useful for generating test data; see <mvgc_demo_nonstationary.html
% |mvgc_demo_nonstationary|>.
%
%% References
%
% [1] L. Barnett and A. K. Seth,
% <http://www.sciencedirect.com/science/article/pii/S0165027013003701 The MVGC
%     Multivariate Granger Causality Toolbox: A New Approach to Granger-causal
% Inference>, _J. Neurosci. Methods_ 223, 2014
% [ <matlab:open('mvgc_preprint.pdf') preprint> ].
%
%% See also
%
% <var_to_tsdata.html |var_to_tsdata|> |
% <genvar_nonstat.html |genvar_nonstat|> |
% <mvgc_demo_nonstationary.html |mvgc_demo_nonstationary|>
%
% (C) Lionel Barnett and Anil K. Seth, 2012. See file license.txt in
% installation directory for licensing terms.
%
%%

function X = var_to_tsdata_nonstat(A,SIG,N)

if nargin < 3 || isempty(N), N = 1; end % single trial

[n,n1,p,m] = size(A);
assert(n1 == n,'VAR coefficients matrices have bad shape');
[n2,n3,m1] = size(SIG);
assert(n2 == n && n3 == n,'residuals covariance matrices do not match VAR coefficients matrices');
statsig = m1 == 1;
if statsig
    [SIG,cholp] = chol(SIG,'lower');
    assert(cholp == 0,'covariance matrix not positive-definite');
else
    assert(m1 == m,'residuals covariance matrices do not match VAR coefficients matrices');
end

X = zeros(n,m,N);

for r = 1:N % for each realization
    X(:,:,r) = genvar_nonstat(A,SIG,n,p,m,statsig);
end

function X = genvar_nonstat(A,SIG,n,p,m,statsig)

% initialise to Gaussian white noise

if statsig
    X = SIG*randn(n,m); % "SIG" is actually Cholesky matrix
else
    X = zeros(n,m);
    for t = 1:m % for each time step
        [C,cholp] = chol(SIG(:,:,t),'lower');
        assert(cholp == 0,'covariance matrix not positive-definite at time step %d',t);
        X(:,t) = C*randn(n,1);
    end
end

% loop through time steps

for t = p+1:m   % for each time step
    for k = 1:p % for each lag
        X(:,t) = X(:,t) + A(:,:,k,t)*X(:,t-k);
    end
end
