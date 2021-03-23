%% var_to_tsdata
%
% Generate random multi-trial Gaussian VAR time series
%
% <matlab:open('var_to_tsdata.m') code>
%
%% Syntax
%
%     [X,E,mtrunc] = var_to_tsdata(A,SIG,m,N,mtrunc,decayfac)
%
%% Arguments
%
% See also <mvgchelp.html#4 Common variable names and data structures>.
%
% _input_
%
%     A          VAR coefficients matrix
%     SIG        residuals covariance matrix
%     m          number of observations per trial
%     N          number of trials (default: 1)
%     mtrunc     number of initial time observations to truncate or (default) empty for automatic calculation
%     decayfac   initial transients decay factor (default: 100)
%
% _output_
%
%     X          multi-trial Gaussian VAR time series
%     E          residuals time series
%     mtrunc     actual number of initial time steps truncated
%
%% Description
%
% Return |N| time series of length |m| sampled from a VAR model with
% coefficients matrix |A|, and iid Gaussian residuals with covariance matrix
% |SIG|:
%
% <<eq_var.png>>
%
% If |mtrunc| is supplied it is taken to be the the number of initial
% (non-stationary transient) observations to truncate; otherwise (default) the
% spectral radius of |A| (see function <var_specrad.html |var_specrad|>) is
% calculated and used to estimate a suitable number |mtrunc| of observations to
% assumed stationarity (roughly till autocovariance decays to near zero). Set
% |decayfac| larger for longer settle time.
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
% <var_specrad.html |var_specrad|> |
% <genvar.html |genvar|>
%
% (C) Lionel Barnett and Anil K. Seth, 2012. See file license.txt in
% installation directory for licensing terms.
%
%%

function [X,E,mtrunc] = var_to_tsdata(A,SIG,m,N,mtrunc,decayfac)

if nargin < 4 || isempty(N), N = 1; end % single trial

if nargin < 5 || isempty(mtrunc)
    if nargin < 6 || isempty(decayfac), decayfac = 100; end % should be more than enough...
    rho = var_specrad(A);
    assert(rho < 1,'unstable VAR');
    mtrunc = round((log(eps)-decayfac)/log(rho)); % enough time for autocovariance to decay to fp accuracy (and then some)
else
    assert(isscalar(mtrunc) && isint(mtrunc) && mtrunc >= 0,'''mtrunc'' parameter must be a non-negative integer');
end

[C,cholp] = chol(SIG,'lower');
assert(cholp == 0,'covariance matrix not positive-definite');

n = size(A,1);

if N > 1 % multi-trial

    X = zeros(n,m,N);
    if nargout > 1
        E = zeros(n,m,N);
        for r = 1:N
            [X(:,:,r),E(:,:,r)] = genvar(A,C*randn(n,m+mtrunc),mtrunc);
        end
    else
        for r = 1:N
            X(:,:,r) = genvar(A,C*randn(n,m+mtrunc),mtrunc);
        end
    end

else

    if nargout > 1
        [X,E] = genvar(A,C*randn(n,m+mtrunc),mtrunc);
    else
        X = genvar(A,C*randn(n,m+mtrunc),mtrunc);
    end

end
