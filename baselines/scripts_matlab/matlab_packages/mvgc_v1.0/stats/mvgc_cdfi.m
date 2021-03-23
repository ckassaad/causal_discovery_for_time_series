%% mvgc_cdfi
%
% Sample MVGC thoretical asymptotic inverse cumulative distribution function
%
% <matlab:open('mvgc_cdfi.m') code>
%
%% Syntax
%
%     x = mvgc_cdfi(P,X,p,m,N,nx,ny,nz,tstat)
%
%% Arguments
%
% See also <mvgchelp.html#4 Common variable names and data structures>.
%
% _input_
%
%     P          vector of MVGC cumulative distribution probabilities
%     X          vector of actual MVGC values
%     p          VAR model order
%     m          number of observations per trial
%     N          number of trials
%     nx         number of target ("to") variables
%     ny         number of source ("from") variables
%     nz         number of conditioning variables (default: 0)
%     tstat      statistic: 'F' or 'chi2' (default: 'F' if nx == 1, else 'chi2')
%
% _output_
%
%     x          vector of MVGC values
%
%% Description
%
% Return theoretical sample MVGC asymptotic inverse cumulative distribution
% function for actual MVGCs in vector |X|, evaluated at probabilities in vector
% P. To calculate the critical MVGC value for a significance level |alpha|,
% assume null hypothesis [H_0]: |X = 0| and set |P = 1-alpha| (see
% <mvgc_cval.html |mvgc_cval|>). For confidence intervals at level |alpha|, set
% |X| to the estimated MVGC and set |P = alpha| for lower bound and |P =
% 1-alpha| for upper bound (see <mvgc_confint.html |mvgc_confint|>).
%
% The theoretical distribution is specified by the |tstat| parameter, which may
% be  'F' for Granger's F-test (default if |nx == 1|) or 'chi2' for Geweke's
% [chi^2] test (default if |nx > 1|). For a multivariate predictee (i.e. |nx >
% 1|) only the [chi^2] test is suitable [1].
%
% For more details see <mvgc_cdf.html |mvgc_cdf|>.
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
% <mvgc_cdf.html |mvgc_cdf|> |
% <mvgc_pval.html |mvgc_pval|> |
% <mvgc_confint.html |mvgc_confint|> |
% <mvgc_cval.html |mvgc_cval|>
%
% (C) Lionel Barnett and Anil K. Seth, 2012. See file license.txt in
% installation directory for licensing terms.
%
%%

function x = mvgc_cdfi(P,X,p,m,N,nx,ny,nz,tstat)

assert(isvector(P),'probabilities must be a vector');
n = length(P);

assert(isvector(X),'MVGC values must be a vector');
assert(all(X >= 0),'MVGC values must be non-negative');
if isscalar(X)
    X = X*ones(n,1);
else
    assert(length(X) == n,'MVGC values must match evaluation values');
end

if nargin < 8 || isempty(nz), nz = 0; end % unconditional

if nargin < 9 || isempty(tstat);
    ftest = nx == 1; % default: use F-test for univariate predictee, chi2 for multivariate
else
    switch lower(tstat)
        case 'f'     % Granger F-test
            assert(nx == 1,'F-distribution is not appropriate for multivariate predictee');
            ftest = true;
        case 'chi2'  % Geweke chi2 test
            ftest = false;
        otherwise
            error('unknown distribution (must be ''chi2'' or ''F'')');
    end
end

x = zeros(n,1);
m = N*(m-p);                      % effective number of observations (p-lag autoregression loses p observations per trial)
if ftest
    if any(X > 0), fprintf(2,'WARNING (mvgc_cdfi): non-central F-test is experimental\n'); end
    d1 = p*ny;                    % #{full model parameters} - #{reduced model parameters}
    d2 = m-p*(1+ny+nz);           % #{observations} - #{full model parameters}
    mm = d2/d1;
    for i = 1:n
        if P(i) < 0 || P(i) > 1
            x(i) = NaN;
        else
            if X(i) > 0           % non-central
                FF = exp(X(i))-1; % Granger form: (RSS_reduced - RSS_full) / RSS_full
                xx = ncfinv(P(i),d1,d2,m*FF)/mm; % NOTE: non-centrality parameter factor might reasonably be m, d2 or d2-2
            else
                xx = finv(P(i),d1,d2)/mm;
            end
            x(i) = log(1+xx);     % Geweke form
        end
    end
else
    d = p*nx*ny;                  % note that d does not depend on the number of conditioning variables
    for i = 1:n
        if P(i) < 0 || P(i) > 1
            x(i) = NaN;
        else
            if X(i) > 0           % non-central
                x(i) = ncx2inv(P(i),d,m*X(i))/m;
            else
                x(i) = chi2inv(P(i),d)/m;
            end
        end
    end
end
