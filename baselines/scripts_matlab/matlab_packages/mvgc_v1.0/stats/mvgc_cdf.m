%% mvgc_cdf
%
% Sample MVGC thoretical asymptotic cumulative distribution function
%
% <matlab:open('mvgc_cdf.m') code>
%
%% Syntax
%
%     P = mvgc_cdf(x,X,p,m,N,nx,ny,nz,tstat)
%
%% Arguments
%
% See also <mvgchelp.html#4 Common variable names and data structures>.
%
% _input_
%
%     x          vector of MVGC values
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
%     P          cumulative distribution probabilities evaluated at x
%
%% Description
%
% Return theoretical sample MVGC asymptotic cumulative distribution function for
% actual MVGCs in vector |X|, evaluated at values in vector |x|. For p-values,
% assume null hypothesis [H_0]: |X = 0| (see <mvgc_pval.html |mvgc_pval|>).
%
% The theoretical distribution is specified by the |tstat| parameter, which may
% be  'F' for Granger's F-test (default if |nx == 1|) or 'chi2' for Geweke's
% [chi^2] test (default if |nx > 1|). For a multivariate predictee (i.e. |nx >
% 1|) only the [chi^2] test is suitable [1].
%
% *_Note 1:_* In theory the F-distribution should be preferable for small
% samples; it has a fatter tail than the corresponding [chi^2] distribution.
% However, for a multivariate predictee (|nx > 1|) it is not appropriate; the
% usual F-test for a nested linear regression demands the Granger form statistic
% |F = (RSS_reduced - RSS_full) / RSS_full|. For a univariate predictee this may
% be derived from a simple transformation of the Geweke form |F =
% log(det(SIGMA_reduced)) / log(det(SIGMA_full))|. However for a multivariate
% predictee the Granger form cannot be derived from the Geweke form; in fact it
% corresponds to the "trace" MVGC variant recommended (for different reasons) in
% [2]; but there are problems with the trace form, in particular lack
% of invariance to certain transformations of variables and filters [2,3] and
% lack of an information-theoretic interpretation [1].
%
% *_Note 2:_* The non-central F-distribution is currently to be treated as
% experimental and a warning is issued if it is invoked. Specifically, the
% problem is: should it be doubly- or singly-noncentral and what should the
% non-centrality parameter(s) be? Currently we assume singly-noncentral and make
% an informed guess at the non-centrality parameter. This appears to work well,
% giving results similar to the [chi^2] test for medium-sized samples, but
% warrants further investigation.
%
%% References
%
% [1] L. Barnett and A. K. Seth,
% <http://www.sciencedirect.com/science/article/pii/S0165027013003701 The MVGC
%     Multivariate Granger Causality Toolbox: A New Approach to Granger-causal
% Inference>, _J. Neurosci. Methods_ 223, 2014
% [ <matlab:open('mvgc_preprint.pdf') preprint> ].
%
% [2] C. Ladroue, S. Guo, K. Kendrick and J. Feng, "Beyond element-wise
% interactions: Identifying complex interactions in biological processes",
% _PLoS ONE_ 4(9), 2009.
%
% [2] A. B. Barrett, L. Barnett and A. K. Seth, "Multivariate Granger causality
% and generalized variance", _Phys. Rev. E_ 81(4), 2010.
%
% [3] L. Barnett and A. K. Seth, "Behaviour of Granger causality under
% filtering: Theoretical invariance and practical application", _J. Neurosci.
% Methods_ 201(2), 2011.
%
%% See also
%
% <mvgc_cdfi.html |mvgc_cdfi|> |
% <mvgc_pval.html |mvgc_pval|> |
% <mvgc_confint.html |mvgc_confint|> |
% <mvgc_cval.html |mvgc_cval|> |
% <mvgc_demo.html |mvgc_demo|>
%
% (C) Lionel Barnett and Anil K. Seth, 2012. See file license.txt in
% installation directory for licensing terms.
%
%%

function P = mvgc_cdf(x,X,p,m,N,nx,ny,nz,tstat)

assert(isvector(x),'evaluation values must be a vector');
n = length(x);

assert(isvector(X),'MVGC values must be a vector');
assert(all(X >= 0),'MVGC values must be non-negative');
if isscalar(X)
    X = X*ones(n,1);
else
    assert(length(X) == n,'MVGC values must match evaluation values');
end

if nargin < 8 || isempty(nz), nz = 0; end % unconditional

if nargin < 9 || isempty(tstat);
    ftest = nx == 1; % default: use F-distribution for univariate predictee, chi2 for multivariate
else
    switch lower(tstat)
        case 'f'     % Granger F-test form
            assert(nx == 1,'F-distribution is not appropriate for multivariate predictee');
            ftest = true;
        case 'chi2'  % Geweke chi2 test form
            ftest = false;
        otherwise
            error('unknown distribution (must be ''chi2'' or ''F'')');
    end
end

P = zeros(n,1);
m = N*(m-p);                  % effective number of observations (p-lag autoregression loses p observations per trial)
if ftest
    if any(X > 0), fprintf(2,'WARNING (mvgc_cdf): non-central F-distribution is experimental\n'); end
    d1 = p*ny;                % #{full model parameters} - #{reduced model parameters}
    d2 = m-p*(1+ny+nz);       % #{observations} - #{full model parameters}
    mm = d2/d1;
    for i = 1:n
        xx = exp(x(i))-1;     % Granger form: (RSS_reduced - RSS_full) / RSS_full
        if X(i) > 0           % non-central
            XX = exp(X(i))-1; % Granger form
            P(i) = ncfcdf(mm*xx,d1,d2,m*XX); % NOTE: non-centrality parameter factor might reasonably be m, d2 or d2-2
        else
            P(i) = fcdf(mm*xx,d1,d2);
        end
    end
else
    d = p*nx*ny;              % note that d does not depend on the number of conditioning variables
    for i = 1:n
        if X(i) > 0           % non-central
            P(i) = ncx2cdf(m*x(i),d,m*X(i));
        else
            P(i) = chi2cdf(m*x(i),d);
        end
    end
end
