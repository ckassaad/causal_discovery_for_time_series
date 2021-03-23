%% GCCA_tsdata_to_pwcgc
%
% Calculate pairwise-conditional time-domain MVGCs (multivariate Granger
% causality) from time series data by "traditional" method (as e.g. in GCCA
% toolbox)
%
% <matlab:open('GCCA_tsdata_to_pwcgc.m') code>
%
%% Syntax
%
%     [F,A,SIG,E] = GCCA_tsdata_to_pwcgc(X,p,regmode)
%
%% Arguments
%
% See also <mvgchelp.html#4 Common variable names and data structures>.
%
% _input_
%
%     X          multi-trial time series data
%     x          vector of indices of target (causee) multi-variable
%     y          vector of indices of source (causal) multi-variable
%     p          model order (number of lags)
%     regmode    regression mode (default as in 'tsdata_to_var')
%
% _output_
%
%     F          Granger causality matrix
%     A          VAR coefficients matrix
%     SIG        residuals covariance matrix
%     E          residuals time series
%
%% Description
%
% Returns the  matrix |F| of pairwise-conditional time-domain MVGCs
%
% <<eq_mvgc_pwc.png>>
%
% (where |[ij]| denotes omission of the |ij|-th variables) between all
% pairs of variables [[ii_inej.png]] represented in |G|, for a stationary VAR
% process with autocovariance sequence |G|. Note that the first index |i|
% of |F| is the target (causee) variable, the second |j| the source
% (causal) variable. See ref. [1] for details. The regression mode is set
% by the |regmode| parameter, which may be |'LWR'| or |'OLS'| (see
% <tsdata_to_var.html |tsdata_to_var|> for details and defaults). The VAR
% parameters |A,SIG| and residuals |E| of the full regression are
% (optionally) returned.
%
% If |p| is a scalar, |p| lags are used for both the full and reduced
% regressions. Otherwise, if |p| is a vector of length 2, then |p(1)| is the
% number of lags for the full regression and |p(2)| the number of lags for the
% reduced regression (*_note_*: this is generally a bad idea! - see the
% accompanying documentation [1] for reasons).
%
% The caller should take note of any warnings issued by this function and test
% results with a call <isbad.html |isbad|>|(F,false)|.
% 
% This routine is included mainly for compatibility with the Granger Causal
% Connectivity Analysis (<http://www.sussex.ac.uk/Users/anils/aks_code.htm
% GCCA>) Toolbox [2]; the MVGCs are calculated by performing separate full
% and reduced regressions - see [1] for details and <mvgc_demo_GCCA.html
% |mvgc_demo_GCCA|> for example usage. Note that the preferred MVGC Toolbox
% method of calculating MVGCs via the autocovariance sequence (see
% <autocov_to_pwcgc.html |autocov_to_pwcgc|>, <mvgc_demo.html
% |mvgc_demo|>) only requires a single regression step and are thus
% generally more accurate.
%
%% References
%
% [1] L. Barnett and A. K. Seth,
% <http://www.sciencedirect.com/science/article/pii/S0165027013003701 The MVGC
%     Multivariate Granger Causality Toolbox: A New Approach to Granger-causal
% Inference>, _J. Neurosci. Methods_ 223, 2014
% [ <matlab:open('mvgc_preprint.pdf') preprint> ].
%
% [2] A. K. Seth, "A MATLAB toolbox for Granger causal connectivity
% analysis", _Journal of Neuroscience Methods_ 186, 2010.
%
%% See also
%
% <GCCA_autocov_to_mvgc.html |GCCA_autocov_to_mvgc|> |
% <autocov_to_pwcgc.html |autocov_to_pwcgc|> |
% <tsdata_to_var.html |tsdata_to_var|> |
% <mvgc_demo_GCCA.html |mvgc_demo_GCCA|> |
% <isbad.html |isbad|>
%
% (C) Lionel Barnett and Anil K. Seth, 2012. See file license.txt in
% installation directory for licensing terms.
%
%%

function [F,A,SIG,E] = GCCA_tsdata_to_pwcgc(X,p,regmode)

if nargin < 3, regmode = []; end % ensure 'tsdata_to_var' default

if isscalar(p)
    pp = p;
else
    assert(isvector(p) && length(p) == 2,'number of lags must be a scalar or a vector of length 2');
    pp = p(2);
    p = p(1);
end

n = size(X,1);
F = nan(n);

% full regression

owstate = warn_supp;
if nargout > 1
    if nargout > 3
        [A,SIG,E] = tsdata_to_var(X,p,regmode);
    else
        [A,SIG]   = tsdata_to_var(X,p,regmode);
    end
else
    [~,SIG]   = tsdata_to_var(X,p,regmode);
end
warn_test(owstate,    'in full regression - data non-stationary or colinear?');
if warn_if(isbad(SIG),'in full regression - regression failed'), return; end % show-stopper!
% we should really check that SIG is positive-definite; we don't, for effciency reasons
LSIG = log(diag(SIG)); % residuals log variances

for j = 1:n;

    % reduced regression

    jo  = [1:j-1 j+1:n]; % omit j
    
    owstate = warn_supp;
    [~,SIGj] = tsdata_to_var(X(jo,:,:),pp,regmode);
    warn_test(owstate,     sprintf('in reduced regression for target node %d - data non-stationary or colinear?',j));
    if warn_if(isbad(SIGj),sprintf('in reduced regression for target node %d - regression failed',j)), continue; end % show-stopper!
    % we should really check that SIGj is positive-definite; we don't, for effciency reasons
    LSIGj = log(diag(SIGj)); % residuals log variances

    % conditional Granger causalities

    for ii=1:n-1;
        i = jo(ii);
        F(i,j) = LSIGj(ii)-LSIG(i);
    end
end
