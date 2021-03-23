%% GCCA_tsdata_to_mvgc
%
% Calculate conditional time-domain MVGC (multivariate Granger causality) from time series data by
% "traditional" method (as e.g. in GCCA toolbox)
%
% <matlab:open('GCCA_tsdata_to_mvgc.m') code>
%
%% Syntax
%
%     F = GCCA_tsdata_to_mvgc(U,x,y,p,regmode)
%
%% Arguments
%
% See also <mvgchelp.html#4 Common variable names and data structures>.
%
% _input_
%
%     U          multi-trial time series data
%     x          vector of indices of target (causee) multi-variable
%     y          vector of indices of source (causal) multi-variable
%     p          model order (number of lags)
%     regmode    regression mode (default as in 'tsdata_to_var')
%
% _output_
%
%     F          Granger causality
%
%% Description
%
% Returns the time-domain MVGC
%
% <<eq_mvgc.png>>
%
% from the variable |Y| (specified by the vector of indices |y|) to the variable
% |X| (specified by the vector of indices |x|), conditional on all other
% variables |Z| in the time series data |U|. The regression mode is set by the
% |regmode| parameter, which may be |'LWR'| or |'OLS'| (see <tsdata_to_var.html
% |tsdata_to_var|> for details and defaults).
%
% If |p| is a vector of length 2, then |p(1)| is the number of lags for the full
% regression and |p(2)| the number of lags for the reduced regression.
% Otherwise, if |p| is a scalar, |p| lags are used for both the full and reduced
% regressions.
%
% This routine is included mainly for compatibility with the Granger Causal
% Connectivity Analysis (<http://www.sussex.ac.uk/Users/anils/aks_code.htm
% GCCA>) Toolbox [2]; the MVGC is calculated by performing separate full and reduced
% regressions - see [1] for details. Note that the preferred MVGC Toolbox method
% of calculating MVGC via the autocovariance sequence (see <autocov_to_mvgc.html
% |autocov_to_mvgc|>, <mvgc_demo.html |mvgc_demo|>) only requires a single
% regression step and is thus generally more accurate.
%
% The caller should take note of any warnings issued by this function and test
% results with a call <isbad.html |isbad|>|(F,false)|.
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
% <tsdata_to_var.html |tsdata_to_var|> |
% <autocov_to_mvgc.html |autocov_to_mvgc|> |
% <isbad.html |isbad|>
%
% (C) Lionel Barnett and Anil K. Seth, 2012. See file license.txt in
% installation directory for licensing terms.
%
%%

function F = GCCA_tsdata_to_mvgc(U,x,y,p,regmode)

if nargin < 5, regmode = []; end % ensure 'tsdata_to_var' default

if isscalar(p)
    pp = p;
else
    assert(isvector(p) && length(p) == 2,'number of lags must be a scalar or a vector of length 2');
    pp = p(2);
    p = p(1);
end

n = size(U,1);

x = x(:)'; % vectorise
y = y(:)'; % vectorise

assert(all(x >=1 & x <= n),'some x indices out of range');
assert(all(y >=1 & y <= n),'some y indices out of range');
assert(isempty(intersect(x,y)),'x and y indices must be distinct');

z = 1:n; z([x y]) = []; % indices of other variables (to condition out)
xz = [x z];
xzy = [xz y];

F = NaN;

% full regression

owstate = warn_supp;
[~,SIG] = tsdata_to_var(U(xzy,:,:),p,regmode);
warn_test(owstate,    'in full regression - data non-stationary or colinear?');
if warn_if(isbad(SIG),'in full regression - regression failed'), return; end % show-stopper!
% we should really check that SIG is positive-definite; we don't, for effciency reasons

% reduced regression

owstate = warn_supp;
[~,SIGR] = tsdata_to_var(U(xz,:,:),pp,regmode);
warn_test(owstate,     'in reduced regression - data non-stationary or colinear?');
if warn_if(isbad(SIGR),'in reduced regression - regression failed'), return; end % show-stopper!
% we should really check that SIGR is positive-definite; we don't, for effciency reasons

x = 1:length(x);
F = log(det(SIGR(x,x)))-log(det(SIG(x,x)));
