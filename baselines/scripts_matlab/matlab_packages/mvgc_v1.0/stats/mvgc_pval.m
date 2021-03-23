%% mvgc_pval
%
% p-values for sample MVGC based on theoretical asymptotic null distribution
%
% <matlab:open('mvgc_pval.m') code>
%
%% Syntax
%
%     pval = mvgc_pval(x,p,m,N,nx,ny,nz,tstat)
%
%% Arguments
%
% See also <mvgchelp.html#4 Common variable names and data structures>.
%
% _input_
%
%     x          matrix of MVGC values
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
%     pval       matrix of p-values
%
%% Description
%
% Returns p-values |pval| for sample MVGCs in |x|, based on theoretical
% (asymptotic) null distribution. |NaN| s are ignored. See <mvgc_cdf.html
% |mvgc_cdf|> for details of other parameters.
%
% *_Important:_* To test p-values for statistical significance you should
% correct for multiple null hypotheses; see routine <significance.html
% |significance|>.
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
% <mvgc_cdfi.html |mvgc_cdfi|> |
% <mvgc_confint.html |mvgc_confint|> |
% <mvgc_cval.html |mvgc_cval|> |
% <mvgc_demo.html |mvgc_demo|> |
% <significance.html |significance|>
%
% (C) Lionel Barnett and Anil K. Seth, 2012. See file license.txt in
% installation directory for licensing terms.
%
%%

function pval = mvgc_pval(x,p,m,N,nx,ny,nz,tstat)

if nargin < 7, nz    = []; end % ensure default
if nargin < 8, tstat = []; end % ensure default

pval = NaN(size(x)); % output p-value matrix is same shape as x matrix
nn   = ~isnan(x);    % indices of non-NaN x values (logical array)
x    = x(nn);        % vectorise non-NaN x values
pval(nn) = 1-mvgc_cdf(x,0,p,m,N,nx,ny,nz,tstat); % assume null hypothesis F = 0
