%% mvgc_cval
%
% Critical values for sample MVGC based on theoretical asymptotic null distribution
%
% <matlab:open('mvgc_cval.m') code>
%
%% Syntax
%
%     x = mvgc_cval(alpha,p,m,N,nx,ny,nz,tstat)
%
%% Arguments
%
% See also <mvgchelp.html#4 Common variable names and data structures>.
%
% _input_
%
%     alpha      vector of significance levels
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
%     x          vector of critical MVGC values
%
%% Description
%
% Returns critical MVGC values |x| at significance levels in |alpha| for sample
% MVGC, based on theoretical (asymptotic) null distribution. See <mvgc_cdfi.html
% |mvgc_cdfi|> for details of other parameters.
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
% <mvgc_pval.html |mvgc_pval|> |
% <mvgc_confint.html |mvgc_confint|> |
% <mvgc_demo_confint.html |mvgc_demo_confint|> |
% <mvgc_demo_nonstationary.html |mvgc_demo_nonstationary|>
%
% (C) Lionel Barnett and Anil K. Seth, 2012. See file license.txt in
% installation directory for licensing terms.
%
%%

function x = mvgc_cval(alpha,p,m,N,nx,ny,nz,tstat)

if nargin < 7, nz    = []; end % ensure default
if nargin < 8, tstat = []; end % ensure default

assert(isvector(alpha),'significance levels must be a scalar or vector');

x = mvgc_cdfi(1-alpha,0,p,m,N,nx,ny,nz,tstat); % assume null hypothesis F = 0
