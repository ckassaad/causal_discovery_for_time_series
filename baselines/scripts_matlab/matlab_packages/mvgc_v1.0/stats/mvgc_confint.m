%% mvgc_confint
%
% Confidence intervals for sample MVGC based on theoretical asymptotic distribution
%
% <matlab:open('mvgc_confint.m') code>
%
%% Syntax
%
%     [xup,xlo] = mvgc_confint(alpha,x,p,m,N,nx,ny,nz,tstat)
%
%% Arguments
%
% See also <mvgchelp.html#4 Common variable names and data structures>.
%
% _input_
%
%     alpha      significance level (scalar)
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
%     xup        matrix of upper confidence bounds
%     xlo        matrix of lower confidence bounds
%
%% Description
%
% Return upper and lower confidence bounds |[xup,xlo]| at significance level
% |alpha| for sample MVGCs in |x|, based on theoretical (asymptotic)
% distributions. |NaN| s are ignored. See <mvgc_cdfi.html |mvgc_cdfi|> for
% details of other parameters.
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
% <mvgc_cval.html |mvgc_cval|> |
% <mvgc_demo_confint.html |mvgc_demo_confint|>
%
% (C) Lionel Barnett and Anil K. Seth, 2012. See file license.txt in
% installation directory for licensing terms.
%
%%

function [xup,xlo] = mvgc_confint(alpha,x,p,m,N,nx,ny,nz,tstat)

if nargin < 8, nz    = []; end % ensure default
if nargin < 9, tstat = []; end % ensure default

assert(isscalar(alpha),'alpha must be a scalar');

s = size(x);
xup = NaN(s);   % output xup matrix is same shape as x matrix
xlo = NaN(s);   % output xlo matrix is same shape as x matrix
nn = ~isnan(x); % indices of non-NaN x values (logical array)
x  = x(nn);     % vectorise non-NaN x values

a = alpha*ones(size(x));
xup(nn) = mvgc_cdfi(1-a,x,p,m,N,nx,ny,nz,tstat);
xlo(nn) = mvgc_cdfi(a,  x,p,m,N,nx,ny,nz,tstat);
