%% empirical_confint
%
% Confidence intervals for sample statistics based on estimated empirical null distribution
%
% <matlab:open('empirical_confint.m') code>
%
%% Syntax
%
%     [xup,xlo] = empirical_confint(alpha,X,ptails,ksmooth)
%
%% Arguments
%
% See also <mvgchelp.html#4 Common variable names and data structures>.
%
% _input_
%
%     alpha      significance level (scalar)
%     X          matrix of sample statistics
%     ptails     Pareto tails lower and upper probabilities (default: no Pareto tails)
%     ksmooth    use kernel smoothing to estimate cdf (default: no smoothing)
%
% _output_
%
%     xup        matrix of upper confidence bounds
%     xlo        matrix of lower confidence bounds
%
%% Description
%
% Return upper and lower confidence bounds |[xup,xlo]| at significance level
% |alpha| based on the empirical distributions in |X| (derived e.g. from a
% bootstrap). The first dimension of |X| must index samples. |NaN| s are ignored.
% See <empirical_cdfi.html |empirical_cdfi|> for details of other parameters.
%
%% See also
%
%
% <empirical_cdf.html |empirical_cdf|> |
% <empirical_cdfi.html |empirical_cdfi|> |
% <empirical_pval.html |empirical_pval|> |
% <empirical_cval.html |empirical_cval|>
% <mvgc_demo_confint.html |mvgc_demo_confint|>
%
% (C) Lionel Barnett and Anil K. Seth, 2012. See file license.txt in
% installation directory for licensing terms.
%
%%

function [xup,xlo] = empirical_confint(alpha,X,ptails,ksmooth)

if nargin < 3, ptails  = []; end % force empirical_cdfi default
if nargin < 4, ksmooth = []; end % force empirical_cdfi default

assert(isscalar(alpha),'alpha must be a scalar');

nn = squeeze(~any(isnan(X),1)); % indices of non-NaNs (logical array)
s = size(nn);
n = nnz(nn);
xup = NaN(s);  % upper bounds matrix
xlo = NaN(s);  % lower bounds matrix
X   = X(:,nn); % vectorise non-NaN X values

xu = zeros(n,1);
xl = zeros(n,1);
for i = 1:n
    xu(i) = empirical_cdfi(1-alpha,X(:,i),ptails,ksmooth);
    xl(i) = empirical_cdfi(alpha,  X(:,i),ptails,ksmooth);
end
xup(nn) = xu; % upper bounds matrix with NaNs in the right place
xlo(nn) = xl; % lower bounds matrix with NaNs in the right place
