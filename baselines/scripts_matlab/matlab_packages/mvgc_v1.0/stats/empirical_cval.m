%% empirical_cval
%
% Critical values for sample statistics based on estimated empirical null distribution
%
% <matlab:open('empirical_cval.m') code>
%
%% Syntax
%
%     x = empirical_cval(alpha,XNULL,ptails,ksmooth)
%
%% Arguments
%
% See also <mvgchelp.html#4 Common variable names and data structures>.
%
% _input_
%
%     alpha      vector of significance levels
%     XNULL      matrix of null sample statistics
%     ptails     Pareto tails lower and upper probabilities (default: no Pareto tails)
%     ksmooth    use kernel smoothing to estimate cdf (default: no smoothing)
%
% _output_
%
%     x          matrix of critical values
%
%% Description
%
% Return critical values |x| at significance levels in |alpha| based on
% empirical null distributions in |XNULL| (derived e.g. from a permutation
% test). The first dimension of |XNULL| must index samples. The leading indices
% of |x| index null distributions, the last significance levels. |NaN| s are
% ignored. See <empirical_cdfi.html |empirical_cdfi|> for details of other
% parameters.
%
%% See also
%
% <empirical_cdf.html |empirical_cdf|> |
% <empirical_cdfi.html |empirical_cdfi|> |
% <empirical_pval.html |empirical_pval|> |
% <empirical_confint.html |empirical_confint|>
%
% (C) Lionel Barnett and Anil K. Seth, 2012. See file license.txt in
% installation directory for licensing terms.
%
%%

function x = empirical_cval(alpha,XNULL,ptails,ksmooth)

if nargin < 3, ptails  = []; end % force empirical_cdfi default
if nargin < 4, ksmooth = []; end % force empirical_cdfi default

assert(~isscalar(XNULL),'null distributions must be a vector or matrix');
assert(isvector(alpha),'significance levels must be a scalar or vector');

s = size(XNULL);
sx = s(2:end);
n = prod(sx);
XNULL = reshape(XNULL,s(1),n);
r = length(alpha);
x = NaN(n,r);
for i = 1:n
    if ~any(isnan(XNULL(:,i)));
        x(i,:) = empirical_cdfi(1-alpha,XNULL(:,i),ptails,ksmooth);
    end
end
x = reshape(x,[sx r]);
