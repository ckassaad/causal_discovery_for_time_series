%% empirical_pval
%
% p-values for sample statistics based on estimated empirical null distribution
%
% <matlab:open('empirical_pval.m') code>
%
%% Syntax
%
%     pval = empirical_pval(x,XNULL,ptails,ksmooth)
%
%% Arguments
%
% See also <mvgchelp.html#4 Common variable names and data structures>.
%
% _input_
%
%     x          matrix of statistic values
%     XNULL      matrix of null sample statistics
%     ptails     Pareto tails lower and upper probabilities (default: no Pareto tails)
%     ksmooth    use kernel smoothing to estimate cdf (default: no smoothing)
%
% _output_
%
%     pval       matrix of p-values
%
%% Description
%
% Returns p-values |pval| for empirical null distribution estimated from sample
% statistics in |XNULL| (derived e.g. from a permutation test), evaluated at
% values in |x|. The first dimension of |XNULL| must index samples, while the
% other dimensions must match the shape of |x|. |NaN| s are ignored. See
% <empirical_cdf.html |empirical_cdf|> for details of other parameters.
%
% *_Important:_* To test p-values for statistical significance you should
% correct for multiple null hypotheses; see routine <significance.html
% |significance|>.
%
%% See also
%
% <empirical_cdf.html |empirical_cdf|> |
% <empirical_cdfi.html |empirical_cdfi|> |
% <empirical_confint.html |empirical_confint|> |
% <empirical_cval.html |empirical_cval|> |
% <significance.html |significance|> |
% <tsdata_to_mvgc_permtest.html |tsdata_to_mvgc_permtest|> |
% <mvgc_demo_permtest.html |mvgc_demo_permtest|> |
% <significance.html |significance|>
%
% (C) Lionel Barnett and Anil K. Seth, 2012. See file license.txt in
% installation directory for licensing terms.
%
%%

function pval = empirical_pval(x,XNULL,ptails,ksmooth)

if nargin < 3, ptails  = []; end % force empirical_cdf default
if nargin < 4, ksmooth = []; end % force empirical_cdf default

s = size(XNULL);
if isvector(x) % including scalars
    goodargs = ndims(XNULL) == 2 && s(2) == length(x);
else
    goodargs = isequal(s(2:end),size(x));
end
assert(goodargs,'empirical null distributions matrix doesn''t match statistics matrix (first index must be samples)');

pval  = NaN(size(x)); % output p-value matrix is same shape as x matrix
nn    = ~isnan(x);    % indices of non-NaN x values (logical array)
x     = x(nn);        % vectorise non-NaN x values
XNULL = XNULL(:,nn);  % vectorise the corrsponding XNULL samples
n     = length(x);    % number of non-NaN x values

pv = zeros(n,1);
for i = 1:n
    pv(i) = 1-empirical_cdf(x(i),XNULL(:,i),ptails,ksmooth);
end
pval(nn) = pv;        % NaNs will be in same positions as they were in original x matrix
