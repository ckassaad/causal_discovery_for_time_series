%% empirical_cdfi
%
% Empirical inverse cumulative distribution function
%
% <matlab:open('empirical_cdfi.m') code>
%
%% Syntax
%
%     x = empirical_cdfi(P,X,ptails,ksmooth)
%
%% Arguments
%
% See also <mvgchelp.html#4 Common variable names and data structures>.
%
% _input_
%
%     P          vector of cumulative distribution probabilities
%     X          vector of sample statistics
%     ptails     Pareto tails lower and upper probabilities (default: no Pareto tails)
%     ksmooth    use kernel smoothing to estimate cdf (default: no smoothing)
%
% _output_
%
%     x          vector of statistic values
%
%% Description
%
% Returns the empirical inverse cumulative distribution function estimated from
% sample statistics in vector |X|, evaluated at probabilities in the vector |P|.
% To calculate the critical value for a significance level |alpha|, |X| should be
% drawn from the appropriate null distribution and |P| set to |1-alpha| (see
% <empirical_cval.html |empirical_cval|>). For confidence intervals at level
% |alpha|, |X| should be drawn from the appropriate sample distribution and |P|
% set to |alpha| for the lower bound and |1-alpha| for the upper bound (see
% <empirical_confint.html |empirical_confint|>).
%
% Uses the Matlab Statistics Toolbox <matlab:doc('paretotails') |paretotails|>
% class to estimate the inverse cdf. Pareto tails for the corresponding
% cumulative distribution may be specified by the 2-vector |ptails|, with lower
% and upper cumulative probabilities in |ptails(1)| and |ptails(2)| respectively
% [where |0 <= ptails(1) < ptails(2) <= 1|]. If the |ksmooth| flag is set, then
% kernel smoothing is used to estimate the cdf, otherwise simple midpoint
% interpolation; see <matlab:doc('paretotails') |paretotails|> for details.
%
%% See also
%
% <empirical_cdf.html |empirical_cdf|> |
% <empirical_pval.html |empirical_pval|> |
% <empirical_confint.html |empirical_confint|> |
% <empirical_cval.html |empirical_cval|> |
% <matlab:doc('paretotails') |paretotails|>
%
% (C) Lionel Barnett and Anil K. Seth, 2012. See file license.txt in
% installation directory for licensing terms.
%
%%

function x = empirical_cdfi(P,X,ptails,ksmooth)

if nargin < 3 || isempty(ptails),  ptails  = [0,1]; end % default: no Pareto tails
if nargin < 4 || isempty(ksmooth), ksmooth = false; end % default: don't use kernel smoothing

if ksmooth, cdffun = 'kernel'; else cdffun = 'ecdf'; end

assert(isvector(X),'samples must be a vector');
assert(isvector(P),'probabilities must be a vector');
assert(isvector(ptails) && length(ptails) == 2 ...
    && ptails(1) >= 0 && ptails(1) <= 1 && ptails(2) >= 0 && ptails(2) <= 1 && ptails(1) < ptails(2), ...
    'Pareto tails must be an ascending 2-vector in [0,1]');

if any(P < 0 || P > 1)
    n = length(P);
    for i = 1:n
        if P(i) < 0, P(i) = 0; elseif P(i) > 1, P(i) = 1; end
    end
end

D = paretotails(X,ptails(1),ptails(2),cdffun);
x = icdf(D,P);
