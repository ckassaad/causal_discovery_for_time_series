%% empirical_cdf
%
% Empirical cumulative distribution function
%
% <matlab:open('empirical_cdf.m') code>
%
%% Syntax
%
%     P = empirical_cdf(x,X,ptails,ksmooth)
%
%% Arguments
%
% See also <mvgchelp.html#4 Common variable names and data structures>.
%
% _input_
%
%     x          vector of statistic values
%     X          vector of sample statistics
%     ptails     Pareto tails lower and upper probabilities (default: no Pareto tails)
%     ksmooth    use kernel smoothing to estimate cdf (default: no smoothing)
%
% _output_
%
%     P          cumulative distribution probabilities evaluated at x
%
%% Description
%
% Returns the empirical cumulative distribution function estimated from sample
% statistics in vector |X|, evaluated at values in vector |x|. For p-values, |X|
% should be drawn from the appropriate null distribution (see
% <empirical_pval.html |empirical_pval|>).
%
% Uses the Matlab Statistics Toolbox <matlab:doc('paretotails') |paretotails|>
% class to estimate the cdf. Pareto tails for the cumulative distribution may be
% specified by the 2-vector |ptails|, with lower and upper cumulative
% probabilities in |ptails(1)| and |ptails(2)| respectively [where |0 <=
% ptails(1) < ptails(2) <= 1|]. If the |ksmooth| flag is set, then kernel
% smoothing is used to estimate the cdf, otherwise simple midpoint
% interpolation; see <matlab:doc('paretotails') |paretotails|> for details.
%
%% See also
%
% <empirical_cdfi.html |empirical_cdfi|> |
% <empirical_pval.html |empirical_pval|> |
% <empirical_confint.html |empirical_confint|> |
% <empirical_cval.html |empirical_cval|> |
% <mvgc_demo_permtest.html |mvgc_demo_permtest|> |
% <matlab:doc('paretotails') |paretotails|>
%
% (C) Lionel Barnett and Anil K. Seth, 2012. See file license.txt in
% installation directory for licensing terms.
%
%%

function P = empirical_cdf(x,X,ptails,ksmooth)

if nargin < 3 || isempty(ptails),  ptails  = [0,1]; end % default: no Pareto tails
if nargin < 4 || isempty(ksmooth), ksmooth = false; end % default: don't use kernel smoothing

if ksmooth, cdffun = 'kernel'; else cdffun = 'ecdf'; end

assert(isvector(X),'samples must be a vector');
assert(isvector(x),'evaluation values must be a vector');
assert(isvector(ptails) && length(ptails) == 2 ...
    && ptails(1) >= 0 && ptails(1) <= 1 && ptails(2) >= 0 && ptails(2) <= 1 && ptails(1) < ptails(2), ...
    'Pareto tails must be an ascending 2-vector in [0,1]');

D = paretotails(X,ptails(1),ptails(2),cdffun);
P = cdf(D,x);

if any(isnan(P))
    n = length(x);
    xmin = min(X);
    xmax = max(X);
    for i = 1:n
        if isnan(P(i))
            if     x(i) <= xmin
                P(i) = 0;
            elseif x(i) >= xmax
                P(i) = 1;
            else
                fprintf(2,'WARNING: empirical_cdf failed for value %d = %f\n',i,x(i));
            end
        end
    end
end
