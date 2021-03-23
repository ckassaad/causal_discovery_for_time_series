%% significance
%
% Statistical significance adjusted for multiple hypotheses
%
% <matlab:open('significance.m') code>
%
%% Syntax
%
%     sig = significance(pval,alpha,correction)
%
%% Arguments
%
% See also <mvgchelp.html#4 Common variable names and data structures>.
%
% _input_
%
%     pval         multi-trial time series data
%     alpha        significance level
%     correction   multiple hypotheses correction (see Description)
%
% _output_
%
%     sig          significance (0 or 1)
%
%% Description
%
% Returns significance (0 or 1) of statistics based on p-values in |pval|,
% which may be a scalar, vector or matrix. NaNs are ignored. The
% |correction| parameter specifies a multiple hypotheses test adjustment,
% and may be one of: |'None'|, |'Bonferroni'|, |'Sidak'|,
% |'Holm'|, |'FDR'| (false discovery rate, independent or
% positively correlated hypotheses [1]) or |'FDRD'| (false discovery rate,
% arbitrary dependencies [2]).
%
% *_Note:_* |correction = 'None'| is not recommended for multiple
% hypotheses, so is _not_ the default! |'FDR'| is generally a good choice.
%
%% References
%
% [1] Y. Benjamini and Y. Hochberg, "Controlling the
% false discovery rate: a practical and powerful approach to multiple
% testing", _J. Royal Stat. Soc. B_, 57(1), 1995.
%
% [2] Y. Benjamini and D. Yekutieli, "The control of the false discovery
% rate in multiple testing under dependency", _Ann. Stat_, 29(4), 2001.
%
%% See also
%
% <mvgc_pval.html |mvgc_pval|> |
% <empirical_pval.html |empirical_pval|>
%
% (C) Lionel Barnett and Anil K. Seth, 2012. See file license.txt in
% installation directory for licensing terms.
%
%%

function sig = significance(pval,alpha,correction)

sig = NaN(size(pval)); % same shape as p-value array
nn  = ~isnan(pval);    % non-NaNs (logical array)
p   = pval(nn);        % vectorise the non-NaN p-values
n   = numel(p);        % number of p-values being tested

switch upper(correction)

    case 'NONE';

        signn = (p < alpha);

    case 'BONFERRONI'

        signn = (p < alpha/n);

    case 'SIDAK'

        signn = (p < 1-realpow(1-alpha,1/n));

    case 'HOLM' % assumes independence of hypothesis tests

        signn = false(1,n);
        [psorted,sortidx] = sort(p);
        for v=1:n
            if psorted(v) < alpha/(n-v+1)
                signn(sortidx(v)) = true;
            else
                break; % remaining null hypotheses accepted
            end
        end

    case 'FDR'   % assumes independence (or positive correlation) of hypothesis tests (more powerful)

        signn  = fdr_bh(p,alpha,true);

    case 'FDRD' %  possible dependencies - no correlation assumptions

        signn  = fdr_bh(p,alpha,false);

    otherwise; error('unknown correction method');
end

sig(nn) = signn; % NaNs will be in same positions as they were in original p-value array

function [h crit_p]=fdr_bh(pvals,q,pdep)

% fdr_bh() - Executes the Benjamini & Hochberg (1995) procedure for
%            controlling the false discovery rate (FDR) of a family of
%            hypothesis tests. FDR is the expected proportion of rejected
%            hypotheses that are mistakenly rejected (i.e., the null
%            hypothesis is actually true for those tests). FDR is a
%            somewhat less conservative/more powerful method for correcting
%            for multiple comparisons than methods like Bonferroni
%            correction that provide strong control of the family-wise
%            error rate (i.e., the probability that one or more null
%            hypotheses are mistakenly rejected).
%
% Usage:
%  >> [h, crit_p]=fdr_bh(pvals,q,method,report);
%
% Required Input:
%   pvals - A vector or matrix (two dimensions or more) containing the
%           p-value of each individual test in a family of tests.
%
% Optional Inputs:
%   q       - The desired false discovery rate. {default: 0.05}
%   method  - ['pdep' or 'dep'] If 'pdep,' the original Bejnamini & Hochberg
%             FDR procedure is used, which is guaranteed to be accurate if
%             the individual tests are independent or positively dependent
%             (e.g., positively correlated).  If 'dep,' the FDR procedure
%             described in Benjamini & Yekutieli (2001) that is guaranteed
%             to be accurate for any test dependency structure (e.g.,
%             positively and/or negatively correlated tests) is used. 'dep'
%             is always appropriate to use but is less powerful than 'pdep.'
%             {default: 'pdep'}
%   report  - ['yes' or 'no'] If 'yes', a brief summary of FDR results are
%             output to the MATLAB command line {default: 'no'}
%
%
% Outputs:
%   h       - A binary vector or matrix of the same size as the input "pvals."
%             If the ith element of h is 1, then the test that produced the
%             ith p-value in pvals is significant (i.e., the null hypothesis
%             of the test is rejected).
%   crit_p  - All p-values less than or equal to crit_p are significant
%             (i.e., their null hypotheses are rejected).  If no p-values are
%             significant, crit_p=0.
%
%
% References:
%   Benjamini, Y. & Hochberg, Y. (1995) Controlling the false discovery
%     rate: A practical and powerful approach to multiple testing. Journal
%     of the Royal Statistical Society, Series B (Methodological). 57(1),
%     289-300.
%
%   Benjamini, Y. & Yekutieli, D. (2001) The control of the false discovery
%     rate in multiple testing under dependency. The Annals of Statistics.
%     29(4), 1165-1188.
%
% Example:
%   [dummy p_null]=ttest(randn(12,15)); %15 tests where the null hypothesis
%                                       %is true
%   [dummy p_effect]=ttest(randn(12,5)+1); %5 tests where the null
%                                          %hypothesis is false
%   [h crit_p]=fdr_bh([p_null p_effect],.05,'pdep','yes');
%
%
% Author:
% David M. Groppe
% Kutaslab
% Dept. of Cognitive Science
% University of California, San Diego
% March 24, 2010

s=size(pvals);
if length(s) > 2 || s(1) > 1
    p_sorted = sort(reshape(pvals,1,prod(s)));
else % p-values are already a row vector
    p_sorted = sort(pvals);
end
m = length(p_sorted); % number of tests

if pdep % BH procedure for independence or positive dependence
    thresh = (1:m)*q/m;
else    % BH procedure for any dependency structure
    thresh = (1:m)*q/(m*sum(1./(1:m)));
end

rej = p_sorted <= thresh;
max_id = find(rej,1,'last'); % find greatest significant pvalue
if isempty(max_id),
    crit_p = 0;
    h = logical(pvals*0);
else
    crit_p = p_sorted(max_id);
    h = pvals <= crit_p;
end
