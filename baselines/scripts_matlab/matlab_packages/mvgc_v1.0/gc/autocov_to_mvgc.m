%% autocov_to_mvgc
%
% Calculate conditional time-domain MVGC (multivariate Granger causality)
%
% <matlab:open('autocov_to_mvgc.m') code>
%
%% Syntax
%
%     F = autocov_to_mvgc(G,x,y)
%
%% Arguments
%
% See also <mvgchelp.html#4 Common variable names and data structures>.
%
% _input_
%
%     G          autocovariance sequence
%     x          vector of indices of target (causee) multi-variable
%     y          vector of indices of source (causal) multi-variable
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
% from the variable |Y| (specified by the vector of indices |y|) to the
% variable |X| (specified by the vector of indices |x|), conditional on all
% other variables |Z| represented in |G|, for a stationary VAR process with
% autocovariance sequence |G|. See ref. [1] for details.
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
%% See also
%
% <autocov_to_var.html |autocov_to_var|> |
% <isbad.html |isbad|>
%
% (C) Lionel Barnett and Anil K. Seth, 2012. See file license.txt in
% installation directory for licensing terms.
%
%%

function F = autocov_to_mvgc(G,x,y)

n = size(G,1);

x = x(:)'; % vectorise
y = y(:)'; % vectorise

assert(all(x >=1 & x <= n),'some x indices out of range');
assert(all(y >=1 & y <= n),'some y indices out of range');
assert(isempty(intersect(x,y)),'x and y indices must be distinct');

z   = 1:n; z([x y]) = []; % indices of other variables (to condition out)
xz  = [x z];
xzy = [xz y];

F = NaN;

% full regression

owstate = warn_supp;
[~,SIG ] = autocov_to_var(G(xzy,xzy,:));
warn_test(owstate,    'in full regression - bad autocovariance matrix? Check output of ''var_info''');
if warn_if(isbad(SIG),'in full regression - regression failed'), return; end % show-stopper!

% reduced regression

owstate = warn_supp;
[~,SIGR] = autocov_to_var(G(xz,xz,:));    % reduced regression
warn_test(owstate,     'in reduced regression - bad autocovariance matrix? Check output of ''var_info''');
if warn_if(isbad(SIGR),'in reduced regression - regression failed'), return; end % show-stopper!

x = 1:length(x);
F = log(det(SIGR(x,x)))-log(det(SIG(x,x)));
