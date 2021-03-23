%% isbad
%
% Determine whether array  is "bad"
%
% <matlab:open('isbad.m') code>
%
%% Syntax
%
%     b = isbad(x)
%
%% Arguments
%
% _input_
%
%     x                 an array
%     demand_allfinite  true (default) if "bad" means at least one NaN or Inf
%
% _output_
%
%     b                 logical scalar
%
%% Description
%
% Simple routine that returns true/false according as to whether the input |x|
% is "bad", meaning (i) if |demand_allfinite| flag set (default) it contains at
% least one |NaN| or |Inf| or (ii) if |demand_allfinite| flag unset it contains
% all |NaN| s or |Inf| s.  This is useful to check results of calculations that
% result in |NaN| s or |Inf| s when they shouldn't (in particular
% <tsdata_to_var.html |tsdata_to_var|>).
%
%% See also
%
% <tsdata_to_var.html |tsdata_to_var|> |
% <matlab:doc('isfinite') |isfinite|>
%
% (C) Lionel Barnett and Anil K. Seth, 2012. See file license.txt in
% installation directory for licensing terms.
%
%%

function b = isbad(x,demand_allfinite)

if nargin < 2 || isempty(demand_allfinite), demand_allfinite = true; end

if demand_allfinite
    b = ~all(isfinite(x(:)));
else
    b = ~any(isfinite(x(:)));
end
