%% mvdiff
%
% Multivariate differencing
%
% <matlab:open('mvdiff.m') code>
%
%% Syntax
%
%     Y = mvdiff(X,dff)
%
%% Arguments
%
% See also <mvgchelp.html#4 Common variable names and data structures>.
%
% _input_
%
%     X          multi-trial time series data
%     dff        vector of number of differencing steps (default: 1 = 1st difference)
%
% _output_
%
%     Y          differenced time series
%
%% Description
%
% Y = mvdiff(X,dff)
%
% Multivariate differencing of time series data in |X|, which may be
% multi-trial. Differences |dff| are supplied as a scalar (applied to each time
% series) or as a vector matching the number of variables. Default is 1st
% difference. The individual time series in the differenced multivariate series
% |Y| are correctly synchronised.
%
%% See also
%
% <matlab:doc('diff') |diff|>
%
% (C) Lionel Barnett and Anil K. Seth, 2012. See file license.txt in
% installation directory for licensing terms.
%
%%

function Y = mvdiff(X,dff)

if nargin < 2 || isempty(dff), dff = 1; end % default

[n,m,N] = size(X);

if isscalar(dff)
    dff = dff*ones(1,n);
else
    assert(isvector(dff) && length(dff) == n,'differences must be a scalar or a vector of same length as number of variables');
end
assert(all(isint(dff)) && all(dff >= 0),'differences must be non-negative integers');

Y = zeros(size(X));
for r = 1:N
    for i = 1:n
        d = dff(i);
        if d > 0
            Y(i,d+1:m,r) = diff(X(i,:,r),d,2);
        else
            Y(i,:,r) = X(i,:,r);
        end
    end
end
Y = Y(:,max(dff)+1:m,:);
