%% demean
%
% Temporal demean of time series data
%
% <matlab:open('demean.m') code>
%
%% Syntax
%
%     Y = demean(X,normalise)
%
%% Arguments
%
% See also <mvgchelp.html#4 Common variable names and data structures>.
%
% _input_
%
%     X          multi-trial time series data
%     normalise  normalise (temporal) variance of each variable to 1 (default: false)
%
% _output_
%
%     Y          demeaned time series data
%
%% Description
%
% Temporally demean time series data |X|, which may be single- or
% multi-trial. If the |normalise| flag is set, data is normalised so that
% the (temporal) variance of each series is 1.
%
% *_Note:_* For multi-trial data we don't demean on a "per-trial" basis,
% since this really doesn't make sense... trials in multi-trial data are
% assumed to be multiple realisations of the same process. In particular,
% demeaning trials separately can introduce large bias in VAR model
% estimation (_cf._ <tsdata_to_var.html |tsdata_to_var|>). The mean
% calculated here is thus the temporal mean of ensemble means. If you feel
% you absolutely have to demean per-trial, then call this function for each
% trial series |X(:,:,r)| _and then_ call it with |X|.
%
%% See also
%
% <tsdata_to_var.html |tsdata_to_var|>
%
% (C) Lionel Barnett and Anil K. Seth, 2012. See file license.txt in
% installation directory for licensing terms.
%
%%

function Y = demean(X,normalise)

if nargin < 2 || isempty(normalise), normalise = false; end

[n,m,N] = size(X);

U = ones(1,N*m);
Y = X(:,:);
Y = Y-mean(Y,2)*U;
if normalise
    Y = Y./(std(Y,[],2)*U);
end
Y = reshape(Y,n,m,N);
