%% mvdetrend
%
% Multivariate polynomial detrend of time series data
%
% <matlab:open('mvdetrend.m') code>
%
%% Syntax
%
%     [Y,P,p,x] = mvdetrend(X,pdeg,x)
%
%% Arguments
%
% See also <mvgchelp.html#4 Common variable names and data structures>.
%
% _input_
%
%     X          multi-trial time series data
%     pdeg       vector of polynomial degrees (default: 1 = linear detrend)
%     x          vector of polynomial evaluation points (default: evenly spaced)
%
% _output_
%
%     Y          detrended time series
%     P          polynomial fit
%     p          polynomial coefficients (cell array)
%     x          normalised evaluation points
%
%% Description
%
% Multivariate polynomial detrend of time series data in |X|, which may be
% multi-trial. Adapted from Matlab <matlab:doc('polyfit') |polyfit|> and
% <matlab:doc('polyval') |polyval|>. Polynomial degrees |pdeg| are supplied as a
% scalar (applied to each time series) or as a vector matching the number of
% variables; the default is a linear detrend. A vector |x| of points at which the
% polynomials are evaluated may be supplied; the default is evenly spaced (i.e.
% at each time step).
%
% The de-trended series are returned in |Y|, the polynomial fits in |P|, the
% polynomial coefficients in the cell array |p| and the normalised evaluation
% points in |x|. *_Note:_* This routine effectively demeans as well.
%
%% See also
%
% <matlab:doc('polyfit') |polyfit|> |
% <matlab:doc('polyval') |polyval|>
%
% (C) Lionel Barnett and Anil K. Seth, 2012. See file license.txt in
% installation directory for licensing terms.
%
%%

function [Y,P,p,x] = mvdetrend(X,pdeg,x)

if nargin < 2 || isempty(pdeg), pdeg = 1; end % default is linear detrend (note that pdeg = 0 simply demeans)

[n,m,N] = size(X);

if nargin < 3 || isempty(x), x = 1:m; end % evenly spaced

assert(isvector(x) && length(x) == m,'evaluation points must be a vector of same length as time series');
x = x(:)'; % ensure row vector

if isscalar(pdeg)
    pdeg = pdeg*ones(1,n);
else
    assert(isvector(pdeg) && length(pdeg) == n,'polynomial degrees must be a scalar or a vector of same length as number of variables');
end
assert(all(isint(pdeg)) && all(pdeg >= 0),'polynomial degrees must be non-negative integers');

mu  = mean(x,2);
sig = std(x,[],2);
x = (x - mu)/sig; % normalise

P = zeros(n,m);

p = cell(n,1);
for i = 1:n
    d = pdeg(i);
    d1 = d+1;
    V = zeros(d1,m);  % Vandermonde matrix
    V(d1,:) = ones(1,m);
    for j = d:-1:1
        V(j,:) = x.*V(j+1,:);
    end
    p{i} = mean(X(i,:,:),3)/V; % mean over trials
    P(i,:) = p{i}(1)*ones(1,m);
    for j = 2:d1
        P(i,:) = x.*P(i,:) + p{i}(j);
    end
end

Y = zeros(size(X));
for r = 1:N
    Y(:,:,r) = X(:,:,r)-P;
end
