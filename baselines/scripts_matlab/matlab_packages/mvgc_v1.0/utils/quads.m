%% quads
%
% Trapezoidal rule quadrature (numerical integration)
%
% <matlab:open('quads.m') code>
%
%% Syntax
%
%     q = quads(x,Y)
%
%% Arguments
%
% _input_
%
%     x          vector of evaluation points
%     Y          matrix of values corresponding to x
%
% _output_
%
%     q          vector of quadratures
%
%% Description
%
% A simple trapezoidal rule quadrature routine (_cf._ <matlab:doc('trapz')
% |trapz|>). Integrate each column of matrix |Y| assumed evaluated at
% points in vector |x|.
%
%
%% See also
%
% <matlab:doc('trapz') |trapz|> |
% <quadsr.html |quadsr|> |
% <smvgc_to_mvgc.html |smvgc_to_mvgc|> 
%
% (C) Lionel Barnett and Anil K. Seth, 2012. See file license.txt in
% installation directory for licensing terms.
%
%%

function q = quads(x,Y)

assert(isvector(x),'evaluation points must be a vector');
assert(size(Y,1) == length(x),'values and evaluation points must have same length');

x = x(:); % force column vector
xx = x(:,ones(1,size(Y,2)));
q  = sum((xx(2:end,:)-xx(1:end-1,:)).*((Y(2:end,:)+Y(1:end-1,:))))/2;
