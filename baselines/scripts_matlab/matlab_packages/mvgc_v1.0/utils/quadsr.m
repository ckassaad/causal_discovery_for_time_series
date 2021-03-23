%% quadsr
%
% Trapezoidal rule quadrature (numerical integration) with sub-range specification
%
% <matlab:open('quadsr.m') code>
%
%% Syntax
%
%     [q,xrsupp] = quadsr(x,y,xrange)
%
%% Arguments
%
% _input_
%
%     x          vector of evaluation points
%     y          vector of values corresponding to x
%     xrange     x sub-range(s) specification vector
%
% _output_
%
%     q          quadrature of y over specified sub-range(s) of x values
%     xrsupp     the "support" (total length) of the  x sub-range(s)
%
%% Description
%
% A trapezoidal rule quadrature routine (_cf._ <matlab:doc('trapz')
% |trapz|>) where a series of sub-ranges of integration may be specified in
% |xrange|. |xrange| must be an ascending vector of even length, specifying
% start and end values of a series of sub-ranges. A |NaN| at the beginning
% (resp. end) of xrange specifies the beginning (resp. end) of the range of
% |x|. The integral (quadrature) of |y| over the sub-range(s) is returned
% in |q|. The "support" - i.e. total length of sub-range(s) - is returned
% in |xrsupp|.
%
%% See also
%
% <matlab:doc('trapz') |trapz|> |
% <quads.html |quads|> |
% <smvgc_to_mvgc.html |smvgc_to_mvgc|> 
%
% (C) Lionel Barnett and Anil K. Seth, 2012. See file license.txt in
% installation directory for licensing terms.
%
%%

function [q,xrsupp] = quadsr(x,y,xrange)

assert(isvector(x),'evaluation points must be a vector');
assert(isvector(y),'values must be a vector');
assert(isvector(xrange),'sub-ranges must be a vector');

x = x(:);
y = y(:);

r2 = length(xrange);
r = floor(r2/2);
assert(2*r == r2, 'sub-range vector must have an even number of points');

if isnan(xrange(1));   xrange(1)   = x(1);   end
if isnan(xrange(end)); xrange(end) = x(end); end

assert(all(xrange(2:end) > xrange(1:end-1)), 'sub-range vector must be ascending');

assert(xrange(1) >= x(1) && xrange(end) <= x(end), 'sub-ranges must lie within x range');

n = length(x);
idx = (1:n)';

q = 0;

xrsupp = 0;

for s=1:r

    xl = xrange(2*s-1);
    il = idx(x == xl);
    interpl = isempty(il);
    if interpl
        il = nnz(x < xl);
        dl = x(il+1)-x(il);
        yl = ((xl-x(il))/dl)*y(il+1)+((x(il+1)-xl)/dl)*y(il);
    end

    xu = xrange(2*s);
    iu = idx(x == xu);
    interpu = isempty(iu);
    if interpu
        iu = nnz(x < xu);
        du = x(iu+1)-x(iu);
        yu = ((xu-x(iu))/du)*y(iu+1)+((x(iu+1)-xu)/du)*y(iu);
    end

    if interpl
        if interpu
            xx = [xl; x(il+1:iu); xu];
            yy = [yl; y(il+1:iu); yu];
        else
            xx = [xl; x(il+1:iu)];
            yy = [yl; y(il+1:iu)];
        end
    else
        if interpu
            xx = [x(il:iu); xu];
            yy = [y(il:iu); yu];
        else
            xx = x(il:iu);
            yy = y(il:iu);
        end
    end

    q = q + quads(xx,yy);

    xrsupp = xrsupp + (xu-xl);

end
