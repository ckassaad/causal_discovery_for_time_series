function r = randi(imax,n)

%fprintf(2,'legacy imax!\n');
assert(nargin > 0,'too few arguments');
if     nargin == 1
    r = 1+floor(imax*rand);
elseif nargin == 2
    r = 1+floor(imax*rand(n));
else
    error('please use array form for matrix dimension argument');
end
