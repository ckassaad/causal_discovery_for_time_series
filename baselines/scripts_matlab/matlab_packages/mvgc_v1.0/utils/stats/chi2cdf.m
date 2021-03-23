function P = chi2cdf(X,V)

assert((isequal(size(X),size(V))),'X and V must be of common size or scalars');

P = gamcdf(X,V/2,2);
