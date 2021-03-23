function X = chi2inv(P,V)

assert((isequal(size(P),size(V))),'P and V must be of common size or scalars');

X = gaminv(P,V/2,2);
