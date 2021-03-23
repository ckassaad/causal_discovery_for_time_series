function P = gamcdf(X,A,B)

assert((isequal(size(A),size(X)) && isequal(size(B),size(X))),'X, A and B must be of common size or scalars');

sz = size(X);
P = zeros(sz);

k = find(~(A > 0) | ~(B > 0) | isnan(X));
if any(k)
    P(k) = NaN;
end

k = find((X > 0) & (A > 0) & (B > 0));
if any(k)
    if isscalar(A) && isscalar(B)
        P(k) = gammainc(X(k)./B,A);
    else
        P(k) = gammainc(X(k)./B(k),A(k));
    end
end
