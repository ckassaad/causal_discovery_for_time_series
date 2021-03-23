function Y = gampdf(X,A,B)

assert((isequal(size(A),size(X)) && isequal(size(B),size(X))),'X, A and B must be of common size or scalars');

sz = size(X);
Y = zeros(sz);

k = find(~(A > 0) | ~(B > 0) | isnan(X));
if any(k)
    Y(k) = NaN;
end

k = find((X > 0) & (A > 0) & (A <= 1) & (B > 0));
if any(k)
    if isscalar(A) && isscalar(B)
        Y(k) = (X(k).^(A - 1)).*exp(-X(k)./B)./gamma(A)./(B .^ A);
    else
        Y(k) = (X(k).^(A(k)-1)).*exp(-X(k)./B(k))./gamma(A(k))./(B(k).^A(k));
    end
end

k = find((X > 0) & (A > 1) & (B > 0));
if any(k)
    if isscalar(A) && isscalar(B)
        Y(k) = exp(-A.*log(B)+(A-1).*log(X(k))-X(k)./B-gammaln(A));
    else
        Y(k) = exp(-A(k).*log(B(k))+(A(k)-1).*log(X(k))-X(k)./B(k)-gammaln(A(k)));
    end
end
