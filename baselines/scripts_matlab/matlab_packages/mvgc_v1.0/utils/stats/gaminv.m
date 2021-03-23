function X = gaminv(P,A,B)

assert((isequal(size(A),size(P)) && isequal(size(B),size(P))),'P, A and B must be of common size or scalars');

sz = size(P);
X = zeros(sz);

k = find((P < 0) | (P > 1) | isnan(P) | ~(A > 0) | ~(B > 0));
if (any(k))
    X(k) = NaN;
end

k = find((P == 1) & (A > 0) & (B > 0));
if any(k)
    X(k) = Inf;
end

k = find((P > 0) & (P < 1) & (A > 0) & (B > 0));
if any(k)
    if ~isscalar(A) || ~isscalar(B)
        A = A(k);
        B = B(k);
        y = A.*B;
    else
        y = A*B*ones(size(k));
    end
    P = P(k);

    if isa(P, 'single')
        myeps = eps('single');
    else
        myeps = eps;
    end

    l = find(P < myeps);
    if any(l)
        y(l) = sqrt(myeps)*ones(length(l), 1);
    end

    y_old = y;
    for i = 1 : 100
        h     = (gamcdf(y_old,A,B)-P)./gampdf(y_old,A,B);
        y_new = y_old-h;
        ind   = find(y_new <= myeps);
        if any(ind)
            y_new(ind) = y_old(ind)/10;
            h = y_old-y_new;
        end
        if max(abs(h)) < sqrt(myeps)
            break;
        end
        y_old = y_new;
    end

    X(k) = y_new;
end
