function fVal = istcBound(w, qlzi, qlxi, loglik)
% The negative lower bound function value

T = size(qlzi,3);
posqi = cat(1, qlzi, qlxi);
% marginal posterior of q

fvq = 0;

for t = 1:T
        posqii = posqi(:,:,t);
        fvq = fvq - sum(log(w(:)).*posqii(:));
end

fvH = 0;

for t = 1:T
    posqii = posqi(:,:,t);
    posqii = posqii(:);
    posqii = posqii(posqii~=0);
    fvH = fvH + sum(posqii.*log(posqii));
end

fVal = fvq - loglik + fvH;
end
