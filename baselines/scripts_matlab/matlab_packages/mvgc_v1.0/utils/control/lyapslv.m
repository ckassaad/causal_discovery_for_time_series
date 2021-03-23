function X = lyapslv(~,A,~,Q)

% called from 'dlyap' as: X = lyapslv('D',A,[],Q);

n = size(A,1);
assert(size(A,2) == n,'matrix A not square');
assert(isequal(size(Q),size(A)),'matrix Q does not match matrix A');

% Schur factorisation

[U,T] = schur(A);
Q = -U'*Q*U;

% solve the equation column ~by column

X = zeros(n);
j = n;
while j > 0
    j1 = j;

    % check Schur block size

    if j == 1
        bsiz = 1;
    elseif T(j,j-1) ~= 0
        bsiz = 2;
        j = j-1;
    else
        bsiz = 1;
    end
    bsizn = bsiz*n;

    Ajj = kron(T(j:j1,j:j1),T)-eye(bsizn);

    rhs = reshape(Q(:,j:j1),bsizn,1);

    if j1 < n
        rhs = rhs + reshape(T*(X(:,(j1+1):n)*T(j:j1,(j1+1):n)'),bsizn,1);
    end

    v = -Ajj\rhs;
    X(:,j) = v(1:n);

    if bsiz == 2
        X(:,j1) = v((n+1):bsizn);
    end

    j = j-1;
end

% transform back to original coordinates

X = U*X*U';
