%% dlyap_aitr
%
% Solve discrete-time Lyapunov equation by Smith's accelerated iterative method
%
% <matlab:open('dlyap_itr.m') code>
%
%% Syntax
%
%     [X,iters] = dlyap_itr(A,Q,maxiters,maxrelerr)
%
%% Arguments
%
% See also <mvgchelp.html#4 Common variable names and data structures>.
%
% _input_
%
%     A          square matrix with spectral radius < 1
%     Q          symmetric positive-definite matrix
%     maxiters   maximum iterations (default: 100)
%     maxrelerr  maximum relative error (default: 1e-8)
%
% _output_
%
%     X          solution of discrete-time Lyapunov equation X = A*X*A'+Q
%     iters      number of iterations performed
%
%% Description
%
% Solves the discrete-time Lyapunov equation (_cf._ <matlab:doc('dlyap')
% dlyap> in the Matlab Control System Toolbox)
%
% <<eq_dlyap.png>>
%
% using Smith's accelerated iterative method [1]:
%
% <<eq_aitr.png>>
%
% |maxiters| is the maximum iterations, |maxrelerr| the maximum relative error
% of [[ii_dlyapre.png]] with respect to |Q| before a convergence test kicks in.
% The actual number of iterations performed is returned in |iters|. If the
% algorithm fails to converge within |maxiters| iterations, an exception
% |MVGC:XMaxItrs| is thrown (this is for consistency with <matlab:doc('dlyap')
% dlyap>). The convergence rate depends on the problem dimensions and maximum
% absolute eigenvalue (spectral radius) of |A|. In general this algorithm
% actually appears to be faster than the Matlab Control System Toolbox
% <matlab:doc('dlyap') dlyap>.
%
% The caller should test that |Q| is positive-definite and that the spectral
% radius of |A| is |< 1| before calling.
%
%% References
%
% [1] R. A. Smith, "Matrix Equation XA + BX = C", _SIAM J. Appl. Math._ 16(1), 1968.
%
%% See also
%
% <matlab:doc('dlyap') dlyap> |
% <dlyap_schur.html |dlyap_schur|> |
% <var_to_autocov.html |var_to_autocov|>
%
% (C) Lionel Barnett and Anil K. Seth, 2012. See file license.txt in
% installation directory for licensing terms.
%
%%

function [X,iters] = dlyap_aitr(A,Q,maxiters,maxrelerr)

if nargin < 3 || isempty(maxiters),  maxiters  = 100;  end
if nargin < 4 || isempty(maxrelerr), maxrelerr = 1e-8; end

assert(size(A,2) == size(A,1),'matrix A not square');
assert(isequal(size(Q),size(A)),'matrix Q does not match matrix A');

X  = Q;
AA = A;
snorm = norm(Q,'fro');
minrelerr = realmax;
for iters = 1:maxiters+1
    relerr = norm(X-A*X*A'-Q,'fro')/snorm;
    if relerr < maxrelerr                  % only start convergence test after max rel error threshold reached
        if relerr >= minrelerr, break; end % deemed converged
    end
    if relerr < minrelerr, minrelerr = relerr; end
    X = AA*X*AA'+X;
    AA = AA*AA;
end

if iters > maxiters
    throw(MException('MVGC:XMaxItrs','exceeded maximum iterations (max. rel. error = %e)',relerr));
end
