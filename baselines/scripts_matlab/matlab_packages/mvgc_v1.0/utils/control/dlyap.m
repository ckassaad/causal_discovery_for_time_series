%% dlyap
%
% Solve discrete-time Lyapunov equation by Schur decomposition
%
% <matlab:open('dlyap.m') code>
%
%% Syntax
%
%     X = dlyap(A,Q)
%
%% Arguments
%
% See also <mvgchelp.html#4 Common variable names and data structures>.
%
% _input_
%
%     A          square matrix with spectral radius < 1
%     Q          symmetric positive-definite matrix
%
% _output_
%
%     X          solution of discrete-time Lyapunov equation X = A*X*A'+Q
%
%% Description
%
% Solves the discrete-time Lyapunov equation
%
% <<eq_dlyap.png>>
%
% using Schur decomposition and column-by-column solution [1,2] (this is quite
% likely the same algorithm as in the Matlab Control System Toolbox
% <matlab:doc('dlyap') dlyap>, except that no balancing is performed). This
% function is substantially slower than <matlab:doc('dlyap') dlyap>; if the
% Matlab Control System Toolbox is not available, the user might like to
% consider using <dlyap_aitr.html |dlyap_aitr|> instead.
%
% The caller should test that |Q| is positive-definite and that the spectral
% radius of |A| is |< 1| before calling.
%
% Adapted from an early scripted (i.e. non- <http://www.slicot.org/ SLICOT>)
% version of <http://www.gnu.org/software/octave/ Octave>'s
% <http://octave.sourceforge.net/control/function/dlyap.html |dlyap|> function.
%
%% References
%
% [1] X. Kitagawa, "An algorithm for solving the matrix equation X = F X F'+S",
%     _Internat. J. Control_ 25(5), 1977.
%
% [2] S. Hammarling, "Numerical solution of the stable, non-negative definite
%    Lyapunov equation", _IMA J. Numer. Anal._ 2, 1982.
%
%% See also
%
% <dlyap_aitr.html |dlyap_aitr|> |
% <var_to_autocov.html |var_to_autocov|>
%
% (C) Lionel Barnett and Anil K. Seth, 2012. See file license.txt in
% installation directory for licensing terms.
%
%%

function X = dlyap(A,Q)

X = lyapslv('D',A,[],-Q);
