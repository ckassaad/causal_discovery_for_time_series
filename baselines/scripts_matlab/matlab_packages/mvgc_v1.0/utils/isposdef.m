%% isposdef
%
% Determine whether (symmetric) matrix is positive-definite
%
% <matlab:open('isposdef.m') code>
%
%% Syntax
%
%     pd = isposdef(A)
%
%% Arguments
%
% _input_
%
%     A          a symmetric matrix
%
% _output_
%
%     pd         logical true if A is positive-definite
%
%% Description
%
% Returns true if A is positive-definite. Uses the Matlab <matlab:doc('chol')
% |chol|> function, and like that function the lower triangle of A is assumed to
% be the (complex conjugate) transpose of the upper triangle; i.e. this routine
% does _not_ check for symmetry!
%
%% See also
%
% <matlab:doc('chol') |chol|>
%
% (C) Lionel Barnett and Anil K. Seth, 2012. See file license.txt in
% installation directory for licensing terms.
%
%%

function pd = isposdef(A)

[~,p] = chol(A);
pd = ~(p > 0);
