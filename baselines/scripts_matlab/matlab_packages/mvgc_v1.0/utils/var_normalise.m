%% var_normalise
%
% Normalise VAR coefficients
%
% <matlab:open('var_normalise.m') code>
%
%% Syntax
%
%     [A,SIG] = var_normalise(A,SIG)
%
%% Arguments
%
% See also <mvgchelp.html#4 Common variable names and data structures>.
%
% _input_
%
%     A          VAR coefficients matrix
%     SIG        residuals covariance matrix
%
% _output_
%
%     A          normalised VAR coefficients matrix
%     SIG        normalised residuals covariance matrix
%
%% Description
%
% Normalise VAR coefficients matrix |A| and residuals covariance matrix |SIG| so
% that residuals have unit variance (are "Studentised").
%
% (C) Lionel Barnett and Anil K. Seth, 2012. See file license.txt in
% installation directory for licensing terms.
%
%%

function [A,SIG] = var_normalise(A,SIG)

[n,n1,p] = size(A);
assert(n1 == n,'VAR coefficients matrix has bad shape');
[nn1,nn2] = size(SIG);
assert(nn1 == nn2,'residuals covariance matrix not square');
assert(nn1 == n  ,'residuals covariance matrix doesn''t match VAR coefficients matrix');

S = sqrt(diag(SIG));
V = diag(S);
VI = diag(1./S);

for k = 1:p
    A(:,:,k) = VI*A(:,:,k)*V;
end
SIG = VI*SIG*VI';
