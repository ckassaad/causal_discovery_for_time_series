%% var_specrad
%
% Calculate VAR spectral radius
%
% <matlab:open('var_specrad.m') code>
%
%% Syntax
%
%     rho    = var_specrad(A)
%     [A,rho]= var_specrad(A,newrho)
%
%% Arguments
%
% See also <mvgchelp.html#4 Common variable names and data structures>.
%
% _input_
%
%     A          VAR coefficients matrix
%     newrho     new VAR spectral radius
%
% _output_
%
%     A          VAR coefficients matrix
%     rho        VAR spectral radius
%
%% Description
%
% _First form:_ return the spectral radius |rho| for a VAR process with
% coefficients matrix |A|. May be used for unit root test (i.e. need |rho < 1|).
%
% _Second form:_, a new value |newrho| for the spectral radius is supplied
% and the coefficients |A| are "decayed" (see function <var_decay.html
% |var_decay|>) so that the new value of the spectral radius becomes
% |newrho|. The new coefficients and old value of the spectral radius are
% returned. Note: |newrho| may be negative, but |newrho < 1| is required
% for the new coefficients to specify a stable process.
%
%% See also
%
% <var_decay.html |var_decay|>
%
% (C) Lionel Barnett and Anil K. Seth, 2012. See file license.txt in
% installation directory for licensing terms.
%
%%

function [out1,out2] = var_specrad(A,newrho)

[n,n1,p] = size(A);
assert(n1 == n,'VAR coefficients matrix has bad shape');
pn1 = (p-1)*n;

% construct VAR coefficients for 1-lag problem

A1 = [reshape(A,n,p*n); eye(pn1) zeros(pn1,n)];

% calculate spectral radius

rho = max(abs(eig(A1)));

if nargin < 2 || isempty(newrho)
    assert(nargout <= 1,'too many output parameters');
    out1 = rho;                     % spectral radius
else
    out1 = var_decay(A,newrho/rho); % adjusted coefficients
    out2 = rho;                     % previous value of spectral radius
end

