function res = multinom(n,k)

% MULTINOM computes the multinomial coefficient.
% --------------------
% res = multinom(n, k)
% --------------------
% Description: computes the multinomial coefficient, defined as:
%              multinom(n,[k1, k2, ..., kl]) = n! / ( k1! k2! ... kl!).
%              The computation is based on the formula
%                   multinom(n,[k1,...,kl]) = binom(n,k1) *
%                             binom(n-k1,k2) * binom(n-k1-k2,k3) * ... *
%                             binom(n-k1-...-k_{l-1},kl).
% Input:       {n} scalar integer.
%              {k} vector of integers summing up to {n}.
% Output:      {res} the result.

% © Liran Carmel
% Classification: Combinatorics
% Last revision date: 19-Aug-2006

% check input argument
% error(nargchk(2,2,nargin));
% error(chkvar(n,'integer','scalar',{mfilename,inputname(1),1}));
% error(chkvar(k,'integer',{'vector',{'sumto',n}},{mfilename,inputname(2),2}));

% compute the first term in the multiplication
res = nchoosek(n,k(1));

% compute the rest of the terms
for ii = 2:length(k)
    n = n - k(ii-1);
    res = res * nchoosek(n,k(ii));
end