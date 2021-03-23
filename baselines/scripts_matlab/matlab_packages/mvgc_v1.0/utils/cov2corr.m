%% cov2corr
%
% Convert (auto)covariance to (auto)correlation
%
% <matlab:open('cov2corr.m') code>
%
%% Syntax
%
%     R = cov2corr(G)
%
%% Arguments
%
% See also <mvgchelp.html#4 Common variable names and data structures>.
%
% _input_
%
%     G          a covariance matrix or autocovariance sequence
%
% _output_
%
%     R          a correlation matrix or autocorrelation sequence
%
%% Description
%
% Calculates (auto)correlation |R| from (auto)covariance |G|.
%
% (C) Lionel Barnett and Anil K. Seth, 2012. See file license.txt in
% installation directory for licensing terms.
%
%%

function R = cov2corr(G)

assert(ndims(G) == 2 || ndims(G) == 3,'(auto)covariance matrix must have 2 or 3 dimensions');
[n,n1,q1] = size(G);
assert(n1 == n,'(auto)covariance matrix has bad shape');

R = zeros(size(G));
D = diag(1./sqrt(diag(G(:,:,1))));
for k = 1:q1
    R(:,:,k) = D*G(:,:,k)*D;
end


