%% trfun2var
%
% Calculate coefficients from VAR transfer function
%
% <matlab:open('trfun2var.m') code>
%
%% Syntax
%
%     [A,p] = trfun2var(H,p)
%
%% Arguments
%
% See also <mvgchelp.html#4 Common variable names and data structures>.
%
% _input_
%
%     H          VAR transfer function matrix
%     p          VAR model order = number of lags (default: calculate automatically)
%
% _output_
%
%     A          VAR coefficients matrix
%     p          VAR model order = number of lags
%
%% Description
%
% Return coefficients |A| for VAR with transfer function |H|. If a VAR
% model order |p| is not specified (default), then it is set automatically
% according to the frequency resolution of the transfer function.
%
%% References
%
% [1] L. Barnett and A. K. Seth,
% <http://www.sciencedirect.com/science/article/pii/S0165027013003701 The MVGC
%     Multivariate Granger Causality Toolbox: A New Approach to Granger-causal
% Inference>, _J. Neurosci. Methods_ 223, 2014
% [ <matlab:open('mvgc_preprint.pdf') preprint> ].
%
%% See also
%
% <var2trfun.html |var2trfun|> |
% <cpsd_to_var.html |cpsd_to_var|>
%
% (C) Lionel Barnett and Anil K. Seth, 2012. See file license.txt in
% installation directory for licensing terms.
%
%%

function [A,p] = trfun2var(H,p)

[n,~,h] = size(H);
fres = h-1;
if nargin < 2 || isempty(p), p = 2*fres-1; end
assert(p < 2*fres,'too many lags');

I = eye(n);
AF = zeros(n,n,2*fres); % over [0,2*pi)
for k = 1:h             % over [0,pi]
    AF(:,:,k) = I/H(:,:,k);
end
for k = 1:h-2           % over (pi,2*pi)
    AF(:,:,h+k) = conj(AF(:,:,h-k));
end
A = real(bifft(AF));
A = -A(:,:,2:p+1);
