%% cpsd_to_autocov
%
% Calculate autocovariance sequence from cross-power spectral density
%
% <matlab:open('cpsd_to_autocov.m') code>
%
%% Syntax
%
%     [G,q] = cpsd_to_autocov(S,q)
%
%% Arguments
%
% See also <mvgchelp.html#4 Common variable names and data structures>.
%
% _input_
%
%     S          cross-power spectral density (cpsd) matrix
%     q          number of autocovariance lags to calculate (default: automatic)
%
% _output_
%
%     G          autocovariance sequence
%     q          number of autocovariance lags actually calculated
%
%% Description
%
% Calculates the autocovariance sequence |G| defined as [[ii_acseq.png]] to |q|
% lags from the cross-power spectral density (cpsd) |S|. This is essentially an
% inverse Fourier transform
%
% <<eq_cpsd2ac.png>>
%
% implemented as an (discrete) inverse fast Fourier transform. If a number of
% lags |q| is not supplied, then the default is to set it to the frequency
% resolution of the cpsd. The actual number of lags calculated is returned in
% |q|.
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
% <autocov_to_cpsd.html |autocov_to_cpsd|>
%
% (C) Lionel Barnett and Anil K. Seth, 2012. See file license.txt in
% installation directory for licensing terms.
%
%%

function [G,q] = cpsd_to_autocov(S,q)

[n,~,h] = size(S);
fres = h-1;
if nargin < 2 || isempty(q), q = fres-1; end
assert(q < 2*fres,'too many lags');
q1 = q+1;

%G = bifft(cat(3,flipdim(permute(S(:,:,2:fres+1),[2 1 3]),3),S(:,:,1:fres)),2*fres); % inverse transform of "circular shifted" spectral density
G = bifft(cat(3,flipdim(conj(S(:,:,2:fres+1)),3),S(:,:,1:fres)),2*fres); % inverse transform of "circular shifted" spectral density

r = ones(1,ceil(q1/2));
sgn = [r; -r];
sgn = sgn(1:q1); % sgn = [1 -1 1 -1 ...]

G = real(reshape(sgn(ones(1,n*n),:).*reshape(G(:,:,1:q1),n*n,q1),n,n,q1));

% note 1: G truncated to q1 lags
% note 2: G should be real, but ...
