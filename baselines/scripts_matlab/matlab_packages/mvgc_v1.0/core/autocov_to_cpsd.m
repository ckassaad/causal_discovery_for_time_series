%% autocov_to_cpsd
%
% Calculate cross-power spectral density from autocovariance sequence
%
% <matlab:open('autocov_to_cpsd.m') code>
%
%% Syntax
%
%     [S,fres] = autocov_to_cpsd(G,fres)
%
%% Arguments
%
% See also <mvgchelp.html#4 Common variable names and data structures>.
%
% _input_
%
%     G          autocovariance sequence
%     fres       frequency resolution to calculate (default: automatic)
%
% _output_
%
%     S          cross-power spectral density (cpsd) matrix
%     fres       frequency resolution actually calculated
%
%% Description
%
% Calculates the cross-power spectral density (cpsd) |S| at frequency resolution
% |fres| from the autocovariance sequence |G| defined as [[ii_acseq.png]]. This
% is essentially a Fourier transform:
%
% <<eq_ac2cpsd.png>>
%
% implemented as a discrete Fast Fourier Transform (FFT). If a frequency resolution
% |fres| is not supplied, then the default is to set it to the number of
% autocovariance lags. The actual frequency resolution calculated is returned in
% |fres|. Call |freqs = <sfreqs.html sfreqs>(fres,fs)|, where |fs| is the
% sampling rate, to get a corresponding vector |freqs| of frequencies on
% |[0,fs/2]|.
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
% <cpsd_to_autocov.html |cpsd_to_autocov|> |
% <sfreqs.html |sfreqs|>
%
% (C) Lionel Barnett and Anil K. Seth, 2012. See file license.txt in
% installation directory for licensing terms.
%
%%

function [S,fres] = autocov_to_cpsd(G,fres)

if nargin < 2 || isempty(fres), fres = size(G,3); end
h = fres+1;

G0 = G(:,:,1);

S = bfft(G,2*fres); % over [0,2*pi)
S = S(:,:,1:h);     % over [0,pi]
S = S + conj(permute(S,[2 1 3])) - G0(:,:,ones(1,1,h)); % last term is repmat(G0,[1 1 h])
