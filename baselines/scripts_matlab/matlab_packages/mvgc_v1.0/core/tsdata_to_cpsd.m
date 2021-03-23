%% tsdata_to_cpsd
%
% Estimate cross-power spectral density from time series data
%
% <matlab:open('tsdata_to_cpsd.m') code>
%
%% Syntax
%
%     S = tsdata_to_cpsd(X,fres,method,window,noverlap,nw,ntapers)
%
%% Arguments
%
% See also <mvgchelp.html#4 Common variable names and data structures>.
%
% _input_
%
%     X          multi-trial time series data
%     fres       frequency resolution
%     method     estimation method: 'WELCH' (default) or 'MT'
%     window     window length (default: min time series trial length and 2*fres)
%     noverlap   window overlap size (default: window/2)
%     nw         [multi-taper only] multi-taper time bandwidth parameter (default: 3)
%     ntapers    [multi-taper only] number of tapers (default: 2*nw-1)
%
% _output_
%
%     S          estimated cross-power spectral density (cpsd) matrix
%
%% Description
%
% Returns an estimate of the cross-power spectral density of a multivariate
% process from time series data |X|, which may be single- or multi-trial. |fres|
% specifies the frequency resolution; call |freqs = <sfreqs.html
% sfreqs>(fres,fs)|, where |fs| is the sampling rate, to get a corresponding
% vector |freqs| of frequencies on |[0,fs/2]|.
%
% NOTE: This function requires the Matlab(R) Signal Processing Toolbox(TM).
%
% The |window| and |noverlap| parameters specify the window length and window
% overlap size, respectively. The estimation |method| may be |'WELCH'| for
% Welch's averaged modified periodogram method (default: see functions
% <matlab:doc('pwelch') |pwelch|>, <matlab:doc('cpsd') |cpsd|> in the Matlab(R)
% Signal Processing Toolbox(TM)) or |'MT'| for a multi-taper estimation
% procedure. For the multi-taper method, |nw| is the multi-taper time bandwidth
% parameter and |ntaper| the number of tapers. Discrete prolate spheroidal
% (Slepian) sequences for the multi-taper method are calculated using the
% function <matlab:doc('dpss') |dpss|> from the Matlab(R) Signal Processing
% Toolbox(TM).
%
% Multi-taper routine adapted from <http://www.chronux.org/ |Chronux|> code [2]
% authored by <http://mitralab.org/ Partha Mitra> and Kenneth D. Harris.
%
%% References
%
% [1] L. Barnett and A. K. Seth,
% <http://www.sciencedirect.com/science/article/pii/S0165027013003701 The MVGC
%     Multivariate Granger Causality Toolbox: A New Approach to Granger-causal
% Inference>, _J. Neurosci. Methods_ 223, 2014
% [ <matlab:open('mvgc_preprint.pdf') preprint> ].
%
% [2] P. Mitra and H. Bokil, "Observed Brain Dynamics", _Oxford University Press_, New York, 2008. 
%
%% See also
%
% <matlab:doc('pwelch') |pwelch|> |
% <matlab:doc('cpsd') |cpsd|> |
% <matlab:doc('dpss') |dpss|>
%
% (C) Lionel Barnett and Anil K. Seth, 2012. See file license.txt in
% installation directory for licensing terms.
%
%%

function S = tsdata_to_cpsd(X,fres,method,window,noverlap,nw,ntapers)

[n,m,N] = size(X);
X = demean(X);
X = permute(X,[2 1 3]); % transpose row, col (works for single-trial data too)

nfft = 2*fres;

if nargin < 3 || isempty(method)
    method = 'WELCH'; % default is WELCH
   %method = 'MT';    % default is multi-taper
end

if nargin < 4 || isempty(window)
    window = min(m,nfft); % according to Chronux ... by default Matlab 'pwelch' splits data into 8 overlapping segments
end
assert(window <= m,'window cannot be longer than data');

if nargin < 5 || isempty(noverlap)
    noverlap = round(window/2);
end
assert(noverlap < window,'overlap must be shorter than window');

if strcmpi(method,'MT')

    nchunks = floor(((m-noverlap)/(window-noverlap))); % FFT chunks per channel

    if nargin < 6 || isempty(nw)
        nw = 3;
    end

    if nargin < 7 || isempty(ntapers)
        ntapers = 2*nw -1;
    end

    tapers   = dpss(window,nw,ntapers,'calc'); % Slepian sequences: tapers is a matrix of size window x ntapers
    taparray = tapers(:,:,ones(1,n));

    S = 0;
    for r = 1:N % works for single-trial too
        S = S + cpsd_mt(X(:,:,r),n,fres+1,window,noverlap,nchunks,taparray);
    end
    S = permute(S,[2 3 1])/N;

elseif strcmpi(method,'WELCH')

    S = 0;
    for r = 1:N % works for single-trial too
        S = S + cpsd_welch(X(:,:,r),n,fres+1,window,noverlap);
    end
    S = pi*S/N; % the 'pi' is for compatibility with 'autocov_to_cpsd' routine

else
    error('unknown method ''%s''',method);
end

% now fill other half of cpsd matrix with complex conjugates

for i = 1:n
    for j = i+1:n
        S(j,i,:) = conj(S(i,j,:));
    end
end

function S = cpsd_mt(X,n,h,window,noverlap,nchunks,taparray)

nfft = 2*(h-1);

S = complex(zeros(h,n,n)); % output array

winstep = window-noverlap;
ntapers = size(taparray,2);

% compute tapered periodogram with FFT

for k = 1:nchunks

    XSEG = X((1:window) + (k-1)*winstep,:);

    % compute periodogram

    P = fft(taparray.*permute(XSEG(:,:,ones(1,ntapers)),[1 3 2]),nfft);
    P = P(1:h,:,:);

    % now make cross-products of them to fill cross-spectrum matrix

    for i = 1:n
        for j = i:n % so we don't compute cross-spectra twice
            S(:,i,j) = S(:,i,j) + mean(P(:,:,i) .* conj(P(:,:,j)),2);
        end
    end

end

S = S/nchunks;

function S = cpsd_welch(X,n,h,window,noverlap)

nfft = 2*(h-1);

S = complex(zeros(n,n,h));

for i = 1:n
    S(i,i,:) = pwelch(X(:,i),window,noverlap,nfft);          % auto-spectra
    for j = i+1:n % so we don't compute cross-spectra twice
        S(i,j,:) = cpsd(X(:,i),X(:,j),window,noverlap,nfft); % cross-spectra
    end
end
