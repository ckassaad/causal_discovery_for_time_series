%% sfreqs
%
% Return vector of frequencies
%
% <matlab:open('sfreqs.m') code>
%
%% Syntax
%
%     freqs = sfreqs(fres,fs)
%
%% Arguments
%
% See also <mvgchelp.html#4 Common variable names and data structures>.
%
% _input_
%
%     fres       frequency resolution
%     fs         sample rate (default: 2*pi)
%
% _output_
%
%     freqs      vector of frequencies
%
%% Description
%
% Returns (column) vector |freqs| of |fres+1| equally spaced frequencies on
% |[0,fs/2]|, where |fs| is a sample rate (so |fs/2| is the Nyqvist frequency).
% If a sample rate is not supplied (default), frequencies are returned over the
% normalised range |[0 pi]|.
%
%% See also
%
% All spectral routines
%
% (C) Lionel Barnett and Anil K. Seth, 2012. See file license.txt in
% installation directory for licensing terms.
%
%%

function freqs = sfreqs(fres,fs)

if nargin < 2 || isempty(fs);
    fs = 2*pi; % normalised
end

freqs = linspace(0,fs/2,fres+1)';
