%% bfft
%
% "Block" Fast Fourier Transform (FFT)
%
% <matlab:open('bfft.m') code>
%
%% Syntax
%
%     AF = bfft(A,q)
%
%% Arguments
%
% See also <mvgchelp.html#4 Common variable names and data structures>.
%
% _input_
%
%     A          n1 x n2 x p matrix 
%     q          frequency resolution (default: p)
%
% _output_
%
%     AF         description
%
%% Description
%
% Returns the |q|-point discrete Fourier transform |AF| of each component
% vector of the |n1 x n2 x p| matrix |A|, where the last index is the vector
% index. If |q < p| a warning is issued that the FFT will truncate.
%
%% See also
%
% <bifft.html |bifft|> |
% <matlab:doc('fft') fft>
%
% (C) Lionel Barnett and Anil K. Seth, 2012. See file license.txt in
% installation directory for licensing terms.
%
%%

function AF = bfft(A,q)

[n1,n2,p] = size(A);

if nargin < 2 || isempty(q)
    q = p;
else
    if q < p
        fprintf(2,'WARNING: frequency resolution too low - fft will truncate!\n');
    end
end

AF = reshape(fft(reshape(A,n1*n2,p).',q).',n1,n2,q); % q-point Fourier transform on [0, 2*pi)
