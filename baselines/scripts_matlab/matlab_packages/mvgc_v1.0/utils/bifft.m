%% bifft
%
% "Block" Inverse Fast Fourier Transform (IFFT)
%
% <matlab:open('bifft.m') code>
%
%% Syntax
%
%     AF = bifft(A,q)
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
% Returns the |q|-point discrete inverse Fourier transform |AF| of each
% component vector of the |n1 x n2 x p| matrix |A|, where the last index is the
% vector index. If |q < p| a warning is issued that the IFFT will truncate.
%
%% See also
%
% <bfft.html |bfft|> |
% <matlab:doc('ifft') ifft>
%
% (C) Lionel Barnett and Anil K. Seth, 2012. See file license.txt in
% installation directory for licensing terms.
%
%%

function AF = bifft(A,q)

[n1,n2,p] = size(A);

if nargin < 2 || isempty(q)
    q = p;
else
    if q < p
        fprintf(2,'WARNING: frequency resolution too low - ifft will truncate!\n');
    end
end

AF = reshape(ifft(reshape(A,n1*n2,p).',q).',n1,n2,q); % q-point inverse Fourier transform on [0, 2*pi)
