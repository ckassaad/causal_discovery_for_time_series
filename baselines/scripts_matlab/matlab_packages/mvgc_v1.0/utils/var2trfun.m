%% var2trfun
%
% Calculate VAR transfer function from VAR coefficients
%
% <matlab:open('var2trfun.m') code>
%
%% Syntax
%
%     H = var2trfun(A,fres)
%
%% Arguments
%
% See also <mvgchelp.html#4 Common variable names and data structures>.
%
% _input_
%
%     A          VAR coefficients matrix
%     fres       frequency resolution
%
% _output_
%
%     H          VAR transfer function matrix
%
%% Description
%
% Return transfer function |H| for VAR with coefficients |A|. |fres|
% specifies the frequency resolution. Call |freqs = <sfreqs.html
% sfreqs>(fres,fs)|, where |fs| is the sampling rate, to get a corresponding
% vector |freqs| of frequencies on |[0,fs/2]|.
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
% <trfun2var.html |trfun2var|> |
% <var_to_cpsd.html |var_to_cpsd|> |
% <autocov_to_smvgc.html |autocov_to_smvgc|> |
% <autocov_to_spwcgc.html |autocov_to_spwcgc|>
%
% (C) Lionel Barnett and Anil K. Seth, 2012. See file license.txt in
% installation directory for licensing terms.
%
%%

function H = var2trfun(A,fres)

n = size(A,1);
I = eye(n);
AF = bfft(cat(3,I,-A),2*fres); % over [0,2*pi)
h = fres+1;
H = zeros(n,n,h);
for k = 1:h % over [0,pi] only
    H(:,:,k) = I/AF(:,:,k);
end
