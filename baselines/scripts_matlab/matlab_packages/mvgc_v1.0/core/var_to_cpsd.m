%% var_to_cpsd
%
% Calculate cross-power spectral density and transfer function from VAR
% parameters
%
% <matlab:open('var_to_cpsd.m') code>
%
%% Syntax
%
%     [S,H] = var_to_cpsd(A,SIG,fres)
%
%% Arguments
%
% See also <mvgchelp.html#4 Common variable names and data structures>.
%
% _input_
%
%     A          VAR coefficients matrix
%     SIG        residuals covariance matrix
%     fres       frequency resolution
%
% _output_
%
%     S          cross-power spectral density (cpsd) matrix
%     H          VAR transfer function matrix
%
%% Description
%
% Calculates cross-spectral density |S| and transfer function |H| of VAR process
% with coefficients matrix |A| and residuals covariance matrix |SIG|. |fres|
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
% <var2trfun.html |var2trfun|> |
% <sfreqs.html |sfreqs|>
%
% (C) Lionel Barnett and Anil K. Seth, 2012. See file license.txt in
% installation directory for licensing terms.
%
%%

function [S,H] = var_to_cpsd(A,SIG,fres)

n = size(A,1);
h = fres+1;
H = var2trfun(A,fres);
S = zeros(n,n,h);
for k = 1:h
    S(:,:,k) = H(:,:,k)*SIG*H(:,:,k)';
end
