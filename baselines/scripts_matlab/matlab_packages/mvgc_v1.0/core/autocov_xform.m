%% autocov_xform
%
% Transform autocovariance sequence for reduced regression
%
% <matlab:open('autocov_xform.m') code>
%
%% Syntax
%
%     G = autocov_xform(G,AR,SIGR,useFFT)
%
%% Arguments
%
% See also <mvgchelp.html#4 Common variable names and data structures>.
%
% _input_
%
%     G          autocovariance sequence
%     AR         VAR coefficients matrix for reduced regression
%     SIGR       residuals covariance matrix for reduced regression
%     useFFT     use FFT method (default: true)
%
% _output_
%
%     G         transformed autocovariance sequence
%
%% Description
%
% Returns the autocovariance sequence |G| for a new variable defined as the
% residuals of a reduced regression, for a VAR with autocovariance sequence |G|.
% |AR| and |SIGR| are the coefficients matrices and residuals covariance matrix
% respectively of the reduced regression, which is is assumed to correspond to
% the first |size(AR,1)| indices of |G|.
%
% If the |useFFT| flag is set (default), then the autocovariance sequence is
% converted to a cpsd via FFT (see <autocov_to_cpsd.html |autocov_to_cpsd|>),
% the transformation effected on the cpsd (_cf._ <cpsd_xform.html |cpsd_xform|>)
% and the result converted back to an autocovariance sequence via IFFT (see
% <cpsd_to_autocov.html |cpsd_to_autocov|>). Otherwise, the autocovariance
% sequence is transformed by explicit convolution. The FFT method is generally
% more efficient than the convolution method, particularly if the number of
% autocovariance lags is large.
%
% This function is crucial to the calculation of spectral causality in the
% conditional case; see <autocov_to_smvgc.html |autocov_to_smvgc|>,
% <autocov_to_spwcgc.html |autocov_to_spwcgc|>. In theory, if the original
% autocovariance sequence is calculated to |q| lags - under the assumption that
% it may not have decayed sufficiently for |k < q| lags (see
% <var_to_autocov.html |var_to_autocov|>) - then the transformed autocovariance
% sequence should be calculated to |2q| lags. In practice we find that
% calculating to |q| lags is generally sufficient for good accuracy. To
% calculate |G| to higher lags, the simplest option is to reduce the |acdectol|
% parameter in the call to <var_to_autocov.html |var_to_autocov|> (e.g. squaring
% it will effectively double the number of lags |q| to which |G| and hence |G|
% is calculated).
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
% <var_to_autocov.html |var_to_autocov|> |
% <autocov_to_var.html |autocov_to_var|> |
% <autocov_to_smvgc.html |autocov_to_smvgc|> |
% <autocov_to_spwcgc.html |autocov_to_spwcgc|> |
% <cpsd_xform.html |cpsd_xform|> |
% <autocov_to_cpsd.html |autocov_to_cpsd|> |
% <cpsd_to_autocov.html |cpsd_to_autocov|>
%
% (C) Lionel Barnett and Anil K. Seth, 2012. See file license.txt in
% installation directory for licensing terms.
%
%%

function G = autocov_xform(G,AR,SIGR,useFFT)

if nargin < 4 || isempty(useFFT), useFFT = true; end

[n,~,q1] = size(G);
q = q1-1;

[nx,nx1,~] = size(AR);
assert(nx1 == nx,'reduced VAR coefficients matrix has bad shape');
assert(nx <= n,'reduced VAR coefficients matrix appears to be for more variables than autocovariance sequence');

[n1,n2] = size(SIGR);
assert(n1 == n2,'reduced VAR residuals covariance matrix not square');
assert(n1 == nx ,'reduced VAR residuals covariance matrix doesn''t match reduced VAR coefficients matrix');

ny = n-nx;
x = 1:nx;
y = nx+1:n;

% transform autocov by reduced regression

if useFFT % convert to cpsd, transform, convert back to autocov
    
    S = autocov_to_cpsd(G);
    AF = bfft(cat(3,eye(nx),-AR),2*q1); % over [0,2*pi)
    for k = 1:q1+1
        S(x,x,k) = SIGR; % flat spectrum, since residuals serially uncorrelated
        S(x,y,k) = AF(:,:,k)*S(x,y,k);
        S(y,x,k) = S(y,x,k)*AF(:,:,k)';
    end
    G = cpsd_to_autocov(S);
    
else      % explicit convolution

    AR = reshape(cat(3,eye(nx),-AR),nx,q1*nx);
    GF = reshape(G(y,x,:),ny,q1*nx);                             % forward  autocovariance sequence
    GB = reshape(permute(flipdim(G(x,y,:),3),[1 3 2]),q1*nx,ny); % backward autocovariance sequence
    G(x,x,1) = SIGR;                                             % just the reduced residuals covariance,
    G(x,x,2:end) = zeros(nx,nx,q);                               % since residuals serially uncorrelated
    for k = 0:q-1
        G(x,y,k+1) = AR(:,1:(k+1)*nx)*GB((q-k)*nx+1:q1*nx,:) + AR(:,(k+1)*nx+1:q1*nx)*GF(:,nx+1:(q1-k)*nx)';
    end
    G(x,y,q1) = AR*GB;
    for k = 0:q
        G(y,x,k+1) = GF(:,k*nx+1:q1*nx)*AR(:,1:(q1-k)*nx)';
    end
    
end
