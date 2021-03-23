%% cpsd_to_var
%
% Spectral factorisation: calculate VAR parameters from cross-power spectral density
%
% <matlab:open('cpsd_to_var.m') code>
%
%% Syntax
%
%     [H,SIG,iters] = cpsd_to_var(S,G0,maxiters,numtol)
%
%% Arguments
%
% See also <mvgchelp.html#4 Common variable names and data structures>.
%
% _input_
%
%     S          cross-power spectral density (cpsd) matrix
%     G0         covariance matrix of VAR process (default: calculated automatically from cpsd)
%     maxiters   maximum iterations (default: 100)
%     numtol     numerical tolerance (default: 1e-10)
%
% _output_
%
%     H          VAR transfer function matrix
%     SIG        VAR residuals covartiance matrix
%     iters      number of iterations performed
%
%% Description
%
% Calculates the transfer function |H| and residuals covariance matrix |SIG|
% from the cpsd |S| and covariance matrix |G0| of a VAR process, using Wilson's
% spectral factorisation algorithm [2]. If |G0| is not supplied (default), then
% |G| is calculated from |S| (see <cpsd_to_autocov.html |cpsd_to_autocov|>) and
% |G0| set to |G(:,:,1)|. The actual number of iterations performed is returned
% in |iters|. If the algorithm fails to converge to numerical tolerance |numtol|
% within |maxiters| iterations, an exception |MVGC:XMaxItrs| is thrown.
%
% *_Note:_* to calculate the VAR coefficients, use the utility function
% <trfun2var.html |trfun2var|>.
%
% Adapted from original code with the kind permission of
% <http://math.iisc.ernet.in/~rangaraj/ Prof. G. Rangarajan> of the Dept. of
% Mathematics, Indian Institute of Science, Bangalore, India; see [3,4] for
% applications to Granger-causal analysis.
%
%% References
%
% [1] L. Barnett and A. K. Seth,
% <http://www.sciencedirect.com/science/article/pii/S0165027013003701 The MVGC
%     Multivariate Granger Causality Toolbox: A New Approach to Granger-causal
% Inference>, _J. Neurosci. Methods_ 223, 2014
% [ <matlab:open('mvgc_preprint.pdf') preprint> ].
%
% [2] G. T. Wilson, "The factorization of matricial spectral densities", _SIAM
% Journal on Applied Mathematics_, 23(4), 1972.
%
% [3] M. Dhamala, G. Rangarajan and M. Ding, "Estimating Granger causality from
% Fourier and wavelet transforms of time series data", _Phys. Rev. Lett._ 100,
% 2008.
%
% [4] M. Dhamala, G. Rangarajan and M. Ding, "Analyzing information flow in
% brain networks with nonparametric Granger causality", _NeuroImage_ 41, 2008.
%
%% See also
%
% <cpsd_to_autocov.html |cpsd_to_autocov|> |
% <autocov_to_var.html |autocov_to_var|> |
% <trfun2var.html |trfun2var|>
%
%% Copyright notice
%
% [(C)] _Lionel Barnett and Anil K. Seth, 2012. See file
% <matlab:open('license.txt') license.txt> in root directory for licensing
% terms._
%
%%

function [H,SIG,iters] = cpsd_to_var(S,G0,maxiters,numtol)

% defaults

if nargin < 2 || isempty(G0), G = cpsd_to_autocov(S); G0 = G(:,:,1); end % calculate covariance matrix from cpsd
if nargin < 3 || isempty(maxiters), maxiters = 100;   end
if nargin < 4 || isempty(numtol),   numtol   = 1e-10; end

[n,~,h] = size(S);
h2 = 2*(h-1);

% extend spectrum

SX = zeros(n,n,h2);
SX(:,:,1) = S(:,:,1);
for k = 2:h
    SX(:,:,k) = S(:,:,k);
    SX(:,:,2*h-k) = S(:,:,k).';
end

% initialise P

C = chol(G0); % this is psi_1
P = zeros(n,n,h2);
for k = 1:h2,
    P(:,:,k) = C; % initialization for the 1st iteration
end

I = eye(n);
g = zeros(n,n,h2);
d = intmax;
for iters = 1:maxiters,

    % calculate g
    
    for k = 1:h2,
        Pkinv = I/P(:,:,k);
        g(:,:,k) = Pkinv*SX(:,:,k)*Pkinv'+I;
    end
    
    % calculate [g]+ (take positive lags only and half of the zero lag)

    beta = bifft(g);
    Gp = beta;
    beta0 = 0.5*beta(:,:,1); 
    Gp(:,:,1) = triu(beta0);  % this is Stau
    Gp(:,:,h+1:end) = 0;
    gp = bfft(Gp);
    
    % update P
    
    Pold = P;
    dold = d;
    for k = 1:h2,
        P(:,:,k) = P(:,:,k)*gp(:,:,k);
    end
    d = maxabs(P-Pold);
    if abs(d-dold) < numtol,
        break;
    end

end 

% calculate coefficients and residuals covariance

AA = real(bifft(P));
A0 = AA(:,:,1); % this is psi_1
SIG = A0*A0';
A0inv = I/A0;
H = zeros(n,n,h);
for k = 1:h,
    H(:,:,k) = P(:,:,k)*A0inv;
end 

if iters == maxiters
    throw(MException('MVGC:XMaxItrs','tolerance not met at %d iterations (norm = %e, numtol = %e)\n',maxiters,maxabs(d-dold),numtol));
end
