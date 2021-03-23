%% var_to_autocov
%
% Return autocovariance sequence for a VAR model
%
% <matlab:open('var_to_autocov.m') code>
%
%% Syntax
%
%     [G,info] = var_to_autocov(A,SIG,acmaxlags,acdectol,aitr,maxiters,maxrelerr)
%
%% Arguments
%
% See also <mvgchelp.html#4 Common variable names and data structures>.
%
% _input_
%
%     A          VAR coefficients matrix
%     SIG        residuals covariance matrix
%     acmaxlags  maximum autocovariance lags to calculate or (default) zero for automatic calculation
%     acdectol   autocovariance decay tolerance (default: 1e-8)
%     aitr       use "accelerated iterative" Lyapunov solver algorithm, else (default) a Schur algorithm (see Description)
%     maxiters   maximum iterations if using iterative algorithm; see 'dlyap_aitr' for defaults
%     maxrelerr  maximum relative error if using iterative algorithm; see 'dlyap_aitr' for defaults
%
% _output_
%
%     G          autocovariance sequence
%     info       info structure, with fields (some may not be present):
%         error      error number (0 if no error)
%         errmsg     error message string
%         warnings   number of warnings (0 if no warnings)
%         warnmsg    warning mesage strings (cell array)
%         rho        VAR spectral radius
%         iters      number of iterations if using iterative algorithm, else 0
%         acrelerr   relative error of associated 1-lag solution
%         acminlags  minimum lags required to achieve specified autocovariance decay factor
%         aclags     actual number of autocovariance lags calculated
%
%% Description
%
% Returns autocovariance sequence |G| defined as [[ii_acseq.png]]
% for a VAR model with coefficients |A| and (positive-definite) residual
% covariance matrix |SIG|, by "reverse-solving" the Yule-Walker equations
%
% <<eq_yweqs.png>>
%
% (where  [[ii_Sigma.png]] = |SIG|). The algorithm solves the associated
% 1-lag problem - a discrete-time Lyapunov equation - and then calculates
% higher lags recursively [1].
%
% Errors, warnings and diagnostics are returned in the |info| struct. Use the
% routine <var_info.html |var_info|> (which should _always_ be called after this
% function) to display this information. Possible errors are
%
%     info.error     info.errmsg
%     ----------------------------------------------------------------------
%         0          (no error, no message)
%         1          unstable VAR (has unit root)
%         2          residuals covariance matrix not positive-definite
%         3          Lyapunov equation solver failed for some reason
%         4          1-lag residuals covariance matrix not positive-definite
%     ----------------------------------------------------------------------
%
% For a stable VAR the the spectral radius |info.rho|  (see
% <var_specrad.html |var_specrad|>) must be < 1; this may be
% considered a unit root test for stationarity [1]. Then the autocovariance
% sequence decays approximately exponentially, by a factor equal to
% |info.rho|. The minimum number of lags required to achieve the specified
% autocovariance decay tolerance |acdectol| is calculated as
% |info.acminlags|, so that |info.rho^info.acminlags < acdectol|. The actual
% number of lags |info.aclags| to which autocovariance is calculated is then
% set to the minimum of |info.acminlags| and the specified maximum number of
% lags, |acmaxlags| (if |acmaxlags| is not supplied - the recommended
% option - it defaults to |info.acminlags|). A warning is issued if
% |info.aclags < info.acminlags|. In this case there is no guarantee that
% MVGCs - particularly in the spectral domain - will be accurate. However,
% if the spectral radius of the VAR model is close to 1, so that
% |info.acminlags| is unfeasibly large, there may be no alternative [note
% that most Granger causality libraries effectively set |info.aclags| to the
% model order]. The calculated 1-lag autocovariance matrix is also checked
% for positive-definitiveness. If this check fails, it may be an indication
% of an ill-conditioned VAR (possibly because residuals variances are too
% small, and/or the process is borderline stationary). In short, the
% |info.error| field _must_ be checked by the caller: 0 signifies success, >
% 0 signifies an error with corresponding message in |info.errmsg|.
% Attention should also be paid to any warnings in |info.WARNINGn|.
%
% The |aitr| flag specifies a fast (but experimental) "accelerated" iterative
% Lyapunov solver algorithm (see <dlyap_aitr.html |dlyap_aitr|>). If not set
% (default), the Control System Toolbox <matlab:doc('dlyap') |dlyap|> solver
% routine is used if available, else a roughly equivalent (but slower) scripted
% algorithm based on Schur decomposition (see <matlab:open('dlyap.m') |dlyap.m|>
% in the |utils/control| directory).
%
%% References
%
% [1] L. Barnett and A. K. Seth,
% <http://www.sciencedirect.com/science/article/pii/S0165027013003701 The MVGC
%     Multivariate Granger Causality Toolbox: A New Approach to Granger-causal
% Inference>, _J. Neurosci. Methods_ 223, 2014
% [ <matlab:open('mvgc_preprint.pdf') preprint> ].
%
% [2] X. Kitagawa, An algorithm for solving the matrix equation X = F X F'+S,
%     _Internat. J. Control_ 25(5), 1977.
%
% [3] S. Hammarling, "Numerical solution of the stable, non-negative definite
%     Lyapunov equation", _IMA J. Numer. Anal._ 2, 1982.
%
%% See also
%
% <var_specrad.html |var_specrad|> |
% <var_info.html |var_info|> |
% <matlab:doc('dlyap') |dlyap (ControlSystem Toolbox)|> |
% <dlyap.html |dlyap (scripted Schur algorithm)|> |
% <dlyap_aitr.html |dlyap_aitr|>
%
% (C) Lionel Barnett and Anil K. Seth, 2012. See file license.txt in
% installation directory for licensing terms.
%
%%

function [G,info] = var_to_autocov(A,SIG,acmaxlags,acdectol,aitr,maxiters,maxrelerr)

%global have_dlyap;

% default parameters

if nargin < 3 || isempty(acmaxlags), acmaxlags = 0;     end % calculate maximum lags automatically
if nargin < 4 || isempty(acdectol),  acdectol  = 1e-8;  end % autocovariance decay tolerance
if nargin < 5 || isempty(aitr),      aitr      = false; end % use "accelerated" iterative Lyapunov equation solver 

% iterative algorithm only: ensure defaults for utils/dlyap_aitr.m.

if nargin < 6, maxiters  = []; end
if nargin < 7, maxrelerr = []; end

[n,n1,p] = size(A);
assert(n1 == n,'VAR coefficients matrix has bad shape');
pn1 = (p-1)*n;

[nn1,nn2] = size(SIG);
assert(nn1 == nn2,'residuals covariance matrix not square');
assert(nn1 == n  ,'residuals covariance matrix doesn''t match VAR coefficients matrix');

% initialise info struct
info.error     = 0;
info.errmsg    = '';
info.warnings  = 0;
info.warnmsg   = cell(0,1);
info.rho       = NaN;
info.iters     = NaN;
info.acrelerr  = NaN;
info.acminlags = NaN;
info.aclags    = NaN;

G = [];

% construct VAR coefficients for 1-lag problem

A1 = [reshape(A,n,p*n); eye(pn1) zeros(pn1,n)];

% calculate spectral radius

info.rho = max(abs(eig(A1)));

if info.rho >= 1
    info.error = 1;
    info.errmsg = 'unstable VAR (unit root)';
    return
end

% construct residual covariances for 1-lag problem

if ~isposdef(SIG);
    info.error = 2;
    info.errmsg = 'residuals covariance matrix not positive-definite';
    return;
end

SIG1 = [SIG zeros(n,pn1); zeros(pn1,n) zeros(pn1)];

% solve the Lyapunov equation for the 1-lag covariance matrix

try
    if aitr
        [G1,info.iters] = dlyap_aitr(A1,SIG1,maxiters,maxrelerr); % experimental: fast, but needs more testing
    else
%       G1 = dlyap(A1,SIG1);           % dlyap seems to work better here without balancing, which seems to break positive-definitiveness
        G1 = lyapslv('D',A1,[],-SIG1); % sometimes. However lyapslv is not an official interface, so this could conceivably break in future.
    end
catch except
    info.error = 3;
    info.errmsg = ['Lyapunov equation solver failed: ' except.message];
    return
end

info.acrelerr = norm(A1*G1*A1'-G1+SIG1)/norm(SIG1); % this should be small (see below)

maxacrelerr = 1e-8; % probably shouldn't be hard-coded :-/
if info.acrelerr > maxacrelerr
    info.warnings = info.warnings+1;
    info.warnmsg{info.warnings} = sprintf('large relative error = %g (tolerance = %g)',info.acrelerr,maxacrelerr);
end

% estimate number of autocov lags

info.acminlags = ceil(log(acdectol)/log(info.rho)); % minimum lags to achieve specified tolerance

if     acmaxlags < 0  % use exactly -acmaxlags lags (not encouraged, hence undocumented!)
    info.aclags = -acmaxlags;
elseif acmaxlags > 0  % use at most acmaxlags lags
    info.aclags = min(info.acminlags,acmaxlags);
else                  % acmaxlags == 0 - use minimum acceptable lags (recommended)
    info.aclags = info.acminlags;
end

if info.aclags < info.acminlags
    info.warnings = info.warnings+1;
    info.warnmsg{info.warnings} = sprintf('too few autocovariance lags = %d (minimum = %d)',info.aclags,info.acminlags);
end

if ~isposdef(G1);
    info.error = 4;
    info.errmsg = '1-lag covariance matrix not positive-definite';
    return
end

q = info.aclags;
q1 = q+1;

% calculate recursively from 1-lag solution (which supplies up to p-1 lags), from p lags up to q

[n,~,p]  = size(A);
assert(info.aclags >= p,'number of lags is too small'); % lags must be at least number of VAR lags
pn = p*n;
G = cat(3,reshape(G1(1:n,:),n,n,p),zeros(n,n,q1-p));   % autocov forward  sequence
B = [zeros((q1-p)*n,n); G1(:,end-n+1:end)];            % autocov backward sequence
A = reshape(A,n,pn);                                   % coefficients
for k = p:q
    r = q1-k;
    G(:,:,k+1) = A*B(r*n+1:r*n+pn,:);
    B((r-1)*n+1:r*n,:) = G(:,:,k+1);
end
