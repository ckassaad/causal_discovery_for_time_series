%% autocov_to_smvgc
%
% Calculate conditional frequency-domain MVGC (spectral multivariate Granger causality)
%
% <matlab:open('autocov_to_smvgc.m') code>
%
%% Syntax
%
%     [f,fres] = autocov_to_smvgc(G,x,y,fres,useFFT)
%
%% Arguments
%
% See also <mvgchelp.html#4 Common variable names and data structures>.
%
% _input_
%
%     G          autocovariance sequence
%     x          vector of indices of target (causee) multi-variable
%     y          vector of indices of source (causal) multi-variable
%     fres       frequency resolution (default: automatic)
%     useFFT     use FFT method for autocovariance transform (default: as for autocov_xform)
%
% _output_
%
%     f          spectral Granger causality
%
%% Description
%
% Returns the frequency-domain (spectral) MVGC
%
% <<eq_smvgc.png>>
%
% from the variable |Y| (specified by the vector of indices |y|) to the
% variable |X| (specified by the vector of indices |x|), conditional on all
% other variables |Z| represented in |G|, for a stationary VAR process with
% autocovariance sequence G.
%
% Spectral causality is calculated up to the Nyqvist frequency at a
% resolution |fres|. If |fres| is not supplied it is calculated optimally
% as the number of autocovariance lags. Call |freqs =
% <sfreqs.html sfreqs>(fres,fs)|, where |fs| is the sampling
% rate, to get a corresponding vector |freqs| of frequencies on |[0,fs/2]|.
%
% The |useFFT| flag specifies the algorithm used to transform the
% autocovariance sequence; see <autocov_xform.html |autocov_xform|> for
% details.
%
% In the conditional case, the algorithm works by transforming the
% autocovariance sequence for the full regression (see
% <autocov_xform.html |autocov_xform|>) to an autocovariance
% sequence for new |X,Z| variables defined as residuals of the reduced
% regression; thus a separate estimation step for the reduced regression,
% which is known to be problematic [2,*], is unnecessary, resulting in high
% efficiency and accuracy. See [1] for details.
%
% The caller should take note of any warnings issued by this function and test
% results with a call <isbad.html |isbad|>|(f,false)|.
%
%% References
%
% [1] L. Barnett and A. K. Seth,
% <http://www.sciencedirect.com/science/article/pii/S0165027013003701 The MVGC
%     Multivariate Granger Causality Toolbox: A New Approach to Granger-causal
% Inference>, _J. Neurosci. Methods_ 223, 2014
% [ <matlab:open('mvgc_preprint.pdf') preprint> ].
%
% [2] Y. Chen, S. L. Bressler and M. Ding, "Frequency decomposition of
% conditional Granger causality and application to multivariate neural
% field potential data", _J. Neurosci. Methods_, 150, 2006. 
%
% [*] In our experience the "partition matrix" method in ref. [2] appears to be
% unsound, producing inaccurate results; hence we do not use it here.
%
%% See also
%
% <autocov_to_var.html |autocov_to_var|> |
% <var_to_cpsd.html |var_to_cpsd|> |
% <autocov_xform.html |autocov_xform|> |
% <var2trfun.html |var2trfun|> |
% <sfreqs.html |sfreqs|> |
% <isbad.html |isbad|>
%
% (C) Lionel Barnett and Anil K. Seth, 2012. See file license.txt in
% installation directory for licensing terms.
%
%%

function [f,fres] = autocov_to_smvgc(G,x,y,fres,useFFT)

[n,~,q1] = size(G);

x = x(:)'; % vectorise
y = y(:)'; % vectorise

assert(all(x >=1 & x <= n),'some x indices out of range');
assert(all(y >=1 & y <= n),'some y indices out of range');
assert(isempty(intersect(x,y)),'x and y indices must be distinct');

z = 1:n; z([x y]) = []; % indices of other variables (to condition out)

if nargin < 4 || isempty(fres), fres = q1; end

if nargin < 5, useFFT = []; end % force autocov_xform default

h = fres+1;
f = nan(1,h);

if isempty(z) % unconditional

    xy = [x y];
    G = G(xy,xy,:); % rearrange
    nx = length(x);
    x = 1:nx;
    y = nx+1:n;

    owstate = warn_supp;
    [A,SIG] = autocov_to_var(G);                  % full regression
    warn_test(owstate,   'in full regression - bad autocovariance matrix? Check output of ''var_info''');
    if warn_if(isbad(A), 'in full regression - regression failed'), return; end % show-stopper!
    
    owstate = warn_supp;
    [S,H] = var_to_cpsd(A,SIG,fres);              % spectrum & transfer function
    warn_test(owstate,   'in spectral calculation');
    if warn_if(isbad(S), 'in spectral calculation - calculation failed'), return; end % show-stopper!
    
    S = S(x,x,:);                                 % only need x part of cpsd
    H = H(x,y,:);                                 % only need x,y part of transfer function
    PSIG = SIG(y,y)-(SIG(y,x)/SIG(x,x))*SIG(x,y); % partial covariance
    for k = 1:h
        f(k) = log(real(det(S(:,:,k)))) - log(real(det(S(:,:,k)-H(:,:,k)*PSIG*H(:,:,k)')));
    end

else % conditional
    
    % transform autocov by reduced regression of x,z

    xzy = [x z y];
    G = G(xzy,xzy,:); % rearrange
    nx = length(x);
    nz = length(z);
    x = 1:nx;
    xz = 1:(nx+nz);
    zy = nx+1:n;
    
    owstate = warn_supp;
    [AR,SIGR] = autocov_to_var(G(xz,xz,:));           % reduced regression
    warn_test(owstate,    'in reduced regression - bad autocovariance matrix? Check output of ''var_info''');
    if warn_if(isbad(AR), 'in reduced regression - regression failed'), return; end % show-stopper!
    
    G = autocov_xform(G,AR,SIGR,useFFT);              % transform autocov

    % now do unconditional with transformed autocov

    owstate = warn_supp;
    [A,SIG] = autocov_to_var(G);                      % full regression
    warn_test(owstate,   'in full regression - bad autocovariance matrix? Check output of ''var_info''');
    if warn_if(isbad(A), 'in full regression - regression failed'), return; end % show-stopper!
    
    S = SIGR(x,x);                                    % only need x part of cpsd - which is flat spectrum!
    
    owstate = warn_supp;
    H = var2trfun(A,fres);                            % transfer function
    warn_test(owstate,   'in transfer function calculation');
    if warn_if(isbad(H), 'in transfer function calculation - calculation failed'), return; end % show-stopper!
    
    H = H(x,zy,:);                                    % only need x,zy part of transfer function
    PSIG = SIG(zy,zy)-(SIG(zy,x)/SIG(x,x))*SIG(x,zy); % partial covariance
    LDS = log(det(S));
    for k = 1:h
        f(k) = LDS-log(real(det(S-H(:,:,k)*PSIG*H(:,:,k)')));
    end

end
