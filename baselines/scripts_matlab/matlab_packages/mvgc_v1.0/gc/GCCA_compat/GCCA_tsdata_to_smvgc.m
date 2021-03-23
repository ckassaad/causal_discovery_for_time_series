%% GCCA_tsdata_to_smvgc
%
% Calculate _unconditional_ frequency-domain MVGC (spectral multivariate
% Granger causality) from time series data by "traditional" method (as e.g.
% in GCCA toolbox)
%
% <matlab:open('GCCA_tsdata_to_smvgc.m') code>
%
%% Syntax
%
%     f = GCCA_tsdata_to_smvgc(U,x,y,p,fres,regmode)
%
%% Arguments
%
% See also <mvgchelp.html#4 Common variable names and data structures>.
%
% _input_
%
%     U          multi-trial time series data
%     x          vector of indices of target (causee) multi-variable
%     y          vector of indices of source (causal) multi-variable
%     p          model order (number of lags)
%     fres       frequency resolution
%     regmode    regression mode (default as in 'tsdata_to_var')
%
% _output_
%
%     f          unconditional spectral Granger causality
%
%% Description
%
% Returns the _unconditional_ frequency-domain (spectral) MVGC
%
% <<eq_smvgc_uncond.png>>
%
% from the variable |Y| (specified by the vector of indices |y|) to the
% variable |X| (specified by the vector of indices |x|) in the time series
% data |U|, for model order |p|.
%
% *_NOTE:_* _This routine will *not* compute conditional causalities_,
% since the "traditional" full/reduced regression method may well produce
% unacceptably inaccurate (or nonesensical) results in the conditional
% case. For an analysis of the reasons for this, see [2,*]. For the same
% reason, there is no pairwise-conditional |GCCA_tsdata_to_smvgc_pwc|
% routine in this toolbox. Since conditional causalities are generally what
% you want (unconditional causalities may be highly misleading [1]), we
% strongly reccommend computation of conditional spectral MVGC using the
% preferred autocovariance-based approach; see <autocov_to_smvgc.html
% |autocov_to_smvgc|>.
%
% Spectral causality is calculated up to the Nyqvist frequency at a
% resolution |fres|. Call |freqs = <sfreqs.html sfreqs>(fres,fs)|, where
% |fs| is the sampling rate, to get a corresponding vector |freqs| of
% frequencies on |[0,fs/2]|. The regression mode is set by the |regmode|
% parameter, which may be |'LWR'| or |'OLS'| (see <tsdata_to_var.html
% |tsdata_to_var|> for details and defaults).
%
% This routine is included mainly for compatibility with the Granger Causal
% Connectivity Analysis (<http://www.sussex.ac.uk/Users/anils/aks_code.htm
% GCCA>) Toolbox [3]; the preferred MVGC Toolbox method of calculating
% spectral MVGC via the autocovariance sequence (see <autocov_to_smvgc.html
% |autocov_to_smvgc|>) will generally be more accurate and, furthermore,
% will calculate conditional spectral causality.
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
% [3] A. K. Seth, "A MATLAB toolbox for Granger causal connectivity
% analysis", _Journal of Neuroscience Methods_ 186, 2010.
%
% [*] In our experience the "partition matrix" method in ref. [2] appears to be
% unsound, producing inaccurate results; hence we do not use it here.
%
%% See also
%
% <tsdata_to_var.html |tsdata_to_var|> |
% <var_to_cpsd.html |var_to_cpsd|> |
% <sfreqs.html |sfreqs|> |
% <autocov_to_smvgc.html |autocov_to_smvgc|> |
% <isbad.html |isbad|>
%
% (C) Lionel Barnett and Anil K. Seth, 2012. See file license.txt in
% installation directory for licensing terms.
%
%%

function f = GCCA_tsdata_to_smvgc(U,x,y,p,fres,regmode)

if nargin < 6, regmode = []; end % ensure 'tsdata_to_var' default

n = size(U,1);

x = x(:)'; % vectorise
y = y(:)'; % vectorise

assert(all(x >=1 & x <= n),'some x indices out of range');
assert(all(y >=1 & y <= n),'some y indices out of range');
assert(isempty(intersect(x,y)),'x and y indices must be distinct');

z = 1:n; z([x y]) = []; % indices of other variables (to condition out)

h = fres+1;
f = nan(h,1);

assert(isempty(z),'conditional spectral MVGC not available in GCCA mode');

xy = [x y];

U = U(xy,:,:); % extract variables, rearrange

nx = length(x);
n = length(xy);
x = 1:nx;
y = nx+1:n;

owstate = warn_supp;
[A,SIG] = tsdata_to_var(U,p,regmode);        % full regression
warn_test(owstate,  'in full regression - data non-stationary or colinear?');
if warn_if(isbad(A),'in full regression - regression failed'), return; end % show-stopper!
% we should really check that SIG is positive-definite; we don't, for effciency reasons

owstate = warn_supp;
[S,H] = var_to_cpsd(A,SIG,fres);             % spectrum & transfer function
warn_test(owstate,  'in spectral calculation');
if warn_if(isbad(S),'in spectral calculation - calculation failed'), return; end % show-stopper!

S = S(x,x,:);
H = H(x,y,:);
SIG = SIG(y,y)-(SIG(y,x)/SIG(x,x))*SIG(x,y); % partial covariance
for k = 1:h
    f(k) = log(real(det(S(:,:,k)))) - log(real(det(S(:,:,k)-H(:,:,k)*SIG*H(:,:,k)')));
end

