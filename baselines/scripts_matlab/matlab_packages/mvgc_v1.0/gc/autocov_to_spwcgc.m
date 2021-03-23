%% autocov_to_spwcgc
%
% Calculate pairwise-conditional frequency-domain MVGCs (spectral multivariate Granger causalites)
%
% <matlab:open('autocov_to_spwcgc.m') code>
%
%% Syntax
%
%     [f,fres] = autocov_to_spwcgc(G,fres,useFFT)
%
%% Arguments
%
% See also <mvgchelp.html#4 Common variable names and data structures>.
%
% _input_
%
%     G          autocovariance sequence
%     fres       frequency resolution
%     useFFT     use FFT method for autocovariance transform (default: as for autocov_xform)
%
% _output_
%
%     f          spectral Granger causality matrix
%
%% Description
%
% Returns the  matrix |f| of pairwise-conditional frequency-domain
% (spectral) MVGCs
%
% <<eq_smvgc_pwc.png>>
%
% (where |[ij]| denotes omission of the |ij|-th variables) between all
% pairs of variables [[ii_inej.png]] represented in |G|, for a stationary VAR
% process with autocovariance sequence |G|. The first index |i| of
% |f| is the target (causee) variable, the second |j| the source (causal)
% variable and the third indexes the frequency. See ref. [1] for details.
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
% The caller should take note of any warnings issued by this function and test
% results with a call <isbad.html |isbad|>|(f,false)|.
%
% For details of the algorithm, see <autocov_to_smvgc.html |autocov_to_smvgc|> and [1].
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
% <autocov_to_smvgc.html |autocov_to_smvgc|> |
% <autocov_to_pwcgc.html |autocov_to_pwcgc|> |
% <autocov_to_var.html |autocov_to_var|> |
% <var2trfun.html |var2trfun|> |
% <autocov_xform.html |autocov_xform|> |
% <sfreqs.html |sfreqs|> |
% <isbad.html |isbad|>
%
% (C) Lionel Barnett and Anil K. Seth, 2012. See file license.txt in
% installation directory for licensing terms.
%
%%

function [f,fres] = autocov_to_spwcgc(G,fres,useFFT)

if nargin < 3, useFFT = []; end % force autocov_xform default

[n,~,q1] = size(G);
if nargin < 2 || isempty(fres);
    fres = q1;
end

h = fres+1;
f = nan(n,n,h);

for j = 1:n
    jo  = [1:j-1 j+1:n]; % omit j
    joj = [jo j];        % rearrange with omitted j indices first

    owstate = warn_supp;
    [Aj,SIGj] = autocov_to_var(G(jo,jo,:));  % reduced regression
    warn_test(owstate,   sprintf('in reduced regression for target node %d - bad autocovariance matrix? Check output of ''var_info''',j));
    if warn_if(isbad(Aj),sprintf('in reduced regression for target node %d - regression failed',j)), continue; end % show-stopper!
    
    Gj = autocov_xform(G(joj,joj,:),Aj,SIGj,useFFT); % transform autocov
    
    owstate = warn_supp;
    [Ajj,SIGjj] = autocov_to_var(Gj);        % transformed full regression
    warn_test(owstate,    sprintf('in transformed regression for target node %d - bad autocovariance matrix? Check output of ''var_info''',j));
    if warn_if(isbad(Ajj),sprintf('in transformed regression for target node %d - regression failed',j)), continue; end % show-stopper!
    
    owstate = warn_supp;
    Hjj = var2trfun(Ajj,fres);               % transfer function
    warn_test(owstate,    sprintf('in transfer function for target node %d',j));
    if warn_if(isbad(Hjj),sprintf('in transfer function for target node %d - calculation failed',j)), continue; end % show-stopper!

    for ii=1:n-1;
        i  = jo(ii);           % i index in omitted j indices
        io = [1:ii-1 ii+1:n];  % omit i

        SIGji = SIGjj(io,io)-(SIGjj(io,ii)*SIGjj(ii,io))/SIGjj(ii,ii); % partial covariance
        Hji = Hjj(ii,io,:);                                            % transfer function
        Sji = SIGj(ii,ii);                                             % i part of spectrum is flat!
        
        LSji = log(Sji);
        for k = 1:h
            f(i,j,k) = LSji - log(real(Sji-Hji(:,:,k)*SIGji*Hji(:,:,k)'));
        end
    end
end
