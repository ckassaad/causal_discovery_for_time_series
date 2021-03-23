%% empirical_var_to_spwcgc
%
% Calculate sampling distribution for pairwise-conditional frequency-domain
% MVGCs from generated time series data for a specified VAR model
%
% <matlab:open('empirical_var_to_spwcgc.m') code>
%
%% Syntax
%
%     fE = empirical_var_to_spwcgc(A,SIG,m,N,fres,H0,nsamps,mtrunc,decayfac,regmode,acmaxlags,acdectol)
%
%% Arguments
%
% See also <mvgchelp.html#4 Common variable names and data structures>.
%
% _input_
%
%     A          VAR coefficients matrix
%     SIG        residuals covariance matrix
%     m          number of observations per trial
%     N          number of trials (default: 1)
%     fres       frequency resolution (default: automatic)
%     H0         flag: impose null hypotheses of zero connectivity?
%     nsamps     number of bootstrap samples
%     mtrunc     number of initial time observations to truncate  (default as for 'var_to_tsdata')
%     decayfac   initial transients decay factor (default as for 'var_to_tsdata')
%     regmode    regression mode (default as for 'tsdata_to_var')
%     acmaxlags  maximum autocovariance lags (default as for 'var_to_autocov')
%     acdectol   autocovariance decay tolerance (default as for 'var_to_autocov')
%
% _output_
%
%     fE         empirical spectral Granger causalities distributions
%
%% Description
%
% Returns |nsamps| samples from the empirical sampling distribution of the
% pairwise-conditional frequency-domain MVGCs for time series data generated from the
% VAR model specified by the coefficients matrix |A| and residuals covariance
% matrix |SIG|. The time series comprise |N| trials of length |m|, and residuals
% are iid Gaussian with covariance matrix |SIG|. If the flag |H0| is set then
% data is generated for the nested null models with zero connectivity; i.e. for
% each target index |i| and source index |j|, |A(i,j,k)| is set to zero for all
% lags |k|. For other parameters see <var_to_tsdata.html |var_to_tsdata|>,
% <tsdata_to_var.html |tsdata_to_var|> and <var_to_autocov.html
% |var_to_autocov|>.
%
% The first dimension of the returned matrix |fE| indexes samples, the second
% indexes the target (causee) variable, the third the source (causal)
% variable and the fourth frequency.
%
% Spectral causality is calculated up to the Nyqvist frequency at a
% resolution |fres|. If |fres| is not supplied it is calculated optimally
% as the number of autocovariance lags. Call |freqs =
% <sfreqs.html sfreqs>(fres,fs)|, where |fs| is the sampling
% rate, to get a corresponding vector |freqs| of frequencies on |[0,fs/2]|.
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
% <empirical_var_to_mvgc.html |empirical_var_to_mvgc|> |
% <empirical_var_to_pwcgc.html |empirical_var_to_pwcgc|> |
% <empirical_var_to_smvgc.html |empirical_var_to_smvgc|> |
% <var_to_tsdata.html |var_to_tsdata|> |
% <tsdata_to_var.html |tsdata_to_var|> |
% <var_to_autocov.html |var_to_autocov|> |
% <autocov_to_spwcgc.html |autocov_to_spwcgc|> |
% <sfreqs.html |sfreqs|>.
%
% (C) Lionel Barnett and Anil K. Seth, 2012. See file license.txt in
% installation directory for licensing terms.
%
%%

function fE = empirical_var_to_spwcgc(A,SIG,m,N,fres,H0,nsamps,mtrunc,decayfac,regmode,acmaxlags,acdectol)

if nargin <  8, mtrunc    = []; end % ensure default
if nargin <  9, decayfac  = []; end % ensure default
if nargin < 10, regmode   = []; end % ensure default
if nargin < 11, acmaxlags = []; end % ensure default
if nargin < 12, acdectol  = []; end % ensure default

[n,n1,p] = size(A);
assert(n1 == n,'VAR coefficients matrix has bad shape');

[nn1,nn2] = size(SIG);
assert(nn1 == nn2,'residuals covariance matrix not square');
assert(nn1 == n  ,'residuals covariance matrix doesn''t match VAR coefficients matrix');

h = fres+1;

fE = nan(nsamps,n,n,h);

if H0 % null hypothesis: no connectivity for each j -> i
        
    for j = 1:n
        for i = 1:n
            if i ~= j
                
                % enforce null hypothesis: no connectivity for j -> i

                A0 = A;
                A0(i,j,:) = 0;
               
                for s = 1:nsamps
                    fprintf('spectral PWCGC from node %d to node %d: empirical sample %d of %d',j,i,s,nsamps);
                                       
                    % generate empirical VAR time series

                    UE = var_to_tsdata(A0,SIG,m,N,mtrunc,decayfac);

                    % estimate empirical VAR parameters

                    [AE,SIGE] = tsdata_to_var(UE,p,regmode);
                    if isbad(AE), fprintf(' *** VAR estimation failed\n'); continue; end % something went badly wrong

                    % calculate empirical PWGC

                    [G,res] = var_to_autocov(AE,SIGE,acmaxlags,acdectol);
                    if res.error, fprintf(' *** bad VAR: %s\n',res.errmsg); continue; end

                    fE(s,i,j,:) = autocov_to_smvgc(G,i,j,fres);
                    
                    fprintf('\n');   
                end
            end
        end
        
    end
    
else
    
    for s = 1:nsamps
        fprintf('empirical spectral PWCGC %d of %d',s,nsamps);

        % generate empirical VAR time series

        UE = var_to_tsdata(A,SIG,m,N,mtrunc,decayfac);

        % estimate empirical VAR parameters

        [AE,SIGE] = tsdata_to_var(UE,p,regmode);
        if isbad(AE), fprintf(' *** VAR estimation failed\n'); continue; end % something went badly wrong

        % calculate empirical PWGC

        [G,res] = var_to_autocov(AE,SIGE,acmaxlags,acdectol);
        if res.error, fprintf(' *** bad VAR: %s\n',res.errmsg); continue; end

        fE(s,:,:,:) = autocov_to_spwcgc(G,fres);

        fprintf('\n');
    end
    
end
