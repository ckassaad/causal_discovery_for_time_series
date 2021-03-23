%% empirical_var_to_pwcgc
%
% Calculate sampling distribution for pairwise-conditional time-domain MVGCs
% from generated time series data for a specified VAR model
%
% <matlab:open('empirical_var_to_pwcgc.m') code>
%
%% Syntax
%
%     FE = empirical_var_to_pwcgc(A,SIG,m,N,H0,nsamps,mtrunc,decayfac,regmode,acmaxlags,acdectol)
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
%     FE         empirical Granger causalities distributions
%
%% Description
%
% Returns |nsamps| samples from the empirical sampling distribution of the
% pairwise-conditional time-domain MVGCs for time series data generated from the
% VAR model specified by the coefficients matrix |A| and residuals covariance
% matrix |SIG|. The time series comprise |N| trials of length |m|, and residuals
% are iid Gaussian with covariance matrix |SIG|. If the flag |H0| is set then
% data is generated for the nested null models with zero connectivity; i.e. for
% each target index |i| and source index |j|, |A(i,j,k)| is set to zero for all
% lags |k|. For other parameters see <var_to_tsdata.html |var_to_tsdata|>,
% <tsdata_to_var.html |tsdata_to_var|> and <var_to_autocov.html
% |var_to_autocov|>.
%
% The first dimension of the returned matrix |FE| indexes samples, the second
% indexes the target (causee) variable and the third the source (causal)
% variable.
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
% <empirical_var_to_smvgc.html |empirical_var_to_smvgc|> |
% <empirical_var_to_spwcgc.html |empirical_var_to_spwcgc|> |
% <var_to_tsdata.html |var_to_tsdata|> |
% <tsdata_to_var.html |tsdata_to_var|> |
% <var_to_autocov.html |var_to_autocov|> |
% <autocov_to_pwcgc.html |autocov_to_pwcgc|>.
%
% (C) Lionel Barnett and Anil K. Seth, 2012. See file license.txt in
% installation directory for licensing terms.
%
%%

function FE = empirical_var_to_pwcgc(A,SIG,m,N,H0,nsamps,mtrunc,decayfac,regmode,acmaxlags,acdectol)

if nargin <  7, mtrunc    = []; end % ensure default
if nargin <  8, decayfac  = []; end % ensure default
if nargin <  9, regmode   = []; end % ensure default
if nargin < 10, acmaxlags = []; end % ensure default
if nargin < 11, acdectol  = []; end % ensure default

[n,n1,p] = size(A);
assert(n1 == n,'VAR coefficients matrix has bad shape');

[nn1,nn2] = size(SIG);
assert(nn1 == nn2,'residuals covariance matrix not square');
assert(nn1 == n  ,'residuals covariance matrix doesn''t match VAR coefficients matrix');

FE = nan(nsamps,n,n);

if H0 % null hypothesis: no connectivity for each j -> i
        
    for j = 1:n
        for i = 1:n
            if i ~= j
                
                % enforce null hypothesis: no connectivity for j -> i

                A0 = A;
                A0(i,j,:) = 0;
               
                for s = 1:nsamps
                    fprintf('PWCGC from node %d to node %d: empirical sample %d of %d',j,i,s,nsamps);
                                       
                    % generate empirical VAR time series

                    UE = var_to_tsdata(A0,SIG,m,N,mtrunc,decayfac);

                    % estimate empirical VAR parameters

                    [AE,SIGE] = tsdata_to_var(UE,p,regmode);
                    if isbad(AE), fprintf(' *** VAR estimation failed\n'); continue; end % something went badly wrong

                    % calculate empirical PWGC

                    [G,res] = var_to_autocov(AE,SIGE,acmaxlags,acdectol);
                    if res.error, fprintf(' *** bad VAR: %s\n',res.errmsg); continue; end

                    FE(s,i,j) = autocov_to_mvgc(G,i,j);
                    
                    fprintf('\n');   
                end
            end
        end
        
    end
    
else
    
    for s = 1:nsamps
        fprintf('empirical PWCGC %d of %d',s,nsamps);

        % generate empirical VAR time series

        UE = var_to_tsdata(A,SIG,m,N,mtrunc,decayfac);

        % estimate empirical VAR parameters

        [AE,SIGE] = tsdata_to_var(UE,p,regmode);
        if isbad(AE), fprintf(' *** VAR estimation failed\n'); continue; end % something went badly wrong

        % calculate empirical PWGC

        [G,res] = var_to_autocov(AE,SIGE,acmaxlags,acdectol);
        if res.error, fprintf(' *** bad VAR: %s\n',res.errmsg); continue; end

        FE(s,:,:) = autocov_to_pwcgc(G);

        fprintf('\n');
    end
    
end
