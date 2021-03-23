%% empirical_var_to_smvgc
%
% Calculate sampling distribution for conditional frequency-domain MVGC from
% generated time series data for a specified VAR model
%
% <matlab:open('empirical_var_to_smvgc.m') code>
%
%% Syntax
%
%     fE = empirical_var_to_smvgc(A,SIG,m,N,x,y,fres,H0,nsamps,mtrunc,decayfac,regmode,acmaxlags,acdectol)
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
%     x          vector of indices of target (causee) multi-variable
%     y          vector of indices of source (causal) multi-variable
%     fres       frequency resolution (default: automatic)
%     H0         flag: impose null hypothesis of zero connectivity?
%     nsamps     number of bootstrap samples
%     mtrunc     number of initial time observations to truncate  (default as for 'var_to_tsdata')
%     decayfac   initial transients decay factor (default as for 'var_to_tsdata')
%     regmode    regression mode (default as for 'tsdata_to_var')
%     acmaxlags  maximum autocovariance lags (default as for 'var_to_autocov')
%     acdectol   autocovariance decay tolerance (default as for 'var_to_autocov')
%
% _output_
%
%     fE         empirical spectral Granger causality distribution
%
%% Description
%
% Returns |nsamps| samples from the empirical sampling distribution of the
% frequency-domain MVGC from the variable |Y| (specified by the vector of indices
% |y|) to the variable |X| (specified by the vector of indices |x|), conditional
% on all other variables in time series data generated from the VAR model
% specified by the coefficients matrix |A| and residuals covariance matrix
% |SIG|. The time series comprise |N| trials of length |m|, and residuals are
% iid Gaussian with covariance matrix |SIG|. If the flag |H0| is set then data
% is generated for the nested null model with zero connectivity from |Y| to |X|;
% i.e. |A(x,y,k)| is set to zero for all lags |k|. For other parameters see
% <var_to_tsdata.html |var_to_tsdata|>, <tsdata_to_var.html |tsdata_to_var|> and
% <var_to_autocov.html |var_to_autocov|>.
%
% The first dimension of the returned matrix |fE| indexes samples, the second
% frequency.
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
% <empirical_var_to_spwcgc.html |empirical_var_to_spwcgc|> |
% <var_to_tsdata.html |var_to_tsdata|> |
% <tsdata_to_var.html |tsdata_to_var|> |
% <var_to_autocov.html |var_to_autocov|> |
% <autocov_to_smvgc.html |autocov_to_smvgc|> |
% <sfreqs.html |sfreqs|>.
%
% (C) Lionel Barnett and Anil K. Seth, 2012. See file license.txt in
% installation directory for licensing terms.
%
%%

function fE = empirical_var_to_smvgc(A,SIG,m,N,x,y,fres,H0,nsamps,mtrunc,decayfac,regmode,acmaxlags,acdectol)

if nargin < 10, mtrunc    = []; end % ensure default
if nargin < 11, decayfac  = []; end % ensure default
if nargin < 12, regmode   = []; end % ensure default
if nargin < 13, acmaxlags = []; end % ensure default
if nargin < 14, acdectol  = []; end % ensure default

[n,n1,p] = size(A);
assert(n1 == n,'VAR coefficients matrix has bad shape');

[nn1,nn2] = size(SIG);
assert(nn1 == nn2,'residuals covariance matrix not square');
assert(nn1 == n  ,'residuals covariance matrix doesn''t match VAR coefficients matrix');

x = x(:)'; % vectorise
y = y(:)'; % vectorise
assert(all(x >=1 & x <= n),     'some x indices out of range');
assert(all(y >=1 & y <= n),     'some y indices out of range');
assert(isempty(intersect(x,y)), 'x and y indices must be distinct');

h = fres+1;

if H0 % null hypothesis: no connectivity y -> x
    A(x,y,:) = 0;
end

fE = nan(nsamps,h);

for s = 1:nsamps
    fprintf('spectral MVGC: empirical sample %d of %d',s,nsamps);
    
    % generate empirical VAR time series
    
    UE = var_to_tsdata(A,SIG,m,N,mtrunc,decayfac);
    
    % estimate empirical VAR parameters
    
    [AE,SIGE] = tsdata_to_var(UE,p,regmode);
    if isbad(AE), fprintf(' *** VAR estimation failed\n'); continue; end % something went badly wrong
    
    % calculate empirical MVGC
    
    [G,res] = var_to_autocov(AE,SIGE,acmaxlags,acdectol);
    if res.error, fprintf(' *** bad VAR: %s\n',res.errmsg); continue; end
    fE(s,:) = autocov_to_smvgc(G,x,y,fres);
    
    fprintf('\n');
end
