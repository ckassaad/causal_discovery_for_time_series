%% MVGC bootstrap demo
%
% Demonstrates confidence interval construction using a nonparametric bootstrap
% on generated VAR data for a 5-node network with known causal structure (see
% <var5_test.html |var5_test|>). Pairwise-conditional Granger causalities are
% estimated and confidence intervals constructed using both the theoretical and
% bootstrap distributions.
%
%% References
%
% [1] L. Barnett and A. K. Seth,
% <http://www.sciencedirect.com/science/article/pii/S0165027013003701 The MVGC
%     Multivariate Granger Causality Toolbox: A New Approach to Granger-causal
% Inference>, _J. Neurosci. Methods_ 223, 2014
% [ <matlab:open('mvgc_preprint.pdf') preprint> ].
%
% [2] D. A. Freedman, Bootstrapping regression models, _Ann. Stats._, 9(6), 1981.
%
% (C) Lionel Barnett and Anil K. Seth, 2012. See file license.txt in
% installation directory for licensing terms.
%
%% Parameters

ntrials   = 10;     % number of trials
nobs      = 100;    % number of observations per trial
nsamps    = 100;    % number of bootstrap samples

regmode   = 'OLS';  % VAR model estimation regression mode ('OLS', 'LWR' or empty for default)

acmaxlags = 1000;   % maximum autocovariance lags (empty for automatic calculation)

tstat     = 'chi2'; % statistical test for MVGC: 'F' for Granger's F-test, 'chi2' for Geweke's chi2 test or leave empty for default
alpha     = 0.05;   % significance level for all statistical tests

seed      = 0;      % random seed (0 for unseeded)

%% Generate VAR data

% Seed random number generator.

rng_seed(seed);

% Get VAR coefficients for 5-node test network.

AT = var5_test;
nvars = size(AT,1);  % number of variables

% Residuals covariance matrix.

SIGT = eye(nvars);

% Generate VAR time series data with normally distributed residuals for
% specified coefficients and covariance matrix.

ptic('\n*** var_to_tsdata... ');
X = var_to_tsdata(AT,SIGT,nobs,ntrials);
ptoc;

%% VAR model estimation and autocovariance calculation

morder = size(AT,3); % actual model order; on real data - i.e. with no generative model
                     % available - use information criteria to estimate (see 'mvgc_demo')
                     
% Calculate VAR model

ptic('*** tsdata_to_var... ');
[A,SIG] = tsdata_to_var(X,morder,regmode);
ptoc;

% Check for failed regression

assert(~isbad(A),'VAR estimation failed');

% Now calculate autocovariance according to the VAR model, to as many lags
% as it takes to decay to below the numerical tolerance level, or to acmaxlags
% lags if specified (i.e. non-empty).

ptic('*** var_to_autocov... ');
[G,res] = var_to_autocov(A,SIG,acmaxlags);
ptoc;

% Report and check for errors.

fprintf('\nVAR check:\n'); disp(res); % report results...
assert(~res.error,'bad VAR');         % ...and bail out if there are errors

%% Granger causality estimation

% Time-domain pairwise conditional causalities

ptic('*** autocov_to_pwcgc... ');
F = autocov_to_pwcgc(G);
ptoc;

% Theoretical confidence intervals.

[FTUP,FTLO] = mvgc_confint(alpha,F,morder,nobs,ntrials,1,1,nvars-2,tstat);

% Critical GC value.

FTCRIT = mvgc_cval(alpha,morder,nobs,ntrials,1,1,nvars-2,tstat);

%% Bootstrap
   
ptic('\n*** bootstrap_tsdata_to_pwcgc\n');
FSAMP = bootstrap_tsdata_to_pwcgc(X,morder,nsamps);
ptoc('*** bootstrap_tsdata_to_pwcgc took ',[],1);

% (We should really check for failed bootstrap estimates here.)

% Bootstrap (empirical) confidence intervals.

[FSUP,FSLO] = empirical_confint(alpha,FSAMP);

% Note: we haven't calculated a bootstrap critical GC value; for this we would
% require an (empirical) null GC distribution, which we could obtain by running
% a permutation test (see e.g. 'permtest_tsdata_to_pwcgc').

%% Plot PWCGC estimates with theoretical and bootstrap confidence intervals

figure(1); clf
subplot(2,1,1);
plot_confints(F,FTUP,FTLO,FTCRIT);
title(sprintf('Theoretical distribution\nconfidence intervals at alpha = %g',alpha));
subplot(2,1,2);
plot_confints(F,FSUP,FSLO);
title(sprintf('Bootstrap distribution\nconfidence intervals at alpha = %g',alpha));

%%
% <mvgc_demo_bootstrap.html back to top>
