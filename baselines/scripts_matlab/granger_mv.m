%% MVGC
%
% Demonstrates typical usage of the MVGC toolbox on generated VAR data for a
% 5-node network with known causal structure (see <var5_test.html |var5_test|>).
% Estimates a VAR model and calculates time- and frequency-domain
% pairwise-conditional Granger causalities (also known as the "causal graph").
% Also calculates Seth's causal density measure [2].
%
%
%% References
%
% [1] L. Barnett and A. K. Seth,
% <http://www.sciencedirect.com/science/article/pii/S0165027013003701 The MVGC
%     Multivariate Granger Causality Toolbox: A New Approach to Granger-causal
% Inference>, _J. Neurosci. Methods_ 223, 2014
% [ <matlab:open('mvgc_preprint.pdf') preprint> ].
%
%% Parameters
function main(path_data, path_sig_level, path_n_lags, dir_path)
plot_bool = false;
%path_data = strcat('./args/', data_name,'.csv');
%path_n_lags = strcat('./args/', n_lags_name,'.txt');
%path_sig_level = strcat('./args/', sig_level_name,'.txt');

tmp = fileread(path_n_lags);
nlags = str2double(tmp);
tmp = fileread(path_sig_level);
sig_level = str2double(tmp);

regmode   = 'OLS';  % VAR model estimation regression mode ('OLS', 'LWR' or empty for default)
icregmode = 'LWR';  % information criteria regression mode ('OLS', 'LWR' or empty for default)

morder    = 'AIC';  % model order to use ('actual', 'AIC', 'BIC' or supplied numerical value)
momax     = nlags;     % maximum model order for model order estimation

acmaxlags = 1000;   % maximum autocovariance lags (empty for automatic calculation)

tstat     = '';     % statistical test for MVGC:  'F' for Granger's F-test (default) or 'chi2' for Geweke's chi2 test
alpha     = sig_level;   % significance level for significance test
mhtc      = 'FDR';  % multiple hypothesis test correction (see routine 'significance')

fs        = 200;    % sample rate (Hz)

seed      = 0;      % random seed (0 for unseeded)

%% Generate VAR test data (<mvgc_schema.html#3 |A3|>)
%
% _*Note:*_ This is where you would read in your own time series data; it should
% be assigned to the variable |X| (see below and <mvgchelp.html#4 Common
% variable names and data structures>).

% Startup: call needed packages

run( './matlab_packages/mvgc_v1.0/startup');

% Seed random number generator.

rng_seed(seed);

% get data
X = readtable(path_data);
names = X.Properties.VariableNames;
X = X{:, :};
X = transpose(X);
%X = fopen(path,'w');
nvars = size(X,1); % number of variables
ntrials   = 1;     % number of trials
nobs      = 1000;   % number of observations per trial

%% Model order estimation (<mvgc_schema.html#3 |A2|>)

% Calculate information criteria up to specified maximum model order.

ptic('\n*** tsdata_to_infocrit\n');
[AIC,BIC,moAIC,moBIC] = tsdata_to_infocrit(X,momax,icregmode);
ptoc('*** tsdata_to_infocrit took ');

% Plot information criteria.
if plot_bool
    figure(1); clf;
    plot_tsdata([AIC BIC]',{'AIC','BIC'},1/fs);
    title('Model order estimation');
end

fprintf('\nbest model order (AIC) = %d\n',moAIC);
fprintf('best model order (BIC) = %d\n',moBIC);

% Select model order.

if strcmpi(morder,'AIC')
    morder = moAIC;
    fprintf('\nusing AIC best model order = %d\n',morder);
elseif strcmpi(morder,'BIC')
    morder = moBIC;
    fprintf('\nusing BIC best model order = %d\n',morder);
else
    fprintf('\nusing specified model order = %d\n',morder);
end

%% VAR model estimation (<mvgc_schema.html#3 |A2|>)

% Estimate VAR model of selected order from data.

ptic('\n*** tsdata_to_var... ');
[A,SIG] = tsdata_to_var(X,morder,regmode);
ptoc;

% Check for failed regression

assert(~isbad(A),'VAR estimation failed');

% NOTE: at this point we have a model and are finished with the data! - all
% subsequent calculations work from the estimated VAR parameters A and SIG.

%% Autocovariance calculation (<mvgc_schema.html#3 |A5|>)

% The autocovariance sequence drives many Granger causality calculations (see
% next section). Now we calculate the autocovariance sequence G according to the
% VAR model, to as many lags as it takes to decay to below the numerical
% tolerance level, or to acmaxlags lags if specified (i.e. non-empty).

ptic('*** var_to_autocov... ');
[G,info] = var_to_autocov(A,SIG,acmaxlags);
ptoc;

% The above routine does a LOT of error checking and issues useful diagnostics.
% If there are problems with your data (e.g. non-stationarity, colinearity,
% etc.) there's a good chance it'll show up at this point - and the diagnostics
% may supply useful information as to what went wrong. It is thus essential to
% report and check for errors here.

var_info(info,true); % report results (and bail out on error)

%% Granger causality calculation: time domain  (<mvgc_schema.html#3 |A13|>)

% Calculate time-domain pairwise-conditional causalities - this just requires
% the autocovariance sequence.

ptic('*** autocov_to_pwcgc... ');
F = autocov_to_pwcgc(G);
ptoc;

% Check for failed GC calculation

assert(~isbad(F,false),'GC calculation failed');

% Significance test using theoretical null distribution, adjusting for multiple
% hypotheses.

pval = mvgc_pval(F,morder,nobs,ntrials,1,1,nvars-2,tstat); % take careful note of arguments!
sig  = significance(pval,alpha,mhtc);

% Plot time-domain causal graph, p-values and significance.

if plot_bool
    figure(2); clf;
    subplot(1,3,1);
    plot_pw(F);
    title('Pairwise-conditional GC');
    subplot(1,3,2);
    plot_pw(pval);
    title('p-values');
    subplot(1,3,3);
    plot_pw(sig);
    title(['Significant at p = ' num2str(alpha)])
end

sig = transpose(sig);
[row,col] = find(sig==1);
for i = 1:length(row)
    sig(row(i), col(i)) = 2;
    sig(col(i), row(i)) = 1;
end
sig(isnan(sig))=1;
sig = array2table(sig,'VariableNames',names);
fprintf(dir_path+'/results/result.txt');
writetable(sig, './results/result.txt');
end