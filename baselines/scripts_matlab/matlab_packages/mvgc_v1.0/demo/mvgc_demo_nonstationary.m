%% MVGC non-stationarity demo
%
% Demonstrates usage of the MVGC toolbox for non-stationary time series by
% "vertical" (windowed) regression on multi-trial data, for a minimal 2-variable
% VAR with linear trend and sinusoidally varying causal coefficient. Data is
% generated using the <var_to_tsdata_nonstat.html |var_to_tsdata_nonstat|>
% routine.
%
%% References
%
% [1] L. Barnett and A. K. Seth,
% <http://www.sciencedirect.com/science/article/pii/S0165027013003701 The MVGC
%     Multivariate Granger Causality Toolbox: A New Approach to Granger-causal
% Inference>, _J. Neurosci. Methods_ 223, 2014
% [ <matlab:open('mvgc_preprint.pdf') preprint> ].
%
% (C) Lionel Barnett and Anil K. Seth, 2012. See file license.txt in
% installation directory for licensing terms.
%
%% Parameters

ntrials   = 10;      % number of trials
nobs      = 1000;    % number of observations per trial
dnobs     = 200;     % initial observations to discard per trial

regmode   = 'OLS';   % VAR model estimation regression mode ('OLS', 'LWR' or empty for default

nlags     = 1;       % number of lags in VAR model
a         =  0.3;    % var 1 -> var 1 coefficient
b         = -0.6;    % var 2 -> var 2 coefficient
cmax      =  2.0;    % var 2 -> var 1 (causal) coefficient maximum

fs        = 1000;    % sample frequency
tr        = [5 -3];  % linear trend slopes
omega     = 2;       % "minimal VAR" X <- Y (causal) coefficient sinusoid frequency

ev        = 10;      % evaluate GC at every ev-th sample
wind      = 4;       % observation regression window size
morder    = 1;       % model order (for real-world data should be estimated)

alpha     = 0.05;    % significance level for statistical tests
mhtc      = 'SIDAK'; % multiple hypothesis test correction (see routine 'significance')

seed      = 0;       % random seed (0 for unseeded)

%% Construct multiple non-stationary VAR time series

rng_seed(seed);

nvars = 2;

tnobs = nobs+dnobs;                 % total observations per trial for time series generation
k = 1:tnobs;                        % vector of time steps
t = (k-dnobs-1)/fs;                 % vector of times
c = cmax*sin(2*pi*omega*t);         % causal coefficient varies sinusoidally
trend = diag(tr)*repmat(t,nvars,1); % linear trend

% set up time-varying VAR parameters

AT = zeros(nvars,nvars,nlags,tnobs); % "minimal VAR" coefficients
for j = k
    AT(:,:,nlags,j) = [a c(j); 0 b];
end
SIGT = eye(nvars);                   % "minimal VAR" residuals covariance

% generate non-stationary VAR

X = var_to_tsdata_nonstat(AT,SIGT,ntrials);

% add linear trend

X = X + repmat(trend,[1 1 ntrials]);

% discard initial observations

X = X(:,dnobs+1:tnobs,:);
c = c(dnobs+1:tnobs);
trend = trend(:,dnobs+1:tnobs);
k = 1:nobs;

figure(1); clf;
xrange = [0,nobs];

% plot causal coeffcient

subplot(3,1,1);
plot(k,c);
xlim(xrange);
title('Causal coefficient sinusoid');
xlabel('time');

% plot first 10 trial time series and trend for variable 1

subplot(3,1,2)
plot(k,squeeze(X(1,:,1:min(10,nobs)))');
hold on
plot(k,trend(1,:),'k');
hold off
xlim(xrange);
title('Time series (variable 1)');
xlabel('time');

%% "Vertical" regression GC calculation

wnobs = morder+wind;   % number of observations in "vertical slice"
ek    = wnobs:ev:nobs; % GC evaluation points
enobs = length(ek);    % number of GC evaluations

F12 = nan(enobs,1);
F21 = nan(enobs,1);

% loop through evaluation points

for e = 1:enobs
    j = ek(e);
    fprintf('window %d of %d at time = %d',e,enobs,j);

    [A,SIG] = tsdata_to_var(X(:,j-wnobs+1:j,:),morder,regmode);
    if isbad(A)
        fprintf(2,' *** skipping - VAR estimation failed\n');
        continue
    end

    [G,info] = var_to_autocov(A,SIG);
    if info.error
        fprintf(2,' *** skipping - bad VAR (%s)\n',info.errmsg);
        continue
    end
    if info.aclags < info.acminlags % warn if number of autocov lags is too small (not a show-stopper)
        fprintf(2,' *** WARNING: minimum %d lags required (decay factor = %e)',info.acminlags,realpow(info.rho,info.aclags));
    end

    FF = autocov_to_pwcgc(G);
    if isbad(FF,false)
        fprintf(2,' *** skipping - GC calculation failed\n');
        continue
    end

    F12(e) = FF(1,2); % estimated GC 2 -> 1 (significant)
    F21(e) = FF(2,1); % estimated GC 1 -> 2 (non-significant)

    fprintf('\n');
end

% theoretical GC 2 -> 1

D = 1+b^2 + c(ek).^2;
F12T = log((D + sqrt(D.^2 - 4*b^2))/2)';

% critical GC value at significance alpha, corrected for multiple hypotheses

nhyp = 2; % number of hypotheses (i.e. 2 -> 1 and 1 -> 2)
switch upper(mhtc)
    case 'NONE',       alpha1 = alpha;
    case 'BONFERRONI', alpha1 = alpha/nhyp;
    case 'SIDAK',      alpha1 = 1-realpow(1-alpha,1/nhyp);
    otherwise, error('unhandled correction method ''%s''',mhtc);
end

Fc = mvgc_cval(alpha1,morder,wnobs,ntrials,1,1,nvars-2);

% plot GCs

subplot(3,1,3);
plot(ek',[F12T F12]);
xlim(xrange);
line(xrange,[Fc Fc],'Color','r'); % critical significance value
legend('theoretical','estimated','cval');
title('variable 2 -> variable 1 causality');
xlabel('time');

%%
% <mvgc_demo_nonstationary.html back to top>
