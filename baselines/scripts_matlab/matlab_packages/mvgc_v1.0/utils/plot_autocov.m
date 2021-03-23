%% plot_autocov
%
% Autocovariance plotting utility
%
% <matlab:open('plot_autocov.m') code>
%
%% Syntax
%
%     plot_autocov(G,leg,dt,trange,auto,acorr)
%
%% Arguments
%
% See also <mvgchelp.html#4 Common variable names and data structures>.
%
% _input_
%
%     G          matrix of autocovariance matrices
%     leg        cell vector of legend strings matching G (default: 'autocov 1', 'autocov 2', etc.)
%     dt         sample time step in seconds or empty (default) for number of lags
%     trange     time range to plot: empty for all (default) else an ascending 2-vector
%     auto       only plot auto- (not cross-) spectra (default: true)
%     acorr      plot autocorrelation rather than autocovariance (default: false)
%
%% Description
%
% Plots autocovariance sequence |G| against lags for each pair of variables on a
% grid. If the |acorr| flag is nonzero, autocorrelation rather than
% autocovariance is plotted (see <cov2corr.html |cov2corr|>). If |acorr| is
% negative, then zero-lag correlation is omitted (this may be useful since
% lag-zero autocorrelation can be an order of magnitude larger than
% autocorrelation at higher lags - e.g. for a near-white noise process).
%
%% See also
%
% <var_to_autocov.html |autocov_to_pwcgc|> |
% <tsdata_to_autocov.html |tsdata_to_autocov|> |
% <cpsd_to_autocov.html |cpsd_to_autocov|> |
% <cov2corr.html |cov2corr|>
%
% (C) Lionel Barnett and Anil K. Seth, 2012. See file license.txt in
% installation directory for licensing terms.
%
%%

function plot_autocov(G,leg,dt,trange,auto,acorr)

[n,n1,q1,N] = size(G);
assert(n1 == n,'autocovariance matrix has bad shape');

if nargin < 2 || isempty(leg), for i=1:N, leg{i} = sprintf('autocov %d',i); end; end % default is 'autocov i';
assert(isvector(leg) && iscellstr(leg),'legend must be a cell vector of strings');
assert(length(leg) == N,'legend does not match autocovariance matrix');

plags = nargin < 3 || isempty(dt); % time = number of lags 
if plags, dt = 1; end

if nargin < 4, trange = []; end; % default to all
if ~isempty(trange)
    assert(isvector(trange) && length(trange) == 2 && trange(1) < trange(2),'frequency range must be an ascending 2-vector of times');
end

if nargin < 5 || isempty(auto),  auto  = true;  end % default is autocovariances only (not cross-covariances)
if nargin < 6 || isempty(acorr), acorr = false; end % default is autocovariance, not autocorrelation

q = q1-1;
tvec = (0:q)*dt;
if ~isempty(trange)
    idx = tvec >= trange(1) & tvec <= trange(2);
    tvec = tvec(idx);
    G = G(:,:,idx,:);
end

xlims = [tvec(1) tvec(end)];
if plags, xlab = 'lags'; else xlab = 'lag (secs)'; end

if acorr, for s = 1:N, G(:,:,:,s) = cov2corr(G(:,:,:,s)); end; end

if auto % covariance only

    for i = 1:n
        subplot(n,1,i);
        plot(tvec,squeeze(G(i,i,:,:)));
        if acorr
            ylim([-1 1]);
            ylabel(sprintf('var %d autocorrelation',i));
        else
            ylabel(sprintf('var %d autocovariance',i));
        end
        xlim(xlims);
        xlabel(xlab);
        legend(leg);
    end

else    % auto- and cross-covariances

    k = 0;
    for i = 1:n
        for j = 1:n
            k = k+1;
            subplot(n,n,k);
            plot(tvec,squeeze(G(i,j,:,:)));
            if acorr
                ylim([-1 1]);
                ylabel(sprintf('var %d,%d autocorrelation',i,j));
            else
                ylabel(sprintf('var %d,%d autocovariance',i,j));
            end
            xlim(xlims);
            xlabel(xlab);
            legend(leg);
        end
    end

end
