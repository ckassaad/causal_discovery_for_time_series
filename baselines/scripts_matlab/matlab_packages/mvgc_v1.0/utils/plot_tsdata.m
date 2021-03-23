%% plot_tsdata
%
% Time series data plotting utility
%
% <matlab:open('plot_tsdata.m') code>
%
%% Syntax
%
%     plot_tsdata(X,leg,dt,trange)
%
%% Arguments
%
% See also <mvgchelp.html#4 Common variable names and data structures>.
%
% _input_
%
%     X          multi-trial time series data
%     leg        cell vector of legend strings matching G (default: 'series 1', 'series 2', etc.)
%     dt         sample time step in seconds or empty (default) for number of time steps
%     trange     time range to plot: empty for all (default) else an ascending 2-vector
%
%% Description
%
% Plot time series data |X| in a given time range |trange|.
%
%% See also
%
% <mvgc_demo.html |mvgc_demo|>
%
% (C) Lionel Barnett and Anil K. Seth, 2012. See file license.txt in
% installation directory for licensing terms.
%
%%

function plot_tsdata(X,leg,dt,trange)

[n,m,N] = size(X);

if nargin < 2 || isempty(leg), for i=1:n, leg{i} = sprintf('series %d',i); end; end % default is 'autocov i';
assert(isvector(leg) && iscellstr(leg),'legend must be a cell vector of strings');
assert(length(leg) == n,'legend does not match time series matrix');

pnts = nargin < 3 || isempty(dt); % time = number of time steps 
if pnts, dt = 1; end

if nargin < 4, trange = []; end; % default to all
if ~isempty(trange)
    assert(isvector(trange) && length(trange) == 2 && trange(1) < trange(2),'frequency range must be an ascending 2-vector of times');
end

m1 = m-1;
tvec = (0:m1)*dt;
if ~isempty(trange)
    idx = tvec >= trange(1) & tvec <= trange(2);
    tvec = tvec(idx);
    X = X(:,idx,:);
end

xlims = [tvec(1) tvec(end)];
if pnts, xlab = 'time steps'; else xlab = 'time (secs)'; end

for r = 1:N
    subplot(N,1,r);
    plot(tvec,X(:,:,r));
    xlim(xlims);
    xlabel(xlab);
    if N > 1, ylabel(sprintf('trial %d',r)); end % multi-trial
    legend(leg);
end
