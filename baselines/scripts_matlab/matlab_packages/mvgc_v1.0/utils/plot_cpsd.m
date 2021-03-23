%% plot_cpsd
%
% Cross-power spectral density plotting utility
%
% <matlab:open('plot_cpsd.m') code>
%
%% Syntax
%
%     plot_cpsd(S,leg,fs,frange,auto)
%
%% Arguments
%
% See also <mvgchelp.html#4 Common variable names and data structures>.
%
% _input_
%
%     S          matrix of cpsd matrices
%     leg        cell vector of legend strings matching S (default: 'cpsd 1', 'cpsd 2', etc.)
%     fs         sample rate in Hz (default: normalised freq as per routine 'sfreqs')
%     frange     frequency range to plot: empty for all (default) else an ascending 2-vector
%     auto       only plot auto- (not cross-) spectra (default: true)
%
%% Description
%
% Plots multiple cpsds vs. frequency on a grid. (|Sn|, |lstrn|) specify (cpsd,
% legend) pairs. If the |auto| flag is set only auto-spectra are plotted.
%
%% See also
%
% <var_to_cpsd.html |var_to_cpsd|> |
% <autocov_to_cpsd.html |autocov_to_cpsd|> |
% <tsdata_to_cpsd.html |tsdata_to_cpsd|> |
% <mvgc_demo_stats.html |mvgc_demo_stats|> |
% <sfreqs.html |sfreqs|>
%
% (C) Lionel Barnett and Anil K. Seth, 2012. See file license.txt in
% installation directory for licensing terms.
%
%%

function plot_cpsd(S,leg,fs,frange,auto)

[n,n1,h,N] = size(S);
assert(n1 == n,'cpsd matrix has bad shape');

if nargin < 2 || isempty(leg),  for i=1:N, leg{i} = sprintf('cpsd %d',i); end; end % default is 'cpsd i';
assert(isvector(leg) && iscellstr(leg),'legend must be a cell vector of strings');
assert(length(leg) == N,'legend does not match cpsd matrix');

if nargin < 3, fs = []; end % default to normalised frequency as per 'sfreqs'

if nargin < 4, frange = []; end; % default to all
if ~isempty(frange)
    assert(isvector(frange) && length(frange) == 2 && frange(1) < frange(2),'frequency range must be an ascending 2-vector of frequencies');
end

if nargin < 5 || isempty(auto), auto = true; end % default is auto-spectra only

fres = h-1;
lam = sfreqs(fres,fs)';
if ~isempty(frange)
    idx = lam >= frange(1) & lam <= frange(2);
    lam = lam(idx);
    S = S(:,:,idx,:);
end
xlims = [lam(1) lam(end)];
if isempty(fs), xlab = 'normalised frequency'; else xlab = 'frequency (Hz)'; end

if auto % auto-spectra only

    for i = 1:n
        subplot(n,1,i);
        plot(lam,squeeze(S(i,i,:,:)));
        xlim(xlims);
        xlabel(xlab);
        ylabel(sprintf('var %d auto-power',i));
        legend(leg);
    end

else    % auto- and cross-spectra

    k = 0;
    for i = 1:n
        for j = 1:n
            k = k+1;
            if i == j
                subplot(n,n,k);
                plot(lam,squeeze(S(i,i,:,:)));
                xlim(xlims);
                xlabel(xlab);
                ylabel(sprintf('var %d auto-power',i));
                legend(leg);
            elseif i < j
                subplot(n,n,k);
                plot(lam,abs(squeeze(S(i,j,:,:))));
                xlim(xlims);
                xlabel(xlab);
                ylabel(sprintf('var %d,%d cross-power',i,j));
                legend(leg);
            end
        end
    end

end
