%% plot_confints
%
% Confidence intervals plotting utility for pairwise-conditional causalities
%
% <matlab:open('plot_confints.m') code>
%
%% Syntax
%
%     plot_confints(F,FUP,FLO,FCRIT)
%
%% Arguments
%
% See also <mvgchelp.html#4 Common variable names and data structures>.
%
% _input_
%
%     F          matrix of pairwise-conditional Granger causalities
%     FUP        matrix of upper confidence bounds
%     FLO        matrix of lower confidence bounds
%     FCRIT      critical GC value (default: none)
%
%% Description
%
% Plot bar chart of pairwise-conditional time-domain causalities |F| (_cf._
% <autocov_to_pwcgc.html |autocov_to_pwcgc|>) with upper and lower
% confidence intervals |FUP|, |FLO| (_cf._  <mvgc_confint.html |mvgc_confint|>,
% <empirical_confint.html |empirical_confint|>). Optionally also plot a critical
% value |FCRIT| (_cf._  <mvgc_cval.html |mvgc_cval|>, <empirical_cval.html
% |empirical_cval|>).
%
%% See also
%
% <autocov_to_pwcgc.html |autocov_to_pwcgc|> |
% <mvgc_confint.html |mvgc_confint|> |
% <empirical_confint.html |empirical_confint|> |
% <mvgc_cval.html |mvgc_cval|> |
% <empirical_cval.html |empirical_cval|> |
% <mvgc_demo_confint.html |mvgc_demo_confint|>
%
% (C) Lionel Barnett and Anil K. Seth, 2012. See file license.txt in
% installation directory for licensing terms.
%
%%

function plot_confints(F,FUP,FLO,FCRIT)

if nargin < 4, FCRIT = []; end

n = size(F,1);
assert(size(F,2) == n,'must be square matrix');
assert(isequal(size(FUP),size(F)),'upper bounds matrix doesn''t match GC matrix');
assert(isequal(size(FLO),size(F)),'upper bounds matrix doesn''t match GC matrix');

r = n*(n-1);
o = ~isnan(F);
xtix = (1:r)';
xrange = [xtix(1)-1 xtix(end)+1];
yrange = [0 1.1*max(FUP(:))];
xtlabs = cell(r,1);
k = 0;
for i = 1:n
    for j = 1:n
        if i ~= j
            k = k+1;
            xtlabs{k} = sprintf('%d -> %d',i,j);
        end
    end
end

xlim(xrange);
ylim(yrange);
hold on
bar(xtix,F(o),'y');
errorbar(xtix,F(o),F(o)-FLO(o),FUP(o)-F(o),'.k');
if ~isempty(FCRIT)
    line(xlim,[FCRIT FCRIT],'Color','r');
end
hold off
set(gca,'TickLength',[0 0]);
set(gca,'XTick',xtix);
set(gca,'XTickLabel',xtlabs);
xlabel('connection')
if ~isempty(FCRIT)
    legend('pairwise-conditional GC','confidence interval','critical significance value');
else
    legend('pairwise-conditional GC','confidence interval');
end
