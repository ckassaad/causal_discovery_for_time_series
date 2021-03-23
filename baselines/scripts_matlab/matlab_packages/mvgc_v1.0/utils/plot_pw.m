%% plot_pw
%
% Plot pairwise quantities on a colourmapped grid
%
% <matlab:open('plot_pw.m') code>
%
%% Syntax
%
%     plot_pw(P,cm)
%
%% Arguments
%
% _input_
%
%     P          square matrix of pairwise quantities
%     cm         colour map (default: something soothing)
%
%% Description
%
% Plot pairwise quantities in |P|, a 2-dim square numerical matrix with
% first index representing target ("to") and second index source ("from")
% quantities, typically causalities, p-values, significances, etc. (see
% e.g. <autocov_to_pwcgc.html |autocov_to_pwcgc|>). Diagonal entries
% are ignored. A colormap |cm| may be supplied.
%
%% See also
%
% <autocov_to_pwcgc.html |autocov_to_pwcgc|> |
% <mvgc_demo.html |mvgc_demo|>
%
% (C) Lionel Barnett and Anil K. Seth, 2012. See file license.txt in
% installation directory for licensing terms.
%
%%

function plot_pw(P,cm)

if nargin < 2 || isempty(cm), cm =flipud(bone); end;

n = size(P,1);
assert(ismatrix(P) && size(P,2) == n,'input must be a square 2D matrix');

colormap(cm);
maxP = max(P(:));
if maxP == 0, maxP = 1; end
imagesc(P,[0 maxP]);
axis('square');
xlabel('from');
ylabel('to');
set(gca,'XTick',1:n);
set(gca,'XTickLabel',1:n);
set(gca,'YTick',1:n);
set(gca,'YTickLabel',1:n);
