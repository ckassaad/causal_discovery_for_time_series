%% plot_varcoeffs
%
% VAR coefficients plotting utility
%
% <matlab:open('plot_varcoeffs.m') code>
%
%% Syntax
%
%     plot_varcoeffs(A)
%
%% Arguments
%
% See also <mvgchelp.html#4 Common variable names and data structures>.
%
% _input_
%
%     A          VAR coefficients matrix
%
%% Description
%
% Plots VAR coefficients |A| against lags for each pair of variables on a grid.
%
%% See also
%
% <tsdata_to_var.html |tsdata_to_var|>
%
% (C) Lionel Barnett and Anil K. Seth, 2012. See file license.txt in
% installation directory for licensing terms.
%
%%

function plot_varcoeffs(A)

[n,n1,q] = size(A);
assert(n1 == n,'VAR coefficients matrix has bad shape');

t = (1:q)';
xlims = [1 q];
ylims = [min(A(:)) max(A(:))];

k = 0;
for i = 1:n
    for j = 1:n
        k = k+1;
        subplot(n,n,k);
        plot(t,squeeze(A(i,j,:)));
        grid on;
        xlabel('lags');
        ylabel(sprintf('var %d,%d',i,j));
        xlim(xlims);
        ylim(ylims);
    end
end
