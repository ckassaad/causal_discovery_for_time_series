%% tsdata_to_autocov_debias (EXPERIMENTAL)
%
% Calculate sample autocovariance sequence from time series data with debias
%
% <matlab:open('tsdata_to_autocov_debias.m') code>
%
%% Syntax
%
%     G = tsdata_to_autocov_debias(X,q)
%
%% Arguments
%
% See also <mvgchelp.html#4 Common variable names and data structures>.
%
% _input_
%
%     X          multi-trial time series data
%     q          number of lags
%
% _output_
%
%     G          sample autocovariance sequence
%
%% Description
%
% Returns |q|-lag sample autocovariance sequence |G| defined as
% [[ii_acseq.png]] for the (presumed stationary) multivariate process |X|.
% |X| may contain single- or multi-trial time series data. See
% <tsdata_to_autocov.html |tsdata_to_autocov|> for more details.
%
% This (experimental) algorithm [2] attempts to reduce bias due to mean estimation
% for small samples. 
%
%% References
%
% [1] L. Barnett and A. K. Seth,
% <http://www.sciencedirect.com/science/article/pii/S0165027013003701 The MVGC
%     Multivariate Granger Causality Toolbox: A New Approach to Granger-causal
% Inference>, _J. Neurosci. Methods_ 223, 2014
% [ <matlab:open('mvgc_preprint.pdf') preprint> ].
%
% [2] Y. Shkolnisky, F. J. Sigworth and A. Singer,
% <http://www.journalogy.org/Publication/5978390/a-note-on-estimating-autocovariance-from-short-time-observations
% A note on estimating autocovariance from short-time observations>, 2008.
%
%% See also
%
% <demean.html |demean|> |
% <tsdata_to_var.html |tsdata_to_var|> |
% <var_to_autocov.html |var_to_autocov|> |
% <tsdata_to_autocov.html |tsdata_to_autocov|>
%
% (C) Lionel Barnett and Anil K. Seth, 2012. See file license.txt in
% installation directory for licensing terms.
%
%%

function G = tsdata_to_autocov_debias(X,q)

[n,m,N] = size(X);

assert(q < m-1,'too many lags');

if N > 1 % multi-trial

    Nm = N*m;
    M = sum(X(:,:),2);
    B = (M*M'-X(:,:)*X(:,:)')/(Nm*(Nm-1)); % debiasing term
    G = zeros(n,n,q+1);
    for k=0:q
        Nmk = Nm-N*k;
        for r = 1:N
            G(:,:,k+1) = G(:,:,k+1) + (X(:,k+1:m,r)*X(:,1:m-k,r)')/Nmk;
        end
        G(:,:,k+1) = G(:,:,k+1)-B;
    end

else

    M = sum(X,2);
    B = (M*M'-X*X')/(m*(m-1)); % debiasing term
    G = zeros(n,n,q+1);
    for k=0:q
        G(:,:,k+1) = (X(:,k+1:m)*X(:,1:m-k)')/(m-k)-B;
    end

end
