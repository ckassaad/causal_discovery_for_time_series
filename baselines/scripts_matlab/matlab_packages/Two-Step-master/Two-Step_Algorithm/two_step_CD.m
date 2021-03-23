function [B,W_m] = two_step_CD(X,lambda)
% function [B,W_m] = two_step_CD(X)
% Two-step method for linear causal discovery that allows cycles and
% confounders
% Input: 
%   Data matrix X (variable number * sample size).
% Output: 
%   B: the causal influence matrix X = BX + E;
%   in B causal direction goes from column to row, the entry Bij, means
%   Xj -> Xi
%   W_m: the ICA de-mixing matrix.
% by Kun Zhang 2016,2017 Carnegie Mellon University

[N,T] = size(X);

% estimate the mask using adaptative lasso
Mask = zeros(N,N);
for i=1:N
    if T<4*N % sample size too small, so preselect the features
        tmp1 = X([1:i-1 i+1:N],:);   %get X\xi
        [tmp2, Ind_t] = sort(abs(corr( tmp1', X(i,:)' )), 'descend');  %compute the correlation of xi with X\xi and sort from larger in absolute value, get values and indices
        X_sel = tmp1(Ind_t(1:floor(N/4)),:); % pre-select N/4 features, the ones more correlated to xi
        [beta_alt, beta_new_nt, beta2_alt, beta2_new_nt] = betaAlasso_grad_2step(X_sel, X(i,:), 0.65^2*var(X(i,:)), log(T)/2); 
        beta2_al = zeros(N-1,1);
        beta2_al(Ind_t(1:floor(N/4))) = beta2_alt;
    else
        [beta_al, beta_new_n, beta2_al, beta2_new_n] = betaAlasso_grad_2step(X([1:i-1 i+1:N],:), X(i,:), 0.65^2*var(X(i,:)), log(T)/2); 
    end
    Mask(i,1:i-1) = abs(beta2_al(1:i-1)) >0.01;
    Mask(i,i+1:N) = abs(beta2_al(i:N-1)) >0.01;
end
Mask = Mask + Mask';
Mask = (Mask~=0);

% perform constrained_ICA
[y_m, W_m, WW_m, Score] = sparseica_W_adasize_Alasso_mask(log(T)*lambda, Mask , X);
B = eye(N) - W_m;

