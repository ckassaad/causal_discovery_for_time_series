function [beta_al, beta_new_n, beta2_al, beta2_new_n] = betaAlasso_grad_2step(x,y,var_noise,lambda)
% function [beta_al, beta_new_n, beta2_al, beta2_new_n] =...
% betaAlasso_2step(x,y,var_noise,lambda)
% Aim: 
%       to find the solution of adaptive Lasso with the given lambda
% Inputs: 
%       x is p*n, and y is 1*n. var_noise is the variance of the noise.
% Outputs: 
%       beta_al contains the obtained beta after step 1 (conventioanl
% ALasso). After convergence of step 1, we repeat upadte \hat{beta} and
% repeat the adaptive Lasso procedure to refine the result.
%       beta_new_n contains the penalty term (beta_al/\hat{beta}) after
%       step 1.
%       beta2_al and beta2_new_n contain the results after step 2.
% Note that in the first step, a combination of gradient-based methods and 
% Newton method is used.
% by Kun Zhang 2016,2017 Carnegie Mellon University

var_noise_back = var_noise;

Trad1 = 0.2; Trad2 = 1 - Trad1;
[N,T] = size(x);
tol = 1E-2; % 1E-10; %%% temp: Aug 5, 17:27
beta_min = 1E-12;
beta_min2 = 1E-2;
sum_adjust_beta = [];
pl = [];

beta_hat = y*x' * inv(x*x');
if var_noise == 0
    var_noise = var(y - beta_hat*x);
end
x_new = diag(beta_hat) * x;
Error = 1;
beta_new_o = ones(N,1);
% beta_new_o = 1E-5 * ones(N,1);

% store for curve plotting
sum_adjust_beta = [sum_adjust_beta sum(abs(beta_new_o))];
pl = [pl (y-beta_new_o'*x_new)*(y-beta_new_o'*x_new)'/2/var_noise + lambda * sum(abs(beta_new_o))];

while Error > tol
    Sigma = diag(1./abs(beta_new_o));
    %     Sigma = diag(1./beta_new_o); % this is wrong?
    %     beta_new_n = inv(x_new*x_new' + var_noise*lambda * Sigma) * (x_new*y');
    % % with gradient trad-off!
    beta_new_n = inv(x_new*x_new' + var_noise*lambda * Sigma) * (x_new*y') * Trad1 + beta_new_o * Trad2 ;
    beta_new_n = sign(beta_new_n) .* max(abs(beta_new_n),beta_min);
    Error = norm(beta_new_n - beta_new_o);
    beta_new_o = beta_new_n;
    sum_adjust_beta = [sum_adjust_beta sum(abs(beta_new_n))];
    pl = [pl (y-beta_new_n'*x_new)*(y-beta_new_n'*x_new)'/2/var_noise + lambda * sum(abs(beta_new_n))];
end

Ind = find(abs(beta_new_n)>1E4*beta_min);
beta_new_n = beta_new_n .* (abs(beta_new_n)>1E4*beta_min);

beta_al = beta_new_n .* beta_hat';

% figure, plot(sum_adjust_beta, pl, 'r.-');

%% step 2
N2 = length(Ind);
x2 = x(Ind,:);
beta2_hat = y*x2' * inv(x2*x2');
if var_noise_back == 0
    var_noise = var(y - beta2_hat*x2);
end

x2_new = diag(beta2_hat) * x2;
beta2_new_o = ones(N2,1);
% beta2_new_o = 1E-5 * ones(N,1);
sum_adjust_beta2 = [];
pl2 = [];
sum_adjust_beta2 = [sum_adjust_beta2 sum(abs(beta2_new_o))];
pl2 = [pl2 (y-beta2_new_o'*x2_new)*(y-beta2_new_o'*x2_new)'/2/var_noise + lambda * sum(abs(beta2_new_o))];
Error = 1;
Iter = 1;
while Error > tol
    Sigma = diag(1./abs(beta2_new_o));
%     Sigma = diag(1./beta2_new_o); % this is wrong?
    if det(x2_new*x2_new' + var_noise*lambda * Sigma) < 0.01
        pause;
    end
    beta2_new_n = pdinv(x2_new*x2_new' + var_noise*lambda * Sigma) * (x2_new*y');
    beta2_new_n = sign(beta2_new_n) .* max(abs(beta2_new_n),beta_min);
    Error = norm(beta2_new_n - beta2_new_o);
    beta2_new_o = beta2_new_n;
    sum_adjust_beta2 = [sum_adjust_beta2 sum(abs(beta2_new_n))];
    pl2 = [pl2 (y-beta2_new_n'*x2_new)*(y-beta2_new_n'*x2_new)'/2/var_noise + lambda * sum(abs(beta2_new_n))];
    Iter = Iter + 1;
    if Iter > 100
        break;
    end
end
beta2_new_n = beta2_new_n .* (abs(beta2_new_n)>beta_min2);
beta2_al = zeros(N,1);

beta2_al(Ind) = beta2_new_n .* beta2_hat';


% hold on, plot(sum_adjust_beta2, pl2, 'k.-');