function [eta, z] = adaptive_size(grad_new, grad_old, eta_old, z_old)

alpha = 0; % 0.7
up = 1.05; %1.1 1.05
down = 0.5; % 0.4 0.5 0.34

z = grad_new + alpha * z_old;

etaup = (grad_new .* grad_old) >= 0;

eta = eta_old .* (up * etaup + down * (1 - etaup));
%%%%
eta = min(eta, 0.03); % 0.03