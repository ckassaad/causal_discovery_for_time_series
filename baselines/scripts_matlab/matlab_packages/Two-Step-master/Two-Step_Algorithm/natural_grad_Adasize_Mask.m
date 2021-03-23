function W = natural_grad_Adasize_Mask(x,Mask)

[N,T] = size(x);
mu = 1E-3;
itmax = 6000;
Tol = 1E-4;
Num_edges = sum(sum(Mask));

% [icasig, AA, W] = fastica(x, 'approach', 'symm', 'g', 'tanh');

% initilization of W
WW = eye(N,N);
for i = 1:N
    Ind_i = find(Mask(i,:)~=0);
    WW(i,Ind_i) = -.5*(x(i,:) * x(Ind_i,:)') * pdinv(x(Ind_i,:) * x(Ind_i,:)');
end
W = .5*(WW+WW');

z = zeros(N,N);
eta = mu * ones(size(W));
W_old = W; 
y_psi = [];


for iter = 1:itmax
  %  fprintf('.');
    if ~mod(iter,100)
  %      fprintf('\n');
    end
    y = W * x;
    
    %     % update W: linear ICA with marginal score function estimated from data...
    %     if mod(iter,10) ==1
    %         for i = 1:N
    %             tem = estim_beta_pham(y(i,:));
    %             y_psi(i,:) = tem(1,:);
    %         end
    %     end
    
    % update W: linear ICA with marginal score function estimated from data...
    if mod(iter,12) == 1
        for i = 1:N
            tem = estim_beta_pham(y(i,:));
            y_psi(i,:) = tem(1,:);
            [y0{i} II]  = sort(y(i,:), 'ascend'); % [y0{i} II JJ]  = unique(y(i,:));
            y_psi0{i} = y_psi(i,II);
        end
    else
        for i = 1:N
            %             if length(y0{i})~=length(y(i,:))
            %                 y0{i} = [y0{i} y0{end}+ 1E-3*std(y0{i}) * [1:length(y(i,:))-length(y0{i})] ];
            %                 y_psi0{i} = [y_psi0{i} y_psi0{i}(end) * ones(1, length(y(i,:))-length(y0{i}))];
            %             end
            [tmp, II2] = sort(y(i,:), 'ascend');
            y_psi(i,II2) = y_psi0{i};
            %             y_psi(i,:) = lininterp1(y0{i}', y_psi0{i}', y(i,:)')';
        end
    end
    
    G = y_psi * y'/T;
    yy = y*y'/T;
    I_N = eye(N);
    %      H = G - diag(diag(G)) + I_N - diag(diag(yy));
    %     H = G - diag(diag(G));
    %     W = (I_N + mu*H) * W;
    Grad_W_n = y_psi*x'/T + inv(W') ;
    if iter ==1
        Grad_W_o = Grad_W_n;
    end
    [eta, z] = adaptive_size(Grad_W_n, Grad_W_o, eta, z);
    delta_W = eta.*z;
    W = W + delta_W .* Mask;
    %     H = Grad_W .*Mask;
    %     W = W + mu * H;
    %     delta_H = sum(abs(H)),
    %     if delta_H < Tol
    % delta_H/N,
    %         return;
    %     end
    % sum(sum(Mask)) * sum(sum(abs(W - W_old))),
    % sum(sum(abs(Grad_W_n)))/N^2,
    %     if sum(sum(Mask)) * sum(sum(abs(W - W_old)))<Tol
    %         break;
    %     end
   % sum(sum(abs(Grad_W_n .* Mask)))/Num_edges;
    if sum(sum(abs(Grad_W_n .* Mask)))/Num_edges<Tol
        break;
    end
    Grad_W_o = Grad_W_n;
    W_old = W;
end
