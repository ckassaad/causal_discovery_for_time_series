function beta=estim_beta_pham(x)
% estim_beta_pham(x)

[t1, t2] = size(x);
if t1 > t2
    error('error in eaastim_beta_pham(x): data must be organized in x in a row fashion')
end

beta = zeros(size(x));

t1 = scorecond(x')';
beta(1,:) = -t1(1,:);

t1 = scorecond(flipud(x)')';
beta(2,:) = -t1(1,:);