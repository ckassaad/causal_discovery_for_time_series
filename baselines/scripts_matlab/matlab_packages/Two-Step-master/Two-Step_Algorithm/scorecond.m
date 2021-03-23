function [psi, entropy] = scorecond(data, q, bdwidth, cova)

% Syntaxe [psi, entropy] = scorecond(data, bdwidth, cova)
%
% Estimate the conditional score function defined as minus the
% gradient of the conditional density of of a random variable X_p
% given x_{p-1}, \dots, x_{q+1}. The values of these random variables are
% provided in the n x p array data.
%
% The estimator is based on the partial derivatives of the
% conditional entropy with respect to the data values, the entropy
% being estimated through the kernel density estimate _after a
% prewhitening operation_, with the kernel being the density of the
% sum of 3 independent uniform random variables in [.5,.5]. The kernel
% bandwidth is set to bdwidth*(standard deviation) and the density is
% evaluated at bandwidth apart. bdwidth defaults to
%   2*(11*sqrt(pi)/20)^((p-q)/(p-q+4))*(4/(3*n))^(1/(p-q+4)
% (n = sample size), which is optimal _for estimating a normal density_
%
% If cova (a p x p matrix) is present, it is used as the covariance
% matrix of the data and the mean is assume 0. This prevents
% recomputation of the mean and covariance if the data is centered
% and/or prewhitenned.
%
% The score function is computed at the data points and returned in
% psi.

if nargin < 1
  error('usage: psi = scorecond(data [, bdwidth, cova])');
end

[n, p] = size(data);

if nargin < 2; q = 0; end

if p < q+1, error('Sorry: not enough variables'); end

if nargin < 4
  tmp = sum(data)/n; % mean of data
  data = data - tmp(ones(n,1),:); % centered data
  cova = data'*data/n; % covariance matrix
end % otherwise centered, cova given

% prewhitening

T = chol(cova);
data = data/T;

if q > 0
  data(:,1:q) = []; % delete first q columns
  p = p-q;
end

if nargin < 3
  bdwidth = 2*(11*sqrt(pi)/20)^(p/(p+4))*(4/(3*n))^(1/(p+4));
end

% Grouping the data into cells, idx gives the index of the cell
% containing a datum, r gives its relative distance to the leftmost
% border of the cell

r = data/bdwidth;
idx = floor(r);
r = r - idx;
tmp = min(idx);
idx = idx - tmp(ones(n,1),:) + 1; % 0 <= idx-1

% Compute the probabilities at grid cells
% The kernel function is
%
%        1/2 + (1/2+u)(1/2-u) |u| <= 1/2
% k(u) = (3/2 - |u|)^2/2 1/2 <= |u| <= 3/2
%        0 otherwise
%
% The contribution to the probability at i-1, i, i+1 by a
% data point at distance r (0 <= r < 1) to i, are respectively:
% (1-r)^2/2, 1/2 + r*(1-r), r^2/2
% The derivative of the contribution to the probability at i-1, i, i+1
% by a data point at distance r to i are respectively: r-1, 1-2r, r

% The array ker contains the contributions to the probability of cells
% The array kerp contains the gradient of these contributions
% The array ix contains the indexes of the cells, arranged in
% _lexicographic order_

ker = [(1-r(:,1)).^2/2 .5+r(:,1).*(1-r(:,1)) r(:,1).^2/2];
ix = [idx(:,1) idx(:,1)+1 idx(:,1)+2];
kerp = [1-r(:,1) 2*r(:,1)-1 -r(:,1)];
mx = max(idx) + 2;
M = cumprod(mx);
nr = 1:n;

for i = 2:p
  ii = i(ones(1,3^(i-1))); %= i repeated 3^(i-1) times
  kerp = [[kerp.*(1-r(nr,ii)).^2/2, ...
           kerp.*(.5 + r(nr,ii).*(1-r(nr,ii))), ...
   kerp.*r(nr,ii).^2/2]
          [ker.*(1-r(:,ii)), ker.*(2*r(:,ii)-1), -ker.*r(:,ii)]];
  nr = [nr 1:n]; %= 1:n repeated i times
  ker = [ker.*(1-r(:,ii)).^2/2, ...
         ker.*(.5 + r(:,ii).*(1-r(:,ii))), ...
ker.*r(:,ii).^2/2];
  Mi = M(i-1);
  ix = [ix+Mi*(idx(:,ii)-1), ix+Mi*idx(:,ii), ix+Mi*(idx(:,ii)+1)];
end
  
pr = sparse(ix(:),1,ker(:),M(p),1)/n; % joint prob. of cells

logp = sparse(M(p),1); % to contain log(cond. prob.)
if p > 1
  pm = sum(reshape(pr,Mi,mx(p)), 2); % marginal prob. (Mi = M(p-1))
  pm = reshape(pm(:,ones(1,mx(p))),M(p),1);
  logp(find(pr)) = log(pr(find(pr))./pm(find(pr))); % avoid 0
else
  logp(find(pr)) = log(pr(find(pr)));
end

% compute the conditional entropy (if asked)

if nargout > 1 % compute the entropy
  entropy = log(bdwidth*T(end,end)) - pr'*logp;
end

% Compute the conditional score

psi = sum(logp(ix(nr,:)).*kerp,2); %nr = 1:n repeated p times
psi = reshape(psi, n, p)/bdwidth;

tmp = sum(psi)/n;
psi = psi - tmp(ones(n,1),:); % center psi

lam = psi'*data/n; % correction
lam = tril(lam) + tril(lam,-1)';
lam(p,p) = lam(p,p) - 1;

if q > 0
  psi = [zeros(n,q), psi - data*lam]/T';
else
  psi = (psi - data*lam)/T';
end
