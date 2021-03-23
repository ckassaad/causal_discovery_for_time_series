%% cpsd_xform
%
% Transform cross-power spectral density for reduced regression
%
% <matlab:open('cpsd_xform.m') code>
%
%% Syntax
%
%     SR = cpsd_xform(S,AR,SIGR)
%
%% Arguments
%
% See also <mvgchelp.html#4 Common variable names and data structures>.
%
% _input_
%
%     S          cross-power spectral density (cpsd) matrix
%     AR         VAR coefficients matrix for reduced regression
%     SIGR       residuals covariance matrix for reduced regression
%
% _output_
%
%     SR         transformed cpsd
%
%% Description
%
% Returns the cpsd |SR| for a new variable defined as the residuals of a reduced
% regression, for a VAR with cpsd |S|. |AR| and |SIGR| are the coefficients
% matrices and residuals covariance matrix respectively of the reduced
% regression, which is is assumed to correspond to the first |size(AR,1)|
% indices of |G|.
%
%% References
%
% [1] L. Barnett and A. K. Seth,
% <http://www.sciencedirect.com/science/article/pii/S0165027013003701 The MVGC
%     Multivariate Granger Causality Toolbox: A New Approach to Granger-causal
% Inference>, _J. Neurosci. Methods_ 223, 2014
% [ <matlab:open('mvgc_preprint.pdf') preprint> ].
%
%% See also
%
% <autocov_xform.html |autocov_xform|>
%
% (C) Lionel Barnett and Anil K. Seth, 2012. See file license.txt in
% installation directory for licensing terms.
%
%%

function S = cpsd_xform(S,AR,SIGR)

[n,~,h] = size(S);
fres = h-1;

[nx,nx1,~] = size(AR);
assert(nx1 == nx,'reduced VAR coefficients matrix has bad shape');
assert(nx <= n,'reduced VAR coefficients matrix appears to be for more variables than autocovariance sequence');

[n1,n2] = size(SIGR);
assert(n1 == n2,'reduced VAR residuals covariance matrix not square');
assert(n1 == nx,'reduced VAR residuals covariance matrix doesn''t match reduced VAR coefficients matrix');

x = 1:nx;
y = nx+1:n;

AF = bfft(cat(3,eye(nx),-AR),2*fres); % over [0,2*pi)

for k = 1:h
    S(x,x,k) = SIGR; % flat spectrum, since residuals serially uncorrelated
    S(x,y,k) = AF(:,:,k)*S(x,y,k);
    S(y,x,k) = S(y,x,k)*AF(:,:,k)';
end
