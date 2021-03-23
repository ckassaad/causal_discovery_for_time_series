%% var_decay
%
% Decay VAR coefficients
%
% <matlab:open('var_decay.m') code>
%
%% Syntax
%
%     A = var_decay(A,dfac)
%
%% Arguments
%
% See also <mvgchelp.html#4 Common variable names and data structures>.
%
% _input_
%
%     A          VAR coefficients matrix
%     dfac       decay factor
%
% _output_
%
%     A          decayed VAR coefficients matrix
%
%% Description
%
% Exponentially decay VAR coefficients |A| by factor |dfac|, which controls
% how fast coefficients decay with lag. This will, in particular,  affect
% the _spectral radius_ (see <var_specrad.html |var_specrad|>) of the
% coefficients.
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
% <var_specrad.html |var_specrad|>
%
% (C) Lionel Barnett and Anil K. Seth, 2012. See file license.txt in
% installation directory for licensing terms.
%
%%

function A = var_decay(A,dfac)

p = size(A,3);
f = dfac;
for k=1:p
    A(:,:,k) = f*A(:,:,k);
    f = dfac*f;
end
