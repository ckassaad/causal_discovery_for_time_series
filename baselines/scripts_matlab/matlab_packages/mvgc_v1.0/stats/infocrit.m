%% infocrit
%
% Calculate Akaike and Bayesian information criteria
%
% <matlab:open('infocrit.m') code>
%
%% Syntax
%
%     [aic,bic] = infocrit(L,k,m)
%
%% Arguments
%
% See also <mvgchelp.html#4 Common variable names and data structures>.
%
% _input_
%
%     L          maximum log-likelihood
%     k          number of free parameters
%     m          number of observations
%
% _output_
%
%     aic        AIC value
%     bic        BIC value
%
%% Description
%
% Calculates the Akaike and Bayesian information criteria from the
% maximised log-likelihood of a model [1]. For the Akaike information
% criterion, the finite sample bias-corrected form ("AICc") of Hurvich and
% Tsai [2] is calculated.
%
%% References
%
% [1] K. P. Burnham and D. R. Anderson, "Model Selection and Multimodel
% Inference: A Practical Information-Theoretic Approach", _2nd ed._,
% Springer-Verlag, 2002.
%
% [2] C. M. Hurvich and C.-L. Tsai, "Regression and time series model selection
% in small samples", _Biometrika_, 76, 1989.
%
%% See also
%
% <tsdata_to_infocrit.html |tsdata_to_infocrit|>
%
% (C) Lionel Barnett and Anil K. Seth, 2012. See file license.txt in
% installation directory for licensing terms.
%
%%

function [aic,bic] = infocrit(L,k,m)

if m-k-1 <= 0
    aic = NaN;
else
    aic = -2*L + 2*k*(m/(m-k-1)); % Note AIC without correction = -2*L + 2*k
end
bic = -2*L + k*log(m);
