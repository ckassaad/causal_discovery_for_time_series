function [y_predicted] = test_linear_dxd(X, w)
% [y_predicted] = test_linear_dxd(X, w)
%
% This function test a linear machine on input patterns.
%
% Input: X: matrix of the input patterns;
%           w: column vector of the linear predictor y = w'x.
%
% Output: y_predicted: predicted values of the linear machine y = w'x.

y_predicted = X' * w;
