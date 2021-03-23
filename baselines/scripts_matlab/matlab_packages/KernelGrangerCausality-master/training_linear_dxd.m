function [w, H] = training_linear_dxd(X, y, d, l, lambda)
% [w, H] = training_linear_dxd(X, y, d, l, lambda)
%
% This function trains a linear machine by solving a linear system dxd large.
%
% Input: X: matrix dxl having the examples as its columns;
%           y: column vector having the output values for each input example;
%           d: number of component of each examples;
%           l: number of examples;
%           lambda: regularization parameter;
%
% Output: w: column vector of the linear predictor y = w'x.
%              H: matrix inv(X*X' + lambda*eye(d)) with size dxd.

H = pinv(X*X' + lambda*eye(d));
w = H*X*y;
