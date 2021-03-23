function [mse] = ms_error(y, y_predicted)
% [mse] = ms_error(y, y_predicted)
%
% This function computes the mean square error between the true y values and the predicted ones.
%
% Input: y: column vector having the true output values for each input example;
%           y_predicted: column vector having the predicted output values for each input example;
%
% Output: mse: mean square error.

z = y - y_predicted;

mse = (z' * z) / length(y);
