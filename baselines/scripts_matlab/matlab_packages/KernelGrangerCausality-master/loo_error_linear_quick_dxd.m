function [loo_error y_predicted_loo] = loo_error_linear_quick_dxd(X, y, l, w, H, problem_type)
% [loo_error] = loo_error_linear_quick_dxd(X, y, l, w, H, problem_type)
%
% This function computes the loo-error of a linear machine training one machine only by using all the examples.
%
% Input: X: matrix dxl having the examples as its columns;
%           y: column vector having the output values for each input example;
%           l: number of examples;
%           w: column vector of the linear predictor y = w'x.
%           H: matrix inv(X*X' + lambda*eye(d)) with size dxd.
%           problem_type: define the type of problem related to your data: regression or classification.
%
% Output: loo_error: leave one out error.

% Compute the loo-error without training l linear machines.

KG = X' * H * X;

for j=1:l % scan all the examples;
    delta(j) = (y(j) - test_linear_dxd(X(:, j), w)) / (1 - KG(j, j));
end

switch problem_type
    
case 'regression'
    
    % disp('regression')
    
    % compute the predicted values of the machines trained leaving one example out.
    y_predicted_loo = y - delta';
    
    % compute the loo-error.
    loo_error = ms_error(y, y_predicted_loo);
    
case 'classification'
    
    % disp('classification')
    
    % compute the predicted values of the machines trained leaving one example out.
    y_predicted_loo = sign(y - delta');
    
    % compute the loo-error.
    loo_error = count_misclassified_patterns(y, y_predicted_loo);
    
otherwise
    disp('Unknown problem type.')
    
end
