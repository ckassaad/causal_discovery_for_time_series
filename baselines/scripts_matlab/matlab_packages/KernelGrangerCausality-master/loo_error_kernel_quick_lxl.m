function [loo_error y_predicted_loo] = loo_error_kernel_quick_lxl(X, y, l, kernel_function_type, kernel_function_parameter, c, K, G, problem_type)
% [loo_error] = loo_error_kernel_quick_lxl(X, y, l, kernel_function_type, kernel_function_parameter, c, K, G, problem_type)
%
% This function computes the loo-error of a kernel machine training one machine only by using all the examples.
%
% Input: X: matrix dxl having the examples as its columns;
%           y: column vector having the output values for each input example;
%           l: number of examples;
%           kernel_function_type: type of kernel function to use: 1 linear, 2 polynomial, 3 gaussian;
%           kernel_function_parameter: parameter of the kernel function;
%           c: column vector of the coefficients of the kernel machine having l components.
%           K: kernel matrix with size lxl.
%           G: this is the matrix inv(K + lambda*eye(l)) with size lxl.
%           problem_type: define the type of problem related to your data: regression or classification.
%
% Output: loo_error: leave one out error.

% Compute the loo-error without training l kernel machines.

KG = K * G;

for j=1:l % scan all the examples;
    
    delta(j) = (y(j) - test_kernel_lxl(X, l, c, kernel_function_type, kernel_function_parameter, X(:, j))) / (1 - KG(j, j));
    
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
