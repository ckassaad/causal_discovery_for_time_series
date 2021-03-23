function [y_predicted] = test_kernel_lxl(X, l, c, kernel_function_type, kernel_function_parameter, x)
% [y_predicted] = test_kernel_lxl(X, c, kernel_function_type, kernel_function_parameter, x)
%
% This function test a kernel machine on input patterns.
%
% Input: X: matrix of the training examples;
%           l: number of examples;
%           c: column vector of the coefficients of the kernel machine having l components;
%           kernel_function_type: type of kernel function to use: 1 linear, 2 polynomial, 3 gaussian;
%           kernel_function_parameter: parameter of the kernel function;
%           x: column vector of input pattern;
%
% Output: y_predicted: predicted values of the kernel machine.

switch kernel_function_type
    
case {'linear'}
    
    % disp('linear kernel')
    k = X' * x;
    y_predicted = c' * k;
    
case 'polynomial'
    
    % disp('polynomial kernel')
    degree = kernel_function_parameter;
    k = (ones(l, 1) + X' * x) .^ degree;
    y_predicted = c' * k;
    
case 'gaussian'
    
    % disp('gaussian kernel')
    sigma = kernel_function_parameter;
    for i=1:l
        k(i, 1) = exp(-(norm(X(:, i) - x) ^ 2) / (2 * sigma^2));
    end
    y_predicted = c' * k;
    
otherwise
    disp('Unknown kernel.')
    
end
