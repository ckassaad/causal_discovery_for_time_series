function [c, K, G] = training_kernel_lxl(X, y, l, lambda, kernel_function_type, kernel_function_parameter)
% [c, K, G] = training_kernel_lxl(X, y, l, lambda, kernel_function_type, kernel_function_parameter)
%
% This function trains a kernel machine by solving a linear system lxl large.
%
% Input: X: matrix dxl having the examples as its columns;
%           y: column vector having the output values for each input example;
%           l: number of examples;
%           lambda: regularization parameter;
%           kernel_function_type: type of kernel function to use: 1 linear, 2 polynomial, 3 gaussian;
%           kernel_function_parameter: parameter of the kernel function;
%
% Output: c: column vector of the coefficients of the kernel machine having l components.
%              K: kernel matrix with size lxl.
%              G: this is the matrix inv(K + lambda*eye(l)) with size lxl.

switch kernel_function_type
    
case {'linear'}
    
    % disp('linear kernel')
    K = X' * X;
    G = pinv(K + lambda*eye(l));
    c = G * y;
    
case 'polynomial'
    
    % disp('polynomial kernel')
    degree = kernel_function_parameter;
    K = (ones(l) + X' * X) .^ degree;
    G = pinv(K + lambda*eye(l));
    c = G * y;
    
case 'gaussian'
    
    % disp('gaussian kernel')
    
    sigma = kernel_function_parameter;
    A = X' * X;
    B = ones(l, 1) * diag(A)';
%     E = B' - 2*A + B;
    E = B';
    E = E - 2*A;
    E = E + B;
    K = exp(- E / (2 * sigma^2));
    G = pinv(K + lambda*eye(l));
    c = G * y;

%     for i=1:l
%         for j=1:l
%             z = X(:, i) - X(:, j);
%             E(i, j) = z' * z;
%         end
%     end
%     K = exp(- E / (2 * sigma^2));
%     G = inv(K + lambda*eye(l));
%     c = G * y;
    
otherwise
    disp('Unknown kernel.')
    
end
