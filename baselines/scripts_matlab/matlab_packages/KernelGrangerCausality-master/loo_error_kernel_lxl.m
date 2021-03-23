function [loo_error y_predicted_loo] = loo_error_kernel_lxl(X, y, l, lambda, kernel_function_type, kernel_function_parameter, problem_type)
% [loo_error] = loo_error_kernel_lxl(X, y, l, lambda, kernel_function_type, kernel_function_parameter, problem_type)
%
% This function computes the loo-error of a kernel machine training one machine only by using all the examples.
%
% Input: X: matrix dxl having the examples as its columns;
%           y: column vector having the output values for each input example;
%           l: number of examples;
%           kernel_function_type: type of kernel function to use: 1 linear, 2 polynomial, 3 gaussian;
%           kernel_function_parameter: parameter of the kernel function;
%           problem_type: define the type of problem related to your data: regression or classification.
%
% Output: loo_error: leave one out error.

% Compute the loo-error in the classic way training l kernel machines.

loo_error = 0;

for j=1:l % the column j is the example to leave out;
    
    j
    
    % leave one example out;
    [X_train, y_train, x_test, y_test] = loo(X, y, j);
    
    % train a kernel machine by solving a linear system (l-1)x(l-1) large.
    [c, K, G] = training_kernel_lxl(X_train, y_train, l-1, lambda, kernel_function_type, kernel_function_parameter);
    
    switch problem_type
        
    case 'regression'
        
        % disp('regression')
        
        % test a kernel machine on 1 input pattern with d components.
        y_predicted_loo = test_kernel_lxl(X_train, l-1, c, kernel_function_type, kernel_function_parameter, x_test);
        
        % compute the loo-error.
        loo_error = loo_error + ms_error(y_test, y_predicted_loo);
        
    case 'classification'
        
        % disp('--------->classification')
        
        % test a kernel machine on 1 input pattern with d components.
        y_predicted_loo = sign(test_kernel_lxl(X_train, l-1, c, kernel_function_type, kernel_function_parameter, x_test));
        
        % compute the loo-error.
        loo_error = loo_error + count_misclassified_patterns(y_test, y_predicted_loo);
        
    otherwise
        disp('Unknown problem type.')
        
    end
    
end

if problem_type == 'regression'
    
    loo_error = loo_error / l;
    
end
