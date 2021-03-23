function [empirical_risk_kernel,q_loo_error_kernel,y_predicted_kernel]=compute_errors(lambda,kernel_function_type,kernel_function_parameter,m,L,X,y)


stdout = 1;
stderr = 2;


problem_type = 'regression'; % define the type of problem related to your data.

% -----------------> Predict s using knowledge about x.



%[X,Y,xt,yt,L]=prepro(m);

number_of_example = L;
number_of_component = m;

% Compare the empirical risk with the loo-error.



  

    % train a kernel machine by solving a linear system number_of_example^2 large.
    [c, K, G] = training_kernel_lxl(X, y, number_of_example, lambda, kernel_function_type, kernel_function_parameter);
    
    % test a kernel machine on l input patterns with d components.
    for i=1:number_of_example
        y_predicted_kernel(i, 1) = test_kernel_lxl(X, number_of_example, c, kernel_function_type, kernel_function_parameter, X(:, i));
    end
    
    % compute the empirical risk.
    empirical_risk_kernel = ms_error(y, y_predicted_kernel);
    
    % Compute the loo-error without training number_of_example kernel machines.
    q_loo_error_kernel = loo_error_kernel_quick_lxl(X, y, number_of_example, kernel_function_type, kernel_function_parameter, c, K, G, problem_type);
    
    % classic leave one out procedure for estimating the best lambda.
%    loo_error_kernel(k) = loo_error_kernel_lxl(X, y, number_of_example, lambda, kernel_function_type, kernel_function_parameter, problem_type);
   
