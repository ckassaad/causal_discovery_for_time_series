function K = kernel(X,type,par)
% Input:    X   : matrix of the training examples;
%           type: type of kernel function to use: 'l' linear, 'p' polynomial, 'g' gaussian;
%           par : parameter of the kernel function;
% Output:   K: Kernel
l=size(X,2);
A= X' * X;
switch type
case {'linear','l'}
    K1 =A;
case {'polynomial','p'}
    K1 = (ones(l) + A).^par;
case {'homogenous polynomial','h'}
    K1 = A.^par;
 case {'gaussian','g'}
    B = ones(l, 1) * diag(A)';
    K1 = exp((-B'+2*A-B)/(2*par^2));
otherwise
    disp('Unknown kernel.')
    K1=0;
end
%simmetrizzaimo
Kx=ones(l,1)*mean(K1,1); 
K=K1-Kx-Kx'+mean(mean(K1));
