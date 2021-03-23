function Kcol = kernel_column(X,ind,ip,type,par,Kmean,Kmm)
% Input:    X: matrix of the training examples;
%           kernel_function_type: type of kernel function to use: 'l' linear, 'p' polynomial, 'g' gaussian;
%           kernel_function_parameter: parameter of the kernel function;
% Output:   K: Kernel
n=size(X,2);
switch type
    case {'linear','l'}
            Kcol=X(:,ind)'*X(:,ip);
    case {'polynomial','p'}
            Kcol=(1+X(:,ind)'*X(:,ip)).^par;
    case {'gaussian','g'}
            Kcol=exp((-X(:,ip)'*X(:,ip)-sum(X(:,ind).^2,1)'+2*X(:,ind)'*X(:,ip))/(2*par^2));
    otherwise
        disp('Unknown kernel.')
        Kmean=0;
end
%simmetrizziamo
Kcol=Kcol-Kmean(ip)-Kmean(ind)'+Kmm;