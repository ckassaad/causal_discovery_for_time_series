function [Kmean Kmm diagK]=init_kernel(X,type,par)
n=size(X,2);
switch type
    case {'linear','l'}
        for i=1:n
            v=X'*X(:,i);
            Kmean(i)=mean(v);
            diagK(i)=v(i);
        end
    case {'polynomial','p'}
        for i=1:n
            v=(1+X'*X(:,i)).^par;
            Kmean(i)=mean(v);
            diagK(i)=v(i);
        end
    case {'gaussian','g'}
        for i=1:n
            v=exp((-X(:,i)'*X(:,i)-sum(X.^2,1)'+2*X'*X(:,i))/(2*par^2));
            Kmean(i)=mean(v);
            diagK(i)=v(i);
        end
    otherwise
        disp('Unknown kernel.')
        Kmean=0;
        ki=0;
end
Kmm=mean(Kmean);
diagK=(diagK-2*Kmean+Kmm)';