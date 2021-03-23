function [L ifail]=Leval(X,m,p,polycall);
global term
clear L
switch p
    case 1
        for k=1:m
            L(:,k)=X(k,:);
        end
    case 2
        k=0;
        for i=1:m
            k=k+1;
            L(:,k)=X(i,:).^2;
            k=k+1;
            L(:,k)=sqrt(2).*X(i,:);
            for j=i+1:m
                k=k+1;
                L(:,k)=sqrt(2).*X(i,:).*X(j,:);
            end
        end
    case 3
        k=0;
        for i=1:m
            k=k+1;
            L(:,k)=X(i,:).^3;
            k=k+1;
            L(:,k)=sqrt(3).*X(i,:).^2;
            k=k+1;
            L(:,k)=sqrt(3).*X(i,:);
            for j=i+1:m
                k=k+1;
                L(:,k)=sqrt(3).*X(i,:).*X(j,:).^2;
                k=k+1;
                L(:,k)=sqrt(3).*X(i,:).^2.*X(j,:);
                k=k+1;
                L(:,k)=sqrt(6).*X(i,:).*X(j,:);
                for l=j+1:m
                    k=k+1;
                    L(:,k)=sqrt(6).*X(i,:).*X(j,:).*X(l,:);
                end
            end
        end
    otherwise
        if polycall
            [term ifail]=polypower('X',m,p);
        end
        if ifail>0
            L=0;
            return
        end
        k=length(term);
        for j=1:k
            L(:,j)=eval(term{j});
        end
end
ifail=0;
for j=1:k
    L(:,j)=(L(:,j)-mean(L(:,j)));
end
