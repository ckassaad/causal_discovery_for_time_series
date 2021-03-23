function [VN ifail]=vnorma(VT,V,VV)
f=1.e-8;
n=size(V,1)-1;
fast=size(VV,2)<n-2;
ifail=0;
if fast
    A=VV'*VT;
    B=(VV'*V)*(V'*VT);
    KKN=A'*A-B*A-A'*B'+B*B';
    [VVN D]=eig(KKN);
    d=chop(diag(D),1.e-8);
    iscomplex=sum(abs(imag(d)))>0;
    if iscomplex
        ifail=2;
    else
        ind=find(abs(d)>f*max(d) & abs(d) >f);
        VN=VV*VVN(:,ind);
        iscomplex=sum(sum(abs(imag(VN))))>0;
        if iscomplex
            ifail=3;
        end
    end
else
    K=VT*VT';
    P=V*V';
    KT=(K-P*K-K*P+P*K*P);
    [VVN D]=eig(KT);
    d=chop(diag(D),1.e-8);
    iscomplex=sum(abs(imag(d)))>0;
    if iscomplex
        ifail=2;
    else
        ind=find(d>f*max(d));
        VN=VVN(:,ind);
        iscomplex=sum(sum(abs(imag(VN(:,ind)))))>0;
        if iscomplex
            ifail=3;
        end
    end
end