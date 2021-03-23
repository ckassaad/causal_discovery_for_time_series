function [term ifail r]=polypower(varname,m,p)
ifail=0;
sp=num2str(p);
NN='[n1';
str=sp;
command=['k=0;for n1=0:' str ';'];
myend='end;';
for i=2:m
    NN=strcat(NN,sprintf(' n%d',i));
    str=strcat(str,sprintf('-n%d',i-1));
    command=strcat(command,sprintf('for n%d=0:',i),str,';');
    myend=strcat(myend,'end;');
end
NN=strcat(NN,sprintf(' n%d]',m+1));
str=strcat(str,sprintf('-n%d',m));
command=strcat(command,sprintf('n%d=',m+1),str,';','k=k+1;N(k,:)=',NN,';nfc(k)=multinom(',sp,',',NN,');',myend);
disp(length(command));
eval(command,'ifail=4;');
if ifail>0
    term=0;
    r=0;
    return
end
nterm=k-1;
for i=1:nterm
    if nfc(i)>1
        tt=sprintf('sqrt(%d).*',nfc(i));
    else
        tt='';
    end
    for j=2:m
        switch N(i,j)
            case 0
            case 1
                tt=strcat(tt,sprintf('%s(%d,:).*',varname,j-1));
            otherwise
                tt=strcat(tt,sprintf('%s(%d,:).^%d.*',varname,j-1,N(i,j)));
        end
    end
    switch N(i,m+1)
        case 0
            tt=tt(1:length(tt)-2);
        case 1
            tt=strcat(tt,sprintf('%s(%d,:)',varname,m));
        otherwise
            tt=strcat(tt,sprintf('%s(%d,:).^%d',varname,m,N(i,m+1)));
    end
    term{i,1}=tt;
end
r=[nfc' N];

