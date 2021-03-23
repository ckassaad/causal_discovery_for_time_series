function [head,tail] = esitmateHeadAndTail(X,Y)
m=size(X,1);
acf = autocorr(X+Y,50);
smallestACF=find(acf<0.2);
head = smallestACF(1);
if (head > min(75,m) )
    warning('possibly long memory process, the output of test might be FALSE.')
end
head = min(head,50);
tail = m;
if (tail - head)< 100
    warning('using less than 100 points for a bootstrap approximation, stability of the test might be affected')
end
end