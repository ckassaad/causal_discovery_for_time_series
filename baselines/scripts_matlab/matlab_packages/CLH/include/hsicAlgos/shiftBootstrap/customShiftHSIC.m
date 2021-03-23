function [test,nullApproximation] = customShiftHSIC(X,Y,alpha,head,tail,sigx,sigy,nullApproximation)

if (nargin < 7)
    error('customShiftHSIC:argumentCheck', 'Wrong number of input arguments')
end

m=size(X,1);

if m~=size(Y,1)
    error('customShiftHSIC:argumentCheck', 'X and Y must have the same lenght')
end

if ( head>tail || head<0 || tail >m ) 
    error('customShiftHSIC:argumentCheck', 'misspecified head and tail')
end
if ( alpha<0 || alpha >1)
   error('customShiftHSIC:argumentCheck', 'alpha must be such 0<alpha<1')
end
if ( sigx<0 || sigy <0)
   error('customShiftHSIC:argumentCheck', 'weights must be larger than zero')
end


K = rbf_dot(X,X,sigx);
L = rbf_dot(Y,Y,sigy);
H = eye(m)-1/m*ones(m,m);
Kc = H*K*H;
testStat = m*mean2(Kc'.*L);
if nargin==7
    nullApproximation = zeros(tail-head,1);
    for cut=head:tail
        indL = [cut:m 1:(cut-1)];
        nullApproximation(cut-head+1) = m*mean2(Kc'.*L(indL,indL));        
    end
end
test.pval  = nnz(nullApproximation > testStat)/(tail-head);
test.areDependent = test.pval < alpha;
end
