function [G,Pinv,ifail] = cholesky(X,type,par,rmax,f)
% INPUT
% X             : matrix of the training examples;
% type          : type of kernel function to use: 'l' linear, 'p' polynomial, 'g' gaussian;
% par           : parameter of the kernel function;
% rmax            : maximal rank
% OUTPUT
% G             : Cholesky decomposition of kernel matrix K -> K(P,P)=G*G' (P is the permutation matrix)
% Pinv          : index of sorted pemutation matrix P
% ifail         : if true algorithm fails
% Copyright (c) Francis R. Bach, 2005.
% Modified by M.Pellicoro 2007

[Kmean Kmm diagK]=init_kernel(X,type,par);
[d n] = size(X);
G = zeros(n,rmax);    % Cholesky factor
P = 1:n;
kadv=0;               % current index of the look ahead steps
Dadv = diagK;
traceK=sum(diagK);
% performs rmax steps of Cholesky decomposition
ifail=true;
r=rmax;
for i=1:rmax
    kadv = kadv + 1;
    % select best index
    diagmax = Dadv(kadv);
    jast = 1;
    for j=1:n-kadv+1
        if (Dadv(j+kadv-1)>diagmax / .99)
            diagmax = Dadv(j+kadv-1);
            jast = j;
        end
    end
    if diagmax<0.1*f;
        % all pivots are too close to zero, stops 
        % this can only happen if the matrix has rank less than m
        r = kadv - 1;
        ifail=false;
        % reduce dimensions of decomposition
        G=G(:,1:r);
        break;
    else
        jast=jast+kadv-1;
        
        % permute indices 
        P( [kadv jast] )        = P( [jast kadv] );
        Dadv( [kadv jast] )     = Dadv( [jast kadv] );
        G([kadv jast],1:kadv-1) = G([ jast kadv],1:kadv-1);
        
        % compute new Cholesky column
        G(kadv,kadv)=Dadv(kadv);
        G(kadv,kadv)=sqrt(G(kadv,kadv));
        newKcol=kernel_column(X,P(kadv+1:n),P(kadv),type,par,Kmean,Kmm);
        G(kadv+1:n,kadv)=1/G(kadv,kadv)*( newKcol - G(kadv+1:n,1:kadv-1)*(G(kadv,1:kadv-1))');
        
        % update diagonal
        Dadv(kadv+1:n) =  Dadv(kadv+1:n) - G(kadv+1:n,kadv).^2;
        Dadv(kadv) = 0;
    end
end
[po Pinv]=sort(P);