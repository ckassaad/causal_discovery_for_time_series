function [cb, cbiv]=new_partial_conditioning(x,kmax,type,par,m,ntrials,th)
%%% input
%%% x data number of time points x nvar
%%% kmax maximum number of variables to condition
%%% type par m  parameters of the model
%%% ntrials = number of trials
%%% th= theshold for statistical testing
%%% output
%%% cb causality partially conditioned
%%% cbiv bivariate granger
nvar=size(x,2);
cbiv=zeros(nvar,nvar);
for i=1:nvar
    for j=1:nvar
        if i ~= j
            cbiv(i,j)=causality_biv_trials(x(:,i),x(:,j),type,par,m,ntrials,th);
        end
    end
end
cb=zeros(nvar,nvar);
%fisso kmax (max numero di variabili su cui condizionare

for j=1:nvar
    [Y,I] = sort(cbiv(:,j),1,'descend');
    for i=1:nvar
        if i ~= j
            ind=setdiff(I(1:kmax),i);
            cb(i,j)=causality_trials(x(:,i),x(:,ind),x(:,j),type,par,m,ntrials,th);
        end
    end
end
