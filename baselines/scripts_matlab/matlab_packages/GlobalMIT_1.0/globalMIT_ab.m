%GlobalMIT: a toolbox for learning optimal dynamic Bayesian network structure with
%the Mutual Information Test (MIT) scoring metric
%(C) 2010-2011 Nguyen Xuan Vinh   
%Email: vinh.nguyen@monash.edu, vinh.nguyenx@gmail.com
%Reference: 
% [1] Vinh, N. X., Chetty, M., Coppel, R. and Wangikar, P. (2011). A polynomial time algorithm 
%     for learning globally optimal dynamic bayesian network.
%     2011-submitted for publication.
%Usage:
% [best_net,best_score]=globalMIT_ab(a,alpha)
% Input:
%       a,b: data, as preprocessed by multi_time_series_preprocessing.m
%       alpha: significance level for the mutual information test for
%       independance.
%       allowSelfLoop: allow self regulated link or not
% Output:
%       best_net: the globally optimal dynamic Bayesian network 
%                 best_net(i,j)=1 -> there exists an edge from node i->j
%       best_score: the optimal S_MIT score


function [best_net,best_score_ori]=globalMIT_ab(a,b,alpha,allowSelfLoop)

if nargin<3 alpha=0.999;end;
if nargin<4 allowSelfLoop=1;end;

best_net=[];best_score=inf;

[n dim]=size(a);
n_state=max(max(a));

chi=zeros(1,dim); %maximally dim parents
for i=1:dim
   chi(i)= chi2inv(alpha,n_state^(i-1)*(n_state-1)^2);
end

g_score=cumsum(chi);  %the precomputed g_MIT score
g_MIT=g_score;

%a need to be preprocessed, so that all variables take value in the range 1-n_state
n_state=max(max(a));

Ne=n;  %#effective samples, assuming single time-series data
best_score_arr=zeros(1,dim);  %best score for each individual node
HX_arr=zeros(1,dim);
best_net=zeros(dim,dim);

%main loop
for i=1:dim  %search the best parent set for each variable independantly
    
    HX=myEntropy(b(:,i),Ne,n_state);
    HX_arr(i)=HX;
    
    %investigate all set from 1->P*-1 elements
    Pstar=sum(g_score<2*Ne*HX)+1;
   
    fprintf('Processing node %d. P*= %d\n',i,Pstar);
    
    best_Pa=[];%empty set
    best_s_MIT=2*Ne*HX; %score of the empty network

    %co-ordinate coding: [Pai] -> p*-1 elements
    score=zeros(1,dim); %1-d array to store the scores of all parent combination
    new_score_arr=[];  
    
    for p=1:Pstar-1 %investigate all set from 1->P*-1 elements
        if g_MIT(p)>=best_s_MIT fprintf('Node %d, stopped at |Pa|=%d. ',i,p-1);break;end;     
        %otherwise, find the best net in this class of Pa with |Pa|=p
        %first: fill in the score array, only the score of the additional
        %variable
        if allowSelfLoop
            all_Pa = nchoosek([1: dim],p);   %warning: all parent set, this is only practical for small dim
        else
            all_Pa = nchoosek([1:i-1,i+1:dim],p);
        end

        %evaluate at level p and fill-in the score matrix
        [nn ll]=size(all_Pa);
        g_mit=g_score(p);
        new_score_arr=zeros(1,nchoosek(dim,p)); 
        
        for j=1:nn   %loop through all potential parent set
            Pa=all_Pa(j,:);
            nPa=length(Pa);
            if nPa==1  %only canculate the score for the 1st level
                CMI=conditional_MI_DBN_ab(a,b,Pa(1), i, 0,n_state);
                d_MIT=2*Ne*(HX-CMI);
                s_MIT=g_mit+d_MIT;
                if best_s_MIT>s_MIT  %found a better score
                     best_s_MIT=s_MIT;   
                     best_Pa=Pa;
                end
                position=Pa(1);  %position of this score in the score caching array
                score(position)=CMI;
            else  %get score from the cache,  and calculate score only for the last added variable
                position=findLexicalIndex(dim, Pa(1:end-1) ); %position from the cache
                score_i=score(position);
                
                %calculate the last score and store
                CMI=conditional_MI_DBN_ab(a,b,Pa(end), i, Pa(1:end-1),n_state);  %Attention: n-1 instead of n!
                score_i=score_i+CMI;
                
                position= findLexicalIndex(dim, Pa ); %position from the cache
                new_score_arr(position)=score_i; %store
                
                d_MIT=2*Ne*(HX-score_i);
                s_MIT=g_mit+d_MIT;
                if best_s_MIT>s_MIT  %found a better score
                     best_s_MIT=s_MIT;   
                     best_Pa=Pa;
                end
            end
        end %of loop through parent sets
        if p>1 score=new_score_arr;end; %update the MI cache
    end
    fprintf('Best score %f. \n',best_s_MIT);
    best_Pa
    
    best_score_arr(i)=best_s_MIT;
    if ~isempty(best_Pa) best_net(best_Pa,i)=ones(length(best_Pa),1);end;
    
end

best_score=sum(best_score_arr);
best_score_ori=2*Ne*sum(HX_arr)-best_score

end


function e=myEntropy(x,n,n_state)
H=zeros(1,n_state);
for j=1:n
   H(x(j))=H(x(j))+1; 
end
H=H/n;
H(find(H==0))=1;
e=-sum(H.*log(H));
end