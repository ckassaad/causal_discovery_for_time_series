%GlobalMIT: a toolbox for learning optimal dynamic Bayesian network structure with
%the Mutual Information Test (MIT) scoring metric
%(C) 2010-2011 Nguyen Xuan Vinh   
%Email: vinh.nguyen@monash.edu, vinh.nguyenx@gmail.com
%Reference: 
% [1] Vinh, N. X., Chetty, M., Coppel, R. and Wangikar, P. (2011). A polynomial time algorithm 
%     for learning globally optimal dynamic bayesian network.
%     2011-submitted for publication.
%Usage: CMI=conditional_MI_DBN(a,Xi, Xj, Pa,n_state)
%Calculating the conditional mutual information between two nodes Xi->Xj
% Input:
%       a: a data matrix, rows are samples, columns are variables
%          values of a must be discrete, taking interger values >=1
%       Xi, Xj: two nodes, with an edge Xi->Xj
%       Pa: other parents, to condition on. Set Pa=0 for unconditioned MI
%       n_state: number of discrete state
% Output:
%       CMI: the conditional mutual information I(Xi;Xj|Pa)

function CMI=conditional_MI_DBN(a,Xi, Xj, Pa,n_state)

if nargin<5 n_state=max(max(a));end;
[n dim]=size(a);


scanned=zeros(1,n);  %mark the scanned row of a
prob_Pa=[];   %probability of each Pa configuration encountered
MI_arr=[];    %the corresponding MI value

Ne=n-1;  %number of effective observation
nPa=length(Pa);

if Pa==0
   Cont=Contingency(a(1:end-1,Xi),a(2:end,Xj),n_state);
   %calculate the unconditional MI
   CMI=Mutu_info(Cont,n_state);
   return;
end

for i=1:n-1  %scan all rows of a, find a new combination of Pa, and construct the corresponding contingency table
    if scanned(i)==0 %found a new combination of Pa
%         a(i,Pa)
        
        count=1;  %the number of times this combination appears in a
        scanned(i)=1;
        T=zeros(n_state,n_state);
        T(a(i,Xi),a(i+1,Xj))=T(a(i,Xi),a(i+1,Xj))+1;
        for j=i+1:n-1 
            if (scanned(j)==0 && sum(a(i,Pa)==a(j,Pa))==nPa)
                scanned(j)=1;
                T(a(j,Xi),a(j+1,Xj))=T(a(j,Xi),a(j+1,Xj))+1;
                count=count+1;
            end
        end
        MI_arr=[MI_arr Mutu_info(T,n_state)];
        prob_Pa=[prob_Pa count];
    end
end

prob_Pa=prob_Pa/sum(prob_Pa);  %normalize the probability 
CMI=sum(prob_Pa.*MI_arr);


function MI=Mutu_info(T,n_state)
%calculate the MI 
a=sum(T');
b=sum(T);
N=sum(a);

MI=0;
for i=1:n_state
    for j=1:n_state
        if T(i,j)>0 MI=MI+T(i,j)*log(T(i,j)*N/(a(i)*b(j)));end;
    end
end
MI=MI/N;

function Cont=Contingency(Mem1,Mem2,n_state)
Cont=zeros(n_state,n_state);
for i = 1:length(Mem1);
   Cont(Mem1(i),Mem2(i))=Cont(Mem1(i),Mem2(i))+1;
end