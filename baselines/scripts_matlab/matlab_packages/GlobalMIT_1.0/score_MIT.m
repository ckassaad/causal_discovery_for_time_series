%GlobalMIT: a toolbox for learning optimal dynamic Bayesian network structure with
%the Mutual Information Test (MIT) scoring metric
%(C) 2010-2011 Nguyen Xuan Vinh   
%Email: vinh.nguyen@monash.edu, vinh.nguyenx@gmail.com
%Reference: 
% [1] Vinh, N. X., Chetty, M., Coppel, R. and Wangikar, P. (2011). A polynomial time algorithm 
%     for learning globally optimal dynamic bayesian network.
%     2011-submitted for publication.
%Usage: score=score_MIT(a,net,alpha)
%Get the S_MIT score for a dynamic Bayesian network
% Input:
%       a: a data matrix, rows are samples, columns are variables
%          values of a must be discrete, taking interger values >=1
%       net: a dynamic Bayesian network, net(i,j)=1 -> there is an edge
%       from node i->j
%       alpha: significance level for the mutual information test for
%       independance.
% Output:
%       score: the  S_MIT score for the network

function score=score_MIT(a,net,alpha)

n_state=max(max(a));
if nargin<3 alpha=0.999;end;  %significance value

[n dim]=size(a);

chi=zeros(1,dim); %maximally n parents
for i=1:dim
   chi(i)= chi2inv(alpha,n_state^(i-1)*(n_state-1)^2);
end

score=0;
%score function is decomposable, get score for each variable
fprintf('S_MIT score: \n');
for i=1:dim
    %columns of net matrix represent parents of a  node
    nPa=sum(net(:,i));  %get the number of parents of xi
    Pa=1:dim;
    Pa=Pa(boolean(net(:,i))); %parent list
    score_i=0;
    
    if nPa>=1   %if this node has any parent at all the calculate the score
        score_i=2*(n-1)*conditional_MI_DBN(a,Pa(1), i, 0,n_state)-chi(1);
    end
    if nPa>=2
       for j=2:nPa
           score_i=score_i+2*(n-1)*conditional_MI_DBN(a,Pa(j), i, Pa(1:j-1),n_state)-chi(j);  %Attention: n-1 instead of n!
       end
    end
    fprintf('Node %d score %f\n',i,score_i);
score=score+score_i;  
end
fprintf('S_MIT score = %f\n',score);
