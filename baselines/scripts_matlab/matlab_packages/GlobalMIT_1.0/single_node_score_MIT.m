%GlobalMIT: a toolbox for learning optimal dynamic Bayesian network structure with
%the Mutual Information Test (MIT) scoring metric
%(C) 2010-2011 Nguyen Xuan Vinh   
%Email: vinh.nguyen@monash.edu, vinh.nguyenx@gmail.com
%Reference: 
% [1] Vinh, N. X., Chetty, M., Coppel, R. and Wangikar, P. (2011). A polynomial time algorithm 
%     for learning globally optimal dynamic bayesian network.
%     2011-submitted for publication.
%Usage: score_i=single_node_score_MIT(a,Pa,i,alpha)
%Get the S_MIT score for node i in the network
% Input:
%       a: a data matrix, rows are samples, columns are variables
%          values of a must be discrete, taking interger values >=1
%       Pa: set of parents for node i
%       alpha: significance level for the mutual information test for
%       independance.
% Output:
%       score_i: the  S_MIT score for node i


function score_i=single_node_score_MIT(a,Pa,i,alpha)

n_state=max(max(a));
if nargin<4 alpha=0.999;end;

chi=zeros(1,length(Pa));
for j=1:length(Pa)
   chi(j)= chi2inv(alpha,n_state^(j-1)*(n_state-1)^2);
end

[n dim]=size(a);

%columns of net matrix represent parents of a  node
nPa=length(Pa);
score_i=0;

if nPa>=1   %if this node has any parent at all the calculate the score
    %shift xi by 1 time step
    score_i=2*(n-1)*conditional_MI_DBN(a,Pa(1), i, 0,n_state)-chi(1);
end
if nPa>=2
    for j=2:nPa
        score_i=score_i+2*(n-1)*conditional_MI_DBN(a,Pa(j), i, Pa(1:j-1),n_state)-chi(j);  %Attention: n-1 instead of n!
    end
end
    