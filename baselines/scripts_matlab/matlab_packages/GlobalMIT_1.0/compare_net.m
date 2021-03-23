%GlobalMIT: a toolbox for learning optimal dynamic Bayesian network structure with
%the Mutual Information Test (MIT) scoring metric
%(C) 2010-2011 Nguyen Xuan Vinh   
%Email: vinh.nguyen@monash.edu, vinh.nguyenx@gmail.com
%Reference: 
% [1] Vinh, N. X., Chetty, M., Coppel, R. and Wangikar, P. (2011). A polynomial time algorithm 
%     for learning globally optimal dynamic bayesian network.
%     2011-submitted for publication.
%Usage: [best_net,score,time]=globalMIT_exe_ab(a,b,alpha)
%Compare two network and calculate the quality metrics
% Input:
%       net: the discovered network
%       true_net: the ground-truth network
% Output:
%       specificity, sensitivity,imprecision


function [specificity, sensitivity,imprecision]=compare_net(net,true_net,verbose)

if nargin<3 verbose=0;end;

[dim dim]=size(net);

FP=0;
FN=0;
TP=0;
TN=0;
for i=1:dim
   for j=1:dim
       if i==j continue;end; %don't care about the self-edge
       if (net(i,j)<true_net(i,j))
           FN=FN+1;
           if verbose fprintf('Missing edge: %d ---> %d\n',i,j);end;
       elseif (net(i,j)>true_net(i,j))
           FP=FP+1; 
           if verbose fprintf('False edge: %d ---> %d\n',i,j);end;
       elseif (net(i,j)==true_net(i,j) && net(i,j)==1)
           TP=TP+1;
       elseif (net(i,j)==true_net(i,j) && net(i,j)==0)
           TN=TN+1;
       end
           
   end
end

if FP==0 && FN==0 
    if verbose fprintf('Congratulations! You''ve done perfectly\n');end;
end;

specificity=TN/(TN+FP);
sensitivity=TP/(TP+FN);
recall=TP/(TP+FN);  %same as sensitivity
imprecision=FP/(FP+TP);

if verbose 
    fprintf('Specificity: %f ; 1-Spe= %f  Sensitivity: %f   \n',specificity,1-specificity,sensitivity);
    fprintf('Recall: %f ; Imprecision= %f  \n',recall, imprecision);
end;
