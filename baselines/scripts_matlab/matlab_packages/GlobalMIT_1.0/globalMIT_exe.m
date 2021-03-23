%GlobalMIT: a toolbox for learning optimal dynamic Bayesian network structure with
%the Mutual Information Test (MIT) scoring metric
%(C) 2010-2011 Nguyen Xuan Vinh   
%Email: vinh.nguyen@monash.edu, vinh.nguyenx@gmail.com
%Reference: 
% [1] Vinh, N. X., Chetty, M., Coppel, R. and Wangikar, P. (2011). A polynomial time algorithm 
%     for learning globally optimal dynamic bayesian network.
%     2011-submitted for publication.
%Usage: [best_net,score,time]=globalMIT_exe(a,alpha)
%Matlab interface for the GlobalMIT C++ 
% Input:
%       a: a data matrix, rows are samples, columns are variables
%          values of a must be discrete, taking interger values >=1
%       alpha: significance level for the mutual information test for
% Output:
%       best_net: the globally optimal dynamic Bayesian network 
%                 best_net(i,j)=1 -> there exists an edge from node i->j
%       best_score: the optimal S_MIT score
%       time: execution time


function [best_net,score,time]=globalMIT_exe(a,alpha,allowSelfLoop)
[n dim]=size(a);

if nargin<2 alpha=0.999;end;
if nargin<3 allowSelfLoop=1;end;
 
a=myDataMapping(a);
n_state=max(max(a));

writeGlobalMITfile(a,alpha);

tic;dos(['globalMIT.exe ' num2str(alpha) ' ' num2str(allowSelfLoop) ' ./myDBNinput.txt ./myDBNoutput.txt']);time=toc;

fid=fopen('./myDBNoutput.txt','r');
s=fgetl(fid);
score=str2num(s);

best_net=zeros(dim,dim);

for i=1:dim
    s=fgetl(fid);
    best_net(i,:)=str2num(s);
end
fclose(fid);