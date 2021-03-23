%GlobalMIT: a toolbox for learning optimal dynamic Bayesian network structure with
%the Mutual Information Test (MIT) scoring metric
%(C) 2010-2011 Nguyen Xuan Vinh   
%Email: vinh.nguyen@monash.edu, vinh.nguyenx@gmail.com
%Reference: 
% [1] Vinh, N. X., Chetty, M., Coppel, R. and Wangikar, P. (2011). A polynomial time algorithm 
%     for learning globally optimal dynamic bayesian network.
%     2011-submitted for publication.
%Usage:  b=myIntervalDiscretize(a,n_inter)
%Discretize the data into n_inter states, using equal bin size
% Input:
%       a: raw data
% Output:
%       b: discretized data

function b=myIntervalDiscretize(a,n_inter)

[n dim]=size(a);
b=zeros(n,dim);
for i=1:dim
   b(:,i)=doDiscretize(a(:,i),n_inter);
end
b=b+1;
% ----------------------------------------
function y_discretized= doDiscretize(y,n_inter)
% ----------------------------------------
% discretize a vector
ymin= min(y);
ymax= max(y);
d=(-ymin+ymax)/n_inter;
vals=d*ones(1,n_inter);
vals=ymin+cumsum(vals);

vals(end)=inf;

y_discretized= y;
for i=1:length(y)
    y_discretized(i)=sum(vals<y(i));
%     if y_discretized(i)==n_inter
%         y_discretized(i)
%     end
end