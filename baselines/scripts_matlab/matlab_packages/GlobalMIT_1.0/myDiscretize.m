%GlobalMIT: a toolbox for learning optimal dynamic Bayesian network structure with
%the Mutual Information Test (MIT) scoring metric
%(C) 2010-2011 Nguyen Xuan Vinh   
%Email: vinh.nguyen@monash.edu, vinh.nguyenx@gmail.com
%Reference: 
% [1] Vinh, N. X., Chetty, M., Coppel, R. and Wangikar, P. (2011). A polynomial time algorithm 
%     for learning globally optimal dynamic bayesian network.
%     2011-submitted for publication.
%Usage: b=myDiscretize(a)
%Discretize the data into 3 states, using equal bin size
% Input:
%       a: raw data
% Output:
%       b: discretized data
function b=myDiscretize(a)

[n dim]=size(a);
b=zeros(n,dim);
for i=1:dim
   b(:,i)=doDiscretize(a(:,i));
end
b=b+2;

% ----------------------------------------
function y_discretized= doDiscretize(y)
% ----------------------------------------
% discretize a vector
ymin= min(y);
ymax= max(y);
d= ymax-ymin;
thresh_low= d/3;
thresh_high= 2*d/3;
y = y-ymin;
y= (y>thresh_low) + (y>thresh_high);
y= y-1;
y_discretized= y;