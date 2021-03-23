%GlobalMIT: a toolbox for learning optimal dynamic Bayesian network structure with
%the Mutual Information Test (MIT) scoring metric
%(C) 2010-2011 Nguyen Xuan Vinh   
%Email: vinh.nguyen@monash.edu, vinh.nguyenx@gmail.com
%Reference: 
% [1] Vinh, N. X., Chetty, M., Coppel, R. and Wangikar, P. (2011). A polynomial time algorithm 
%     for learning globally optimal dynamic bayesian network.
%     2011-submitted for publication.
%Usage:b=myDataMapping(a)
%Mapping the data to the correct range, omit any missing state
% Input:
%       a: a data matrix, rows are samples, columns are variables
%          values of a must be discrete, taking interger values >=1
% Output:
%       b: the processed data matrix, with values in each column taking all values in the
%          range 1 -> number of states
% Example: a=[1 3 5;3 4 2]   (state no 2 and 4 are missing)
%        =>b=[1 1 2;2 2 1]   (state takes continuous values)



function b=myDataMapping(a)

[n dim]=size(a);

for i=1:dim
  a(:,i)=map(a(:,i));
end
b=a;



function y=map(x)

b = sort(unique(x));  %set of unique value

for j=1:length(b)
    x(find(x==b(j)))=j; %do the mapping for each state    
end
y=x;