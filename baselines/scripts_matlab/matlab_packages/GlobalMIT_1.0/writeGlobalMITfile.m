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
%       a: a data matrix, rows are samples, columns are variables
%          values of a must be discrete, taking interger values >=1
%       alpha: significance level for the mutual information test for
%       independance.
% Output:
%       files necessary for GlobalMIT C++ version

function writeGlobalMITfile(a,alpha,net,filename)

[n dim]=size(a);

if nargin<2 alpha=0.999;end;
if nargin<3 net=round(rand(dim,dim));end;
if nargin<4 filename='./myDBNinput.txt';end;

fid=fopen(filename, 'w');
fprintf(fid,'%d %d\n',n,dim);
for i=1:n
    for j=1:dim
        fprintf(fid,'%d ',a(i,j));
    end
    fprintf(fid,'\n');
end
fclose(fid);

chi=zeros(1,dim); %maximally n parents
n_state=max(max(a));
for i=1:dim
   chi(i)= chi2inv(alpha,n_state^(i-1)*(n_state-1)^2);
end

filename='./myChiValue.txt';
fid=fopen(filename, 'w');
for i=1:dim
   fprintf(fid,'%f \n',chi(i)) ;
end
fclose(fid);

filename='./myChiValue.txt';
fid=fopen(filename, 'w');
for i=1:dim
   fprintf(fid,'%f \n',chi(i)) ;
end
fclose(fid);

filename='./myNet.txt';
fid=fopen(filename, 'w');
for i=1:dim
    for j=1:dim
        fprintf(fid,'%d ',net(i,j)) ;
    end
    fprintf(fid,'\n');
end
fclose(fid);