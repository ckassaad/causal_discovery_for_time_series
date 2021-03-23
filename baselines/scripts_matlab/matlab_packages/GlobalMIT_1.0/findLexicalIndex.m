%GlobalMIT: a toolbox for learning optimal dynamic Bayesian network structure with
%the Mutual Information Test (MIT) scoring metric
%(C) 2010-2011 Nguyen Xuan Vinh   
%Email: vinh.nguyen@monash.edu, vinh.nguyenx@gmail.com
%Reference: 
% [1] Vinh, N. X., Chetty, M., Coppel, R. and Wangikar, P. (2011). A polynomial time algorithm 
%     for learning globally optimal dynamic bayesian network.
%     2011-submitted for publication.
%Usage: [ pos ] = findLexicalIndex(n, Pa )
%Find the lexical order of a combination of parent
% Input:
%       n: total number of variable
%       Pa: a parent combination
% Output:
%       pos: the lexical order of Pa

function [ pos ] = findLexicalIndex(n, Pa )
%Find lexical index of a parent set, variable IDs are from 1->n
%variable in Pa must be in increasing order
p=length(Pa);
if p==1 pos=Pa;return;end;

pos=0;
for i=1:p-1
    if (i==1) last_pos=1; else last_pos=Pa(i-1)+1;end;
    
    for j=last_pos:Pa(i)-1
        pos=pos+nchoosek(n-j,p-i);
    end    
end
pos=pos+Pa(end)-Pa(end-1);

end

