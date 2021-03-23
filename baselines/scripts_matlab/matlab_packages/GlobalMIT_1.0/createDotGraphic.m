%GlobalMIT: a toolbox for learning optimal dynamic Bayesian network structure with
%the Mutual Information Test (MIT) scoring metric
%(C) 2010-2011 Nguyen Xuan Vinh   
%Email: vinh.nguyen@monash.edu, vinh.nguyenx@gmail.com
%Reference: 
% [1] Vinh, N. X., Chetty, M., Coppel, R. and Wangikar, P. (2011). A polynomial time algorithm 
%     for learning globally optimal dynamic bayesian network.
%     2011-submitted for publication.
%Usage: createDotGraphic(net,tit,filename)
%Create dot graphic
% Input:
%       net: a dynamic Bayesian network, net(i,j)=1 -> there is an edge
%       from node i->j
%       tit: title for the graph
%       nodeName: name of the node, optional, use empty variable []  for
%       default name (using number)
% Output:
%       a text file input for graphviz display
%       replace the actual path to Graphviz for display in Matlab

function createDotGraphic(net,nodeName,tit,filename)
if nargin<4 filename='.\default_dot_file.txt';end;
if nargin<3 tit='My graph';end;

fid=fopen(filename,'w');

fprintf(fid,'digraph abstract { \n');
fprintf(fid,['label = "' tit '";\n']); 
fprintf(fid,'labeljust="l"; \n');

[dim dim]=size(net);
for i=1:dim
    for j=1:dim
        if net(i,j)==1 
            if isempty(nodeName)
                fprintf(fid,'%d->%d;\n',i,j);
            else
                fprintf(fid,'%s->%s;\n',nodeName{i},nodeName{j});
            end;        
        end
    end
end
fprintf(fid,'}\n');
fclose(fid);

%%Replace with actual path to Graphviz
dos(['"C:\Program Files\Graphviz2.26.3\bin\dot" -Tjpg  ' filename ' -o .\myDefaultGraph.jpg']);
dos(['"C:\Program Files\Graphviz2.26.3\bin\dot" -Tpdf  ' filename ' -o .\myDefaultGraph.pdf']);

figure;
im=imread('.\myDefaultGraph.jpg');
image(im);
title('Dynamic Bayesian Network');