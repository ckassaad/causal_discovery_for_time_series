%Radial basis function inner product
%Arthur Gretton

%Pattern input format : [pattern1 ; pattern2 ; ...]
%Output : Matrix of RBF values k(x1,x2)
%Deg is kernel size


function [H]=rbf_dot(patterns1,patterns2,deg);

%Note : patterns are transposed for compatibility with C code.

size1=size(patterns1);
size2=size(patterns2);


G = sum((patterns1.*patterns1),2);
H = sum((patterns2.*patterns2),2);

Q = repmat(G,1,size2(1));
R = repmat(H',size1(1),1);

H = Q + R - 2*patterns1*patterns2';


H=exp(-H/2/deg^2);

