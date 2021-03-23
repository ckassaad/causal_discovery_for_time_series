%% var5_test
%
% Create VAR coefficients for 5-node test network
%
% <matlab:open('var5_test.m') code>
%
%% Syntax
%
%     A = var5_test
%
%% Arguments
%
% See also <mvgchelp.html#4 Common variable names and data structures>.
%
% _output_
%
%     A          VAR coefficients matrix for 5-node test network
%
%% Description
%
% Returns coefficients |A| for a VAR(3) based on a 5-node test network with causal
% architecture:
%
% <<var5test.png>>
%
% (C) Lionel Barnett and Anil K. Seth, 2012. See file license.txt in
% installation directory for licensing terms.
%
%%

function A = var5_test

n = 5;
p = 3;

r = sqrt(2);

A = zeros(n,n,p);

A(1,1,1) =  0.95*r;
A(1,1,2) = -0.9025;

A(2,1,2) =  0.5;
A(3,1,3) = -0.4;

A(4,1,2) = -0.5;
A(4,4,1) =  0.25*r;
A(4,5,1) =  0.25*r;

A(5,4,1) = -0.25*r;
A(5,5,1) =  0.25*r;
