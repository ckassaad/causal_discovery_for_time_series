%% isint
%
% Determine whether numerical values are integers
%
% <matlab:open('isint.m') code>
%
%% Syntax
%
%     I = isint(x)
%
%% Arguments
%
% _input_
%
%     x          a numeric array
%
% _output_
%
%     I          a logical array
%
%% Description
%
% Simple routine that returns a logical array with logical |true| wherever the
% corresponding entry in |x| is an integer.
%
% (C) Lionel Barnett and Anil K. Seth, 2012. See file license.txt in
% installation directory for licensing terms.
%
%%

function I = isint(x)

assert(isnumeric(x),'not a numeric array');
I = x == floor(x);
