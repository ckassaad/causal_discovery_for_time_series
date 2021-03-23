%% maxabs
%
% Calculate maximum absolute value of all entries in an array
%
% <matlab:open('maxabs.m') code>
%
%% Syntax
%
%     d = maxabs(X)
%
%% Arguments
%
% _input_
%
%     X          an array
%
% _output_
%
%     d          maximum absolute value of all entries in X
%
%% Description
%
% Returns the maximum absolute value |d| of all entries in |X|. |NaN| s are ignored.
%
% (C) Lionel Barnett and Anil K. Seth, 2012. See file license.txt in
% installation directory for licensing terms.
%
%%

function d = maxabs(X)

d = max(abs(X(:)));
