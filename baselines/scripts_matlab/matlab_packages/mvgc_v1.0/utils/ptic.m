%% ptic
%
% Print message and start timer
%
% <matlab:open('ptic.m') code>
%
%% Syntax
%
%     ptic(s)
%
%% Arguments
%
% _input_
%
%     s          message string
%
%% Description
%
% Simple wrapper for Matlab <matlab:doc('tic') |tic|> timer function. Print
% message string |s| and start a timer. Stop timer with <ptoc.html |ptoc|>.
%
%% See also
%
% <ptoc.html |ptoc|> |
% <matlab:doc('tic') |tic|> |
% <matlab:doc('toc') |toc|>
%
% (C) Lionel Barnett and Anil K. Seth, 2012. See file license.txt in
% installation directory for licensing terms.
%
%%

function ptic(s)

if nargin < 1 || isempty(s); s = 'TIMER: ';   end

fprintf(s);

tic
