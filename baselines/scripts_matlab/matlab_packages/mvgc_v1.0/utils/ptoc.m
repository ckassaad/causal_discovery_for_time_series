%% ptoc
%
% Stop timer and print message
%
% <matlab:open('ptoc.m') code>
%
%% Syntax
%
%     ptoc(s1,s2,inhms)
%
%% Arguments
%
% _input_
%
%     s1         pre-time message string
%     21         post-time message string
%     inhms      use `timestr' time formatting utility to print time
%
%% Description
%
% Simple wrapper for Matlab <matlab:doc('toc') |toc|> timer function. Stop
% timer, print message string |s1|, followed by the elapsed time (see
% <ptic.html |ptic|>), followed by message string |s2|. If the |inhms| flag
% is set, the utility <timestr.html |timestr|> is used to format the
% elapsed time printout, else elapsed time is reported in seconds (see
% <matlab:doc('toc') |toc|>).
%
%% See also
%
% <ptic.html |ptic|> |
% <timestr.html |timestr|> |
% <matlab:doc('tic') |tic|> |
% <matlab:doc('toc') |toc|>
%
% (C) Lionel Barnett and Anil K. Seth, 2012. See file license.txt in
% installation directory for licensing terms.
%
%%

function ptoc(s1,s2,inhms)

t=toc;

if nargin < 1 || isempty(s1); s1 = '';   end
if nargin < 2 || isempty(s2); s2 = '\n'; end

if nargin > 2 && inhms
    fprintf([s1 '%s' s2],timestr(t));
else
    fprintf([s1 '%f secs' s2],t);
end
