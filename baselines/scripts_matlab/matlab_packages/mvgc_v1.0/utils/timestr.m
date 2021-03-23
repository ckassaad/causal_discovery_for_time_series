%% timestr
%
% Format time string in hours/minutes/seconds
%
% <matlab:open('timestr.m') code>
%
%% Syntax
%
%     tstr = timestr(t)
%
%% Arguments
%
% _input_
%
%     t          time (seconds)
%
% _output_
%
%     tstr       formatted time string
%
%% Description
%
% Return formatted string |tstr| of time |t| (assumed in seconds).
%
%% See also
%
% <ptic.html |ptic|> |
% <ptoc.html |ptoc|> |
% <secs2hms.html |secs2hms|>
%
%% Copyright notice
%
% [(C)] _Lionel Barnett and Anil K. Seth, 2012. See file
% <matlab:open('license.txt') license.txt> in root directory for licensing
% terms._
%
%%

function tstr = timestr(t)

[h,m,s] = secs2hms(t);
if h > 0
    tstr = sprintf('%d:%d:%05.2f secs',h,m,s);
else
    if m > 0
        tstr = sprintf('%d:%05.2f secs',m,s);
    else
        tstr = sprintf('%.2f secs',s);
    end
end

