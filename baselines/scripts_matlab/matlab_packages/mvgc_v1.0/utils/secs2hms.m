%% secs2hms
%
% Convert time in seconds to hours/minutes/seconds
%
% <matlab:open('secs2hms.m') code>
%
%% Syntax
%
%     [h,m,s] = secs2hms(t)
%
%% Arguments
%
% _input_
%
%     t          time in seconds
%
% _output_
%
%     h          hours
%     m          minutes
%     s          seconds
%
%% Description
%
% Return hours, minutes, seconds |h|, |m|, |s| of time |t| in seconds.
%
%% See also
%
% <timestr.html |timestr|>
%
%% Copyright notice
%
% [(C)] _Lionel Barnett and Anil K. Seth, 2012. See file
% <matlab:open('license.txt') license.txt> in root directory for licensing
% terms._
%
%%

function [h,m,s] = secs2hms(t)

s = t;
h = fix(s/3600);
s = s - 3600*h;
m = fix(s/60);
s = s - 60*m;

