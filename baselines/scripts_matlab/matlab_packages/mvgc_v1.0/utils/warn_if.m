%% warn_if
%
% Test a condition and warn if true
%
% <matlab:open('warn_if.m') code>
%
%% Syntax
%
%     wcond = warn_if(wcond,msgstr,warnfunc)
%
%% Arguments
%
% _input_
%
%     wcond      condition to test
%     msgstr     string specifying warning message
%     warnfunc   string specifying which function issued the warning
%
% _output_
%
%     wcond      condition to test
%
%% Description
%
% Tests the condition |wcond|; if |false| the warning message |msgstr| is
% printed. A string specifying which function issued the warning may be supplied
% in |warnfunc| - if not present it defaults to the enclosing function.
%
% This utility is included for consistency with <warn_test.html |warn_test|>.
%
%% Example
%
%     A = 10;
%     if warn_if(A > 9, 'A is > 9 !!');
%         return % can't carry on
%     end
%
%% See also
%
% <warn_test.html |warn_test|>
%
% (C) Lionel Barnett and Anil K. Seth, 2012. See file license.txt in
% installation directory for licensing terms.
%
%%

function wcond = warn_if(wcond,msgstr,wfunc)

assert(nargin > 1 && ischar(msgstr),'must specify a consition and a warning message string');
if nargin <  3, wfunc = []; end

if ~wcond, return; end

if isempty(wfunc) % no function specified
    ST = dbstack;
    if length(ST) > 1
        fprintf(2,'WARNING in ''%s'': %s\n',ST(2).name,msgstr);
    else
        fprintf(2,'WARNING in unknown function: %s\n',msgstr);
    end
else
    assert(ischar(wfunc),'bad warning function specification');
    fprintf(2,'WARNING in ''%s'': %s\n',wfunc,msgstr);
end
