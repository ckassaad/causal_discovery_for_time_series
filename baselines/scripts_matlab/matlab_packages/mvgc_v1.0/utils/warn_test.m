%% warn_test
%
% Test for (suppressed) Matlab warning message
%
% <matlab:open('warn_test.m') code>
%
%% Syntax
%
%     [waswarn,msgstr,warnid] = warn_test(oldstate, false)
%     [waswarn,msgstr,warnid] = warn_test(oldstate, warnmsg, warnfunc,rwarnid)
%
%% Arguments
%
% _input_
%
%     oldstate   old warning state to restore (required)
%     warnmsg    string specifying (user-supplied) warning message (or 'false' to suppress warning)
%     warnfunc   string specifying which function issued the warning
%     rwarnid    report Matlab warning identifier
%
% _output_
%
%     waswarn    true if a warning was issued
%     msgstr     warning message
%     warnid     warning Matlab warning identifier
%
%% Description
%
% Tests for Matlab warning suppressed by <warn_supp.html |warn_supp|>; see also
% <matlab:doc('warning') |warning|>. After testing, the previous warning state
% (as returned by <warn_supp.html |warn_supp|>) is restored.
%
% The <warn_supp.html |warn_supp|> ... <warn_test.html |warn_test|> routines
% exist for two reasons: (1) if routines issuing warnings are run in a loop
% (which they frequently are in the MVGC Toolbox), Matlab will by default spew
% confusing warnings, obscuring meaningful ouput; (2) we want to be able to
% issue meaningful context-dependent warning messages.
%
% If |warnmsg| is a logical flag set to |false|, reporting of warning messages
% is suppressed (presumably the user will handle warnings based on the return
% values). Otherwise, |warnmsg| is interpreted as a user-defined warning message
% to report; if empty (default), the Matlab warning message is reported. A
% string specifying which function issued the warning may be supplied in
% |warnfunc| - if not present it defaults to the enclosing function. If the flag
% |rwid| is set, the Matlab warining identifier is reported (default: false).
%
%% Example
%
%     A = magic(6);
%     oldstate = warn_supp('MATLAB:nearlySingularMatrix');
%     B = inv(A);
%     warn_test(oldstate,'badly-conditioned matrix?');
%
%% See also
%
% <matlab:doc('warning') |warning|> |
% <warn_supp.html |warn_supp|> |
% <warn_if.html |warn_if|> |
% <startup.html |startup|>
%
% (C) Lionel Barnett and Anil K. Seth, 2012. See file license.txt in
% installation directory for licensing terms.
%
%%

function [waswarn,msgstr,warnid] = warn_test(oldstate,umsgstr,wfunc,rwid)

assert(nargin > 0,'old warning state must be specified');

[msgstr,warnid] = lastwarn; % get last warning
waswarn = ~isempty(warnid); % was there a warning?
warning(oldstate);          % reset warning state

if ~waswarn, return; end

% report the warning
    
if nargin <  2, umsgstr = []; end
if nargin <  3 || isempty(wfunc), wfunc = true;  end
if nargin <  4 || isempty(rwid),  rwid  = false; end

if ~isempty(umsgstr) % override Matlab message
    assert(ischar(umsgstr),'bad reporting specification');
    msgstr = umsgstr;
end

if islogical(wfunc) % no function specified
    if wfunc % report warning
        ST = dbstack;
        if length(ST) > 1
            if rwid
                fprintf(2,'WARNING in ''%s'' (%s): %s\n',ST(2).name,warnid,msgstr);
            else
                fprintf(2,'WARNING in ''%s'': %s\n',ST(2).name,msgstr);
            end
        else
            if rwid
                fprintf(2,'WARNING in unknown function (%s): %s\n',warnid,msgstr);
            else
                fprintf(2,'WARNING in unknown function: %s\n',msgstr);
            end
        end
    %else    % don't report warning (presumably the user will handle it...)
    end
else
    assert(ischar(wfunc),'bad reporting specification');
    if rwid
        fprintf(2,'WARNING in ''%s'' (%s): %s\n',wfunc,warnid,msgstr);
    else
        fprintf(2,'WARNING in ''%s'': %s\n',wfunc,msgstr);
    end
end
