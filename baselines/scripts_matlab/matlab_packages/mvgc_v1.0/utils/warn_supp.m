%% warn_supp
%
% Suppress printing of Matlab warning message
%
% <matlab:open('warn_supp.m') code>
%
%% Syntax
%
%     oldstate = warn_supp(warnid)
%
%% Arguments
%
% _input_
%
%     warnid     Matlab warning identifier (default: 'all', for all warnings)
%
% _output_
%
%     oldstate   previous state of warning corresponding to warnid
%
%% Description
%
% Suppresses printing of Matlab warning message corresponding to warning ID
% |warnid|; see <matlab:doc('warning') |warning|>. The default is to
% suppress _all_ warnings. The routine <warn_test.html |warn_test|> can be
% used to test whether a condition triggering the given warning was
% actually issued.
%
% The <warn_supp.html |warn_supp|> ... <warn_test.html |warn_test|> routines
% exist for two reasons: (1) if routines issuing warnings are run in a loop
% (which they frequently are in the MVGC Toolbox), Matlab will by default spew
% confusing warnings, obscuring meaningful ouput; (2) we want to be able to
% issue meaningful context-dependent warning messages.
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
% <warn_test.html |warn_test|> |
% <startup.html |startup|>
%
% (C) Lionel Barnett and Anil K. Seth, 2012. See file license.txt in
% installation directory for licensing terms.
%
%%

function oldstate = warn_supp(warnid)

if nargin < 1 || isempty(warnid), warnid = 'all'; end % default: suppress all warnings

oldstate = warning('off', warnid); % suppress specified warning(s)
lastwarn('');                      % clear last warning
