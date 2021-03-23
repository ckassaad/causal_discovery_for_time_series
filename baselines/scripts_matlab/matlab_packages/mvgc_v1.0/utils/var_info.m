%% var_info
%
% Display VAR information as returned by |var_to_autocov|
%
% <matlab:open('var_info.m') code>
%
%% Syntax
%
%     acerr = var_info(info,abort_on_error)
%
%% Arguments
%
% See also <mvgchelp.html#4 Common variable names and data structures>.
%
% _input_
%
%     info            info structure returned by var_to_autocov
%     abort_on_error  if set to true, abort if var_to_autocov reports an error 
%
% _output_
%
%     acerr           flag set to true if var_to_autocov reports an error
%
%% Description
%
% Displays information (errors, warnings, diagnostics, _etc_.) reported by the
% routine <var_to_autocov.html |var_to_autocov|>; see that routine for details.
% This function should _always_ be called after a call to <var_to_autocov.html
% |var_to_autocov|>.
%
%% See also
%
% <var_to_autocov.html |var_to_autocov|>
%
% (C) Lionel Barnett and Anil K. Seth, 2012. See file license.txt in
% installation directory for licensing terms.
%
%%

function acerr = var_info(info,abort_on_error)

if nargin < 2 || isempty(abort_on_error), abort_on_error = true; end

acerr = info.error;

fprintf('\nVAR info:\n');
if info.error
    fprintf(2,'    ERROR: %s\n',info.errmsg);
else
    fprintf('    no errors\n');
end

if info.warnings > 0
    for n = 1:info.warnings
        fprintf(2,'    WARNING: %s\n',info.warnmsg{n});
    end
else
    fprintf('    no warnings\n');
end

if ~isnan(info.rho)
    fprintf('    spectral radius   : %f\n',info.rho);
end

if ~isnan(info.iters)
    fprintf('    iterations        : %d\n',info.iters);
end

if ~isnan(info.acrelerr)
    fprintf('    ac relative error : %g\n',info.acrelerr);
end

if ~isnan(info.acminlags)
    fprintf('    minimum ac lags   : %d\n',info.acminlags);
end

if ~isnan(info.aclags)
    fprintf('    actual  ac lags   : %d\n',info.aclags);
end

fprintf('\n');

assert(~(info.error && abort_on_error),'VAR info: aborting on error.');
