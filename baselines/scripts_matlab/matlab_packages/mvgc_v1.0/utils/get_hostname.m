%% get_hostname
%
% Get system host name
%
% <matlab:open('get_hostname.m') code>
%
%% Syntax
%
%     hostname = get_hostname
%
%% Arguments
%
% _output_
%
%     hostname   string; system short host name
%
%% Description
%
% Returns system short hostname for Unix, Mac and Windows machines. This is
% handy if you run Matlab on different machines and a program needs to know
% which machine it's running on.
%
% On Unix and Mac, it uses the system |hostname -s| call; if this fails, it
% tries the system |uname -n| call. On Windows machines it uses the call
% |getenv('COMPUTERNAME')|. If |get_hostname| fails it issues a warning and
% returns an empty string. If called without a return argument it prints out the
% name of the host name.
%
% (C) Lionel Barnett and Anil K. Seth, 2012. See file license.txt in
% installation directory for licensing terms.
%
%%

function hostname = get_hostname

hostname = '';

if isunix || ismac

    [status,hostname] = system('hostname -s');
    if status ~= 0 || isempty(hostname)
        [status,hostname] = system('uname -n');
        if status ~= 0
            hostname = '';
        end
    end
    if isempty(hostname) % give up
        fprintf(2,'WARNING(get_hostname): failed to get host name\n');
    else
        hostname = hostname(1:end-1); % system calls seem to append a newline - cull it
        if nargout == 0, fprintf('host name is ''%s''\n', hostname); end
    end

elseif ispc

    hostname = getenv('COMPUTERNAME');
    if isempty(hostname) % give up
        fprintf(2,'WARNING(get_hostname): failed to get host name\n');
    else
        if nargout == 0, fprintf('host name is ''%s''\n', hostname); end
    end

else
    fprintf(2,'WARNING(get_hostname): we only know about Unix, Mac or Windows\n');
end
