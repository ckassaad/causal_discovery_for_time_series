%% helpon
%
% Display the <mvgcfuncref.html Function Reference> page for an MVGC function or
% script in the Matlab Help Browser.
%
% <matlab:open('get_hostname.m') code>
%
%% Syntax
%
%     helpon mname
%     helpon(mname)
%
%% Arguments
%
% _input_
%
%     mname      string; the name of an MVGC function or script
%
%% Description
%
% Displays help in the Matlab Help Browser on |mname|, provided an html help
% file called |mname.html| exists in the docs/html subdirectory of the MVGC root
% directory. If |mname.html| is not found, a "not found" page is displayed.
%
%% See also
%
% <mvgchelp.html#2 Help and documentation> |
% <mvgcfuncref.html Function Reference>
%
% (C) Lionel Barnett and Anil K. Seth, 2012. See file license.txt in
% installation directory for licensing terms.
%
%%

function helpon(mname)

global mvgc_root;

if nargin < 1 || isempty(mname), mname = 'mvgchelp'; end

htmlfile = fullfile(mvgc_root,'docs','html',[mname '.html']);
if exist(htmlfile,'file') == 2
    web(htmlfile,'-helpbrowser');
else
    notfound = fullfile(mvgc_root,'docs','html','notfound.html');
    web(notfound,'-helpbrowser');
    fprintf(2,'\nSorry, no help on ''%s''\n\n',mname);
end
