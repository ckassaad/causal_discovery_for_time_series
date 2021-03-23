%% mvgc_makemex
%
% Build MVGC |mex| files
%
% <matlab:open('mvgc_makemex.m') code>
%
%% Syntax
%
%     mvgc_makemex(force_recompile,verbose)
%
%% Arguments
%
% _input_
%
%     force_recompile  forced recompilation flag (default: false)
%     verbose          verbosity flag (default: false)
%
%% Description
%
% Builds and then tests all MVGC |mex| files from |C| source in the |C|
% subdirectory. If a |mex| file for the current platform already exists in the
% |mex| subdirectory it is just tested, unless the |force_recompile| flag is
% set, in which case it is recompiled. A build is not considered successful
% unless it has also been successfully tested.
%
% Assumes a working Matlab-compatible |C| compiler for your platform.
% Sensible defaults are hard-coded, but you may have to (or want to) tweak this
% function for your platform/compiler.
%
% _*Note 1:*_ The toolbox is currently distributed with pre-built and tested
% |mex| files for 64-bit Unix (including Linux), Windows and Mac, as these were
% the only test platforms available to us. If Matlab crashes on you, there is
% a very good chance that a pre-built |mex| is to blame. In this case (assuming
% you have a Matlab-compatible C compiler available) you should try running
% <mvgc_makemex.html |mvgc_makemex|> with the |force_recompile| flag set.
%
% _*Note 2:*_ The pre-built Windows 64-bit |mex| files distributed with the
% toolbox were compiled with Microsoft(R) Visual Studio 2010. Apparently code
% compiled with this compiler requires the Microsoft(R) Visual Studio 2010
% runtime components. There is not much we can do about this; if you do not have
% Microsoft(R) Visual Studio 2010 installed on your 64-bit Windows system you
% can install the required components from
% <http://www.microsoft.com/en-us/download/details.aspx?id=14632 here>, or
% recompile the |mex| files using a different compiler, again by running
% <mvgc_makemex.html |mvgc_makemex|> with the |force_recompile| flag.
%
%% See also
%
% <startup.html |startup|> |
% <matlab:doc('mex') |mex|>
%
% (C) Lionel Barnett and Anil K. Seth, 2012. See file license.txt in
% installation directory for licensing terms.
%
%%

function mvgc_makemex(force_recompile,verbose)

% Build MVGC 'mex' files

if nargin < 1 || isempty(force_recompile), force_recompile = false; end;
if nargin < 2 || isempty(verbose),         verbose         = false; end;

% default mex command

MFLAGS = '-O -largeArrayDims';
if verbose, MFLAGS = ['-v ' MFLAGS]; end % spews lots of details
mexcmd = ['mex ' MFLAGS];

if isunix
    plat = 'Unix';
    CFLAGS = '-Wall -Werror -O3';  % gcc 4.4 ok with these
    mexcmd = [mexcmd ' CFLAGS="\$CFLAGS ' CFLAGS '"'];
elseif ispc
    plat = 'Windows';
    % If you want to override compiler flags, you're on your own...
elseif ismac
    plat = 'Mac';
    % If you want to override compiler flags, you're on your own...
else
    plat = 'Unknown';
    fprintf(2,'\nNOTE: At present ''mvgc_makemex'' has only been tested on Unix (gcc) and Windows (Visual Studio).\nIf you are on a different  platform/compiler and have any useful compilation tips, it would be\nvery helpfull if you could tell the MVGC maintainers about it.\n');
end

cc = mex.getCompilerConfigurations();

fprintf('\nYour platform   appears to be : %s (mex extension: %s)\n',plat,mexext);
fprintf('Your C compiler appears to be : %s\n',cc(1).Name);
fprintf('Will use the mex command      : %s <cfile> -outdir <mexdir>\n',mexcmd);

global mvgc_root;
cdir = fullfile(mvgc_root,'C');
mexdir = fullfile(mvgc_root,'mex');

% genvar

cfroot   = 'genvar';
mfroot   = [cfroot '_mex'];
cfile    = [mfroot '.c'];
mexfile  = [mfroot '.' mexext];
mexinvoc = [mexcmd ' ' cdir filesep cfile ' -outdir ' mexdir];

global have_genvar_mex;
if exist(mfroot,'file') == 3 && ~force_recompile
    fprintf('\nA mex file ''%s'' already exists for your platform. If you want\nto recompile it, then re-run this routine with the ''force_recompile'' flag set.\n\n',mexfile);
    have_genvar_mex = testmex(cfroot);
else
    fprintf('\nGoing to compile ''%s''\n',cfile);
    try
        eval(mexinvoc);
        fprintf('Looks like ''%s'' compiled ok. ',cfile);
        have_genvar_mex = testmex(cfroot);
    catch
        fprintf(2,'Hrmmph. ''%s'' failed to compile. Please tell the maintainer about this.\nMeanwhile don''t panic, an equivalent (but slower) scripted ''%s'' will be used.\n\n',cfile,cfroot);
    end
end

% More compilations and test functions go here

return

function success = genvar_mex_test

A = var9_test;
oldstate = rng_seed(67132);
E = randn(9,1000);
rng_restore(oldstate);
have_genvar_mex = true;
X1 = genvar(A,E);
have_genvar_mex = false;
X2 = genvar(A,E);
maxad = maxabs(X1-X2);
fprintf('\nMaximum absolute difference = %.g\n',maxad);
success = maxad < 1e-10;
return

% Test the mex file

function success = testmex(cfroot)

mfroot  = [cfroot '_mex'];
mexfile = [mfroot '.' mexext];
testfun = [cfroot '_mex_test'];

fprintf('Do you want to test ''%s'' (note: If Matlab crashes the test failed ;)\n',mexfile);
reply = input('Go for test? y/n [y]: ', 's');

success = false;
if isempty(reply) || strcmpi(reply,'y')
    if eval(testfun)
        fprintf('Congratulations, ''%s'' passed the test! It will now be used by default.\n\n',mexfile);
        success = true;
    else
        fprintf(2,'Drat. ''%s'' failed the test. Please tell the maintainer about this.\nMeanwhile don''t panic, an equivalent (but slower) scripted ''%s'' will be used.\n\n',mexfile,cfroot);
    end
else
    fprintf(2,'\nOk, you bottled out. An equivalent (but slower) scripted ''%s'' will be used.\n\n',cfroot);
end

return
