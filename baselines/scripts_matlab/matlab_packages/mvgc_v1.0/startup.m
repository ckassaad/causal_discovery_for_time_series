%% MVGC Toolbox "startup" script
%
% Initialise MVGC Toolbox. This file is run automatically if Matlab is started
% in the toolbox root (installation) directory.
%
% You may have to (or want to) customise this script for your computing
% environment.
%
% (C) Lionel Barnett and Anil K. Seth, 2012. See file license.txt in
% installation directory for licensing terms.
%
%% Set toolbox version

global mvgc_version;
mvgc_version.major = 1;
mvgc_version.minor = 0;

fprintf('[mvgc startup] Initialising MVGC toolbox version %d.%d\n', mvgc_version.major, mvgc_version.minor);

%% Set path

% Add mvgc root directory and appropriate subdirectories to path

global mvgc_root;
mvgc_root = fileparts(mfilename('fullpath')); % directory containing this file

% essentials
addpath(mvgc_root);
addpath(fullfile(mvgc_root,'core'));
addpath(fullfile(mvgc_root,'gc'));
addpath(fullfile(mvgc_root,'gc','GCCA_compat'));
addpath(fullfile(mvgc_root,'gc','subsample'));
addpath(fullfile(mvgc_root,'stats'));
addpath(fullfile(mvgc_root,'utils'));
if ~fexists(@rng) || ~fexists(@randi) % legacy hack
    addpath(fullfile(mvgc_root,'utils','legacy'));
    if ~fexists(@rng),   addpath(fullfile(mvgc_root,'utils','legacy','rng'));   end
    if ~fexists(@randi), addpath(fullfile(mvgc_root,'utils','legacy','randi')); end
end
addpath(fullfile(mvgc_root,'demo'));
addpath(fullfile(mvgc_root,'mex'));
addpath(fullfile(mvgc_root,'experimental'));
addpath(fullfile(mvgc_root,'docs')); % don't add the 'html' subdirectory

% comment out for release
% addpath(fullfile(mvgc_root,'testing'));
% addpath(fullfile(mvgc_root,'maintainer'));

fprintf('[mvgc startup] Added MVGC root directory %s and subdirectories to path\n',mvgc_root);

%% Check |mex| files

% Check for |mex| files and set flags appropriately

global have_genvar_mex;
have_genvar_mex = exist('genvar_mex','file') == 3;
if ~have_genvar_mex
    fprintf(2,'[mvgc startup] WARNING: no ''genvar'' mex file found; please run the ''mvgc_makemex'' script.\n');
    fprintf(2,'[mvgc startup]          Meanwhile, a slower scripted VAR simulation routine will be used.\n');
else
    fprintf('[mvgc startup] All MVGC ''mex'' files for your platform exist\n');
end

%% Check for dependencies on other Matlab(R) toolboxes

% Check if we have Statistics toolbox - see if ch2cdf is present

if fexists(@chi2cdf)
    fprintf('[mvgc startup] Statistics Toolbox(TM) seems to be present.\n');
else
	addpath(fullfile(mvgc_root,'utils','stats'));
    fprintf(2,'[mvgc startup] WARNING: Statistics Toolbox(TM) does not seem to be present.\n');
    fprintf(2,'[mvgc startup]          Will use slower scripted routines (see utils/stats directory).\n');
end

% Check if we have Signal Processing toolbox - see if pwelch is present

if fexists(@pwelch)
    fprintf('[mvgc startup] Signal Processing Toolbox(TM) seems to be present.\n');
else
    fprintf(2,'[mvgc startup] WARNING: Signal Processing Toolbox(TM) does not seem to be present.\n');
    fprintf(2,'[mvgc startup]          Some spectral estimation routines may not work.\n');
end

% Check if we have 'dlyap' from the Control System toolbox

if fexists(@dlyap)
    fprintf('[mvgc startup] Control System Toolbox(TM) seems to be present.\n');
else
	addpath(fullfile(mvgc_root,'utils','control'));
    fprintf(2,'[mvgc startup] WARNING: Control System Toolbox(TM) does not seem to be present.\n');
    fprintf(2,'[mvgc startup]          Will use slower scripted routines (see utils/control directory).\n');
end

%% Initialise default random number stream

% Have we got global rng control? Otherwise we're in annoying legacy territory

% Initialise rng to avoid predictability of sessions

rng_seed(-1); % seed from /dev/urandom (Unix/Mac) else from clock

fprintf('[mvgc startup] Random number generator initialised\n');

%% Enable all warnings

warning on all
fprintf('[mvgc startup] All warnings enabled\n');

% Done

fprintf('[mvgc startup] Initialisation complete (you may re-run ''startup'' at any time)\n');

fprintf('[mvgc startup] Type ''helpon'' to get started\n');

%%
% <startup.html back to top>
