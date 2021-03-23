function mdir = make_legacy(tdir)

global mvgc_root;

if nargin < 1 || isempty(tdir), tdir = [tempdir 'mvgc_legacy']; end

fprintf('Populating target directory...');
syscmd = ['cp -r ' mvgc_root filesep '* ' tdir];
status = system(syscmd,'-echo');
if status == 0
    fprintf(' done\n');
else
   fprintf(2,' failed\n');
   return
end

mvgctree = genpath(mvgc_root);
while length(mvgctree) > 1
    [mdir, mvgctree] = strtok(mvgctree,pathsep);
    if isempty(strfind(mdir,'/docs')) && isempty(strfind(mdir,'/C')) && isempty(strfind(mdir,'/mex')) 
        make_leg(mdir,tdir)
    end
end

function make_leg(mdir,tdir)

global mvgc_root;

% Replace stuff in m files

tdir = strrep(mdir,mvgc_root,tdir);

fprintf('Target ''%s''\n',tdir);

mcells = struct2cell(dir([mdir filesep '*.m']));
n = size(mcells,2); % 1st two are ',' and '..'
for i = 1:n
    mfile = mcells{1,i};
    fprintf('\tProcessing m file %2d of %d = ''%s''...',i,n,mfile);

    mfilef = fullfile(mdir,mfile);
    fidr = fopen(mfilef, 'r');
    if fidr == -1
        fprintf(2,' failed to open input file ''%s''\n',mfilef);
        continue
    end    
    str = fread(fidr,inf,'*char')';    
    status = fclose(fidr);
    if status ~=  0
        fprintf(2,' failed to close input file ''%s''\n',mfilef);
        continue
    end
    
    % replace '~' in return values with 'dummy'

    expr   = '[~,';
    repstr = '[dummy,';
    str    = strrep(str,expr,repstr);

    expr   = ',~,';
    repstr = ',dummy,';
    str    = strrep(str,expr,repstr);

    expr   = ',~]';
    repstr = ',dummy]';
    str    = strrep(str,expr,repstr);

    expr   = '[~]';
    repstr = 'dummy';
    str    = strrep(str,expr,repstr);

    tfilef = fullfile(tdir,mfile);
    fidw = fopen(tfilef, 'w');
    if fidw == -1
        fprintf(2,' failed to open output file ''%s''\n',tfilef);
        continue
    end
    fwrite(fidw,str);
    status = fclose(fidw);
    if status ~=  0
        fprintf(2,' failed to close output file ''%s''\n',tfilef);
        continue
    end

    fprintf(' done\n');
end
