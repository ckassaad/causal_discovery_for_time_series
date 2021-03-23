%% rng_seed
%
% Seed the Matlab default (global) random number generator
%
% <matlab:open('rng_seed.m') code>
%
%% Syntax
%
%     state = rng_seed(seed)
%
%% Arguments
%
% _input_
%
%     seed       random seed: an integer in the range [0, 2^32 ? 1], or a negative number
%
% _output_
%
%     state      previous rng state
%
%% Description
%
% Seed the Matlab default (global) <matlab:doc('rng') |random number generator|>
% (which is used for calls such as <matlab:doc('rand') |rand|>,
% <matlab:doc('randn') |randn|>, etc.). If |seed| is positive, use that seed. If
% zero, do nothing. If negative, seed with something really unpredictable, using
% either <get_urand.html |get_urand|> (Unix, Mac) or <get_crand.html
% |get_crand|> (Windows); see e.g. the MVGC <startup.html |startup|> script.
% Optionally return the previous rng state in |state| (see <rng_restore.html
% |rng_restore|>).
%
%% Example
%
%     state = rng_seed(seed);
%     X = rand(10,3);
%     rng_restore(state);
%
%% See also
%
% <matlab:doc('rng') |rng|> |
% <rng_save.html |rng_save|> |
% <rng_restore.html |rng_restore|> |
% <get_urand.html |get_urand|> |
% <get_crand.html |get_crand|> |
% <startup.html |startup|>
%
% (C) Lionel Barnett and Anil K. Seth, 2012. See file license.txt in
% installation directory for licensing terms.
%
%%

function state = rng_seed(seed)

if nargout > 0
    if ~isstruct(seed) && isscalar(seed) && seed == 0
        state = 0;        % state = 0 means don't restore  (see rng_restore)
    else
        state = rng_save; % save state (see rng_restore)
    end
end

if seed == 0, return; end % do nothing

if seed < 0               % seed with unpredictable seed
    if isunix || ismac
        seed = get_urand; % from Unix /dev/urandom
    else
        seed = get_crand; % from clock time
    end
end

rng(seed);                % seed rng
