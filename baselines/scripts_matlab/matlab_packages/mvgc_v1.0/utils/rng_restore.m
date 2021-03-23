%% rng_restore
%
% Restore the Matlab default (global) random number generator state
%
% <matlab:open('rng_restore.m') code>
%
%% Syntax
%
%     rng_restore(state)
%
%% Arguments
%
% _input_
%
%     state      previous rng state
%
%% Description
%
% Restore the Matlab default (global) <matlab:doc('rng') |random number
% generator|> state from |state|. If |state == 0| it is ignored. See
% <rng_seed.html |rng_seed|>, <rng_save.html |rng_save|>.
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
% <rng_seed.html |rng_seed|> |
% <rng_save.html |rng_save|>
%
% (C) Lionel Barnett and Anil K. Seth, 2012. See file license.txt in
% installation directory for licensing terms.
%
%%

function rng_restore(state)

if ~isstruct(state) && isscalar(state) && state == 0, return; end % do nothing

rng(state); % restore rng state
