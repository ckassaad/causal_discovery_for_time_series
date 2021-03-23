%% rng_save
%
% Save the Matlab default (global) random number generator state
%
% <matlab:open('rng_save.m') code>
%
%% Syntax
%
%     state = rng_save
%
%% Arguments
%
% _output_
%
%     state      current rng state
%
%% Description
%
% Save the Matlab default (global) <matlab:doc('rng') |random number generator|>
% state. It may be restored later using <rng_restore.html |rng_restore|>.
%
%% See also
%
% <matlab:doc('rng') |rng|> |
% <rng_seed.html |rng_seed|> |
% <rng_restore.html |rng_restore|>
%
% (C) Lionel Barnett and Anil K. Seth, 2012. See file license.txt in
% installation directory for licensing terms.
%
%%

function state = rng_save

state = rng; % save rng state
