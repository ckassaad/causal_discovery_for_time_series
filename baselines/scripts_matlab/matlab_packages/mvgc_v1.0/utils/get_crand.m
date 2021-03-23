%% get_crand
%
% Generate a high-entropy random number from clock time
%
% <matlab:open('get_crand.m') code>
%
%% Syntax
%
%     u = get_crand
%
%% Arguments
%
% _output_
%
%     u          a random number (double)
%
%% Description
%
% Generates a 32-bit integer from clock time via an
% <http://www.mathworks.com/matlabcentral/newsreader/view_thread/122364
% inscrutable munge> and converts it to double. (Nearly always) increases
% sequentially every tenth of a millisecond, and will not repeat for about five
% days. Handy for initial seeding of RNGs; see <rng_seed.html |rng_seed|>,
% <startup.html startup>. On Unix or Mac systems, <get_urand.html |get_urand|>
% is probably to be preferred (for Windows we might implement something based on
% <http://msdn.microsoft.com/en-us/library/sxtz2fa8%28v=vs.80%29.aspx |rand_s|>
% at some point).
%
%% See also
%
% <get_urand.html |get_urand|> |
% <rng_seed.html |rng_seed|> |
% <startup.html startup>
%
% (C) Lionel Barnett and Anil K. Seth, 2012. See file license.txt in
% installation directory for licensing terms.
%
%%

function u = get_crand

t = clock;
u = mod(round(cumprod([10000,60,60,24,31,12])*t(6:-1:1)'),2^32); % don't ask
