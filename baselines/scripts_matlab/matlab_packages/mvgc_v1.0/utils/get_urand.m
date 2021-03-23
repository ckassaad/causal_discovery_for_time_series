%% get_urand
%
% Read high-entropy random numbers from from |/dev/urandom| (Unix and Mac only)
%
% <matlab:open('get_urand.m') code>
%
%% Syntax
%
%     u = get_urand(n)
%
%% Arguments
%
% _input_
%
%     n          number of random numbers (default: 1)
%
% _output_
%
%     u          random numbers (double)
%
%% Description
%
% Reads |n| unsigned 32-bit integers from |/dev/urandom| and converts them to
% double. Handy for initial seeding of RNGs; see <rng_seed.html |rng_seed|>,
% <startup.html startup>. On systems without |/dev/urandom| an alternative is
% <get_crand.html |get_crand|>.
%
%% See also
%
% <get_crand.html |get_crand|> |
% <rng_seed.html |rng_seed|> |
% <startup.html startup>
%
% (C) Lionel Barnett and Anil K. Seth, 2012. See file license.txt in
% installation directory for licensing terms.
%
%%

function u = get_urand(n)

if nargin < 1, n = 1; end

fid = fopen('/dev/urandom');
assert(fid >= 0,'get_urand: failed to open /dev/urandom (does it exist on your system?)');
[u,count] = fread(fid,n,'uint32');
status = fclose(fid);
assert(status == 0,'get_urand: failed to close /dev/urandom');
assert(count == n,'get_urand: read %d of %d unsigned integers',count,n);
