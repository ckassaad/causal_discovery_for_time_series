%% smvgc_to_mvgc
%
% Average (integrate) frequency-domain causality over specified frequency range
%
% <matlab:open('smvgc_to_mvgc.m') code>
%
%% Syntax
%
%     F = smvgc_to_mvgc(f,B)
%
%% Arguments
%
% See also <mvgchelp.html#4 Common variable names and data structures>.
%
% _input_
%
%     f          spectral (frequency-domain) Granger causality
%     B          frequency range; if unspecified (default) entire frequency
%                range is used
% _output_
%
%     F          Granger causality (time domain)
%
%% Description
%
% Calculates (conditional or unconditional) time-domain causality |F| from
% spectral causality |f| by integration (numerical quadrature - see <quads.html
% |quads|>) over the frequency range |B|. If a frequency band |B| is not
% supplied (default), the spectral causality is averaged from zero to the
% Nyqvist frequency. In that case the formula
%
% <<eq_smvgc_int.png>>
%
% (see [1]), where [[ii_nu.png]] is the Nyqvist frequency, should hold (at least
% approximately, numerically; see e.g. <mvgc_demo.html |mvgc_demo|>).
%
% A frequency band |B| is specified by a vector comprising pairs of points in
% ascending order in the range |[0,1]| - corresponding to zero up to the Nyqvist
% frequency. See <quadsr.html |quadsr|> for more details. In this case
% _band-limited_ time-domain causality [1,2] is calculated.
%
%% References
%
% [1] L. Barnett and A. K. Seth,
% <http://www.sciencedirect.com/science/article/pii/S0165027013003701 The MVGC
%     Multivariate Granger Causality Toolbox: A New Approach to Granger-causal
% Inference>, _J. Neurosci. Methods_ 223, 2014
% [ <matlab:open('mvgc_preprint.pdf') preprint> ].
%
% [2] L. Barnett and A. K. Seth, "Behaviour of Granger causality under
% filtering: Theoretical invariance and practical application", _J. Neurosci.
% Methods_ 201(2), 2011.
%
%% See also
%
% <quadsr.html |quadsr|> | <quadsr.html |quadsr|> | <sfreqs.html |sfreqs|> |
% <mvgc_demo.html |mvgc_demo|>
%
% (C) Lionel Barnett and Anil K. Seth, 2012. See file license.txt in
% installation directory for licensing terms.
%
%%

function F = smvgc_to_mvgc(f,B)

% FIXME - for dim = 1 !!!

sz = size(f);
nd = length(sz);
h  = sz(nd);
od = sz(1:nd-1);
vd = prod(od);

lam = linspace(0,1,h)';
f   = reshape(f,vd,h);
F   = nan(vd,1);
for i = 1:vd
    if nargin < 2 || isempty(B) % integrate over whole range
        F(i) = quads(lam,f(i,:)');
    else                        % integrate over sub-ranges in B
        [F(i),L] = quadsr(lam,f(i,:)',B);
        F(i) = F(i)/L;
    end
end
if ~isscalar(od)
    F = reshape(F,od);
end
