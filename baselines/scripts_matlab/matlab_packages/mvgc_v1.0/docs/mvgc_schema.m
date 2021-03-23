%% Overview: MVGC computational pathways
%
% The diagram below, adapted from the MVGC Toolbox reference document> [1]
% (Section 3) illustrates key computational pathways implemented in the toolbox.
%
% <<mvgcschema.png>>
%
% A great deal of effort has gone into identifying useful, accurate and
% numerically efficient computational pathways; these are represented as bold
% arrows in the diagram. Blue arrows represent actual Granger causality
% computation, while the dashed arrow represents simulation for testing
% purposes. The pink shaded circles represent the equivalent VAR representations
% [1], while the |Ann| label computational algorithms implemented in the toolbox
% (see table below).
%
% See the <mvgc_demo.html mvgc_demo> script for a demonstration of how this works out in
% practice.
%
%% Data structures
%
% See also <mvgchelp.html#4 common variable names and data structures> in the MVGC
% help page.
%
% <html>
% <table cellpadding=1>
% <tr><th align=left>name</th><th>description</th></tr>
% <tr valign=top><td>X</td><td>multivariate (possibly multi-trial) time series data</td></tr>
% <tr valign=top><td>A,SIG</i></td><td>VAR parameters (coefficients and residuals covariance)</td></tr>
% <tr valign=top><td>G</i></td><td>autocovariance sequence</td></tr>
% <tr valign=top><td>S</i></td><td>cross-power spectral density (cpsd)</td></tr>
% <tr valign=top><td>F</i></td><td>Granger causality (time domain)</td></tr>
% <tr valign=top><td>f</i></td><td>Granger causality (frequency domain)</td></tr>
% </table>
% </html>
%
%% Algorithms
%
% <html>
% <table cellpadding=1>
% <tr><th align=left>&nbsp;</th><th>description</th><th align=left>implementation</th></tr>
% <tr valign=top>
%       <td><b>A1</b></td>
%       <td>estimate autocovariance sequence from time series data</td>
%       <td><a href="tsdata_to_autocov.html">tsdata_to_autocov</a></td>
% </tr>
% <tr valign=top>
%       <td><b>A2</b></td>
%       <td>estimate VAR model from time series data</td>
%       <td><a href="tsdata_to_var.html">tsdata_to_var</a><br><a href="tsdata_to_infocrit.html">tsdata_to_infocrit</a></td>
% </tr>
% <tr valign=top>
%       <td><b>A3</b></td>
%       <td>simulate (multiple) VAR processess</td>
%       <td><a href="var_to_tsdata.html">var_to_tsdata</a><br><a href="var_to_tsdata_nonstat.html">var_to_tsdata_nonstat</a></td>
% </tr>
% <tr valign=top>
%       <td><b>A4</b></td>
%       <td>estimate cpsd from time series data</td>
%       <td><a href="tsdata_to_cpsd.html">tsdata_to_cpsd</a></td>
% </tr>
% <tr valign=top>
%       <td><b>A5</b></td>
%       <td>calculate autocovariance sequence from VAR parameters</td>
%       <td><a href="var_to_autocov.html">var_to_autocov</a></td>
% </tr>
% <tr valign=top>
%       <td><b>A6</b></td>
%       <td>calculate VAR parameters from autocovariance sequence</td>
%       <td><a href="autocov_to_var.html">autocov_to_var</a></td>
% </tr>
% <tr valign=top>
%       <td><b>A7</b></td>
%       <td>calculate VAR parameters from cpsd</td>
%       <td><a href="cpsd_to_var.html">cpsd_to_var</a></td>
% </tr>
% <tr valign=top>
%       <td><b>A8</b></td>
%       <td>calculate cpsd from VAR parameters</td>
%       <td><a href="var_to_cpsd.html">var_to_cpsd</a></td>
% </tr>
% <tr valign=top>
%       <td><b>A9</b></td>
%       <td>calculate cpsd from autocovariance sequence (fft)</td>
%       <td><a href="autocov_to_cpsd.html">autocov_to_cpsd</a></td>
% </tr>
% <tr valign=top>
%       <td><b>A10</b></td>
%       <td>calculate autocovariance sequence from cpsd (ifft)</td>
%       <td><a href="cpsd_to_autocov.html">cpsd_to_autocov</a></td>
% </tr>
% </tr>
% <tr valign=top>
%       <td><b>A11</b></td>
%       <td>transform autocovariance for reduced regression</td>
%       <td><a href="autocov_xform.html">autocov_xform</a></td>
% </tr>
% <tr valign=top>
%       <td><b>A12</b></td>
%       <td>transform cpsd for reduced regression</td>
%       <td><a href="cpsd_xform.html">cpsd_xform</a></td>
% </tr>
% <tr valign=top>
%       <td><b>A13</b></td>
%       <td>calculate time-domain causality from VAR parameters</td>
%       <td>
%           <a href="var_to_autocov.html">var_to_autocov</a><br>
%           <a href="autocov_to_var.html">autocov_to_var</a><br>
%           <a href="autocov_to_mvgc.html">autocov_to_mvgc</a><br>
%           <a href="autocov_to_pwcgc.html">autocov_to_pwcgc</a>
%       </td>
% </tr>
% <tr valign=top>
%       <td><b>A14</b></td>
%       <td>calculate frequency-domain (spectral) causality from VAR parameters and cpsd</td>
%       <td><a href="var_to_autocov.html">var_to_autocov</a><br>
%           <a href="autocov_to_var.html">autocov_to_var</a><br>
%           <a href="autocov_xform.html">autocov_xform</a><br>
%           <a href="autocov_to_smvgc.html">autocov_to_smvgc</a><br>
%           <a href="autocov_to_spwcgc.html">autocov_to_spwcgc</a>
%       </td>
% </tr>
% <tr valign=top>
%       <td><b>A15</b></td>
%       <td>calculate (band-limited) time-domain causality from spectral causality (integrate)</td>
%       <td><a href="smvgc_to_mvgc.html">smvgc_to_mvgc</a></td>
% </tr>
% </table>
% </html>
%
%% References
%
% [1] L. Barnett and A. K. Seth,
% <http://www.sciencedirect.com/science/article/pii/S0165027013003701 The MVGC
%     Multivariate Granger Causality Toolbox: A New Approach to Granger-causal
% Inference>, _J. Neurosci. Methods_ 223, 2014
% [ <matlab:open('mvgc_preprint.pdf') preprint> ].
%
% (C) Lionel Barnett and Anil K. Seth, 2012. See file license.txt in
% installation directory for licensing terms.
%
%%
