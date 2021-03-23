%% MVGC Multivariate Granger Causality Toolbox - Function Reference
%
%% Toolbox initialisation
%
% * <startup.html startup>
%
%% Demonstration scripts
%
% * <mvgc_demo.html mvgc_demo>
% * <mvgc_demo_bootstrap.html mvgc_demo_bootstrap>
% * <mvgc_demo_GCCA.html mvgc_demo_GCCA>
% * <mvgc_demo_nonstationary.html mvgc_demo_nonstationary>
% * <mvgc_demo_permtest.html mvgc_demo_permtest>
% * <mvgc_demo_stats.html mvgc_demo_stats>
% * <var5_test.html var5_test>
% * <var9_test.html var9_test>
%
%% Core functionality
%
% * <autocov_to_cpsd.html autocov_to_cpsd>
% * <autocov_to_var.html autocov_to_var>
% * <autocov_xform.html autocov_xform>
% * <cpsd_to_autocov.html cpsd_to_autocov>
% * <cpsd_to_var.html cpsd_to_var>
% * <cpsd_xform.html cpsd_xform>
% * <tsdata_to_autocov.html tsdata_to_autocov>
% * <tsdata_to_cpsd.html tsdata_to_cpsd>
% * <tsdata_to_infocrit.html tsdata_to_infocrit>
% * <tsdata_to_var.html tsdata_to_var>
% * <var_to_autocov.html var_to_autocov>
% * <var_to_cpsd.html var_to_cpsd>
% * <var_to_tsdata.html var_to_tsdata>
% * <var_to_tsdata_nonstat.html var_to_tsdata_nonstat>
%
%% Granger causality
%
% * <autocov_to_mvgc.html autocov_to_mvgc>
% * <autocov_to_pwcgc.html autocov_to_pwcgc>
% * <autocov_to_smvgc.html autocov_to_smvgc>
% * <autocov_to_spwcgc.html autocov_to_spwcgc>
% * <smvgc_to_mvgc.html smvgc_to_mvgc>
%
%% Granger causality: GCCA compatibility
%
% * <GCCA_tsdata_to_mvgc.html GCCA_tsdata_to_mvgc>
% * <GCCA_tsdata_to_pwcgc.html GCCA_tsdata_to_pwcgc>
% * <GCCA_tsdata_to_smvgc.html GCCA_tsdata_to_smvgc>
%
%% Granger causality: subsampling
%
% * <bootstrap_tsdata_to_mvgc.html bootstrap_tsdata_to_mvgc>
% * <bootstrap_tsdata_to_pwcgc.html bootstrap_tsdata_to_pwcgc>
% * <bootstrap_tsdata_to_smvgc.html bootstrap_tsdata_to_smvgc>
% * <bootstrap_tsdata_to_spwcgc.html bootstrap_tsdata_to_spwcgc>
% * <empirical_var_to_mvgc.html empirical_var_to_mvgc>
% * <empirical_var_to_pwcgc.html empirical_var_to_pwcgc>
% * <empirical_var_to_smvgc.html empirical_var_to_smvgc>
% * <empirical_var_to_spwcgc.html empirical_var_to_spwcgc>
% * <permtest_tsdata_to_mvgc.html permtest_tsdata_to_mvgc>
% * <permtest_tsdata_to_pwcgc.html permtest_tsdata_to_pwcgc>
% * <permtest_tsdata_to_smvgc.html permtest_tsdata_to_smvgc>
% * <permtest_tsdata_to_spwcgc.html permtest_tsdata_to_spwcgc>
%
%% Statistics
%
% * <consistency.html consistency>
% * <demean.html demean>
% * <empirical_cdf.html empirical_cdf>
% * <empirical_cdfi.html empirical_cdfi>
% * <empirical_confint.html empirical_confint>
% * <empirical_cval.html empirical_cval>
% * <empirical_pval.html empirical_pval>
% * <infocrit.html infocrit>
% * <mvgc_cdf.html mvgc_cdf>
% * <mvgc_cdfi.html mvgc_cdfi>
% * <mvgc_confint.html mvgc_confint>
% * <mvgc_cval.html mvgc_cval>
% * <mvgc_pval.html mvgc_pval>
% * <rsquared.html rsquared>
% * <significance.html significance>
% * <whiteness.html whiteness>
%
%% Utilities
%
% * <bfft.html bfft>
% * <bifft.html bifft>
% * <cov2corr.html cov2corr>
% * <dlyap_aitr.html dlyap_aitr>
% * <genvar.html genvar>
% * <get_crand.html get_crand>
% * <get_hostname.html get_hostname>
% * <get_urand.html get_urand>
% * <helpon.html helpon>
% * <isbad.html idbad>
% * <isint.html isint>
% * <isposdef.html isposdef>
% * <maxabs.html maxabs>
% * <mvdetrend.html mvdetrend>
% * <mvdiff.html mvdiff>
% * <mvgc_makemex.html mvgc_makemex>
% * <plot_autocov.html plot_autocov>
% * <plot_confints.html plot_confints>
% * <plot_cpsd.html plot_cpsd>
% * <plot_pw.html plot_pw>
% * <plot_spw.html plot_spw>
% * <plot_varcoeffs.html plot_varcoeffs>
% * <ptic.html ptic>
% * <ptoc.html ptoc>
% * <quads.html quads>
% * <quadsr.html quadsr>
% * <rng_restore.html rng_restore>
% * <rng_save.html rng_save>
% * <rng_seed.html rng_seed>
% * <secs2hms.html secs2hms>
% * <sfreqs.html sfreqs>
% * <timestr.html timestr>
% * <trfun2var.html trfun2var>
% * <var2trfun.html var2trfun>
% * <var_decay.html var_decay>
% * <var_info.html var_info>
% * <var_normalise.html var_normalise>
% * <var_specrad.html var_specrad>
% * <warn_if.html warn_if>
% * <warn_supp.html warn_supp>
% * <warn_test.html warn_test>
%
%% Experimental
%
% * <mvgc_adf.html mvgc_adf>
% * <mvgc_kpss.html mvgc_kpss>
% * <tsdata_to_autocov_debias.html tsdata_to_autocov_debias>
%
% (C) Lionel Barnett and Anil K. Seth, 2012. See file license.txt in
% installation directory for licensing terms.
%
%%
