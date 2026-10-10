---
name: duckdb-stats-catalog
description: >
  Function catalog for duckdb-stats: all 64 stats_fill.sql macros organized
  by task — descriptives, confidence intervals, power & sample size, tests,
  ANOVA & post-hoc, effect sizes, categorical (OR/RR/NNT), accuracy metrics,
  special functions, weighted stats, concentration (Gini/HHI). Use to pick
  the right macro for a statistical question without reading the repo.
version: 0.1.0
---

# duckdb-stats: the 64 macros by task

All are SQL macros from stats_fill.sql (scipy-verified). Signatures with
return fields: `docs/FUNCTION_REFERENCE.md` in the repo.

## A. Descriptives

`stat_pctl` `stat_median` `stat_geometric_mean` `stat_harmonic_mean`
`stat_cv` (coefficient of variation) `stat_iqr` `stat_mad` `stat_mad_scaled`
`stat_trimmed_mean`

## B. Confidence intervals

`ci_mean` `ci_mean_z` `ci_proportion` `ci_proportion_wald` `ci_variance`
`ci_mean_diff_welch` `ci_mean_diff_pooled` `ci_mean_diff_paired`
`ci_proportion_diff` + `margin_of_error_mean` / `margin_of_error_proportion`

## C. Power & sample size

`power_2samp(n, d)` `power_1samp` — statistical power
`sample_size_2samp` `sample_size_1samp` — plan n for target power
`sample_size_moe_prop` `sample_size_moe_mean` — plan n for target margin

## D. Tests

`z_test` `z_test_proportion` `z_test_proportion_2samp` `f_test_var`
(variance ratio) `kruskal_wallis` (non-parametric multi-group)
(ANOVA itself: `oneway_anova` / `twoway_anova` below; t-tests live in
stats_duck: `t_test` family)

## E. ANOVA & post-hoc

`oneway_anova` `twoway_anova` `tukey_hsd`* `scheffe_p` `dunnett_approx`*
(*[APPROX] z-approximation p-values)

## F. Effect sizes & conversions

`cohens_d` `cohens_d_2samp` `hedges_g` `hedges_g_2samp` `glass_delta`
`t_to_r` `d_to_r` `r_to_d` `eta_squared_from_f` `cohens_h`

## G. Categorical 2×2

`odds_ratio` `risk_ratio` `risk_difference` `fisher_exact_p(a,b,c,d)`
(exact, two-sided)

## H. Time series & accuracy

`acf` / `acf_vector` (autocorrelation) `mae` `mse` `rmse` `mape`

## I. Special functions

`erf` `erfc`

## J. Weighted & concentration

`weighted_mean` `weighted_stddev` `quantile_weighted` `gini` `hhi`
(Herfindahl–Hirschman)

## Picking examples

| Question | Macro |
|---|---|
| Is this A/B lift significant? | `z_test_proportion_2samp` or stats_duck `t_test` |
| How big is the effect? | `cohens_d_2samp` / `hedges_g_2samp` |
| Did I collect enough users? | `power_2samp` → `sample_size_2samp` |
| Which groups differ after ANOVA? | `tukey_hsd` (mind [APPROX]) |
| Market concentration? | `hhi`, `gini` |
