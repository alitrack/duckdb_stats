# 统计宏函数参考

`stats_fill.sql` 共 **64 个 SQL 宏**，加载后即可用。全部经 **scipy 1.18 数值验证**（`tests/verify_fill.py`，0 失败）。

## A. 描述统计补全

### `stat_pctl(l, p)`

### `stat_median(l)`

### `stat_geometric_mean(l)`

### `stat_harmonic_mean(l)`

### `stat_cv(l)`

### `stat_iqr(l)`

### `stat_mad(l)`（TABLE 宏：`SELECT * FROM ...`）

### `stat_mad_scaled(l)`（TABLE 宏：`SELECT * FROM ...`）

### `stat_trimmed_mean(l, trim_ratio)`

## B. 置信区间

### `ci_mean(x, alpha (默认 0.05))`

### `ci_mean_z(x, sigma2, alpha (默认 0.05))`

### `ci_proportion(k, n, alpha (默认 0.05))`

### `ci_proportion_wald(k, n, alpha (默认 0.05))`

### `ci_variance(x, alpha (默认 0.05))`

### `ci_mean_diff_welch(x1, x2, alpha (默认 0.05))`

### `ci_mean_diff_pooled(x1, x2, alpha (默认 0.05))`

### `ci_mean_diff_paired(d, alpha (默认 0.05))`

### `ci_proportion_diff(k1, n1, k2, n2, alpha (默认 0.05))`

### `margin_of_error_mean(sd, n, alpha (默认 0.05))`

### `margin_of_error_proportion(n, alpha (默认 0.05))`

## C. 样本量 / 功效分析 (正态近似, 与 statcpp 口径一致)

### `sample_size_moe_prop(moe, alpha (默认 0.05))`

### `sample_size_moe_mean(moe, sd, alpha (默认 0.05))`

### `sample_size_2samp(d, alpha (默认 0.05), power (默认 0.80))`

### `sample_size_1samp(d, alpha (默认 0.05), power (默认 0.80))`

### `power_2samp(n, d, alpha (默认 0.05))`

### `power_1samp(n, d, alpha (默认 0.05))`

## D. 缺失的假设检验

### `z_test(x, mu, sigma, alpha (默认 0.05))`

### `z_test_proportion(k, n, p0, alpha (默认 0.05))`

### `z_test_proportion_2samp(k1, n1, k2, n2, alpha (默认 0.05))`

### `f_test_var(x1, x2, alpha (默认 0.05))`

### `kruskal_wallis(xl, gl)`（TABLE 宏：`SELECT * FROM ...`）

## E. ANOVA 与事后检验

### `oneway_anova(xl, gl)`（TABLE 宏：`SELECT * FROM ...`）

### `twoway_anova(xl, al, bl)`（TABLE 宏：`SELECT * FROM ...`）

### `tukey_hsd(xl, gl, alpha (默认 0.05))`（TABLE 宏：`SELECT * FROM ...`） `[APPROX]`

### `scheffe_p(mean_diff, n1, n2, k, df_error, mse)`

### `dunnett_approx(diffs, se_each, alpha (默认 0.05))`（TABLE 宏：`SELECT * FROM ...`） `[APPROX]`

## F. 效应量

### `cohens_d(x, mu0 (默认 0.0))`

### `cohens_d_2samp(x1, x2)`

### `hedges_g(x, mu0 (默认 0.0))`

### `hedges_g_2samp(x1, x2)`

### `glass_delta(x1, x2)`

### `t_to_r(t, df)`

### `d_to_r(d, df)`

### `r_to_d(r)`

### `eta_squared_from_f(f, df1, df2)`

### `cohens_h(p1, p2)`

## G. 分类变量 (2x2 表)

### `odds_ratio(a, b, c, d, alpha (默认 0.05))`

### `risk_ratio(a, b, c, d, alpha (默认 0.05))`

### `risk_difference(a, b, c, d, alpha (默认 0.05))`

### `fisher_exact_p(a, b, c, d)`

## H. 时间序列

### `acf(x, lag)`

### `acf_vector(x, max_lag)`（TABLE 宏：`SELECT * FROM ...`）

### `mae(actual, pred)`

### `mse(actual, pred)`

### `rmse(actual, pred)`

### `mape(actual, pred)`

## I. 特殊函数 (statcpp special_functions 里 core/stats_duck 没有的)

### `erf(x)`

### `erfc(x)`

### `norm_sf(x)`

## J. 加权统计 + 集中度 (statcpp weighted/concentration 模块)

### `weighted_mean(v, wt)`（TABLE 宏：`SELECT * FROM ...`）

### `weighted_stddev(v, wt)`（TABLE 宏：`SELECT * FROM ...`）

### `quantile_weighted(v, wt, p (默认 0.5))`（TABLE 宏：`SELECT * FROM ...`）

### `gini(v)`（TABLE 宏：`SELECT * FROM ...`）

### `hhi(v)`（TABLE 宏：`SELECT * FROM ...`）
