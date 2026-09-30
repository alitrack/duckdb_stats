# 统计宏函数参考

`stats_fill.sql` 共 **64 个 SQL 宏**，`LOAD` 后即可用。全部经 **scipy 1.18 数值验证**
（[`tests/verify_fill.py`](../tests/verify_fill.py)，64/64 加载、58 项断言 0 失败）。

## 通用调用约定

| 宏形态 | 写法 |
|---|---|
| 单列标量宏 | `SELECT ci_mean(array_agg(x)) FROM t` |
| 两列标量宏（**必须同序对齐**） | `SELECT cohens_d_2samp(array_agg(a ORDER BY id), array_agg(b ORDER BY id)) FROM t` |
| TABLE 宏（返回多行） | `SELECT * FROM kruskal_wallis((SELECT array_agg(x ORDER BY id) FROM t), (SELECT array_agg(grp ORDER BY id) FROM t))` |

- 入参 `LIST` 由 `array_agg(...)` 聚出；多个 `array_agg` 是独立聚合，**行序不保证**，
  分组宏必须共用同一 `ORDER BY` 键。
- 返回 `STRUCT` 的宏，字段列在下方"返回"列；TABLE 宏每列一字段。
- `[APPROX]` = p 值用 z 近似（DuckDB 无 studentized-range CDF，小样本偏松）。
- 默认 `alpha := 0.05`（置信水平 0.95）。

---

## A. 描述统计补全（9）

返回单个 `DOUBLE`。

| 宏 | 签名 | 说明 |
|---|---|---|
| `stat_pctl` | `(l, p)` | 分位数（线性插值，等价 numpy `percentile`）；`p∈[0,1]` |
| `stat_median` | `(l)` | 中位数 `= stat_pctl(l, 0.5)` |
| `stat_geometric_mean` | `(l)` | 几何均值 `exp(avg(ln x))`；x 须 > 0 |
| `stat_harmonic_mean` | `(l)` | 调和均值 `n / avg(1/x)` |
| `stat_cv` | `(l)` | 变异系数 `stddev_samp / avg` |
| `stat_iqr` | `(l)` | 四分位距 `Q3 − Q1` |
| `stat_mad` | `(l)` `⊞` | 中位数绝对偏差 `median(\|x − median\|)` |
| `stat_mad_scaled` | `(l)` `⊞` | 缩放 MAD `1.4826 · median(\|x − median\|)`（正态下 σ 的无偏估计） |
| `stat_trimmed_mean` | `(l, trim_ratio)` | 截尾均值（两端各切 `trim_ratio` 比例） |

`⊞` = TABLE 宏（`SELECT * FROM stat_mad(array_agg(x)) FROM t`），单值结果。

## B. 置信区间（11）

返回 `STRUCT`。

| 宏 | 签名 | 返回字段 | 说明 |
|---|---|---|---|
| `ci_mean` | `(x, alpha)` | point, se, lower, upper, df, n | 单均值 t 分布 CI |
| `ci_mean_z` | `(x, sigma2, alpha)` | point, se, lower, upper, n | 已知方差 `sigma2` 的 z CI |
| `ci_proportion` | `(k, n, alpha)` | point, lower, upper, n | 二项比例 **Wilson** CI |
| `ci_proportion_wald` | `(k, n, alpha)` | point, lower, upper, n | Wald（正态近似）CI |
| `ci_variance` | `(x, alpha)` | point, lower, upper, df, n | 方差的卡方 CI |
| `ci_mean_diff_welch` | `(x1, x2, alpha)` | diff, se, df, lower, upper | 两独立均值差，**Welch**（不等方差） |
| `ci_mean_diff_pooled` | `(x1, x2, alpha)` | diff, sp2, se, df, lower, upper | 两独立均值差，合并方差（等方差） |
| `ci_mean_diff_paired` | `(d, alpha)` | （同 `ci_mean`） | 配对均值差；`d` = 差值列表 |
| `ci_proportion_diff` | `(k1, n1, k2, n2, alpha)` | diff, se, lower, upper | 两独立比例差 CI |
| `margin_of_error_mean` | `(sd, n, alpha)` | — | 均值误差界限 `qt · sd / √n` |
| `margin_of_error_proportion` | `(n, alpha)` | — | 比例误差界限（最保守 p=0.5） |

## C. 样本量 / 功效分析（6）

正态近似（与 statcpp 口径一致）。返回 `DOUBLE`。

| 宏 | 签名 | 说明 |
|---|---|---|
| `sample_size_moe_prop` | `(moe, alpha)` | 给定误差界限 `moe` 估计比例所需 n |
| `sample_size_moe_mean` | `(moe, sd, alpha)` | 给定 `moe` 与 `sd` 估计均值所需 n（t 迭代收敛） |
| `sample_size_2samp` | `(d, alpha, power)` | 两独立样本、效应量 `d`、目标功效所需**每组** n |
| `sample_size_1samp` | `(d, alpha, power)` | 单样本所需 n |
| `power_2samp` | `(n, d, alpha)` | 已知每组 `n` 的功效 |
| `power_1samp` | `(n, d, alpha)` | 单样本功效 |

## D. 假设检验（5）

| 宏 | 签名 | 返回字段 | 说明 |
|---|---|---|---|
| `z_test` | `(x, mu, sigma, alpha)` | statistic, p_value, alternative | 单样本 z（已知 `sigma`，双侧） |
| `z_test_proportion` | `(k, n, p0, alpha)` | statistic, p_value, p_hat, n | 单比例 z |
| `z_test_proportion_2samp` | `(k1, n1, k2, n2, alpha)` | statistic, p_value | 两比例 z（合并 p̂） |
| `f_test_var` | `(x1, x2, alpha)` | statistic, p_value, df1, df2 | 方差齐性 F 检验（双侧） |
| `kruskal_wallis` | `(xl, gl)` `⊞` | h_statistic, df, p_value | 多组秩和检验（任意组数） |

## E. ANOVA 与事后检验（5）

| 宏 | 签名 | 返回字段 | 说明 |
|---|---|---|---|
| `oneway_anova` | `(xl, gl)` `⊞` | ss_between, ss_within, df_between, df_within, f_stat, p_value, eta_squared, cohen_f | 单因素方差分析 |
| `twoway_anova` | `(xl, al, bl)` `⊞` | ss_a, ss_b, ss_ab, ss_error, df_a, df_b, df_ab, df_error, f_a, p_a, f_b, p_b, f_ab, p_ab | 双因素（不等样本量一般定义） |
| `tukey_hsd` | `(xl, gl, alpha)` `⊞` | group1, group2, mean_diff, se, q_stat, p_z_approx, significant_z95 | Tukey 成对比较 `[APPROX]` |
| `scheffe_p` | `(mean_diff, n1, n2, k, df_error, mse)` | — | Scheffé 事后单对 p 值 |
| `dunnett_approx` | `(diffs, se_each, alpha)` `⊞` | diff, se, z_stat, p_unadjusted, p_bonferroni, significant | 各处理组 vs 对照 `[APPROX]` |

## F. 效应量（10）

返回 `DOUBLE`。

| 宏 | 签名 | 说明 |
|---|---|---|
| `cohens_d` | `(x, mu0)` | `(mean − mu0) / stddev_samp` |
| `cohens_d_2samp` | `(x1, x2)` | 双样本 d（合并标准差） |
| `hedges_g` | `(x, mu0)` | Hedges g（J 校正小样本偏差） |
| `hedges_g_2samp` | `(x1, x2)` | 双样本 Hedges g |
| `glass_delta` | `(x1, x2)` | Glass's Δ（用控制组 sd） |
| `t_to_r` | `(t, df)` | t 统计量 → 相关系数 r |
| `d_to_r` | `(d, df)` | d → r |
| `r_to_d` | `(r)` | r → d |
| `eta_squared_from_f` | `(f, df1, df2)` | F → η² |
| `cohens_h` | `(p1, p2)` | 两比例 Cohen's h |

## G. 分类变量 2×2（4）

列约定 `a`=暴露+事件、`b`=暴露−事件、`c`=对照+事件、`d`=对照−事件。

| 宏 | 签名 | 返回字段 | 说明 |
|---|---|---|---|
| `odds_ratio` | `(a, b, c, d, alpha)` | point, se_log, lower, upper | OR + log 尺度 CI |
| `risk_ratio` | `(a, b, c, d, alpha)` | point, se, lower, upper | RR + CI |
| `risk_difference` | `(a, b, c, d, alpha)` | point, se, lower, upper, nnt | 风险差 + CI + NNT |
| `fisher_exact_p` | `(a, b, c, d)` | — | 超几何**精确**双侧 p |

## H. 时间序列（6）

返回 `DOUBLE`（`acf_vector` 除外）。

| 宏 | 签名 | 说明 |
|---|---|---|
| `acf` | `(x, lag)` | `lag` 阶样本自相关 |
| `acf_vector` | `(x, max_lag)` `⊞` | 1..max_lag 的 ACF 向量 `{lag, r}` |
| `mae` | `(actual, pred)` | 平均绝对误差 |
| `mse` | `(actual, pred)` | 均方误差 |
| `rmse` | `(actual, pred)` | 均方根误差 |
| `mape` | `(actual, pred)` | 平均绝对百分比误差（%） |

## I. 特殊函数（3）

core / stats_duck 都没有、statcpp `special_functions` 里的。返回 `DOUBLE`。

| 宏 | 签名 | 说明 |
|---|---|---|
| `erf` | `(x)` | 误差函数 `= 2·pnorm(x√2) − 1`（验证到 1e-8 vs scipy） |
| `erfc` | `(x)` | 互补误差函数 `= 1 − erf(x)` |
| `norm_sf` | `(x)` | 正态上尾 `1 − pnorm(x)`（大 x 时比 `1−pnorm` 稳） |

> `betainc`/`gammainc` 系列已被 stats_duck 的分布族覆盖：`pbeta(q,a,b)`、`dbeta`、`qbeta`、
> `pgamma`、`qgamma` 等 15+ 族。

## J. 加权统计 + 集中度（5）

| 宏 | 签名 | 返回字段 | 说明 |
|---|---|---|---|
| `weighted_mean` | `(v, wt)` `⊞` | weighted_mean | 加权均值 `Σw·x / Σw` |
| `weighted_stddev` | `(v, wt)` `⊞` | weighted_stddev | 加权样本标准差（`W − 1`） |
| `quantile_weighted` | `(v, wt, p)` `⊞` | quantile | 加权分位数 |
| `gini` | `(v)` `⊞` | gini | Gini 系数（0=完全均等 → 1=极端集中） |
| `hhi` | `(v)` `⊞` | hhi | 赫芬达尔指数（基于占比；1/k=完全竞争 → 1=垄断） |
