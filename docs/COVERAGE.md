# statcpp 386 函数 → DuckDB 三层覆盖矩阵

对照基准：DuckDB **三层叠加** —— `core` + 社区扩展 `stats_duck`(KoliStat) + `duckdb-ml`(34 算法 Rust 扩展)。
每一项都在真实 DuckDB（1.5.x）上实测或读 `duckdb-ml` 源码确认过。

**图例**：✅ 覆盖 · 🟡 部分覆盖/个别函数缺 · ❌ 整块缺（由本仓库 `stats_fill.sql` 或可选 Rust 补）

> 状态 = 三层合计。`补法` 列里"**本仓库已做**"指 `stats_fill.sql` 的 64 宏已实现。

| # | statcpp 模块 | 状态 | 已有覆盖（core / stats_duck / duckdb-ml） | 缺口 | 补法 |
|---|---|:--:|---|---|---|
| 1 | 描述统计 | 🟡 | core: sum/count/avg/median/mode/variance/stddev/pop 变体 | geometric/harmonic/logarithmic/trimmed/weighted_mean, argmin/argmax | `stats_fill.sql` **已做**（A/J 节） |
| 2 | 离散度 | 🟡 | core: variance/stddev 变体, mad | iqr, cv, MAD 缩放, 加权方差/标准差 | `stats_fill.sql` **已做**（A/J 节） |
| 3 | 顺序统计 | 🟡 | core: min/max/quantile/median | 加权分位, 五数概括 | `stats_fill.sql` **已做**（`quantile_weighted`） |
| 4 | 分布形态 | ✅ | core: skewness/kurtosis（pop/sample 变体） | — | — |
| 5 | 相关/协方差 | 🟡 | core: corr/covar_samp/pop；stats_duck: pearson/spearman/kendall_test/corr_matrix | 频率加权协方差 | 低频，SQL 宏可补 |
| 6 | 频数分布 | ✅ | SQL 原生 `GROUP BY` + count/比例/累计窗口 | 形式化函数 | SQL 即可，不必包 |
| 7 | 特殊函数 | 🟡 | core: gamma/lgamma；stats_duck: pgamma/qgamma 等 | erf/erfc, norm_sf | `stats_fill.sql` **已做**（I 节）；betainc 系列已被分布族覆盖 |
| 8 | 随机引擎 | ✅ | core: set_seed/rnd；stats_duck: rnorm/rpois 全 r* 族 | — | — |
| 9 | 连续分布 | ✅ | stats_duck d/p/q/r 四件套：normal/t/chisq/f/gamma/beta/exp/weibull/lognormal/uniform | studentized_range（Tukey p 值用） | 正态近似宏 `[APPROX]`，或 Rust |
| 10 | 离散分布 | ✅ | stats_duck: pois/nbinom/hyper/binom 四件套 | discrete_uniform | SQL `(a+(b-a+1)·rnd)::INT` |
| 11 | 置信区间/估计 | ❌→✅ | 无 | ci_mean/proportion/variance/diff, margin_of_error, sample_size_for_moe | `stats_fill.sql` **已做**（B 节，全标量数学） |
| 12 | 参数检验 | 🟡 | stats_duck: ttest_1samp/2samp/paired, chisq_gof/independence, `adjust_p`(bonferroni/holm/hochberg/fdr/BY) | z_test/z_test_proportion(±2samp), f_test | `stats_fill.sql` **已做**（D 节） |
| 13 | 非参数检验 | 🟡 | stats_duck: shapiro_wilk/anderson_darling/jarque_bera/ks/mann_whitney_u/wilcoxon/sign_test | lilliefors, levene, bartlett, **kruskal_wallis, fisher_exact** | KW + Fisher `stats_fill.sql` **已做**；Levene/Bartlett 可选补 |
| 14 | 效应量 | ❌→✅ | 无（duckdb-ml 也无） | cohens d/hedges g/glass, t/d/r 互转, η², cohens_h, OR/RR | `stats_fill.sql` **已做**（F 节）——与 R 对标最常被问的层 |
| 15 | 重抽样 | 🟡 | stats_duck: bootstrap(aggregate) | permutation_test, bootstrap_bca | 低频；Rust 或 Python 侧 |
| 16 | 功效/样本量 | ❌→✅ | 无 | power/sample_size t-test & proportion | `stats_fill.sql` **已做**（C 节，正态近似） |
| 17 | 线性回归 | 🟡 | stats_duck: lm/lm_fit/lm_summary(HC/cluster 稳健 SE)；duckdb-ml: linear/ridge/lasso/elastic/poly/robust | prediction_interval, VIF, PRESS, LOOCV, 多重共线性 | VIF=corr_matrix 行列式宏；PI 宏；CV 低频 |
| 18 | ANOVA/事后 | ❌→✅ | stats_duck: anova_oneway | two_way, **tukey_hsd**, dunnett, scheffe, ancova, cohens_f | `stats_fill.sql` **已做**（E 节；Tukey/Dunnett `[APPROX]`） |
| 19 | GLM | 🟡 | duckdb-ml: logistic/multinomial/ordinal **训练+预测**（源码确认） | GLM 推断层：SE/z/p/AIC、odds_ratios_ci、过离散、Poisson、伪 R² | Rust 包 `statcpp::glm`（ML 场景 duckdb-ml 够则不做） |
| 20 | 模型选择 | 🟡 | duckdb-ml: ridge/lasso/elastic_net/cv | aic/aicc/bic 标量, press, loocv 独立函数 | SQL 宏 `AIC = k − 2·lnL` 两行 |
| 21 | 距离度量 | 🟡 | core: cosine_distance/similarity | euclidean/manhattan/minkowski/chebyshev(数组), mahalanobis | 欧氏/曼哈顿=SQL 数组展开；mahalanobis 需矩阵逆→Rust |
| 22 | 数值工具 | — | 内部 helper（log1p/expm1/safe_divide） | — | 跳过 |
| 23 | 多元分析 | 🟡 | duckdb-ml: pca/lda；stats_duck: corr_matrix | covariance_matrix 独立函数, mahalanobis, power_iteration | cov_matrix=SQL 宏；mahalanobis 同 #21 |
| 24 | 时间序列 | 🟡 | duckdb-ml: arima 训练+预测 | **acf/pacf**, mae/mse/rmse/mape, moving_average/EMA, diff/lag | `stats_fill.sql` **已做**（H 节 ACF/误差）；pacf 需迭代→Rust |
| 25 | 分类变量 | ❌→✅ | 无 | contingency_table, **OR/RR/NNT + CI** | OR/RR/NNT `stats_fill.sql` **已做**（G 节）；contingency=`PIVOT` 原生 |
| 26 | 生存分析 | 🟡 | duckdb-ml: kaplan_meier/cox 训练 | logrank_test（源码未暴露）, nelson_aalen, median_survival_time | Nelson-Aalen=窗口宏；logrank 需新写 |
| 27 | 稳健统计 | 🟡 | core: mad | mad_scaled, winsorize, cooks_distance, dffits, hodges_lehmann, biweight | MAD 缩放/IQR 离群=宏（`stat_mad_scaled` **已做**）；cooks/dffits 依赖回归杠杆→缓 |
| 28 | 聚类 | ✅ | duckdb-ml: kmeans/dbscan/agglomerative/fuzzy_cmeans/tsne | silhouette_score, kmeans++ 细节 | silhouette=SQL 宏（可选） |
| 29 | 数据操作 | ✅ | SQL 原生：COALESCE/窗口/GROUP BY/UNPIVOT | boxcox, label_encode | SQL 已够；boxcox=宏 |
| 30 | 高级缺失数据 | ❌ | 无 | MCAR 检验, PMM/bootstrap 多重插补, 敏感性分析 | MCAR=卡方宏可做；**PMI 研究级，建议不做**（超出 SQL 扩展合理边界） |

## 一句话结论

statcpp 386 函数 = 分布族(#8-10,✅stats_duck) + 检验(#12-13,✅stats_duck) + ML/聚类/生存(#19,26,28,🟡duckdb-ml)
+ **推断层(#11,14,16,18,25,全部由本仓库 `stats_fill.sql` 补齐)**。30 个模块里 **12 个 ✅、16 个 🟡、2 个
明确不做**（数值工具、高级缺失数据）。`[APPROX]` 仅剩 Tukey/Dunnett 的 p 值。

> 原始数据（CSV 格式，便于程序处理）见 [`statcpp_vs_duckdb.csv`](statcpp_vs_duckdb.csv)。
