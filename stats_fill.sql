-- ============================================================================
-- stats_fill.sql — DuckDB 统计缺口补丁 (SQL 宏层)
-- 补齐 statcpp 有而 DuckDB (core + stats_duck + duckdb-ml) 没有的推断统计层:
--   描述补全 / 置信区间 / 样本量·功效 / 缺失检验 / 效应量 / 分类变量 /
--   Kruskal-Wallis / ANOVA(单·双因素) / 事后检验 / 时间序列
--
-- 依赖:
--   stats_duck 社区扩展:  INSTALL stats_duck FROM community; LOAD stats_duck;
--     (提供 pt/pf/pchisq/qchisq/qnorm/pnorm — 注意: pf/pt 是 CDF, 上尾用 1-*)
--   DuckDB >= 1.4
--
-- 调用约定 (逐条验证过, 2026-09-30, duckdb 1.5.5/1.5.6):
--   * 单列函数收 LIST:      SELECT ci_mean(array_agg(x)) FROM t
--     (无 list(); array_agg 不能嵌套在另一个聚合里 → 需要时包子查询)
--   * 两列函数收两个 LIST:   SELECT ci_mean_diff_welch(array_agg(a), array_agg(b)) FROM t
--   * 分组函数是 TABLE 宏收 LIST, 参数必须是【标量子查询】, 且各 list 必须同序对齐:
--       多个 array_agg 是独立聚合, 顺序不保证! 必须用同一个 ORDER BY 键:
--       SELECT * FROM oneway_anova((SELECT array_agg(x ORDER BY x) FROM t),
--                                  (SELECT array_agg(grp ORDER BY x) FROM t))
--   * 多值结果返回 STRUCT (字段见各宏); 单值返回 DOUBLE
--   * 默认参数用 SQL 宏原生 `alpha := 0.05`
--
-- 标注 [APPROX] 的 p 值用了 z 近似 (DuckDB 无 studentized-range CDF), 小样本偏松;
--   精确值请查统计表或用 scipy/Python.
-- 数值验证: verify_fill.py (scipy 1.18.1 交叉核对)
-- ============================================================================

-- ============ A. 描述统计补全 ============

CREATE OR REPLACE MACRO stat_pctl(l, p) AS (
  WITH s AS (SELECT list_sort(l) AS s, (len(l)-1)*p::DOUBLE AS k),
       i AS (SELECT floor(s.k)::BIGINT AS lo, ceil(s.k)::BIGINT AS hi, s.k - floor(s.k) AS frac FROM s)
  SELECT s.s[i.lo+1] + (s.s[i.hi+1] - s.s[i.lo+1]) * i.frac
  FROM s CROSS JOIN i
);

CREATE OR REPLACE MACRO stat_median(l) AS stat_pctl(l, 0.5);

CREATE OR REPLACE MACRO stat_geometric_mean(l) AS
  exp(list_sum(list_transform(l, v -> ln(v))) / len(l));

CREATE OR REPLACE MACRO stat_harmonic_mean(l) AS
  len(l) / list_sum(list_transform(l, v -> 1.0 / v));

CREATE OR REPLACE MACRO stat_cv(l) AS list_stddev_samp(l) / list_avg(l);

CREATE OR REPLACE MACRO stat_iqr(l) AS stat_pctl(l, 0.75) - stat_pctl(l, 0.25);

-- MAD 用表宏 (lambda 不能引用派生标量; unnest+聚合更直接)
CREATE OR REPLACE MACRO stat_mad(l) AS TABLE
  WITH u AS (SELECT unnest(l) AS x),
       m AS (SELECT median(x) AS m FROM u)
  SELECT median(abs(u.x - m.m)) AS mad FROM u CROSS JOIN m;

CREATE OR REPLACE MACRO stat_mad_scaled(l) AS TABLE
  WITH u AS (SELECT unnest(l) AS x),
       m AS (SELECT median(x) AS m FROM u)
  SELECT 1.4826 * median(abs(u.x - m.m)) AS mad FROM u CROSS JOIN m;

CREATE OR REPLACE MACRO stat_trimmed_mean(l, trim_ratio) AS (
  WITH s AS (SELECT list_sort(l) AS s, len(l) AS n)
  SELECT list_sum(list_slice(s.s, floor(s.n*trim_ratio)::INT+1, s.n - floor(s.n*trim_ratio)::INT))
         / (s.n - 2*floor(s.n*trim_ratio)::INT)
  FROM s
);

-- ============ B. 置信区间 ============

CREATE OR REPLACE MACRO ci_mean(x, alpha := 0.05) AS (
  WITH z AS (SELECT list_avg(x) AS m, list_stddev_samp(x) AS s, len(x) AS n)
  SELECT struct_pack(point := z.m, se := z.s / sqrt(z.n),
    lower := z.m - qt(1 - alpha/2, z.n - 1) * z.s / sqrt(z.n),
    upper := z.m + qt(1 - alpha/2, z.n - 1) * z.s / sqrt(z.n),
    df := z.n - 1, n := z.n)
  FROM z
);

CREATE OR REPLACE MACRO ci_mean_z(x, sigma2, alpha := 0.05) AS (
  WITH z AS (SELECT list_avg(x) AS m, len(x) AS n)
  SELECT struct_pack(point := z.m, se := sqrt(sigma2 / z.n),
    lower := z.m - qnorm(1 - alpha/2) * sqrt(sigma2 / z.n),
    upper := z.m + qnorm(1 - alpha/2) * sqrt(sigma2 / z.n), n := z.n)
  FROM z
);

CREATE OR REPLACE MACRO ci_proportion(k, n, alpha := 0.05) AS (
  WITH z AS (SELECT k::DOUBLE/n AS p, n::DOUBLE AS n, qnorm(1 - alpha/2) AS q)
  SELECT struct_pack(point := z.p,
    lower := (z.p + z.q*z.q/(2*z.n) - z.q*sqrt(z.p*(1-z.p)/z.n + z.q*z.q/(4*z.n*z.n))) / (1 + z.q*z.q/z.n),
    upper := (z.p + z.q*z.q/(2*z.n) + z.q*sqrt(z.p*(1-z.p)/z.n + z.q*z.q/(4*z.n*z.n))) / (1 + z.q*z.q/z.n),
    n := z.n)
  FROM z
);

CREATE OR REPLACE MACRO ci_proportion_wald(k, n, alpha := 0.05) AS (
  WITH z AS (SELECT k::DOUBLE/n AS p, n::DOUBLE AS n, qnorm(1 - alpha/2) AS q)
  SELECT struct_pack(point := z.p,
    lower := z.p - z.q*sqrt(z.p*(1-z.p)/z.n),
    upper := z.p + z.q*sqrt(z.p*(1-z.p)/z.n), n := z.n)
  FROM z
);

CREATE OR REPLACE MACRO ci_variance(x, alpha := 0.05) AS (
  WITH z AS (SELECT list_stddev_samp(x) AS s, len(x) AS n)
  SELECT struct_pack(point := z.s*z.s,
    lower := (z.n-1)*z.s*z.s / qchisq(1 - alpha/2, z.n - 1),
    upper := (z.n-1)*z.s*z.s / qchisq(alpha/2, z.n - 1),
    df := z.n - 1, n := z.n)
  FROM z
);

CREATE OR REPLACE MACRO ci_mean_diff_welch(x1, x2, alpha := 0.05) AS (
  WITH z AS (SELECT list_avg(x1) AS m1, list_stddev_samp(x1) AS s1, len(x1) AS n1,
                    list_avg(x2) AS m2, list_stddev_samp(x2) AS s2, len(x2) AS n2),
       w AS (SELECT (z.m1 - z.m2) AS d,
                   sqrt(z.s1*z.s1/z.n1 + z.s2*z.s2/z.n2) AS se,
                   pow(z.s1*z.s1/z.n1 + z.s2*z.s2/z.n2, 2) /
                   (pow(z.s1*z.s1/z.n1, 2)/(z.n1-1) + pow(z.s2*z.s2/z.n2, 2)/(z.n2-1)) AS df
             FROM z)
  SELECT struct_pack(diff := w.d, se := w.se, df := w.df,
    lower := w.d - qt(1-alpha/2, w.df)*w.se,
    upper := w.d + qt(1-alpha/2, w.df)*w.se)
  FROM w
);

CREATE OR REPLACE MACRO ci_mean_diff_pooled(x1, x2, alpha := 0.05) AS (
  WITH z AS (SELECT list_avg(x1) AS m1, list_stddev_samp(x1) AS s1, len(x1) AS n1,
                    list_avg(x2) AS m2, list_stddev_samp(x2) AS s2, len(x2) AS n2),
       w AS (SELECT (z.m1 - z.m2) AS d,
                   ((z.n1-1)*z.s1*z.s1 + (z.n2-1)*z.s2*z.s2) / (z.n1+z.n2-2) AS sp2,
                   z.n1 + z.n2 - 2 AS df FROM z)
  SELECT struct_pack(diff := w.d, sp2 := w.sp2, se := sqrt(w.sp2*(1.0/z.n1 + 1.0/z.n2)),
    df := w.df,
    lower := w.d - qt(1-alpha/2, w.df)*sqrt(w.sp2*(1.0/z.n1 + 1.0/z.n2)),
    upper := w.d + qt(1-alpha/2, w.df)*sqrt(w.sp2*(1.0/z.n1 + 1.0/z.n2)))
  FROM w
);

CREATE OR REPLACE MACRO ci_mean_diff_paired(d, alpha := 0.05) AS ci_mean(d, alpha);

CREATE OR REPLACE MACRO ci_proportion_diff(k1, n1, k2, n2, alpha := 0.05) AS (
  WITH z AS (SELECT k1::DOUBLE/n1 AS p1, k2::DOUBLE/n2 AS p2, n1::DOUBLE AS n1, n2::DOUBLE AS n2,
                    qnorm(1 - alpha/2) AS q)
  SELECT struct_pack(diff := z.p1 - z.p2,
    se := sqrt(z.p1*(1-z.p1)/z.n1 + z.p2*(1-z.p2)/z.n2),
    lower := z.p1 - z.p2 - z.q*sqrt(z.p1*(1-z.p1)/z.n1 + z.p2*(1-z.p2)/z.n2),
    upper := z.p1 - z.p2 + z.q*sqrt(z.p1*(1-z.p1)/z.n1 + z.p2*(1-z.p2)/z.n2))
  FROM z
);

CREATE OR REPLACE MACRO margin_of_error_mean(sd, n, alpha := 0.05) AS
  qt(1 - alpha/2, n - 1) * sd / sqrt(n);

CREATE OR REPLACE MACRO margin_of_error_proportion(n, alpha := 0.05) AS
  qnorm(1 - alpha/2) * sqrt(0.25 / n);

-- ============ C. 样本量 / 功效分析 (正态近似, 与 statcpp 口径一致) ============

CREATE OR REPLACE MACRO sample_size_moe_prop(moe, alpha := 0.05) AS
  ceil(pow(qnorm(1 - alpha/2), 2) * 0.25 / (moe * moe));

CREATE OR REPLACE MACRO sample_size_moe_mean(moe, sd, alpha := 0.05) AS (
  WITH n1 AS (SELECT ceil(pow(qnorm(1-alpha/2) * sd / moe, 2)) AS n),
       n2 AS (SELECT ceil(pow(qt(1-alpha/2, n - 1) * sd / moe, 2)) AS n FROM n1),
       n3 AS (SELECT ceil(pow(qt(1-alpha/2, n - 1) * sd / moe, 2)) AS n FROM n2)
  SELECT n FROM n3
);

CREATE OR REPLACE MACRO sample_size_2samp(d, alpha := 0.05, power := 0.80) AS
  ceil(2.0 * pow(qnorm(1 - alpha/2) + qnorm(power), 2) / (d * d));

CREATE OR REPLACE MACRO sample_size_1samp(d, alpha := 0.05, power := 0.80) AS
  ceil(pow(qnorm(1 - alpha/2) + qnorm(power), 2) / (d * d));

CREATE OR REPLACE MACRO power_2samp(n, d, alpha := 0.05) AS
  pnorm(d * sqrt(n / 2.0) - qnorm(1 - alpha/2));

CREATE OR REPLACE MACRO power_1samp(n, d, alpha := 0.05) AS
  pnorm(d * sqrt(n) - qnorm(1 - alpha/2));

-- ============ D. 缺失的假设检验 ============

CREATE OR REPLACE MACRO z_test(x, mu, sigma, alpha := 0.05) AS (
  WITH z AS (SELECT (list_avg(x) - mu) / (sigma / sqrt(len(x))) AS zstat)
  SELECT struct_pack(statistic := z.zstat,
    p_value := 2 * (1 - pnorm(abs(z.zstat))), alternative := 'two-sided')
  FROM z
);

CREATE OR REPLACE MACRO z_test_proportion(k, n, p0, alpha := 0.05) AS (
  WITH z AS (SELECT (k::DOUBLE/n - p0) / sqrt(p0*(1-p0)/n) AS zstat)
  SELECT struct_pack(statistic := z.zstat,
    p_value := 2 * (1 - pnorm(abs(z.zstat))), p_hat := k::DOUBLE/n, n := n)
  FROM z
);

CREATE OR REPLACE MACRO z_test_proportion_2samp(k1, n1, k2, n2, alpha := 0.05) AS (
  WITH z AS (SELECT k1::DOUBLE/n1 AS p1, k2::DOUBLE/n2 AS p2, n1::DOUBLE AS n1, n2::DOUBLE AS n2,
                    (k1+k2)::DOUBLE/(n1+n2) AS pp)
  SELECT struct_pack(
    statistic := (z.p1 - z.p2) / sqrt(z.pp*(1-z.pp)*(1.0/z.n1 + 1.0/z.n2)),
    p_value := 2 * (1 - pnorm(abs((z.p1 - z.p2) / sqrt(z.pp*(1-z.pp)*(1.0/z.n1 + 1.0/z.n2))))))
  FROM z
);

CREATE OR REPLACE MACRO f_test_var(x1, x2, alpha := 0.05) AS (
  WITH z AS (SELECT pow(list_stddev_samp(x1),2) AS v1, pow(list_stddev_samp(x2),2) AS v2,
                    len(x1) AS n1, len(x2) AS n2)
  SELECT struct_pack(statistic := z.v1 / z.v2,
    p_value := 2 * least(pf(z.v1/z.v2, z.n1-1, z.n2-1), 1 - pf(z.v1/z.v2, z.n1-1, z.n2-1)),
    df1 := z.n1 - 1, df2 := z.n2 - 1)
  FROM z
);

-- Kruskal-Wallis (TABLE 宏; xl=值列表, gl=同序分组标签列表, 任意组数)
CREATE OR REPLACE MACRO kruskal_wallis(xl, gl) AS TABLE
  WITH base AS (SELECT xl[i] AS x, gl[i] AS grp FROM unnest(range(1, len(xl) + 1)) t(i)),
       r AS (SELECT x, grp, rank() OVER (ORDER BY x) AS rk FROM base),
       g AS (SELECT grp, avg(rk) AS rk_mean, count(*)::DOUBLE AS n FROM r GROUP BY grp),
       tot AS (SELECT avg(rk) AS rk_all, count(*)::DOUBLE AS n_all FROM r),
       hs AS (SELECT sum(pow(rk_mean - tot.rk_all, 2) * n) AS s, count(*) AS k
              FROM g CROSS JOIN tot)
  SELECT 12.0 * hs.s / (tot.n_all*(tot.n_all+1)) AS h_statistic,
         hs.k - 1 AS df,
         1 - pchisq(12.0 * hs.s / (tot.n_all*(tot.n_all+1)), hs.k - 1) AS p_value
  FROM hs CROSS JOIN tot

-- ============ E. ANOVA 与事后检验 ============

CREATE OR REPLACE MACRO oneway_anova(xl, gl) AS TABLE
  WITH base AS (SELECT xl[i] AS x, gl[i] AS grp FROM unnest(range(1, len(xl) + 1)) t(i)),
       gm AS (SELECT grp, avg(x) AS g_mean, count(*)::DOUBLE AS n, var_samp(x) AS v
              FROM base GROUP BY grp),
       grand AS (SELECT avg(x) AS g_all, count(*)::DOUBLE AS n_all FROM base),
       sst AS (SELECT sum(pow(g_mean - grand.g_all, 2)*n) AS ssb, sum(v*(n-1)) AS ssw,
                      count(*) AS k
               FROM gm CROSS JOIN grand)
  SELECT sst.ssb AS ss_between, sst.ssw AS ss_within, sst.k-1 AS df_between,
         grand.n_all - sst.k AS df_within,
         (sst.ssb/(sst.k-1))/(sst.ssw/(grand.n_all - sst.k)) AS f_stat,
         1 - pf((sst.ssb/(sst.k-1))/(sst.ssw/(grand.n_all - sst.k)), sst.k-1, grand.n_all - sst.k) AS p_value,
         sst.ssb/(sst.ssb+sst.ssw) AS eta_squared,
         sqrt(sst.ssb/(sst.ssb+sst.ssw) / (1 - sst.ssb/(sst.ssb+sst.ssw))) AS cohen_f
  FROM sst CROSS JOIN grand

-- 两因素 ANOVA (xl 值, al a因素, bl b因素; 交叉设计)
CREATE OR REPLACE MACRO twoway_anova(xl, al, bl) AS TABLE
  WITH base AS (SELECT xl[i] AS x, al[i] AS a, bl[i] AS b
                FROM unnest(range(1, len(xl) + 1)) t(i)),
       fa AS (SELECT a, avg(x) AS ma, count(*)::DOUBLE AS na FROM base GROUP BY a),
       fb AS (SELECT b, avg(x) AS mb, count(*)::DOUBLE AS nb FROM base GROUP BY b),
       cell AS (SELECT a, b, avg(x) AS m, count(*)::DOUBLE AS n FROM base GROUP BY a, b),
       gm AS (SELECT avg(x) AS mg, count(*)::DOUBLE AS ng FROM base),
       ssa AS (SELECT sum(pow(fa.ma - gm.mg, 2)*fa.na) AS s FROM fa CROSS JOIN gm),
       ssb AS (SELECT sum(pow(fb.mb - gm.mg, 2)*fb.nb) AS s FROM fb CROSS JOIN gm),
       -- 一般定义 (不等样本量也正确): sum (m_ab - m_a. - m_.b + grand)^2 * n_ab
       ssab AS (SELECT sum(pow(cell.m - fa.ma - fb.mb + gm.mg, 2) * cell.n) AS v
                FROM cell JOIN fa ON cell.a = fa.a JOIN fb ON cell.b = fb.b CROSS JOIN gm),
       sse AS (SELECT sum(pow(b.x - c.m, 2)) AS s FROM base b JOIN cell c ON b.a = c.a AND b.b = c.b),
       dim AS (SELECT (SELECT count(*) FROM fa) AS k, (SELECT count(*) FROM fb) AS j,
                     (SELECT ng FROM gm) AS ng),
       f AS (SELECT ssa.s/(dim.k-1) / (sse.s/(dim.ng - dim.k*dim.j)) AS fa,
                    ssb.s/(dim.j-1) / (sse.s/(dim.ng - dim.k*dim.j)) AS fb,
                    ssab.v/((dim.k-1)*(dim.j-1)) / (sse.s/(dim.ng - dim.k*dim.j)) AS fab
             FROM ssa CROSS JOIN ssb CROSS JOIN ssab CROSS JOIN sse CROSS JOIN dim)
  SELECT ssa.s AS ss_a, ssb.s AS ss_b, ssab.v AS ss_ab, sse.s AS ss_error,
         dim.k-1 AS df_a, dim.j-1 AS df_b, (dim.k-1)*(dim.j-1) AS df_ab,
         dim.ng - dim.k*dim.j AS df_error,
         f.fa AS f_a, 1 - pf(f.fa, dim.k-1, dim.ng - dim.k*dim.j) AS p_a,
         f.fb AS f_b, 1 - pf(f.fb, dim.j-1, dim.ng - dim.k*dim.j) AS p_b,
         f.fab AS f_ab, 1 - pf(f.fab, (dim.k-1)*(dim.j-1), dim.ng - dim.k*dim.j) AS p_ab
  FROM ssa CROSS JOIN ssb CROSS JOIN ssab CROSS JOIN sse CROSS JOIN dim CROSS JOIN f

-- Tukey HSD (TABLE 宏; 返回所有组对). q 统计量精确, p 为 z 近似 [APPROX]
CREATE OR REPLACE MACRO tukey_hsd(xl, gl, alpha := 0.05) AS TABLE
  WITH base AS (SELECT xl[i] AS x, gl[i] AS grp FROM unnest(range(1, len(xl) + 1)) t(i)),
       gm AS (SELECT grp, avg(x) AS g_mean, count(*)::DOUBLE AS n, var_samp(x) AS v
              FROM base GROUP BY grp),
       gsum AS (SELECT sum(n) AS n_all, count(*) AS k FROM gm),
       mse AS (SELECT sum(v*(n-1)) / ANY_VALUE(gsum.n_all - gsum.k) AS mse FROM gm CROSS JOIN gsum)
  SELECT g1.grp AS group1, g2.grp AS group2,
         g1.g_mean - g2.g_mean AS mean_diff,
         sqrt(mse.mse/2.0*(1.0/g1.n + 1.0/g2.n)) AS se,
         abs(g1.g_mean - g2.g_mean) / (sqrt(mse.mse/2.0*(1.0/g1.n + 1.0/g2.n))) AS q_stat,
         2*(1-pnorm(abs(g1.g_mean - g2.g_mean) / sqrt(mse.mse/2.0*(1.0/g1.n + 1.0/g2.n))))
           AS p_z_approx,
         abs(g1.g_mean - g2.g_mean) > 1.96*sqrt(mse.mse/2.0*(1.0/g1.n + 1.0/g2.n))
           AS significant_z95
  FROM gm g1, gm g2, mse
  WHERE g1.grp < g2.grp

-- Scheffé 事后 (单值; mean_diff/n1/n2/组数k/df_error/mse)
CREATE OR REPLACE MACRO scheffe_p(mean_diff, n1, n2, k, df_error, mse) AS (
  WITH z AS (SELECT (k-1)*pow(mean_diff/sqrt(mse*(1.0/n1 + 1.0/n2)), 2) AS f,
                    k-1 AS df1, df_error AS df2)
  SELECT 1 - pf(z.f, z.df1, z.df2) AS p_value
  FROM z
);

-- Dunnett vs 对照 (TABLE 宏; diffs=各处理组对对照的均值差, se_each=对应 SE)
CREATE OR REPLACE MACRO dunnett_approx(diffs, se_each, alpha := 0.05) AS TABLE
  WITH d AS (SELECT diffs[i] AS diff, se_each[i] AS se
             FROM unnest(range(1, len(diffs) + 1)) t(i))
  SELECT diff, se, diff/se AS z_stat,
         2*(1-pnorm(abs(diff/se))) AS p_unadjusted,
         2*(1-pnorm(abs(diff/se)))/len(diffs) AS p_bonferroni,
         2*(1-pnorm(abs(diff/se))) < alpha/len(diffs) AS significant
  FROM d

-- ============ F. 效应量 ============

CREATE OR REPLACE MACRO cohens_d(x, mu0 := 0.0) AS
  (list_avg(x) - mu0) / list_stddev_samp(x);

CREATE OR REPLACE MACRO cohens_d_2samp(x1, x2) AS (
  WITH z AS (SELECT list_avg(x1) AS m1, list_stddev_samp(x1) AS s1, len(x1) AS n1,
                    list_avg(x2) AS m2, list_stddev_samp(x2) AS s2, len(x2) AS n2)
  SELECT (z.m1 - z.m2) /
         sqrt(((z.n1-1)*z.s1*z.s1 + (z.n2-1)*z.s2*z.s2)/(z.n1+z.n2-2))
  FROM z
);

CREATE OR REPLACE MACRO hedges_g(x, mu0 := 0.0) AS (
  WITH z AS (SELECT (list_avg(x)-mu0)/list_stddev_samp(x) AS d, len(x) AS n)
  SELECT (1 - 3.0/(4*(z.n-1) - 1)) * z.d
  FROM z
);

CREATE OR REPLACE MACRO hedges_g_2samp(x1, x2) AS (
  WITH z AS (SELECT list_avg(x1) AS m1, list_stddev_samp(x1) AS s1, len(x1) AS n1,
                    list_avg(x2) AS m2, list_stddev_samp(x2) AS s2, len(x2) AS n2)
  SELECT (1 - 3.0/(4*(z.n1+z.n2-2) - 1)) *
         (z.m1 - z.m2) / sqrt(((z.n1-1)*z.s1*z.s1 + (z.n2-1)*z.s2*z.s2)/(z.n1+z.n2-2))
  FROM z
);

CREATE OR REPLACE MACRO glass_delta(x1, x2) AS
  (list_avg(x1) - list_avg(x2)) / list_stddev_samp(x2);

CREATE OR REPLACE MACRO t_to_r(t, df) AS sqrt(t*t / (t*t + df));
CREATE OR REPLACE MACRO d_to_r(d, df) AS d / sqrt(d*d + df);
CREATE OR REPLACE MACRO r_to_d(r) AS r / sqrt(1 - r*r);
CREATE OR REPLACE MACRO eta_squared_from_f(f, df1, df2) AS (f*df1) / (f*df1 + df2);
CREATE OR REPLACE MACRO cohens_h(p1, p2) AS 2*asin(sqrt(p1)) - 2*asin(sqrt(p2));

-- ============ G. 分类变量 (2x2 表) ============

-- a=暴露+事件, b=暴露-事件, c=对照+事件, d=对照-事件
CREATE OR REPLACE MACRO odds_ratio(a, b, c, d, alpha := 0.05) AS (
  WITH z AS (SELECT (a*d)::DOUBLE/(b*c) AS or_, ln((a*d)::DOUBLE/(b*c)) AS lor,
                    sqrt(1.0/a + 1.0/b + 1.0/c + 1.0/d) AS se, qnorm(1 - alpha/2) AS q)
  SELECT struct_pack(point := z.or_, se_log := z.se,
    lower := exp(z.lor - z.q*z.se), upper := exp(z.lor + z.q*z.se))
  FROM z
);

CREATE OR REPLACE MACRO risk_ratio(a, b, c, d, alpha := 0.05) AS (
  WITH z AS (SELECT (a/(a+b))::DOUBLE/(c/(c+d)) AS rr,
                    sqrt(1.0/a - 1.0/(a+b) + 1.0/c - 1.0/(c+d)) AS se,
                    qnorm(1 - alpha/2) AS q)
  SELECT struct_pack(point := z.rr, se := z.se,
    lower := z.rr*exp(-z.q*z.se), upper := z.rr*exp(z.q*z.se))
  FROM z
);

CREATE OR REPLACE MACRO risk_difference(a, b, c, d, alpha := 0.05) AS (
  WITH z AS (SELECT a::DOUBLE/(a+b) - c::DOUBLE/(c+d) AS rd,
                    sqrt(a*(a+b)/pow(a+b,3) + c*(c+d)/pow(c+d,3)) AS se,
                    qnorm(1 - alpha/2) AS q)
  SELECT struct_pack(point := z.rd, se := z.se,
    lower := z.rd - z.q*z.se, upper := z.rd + z.q*z.se, nnt := 1.0/abs(z.rd))
  FROM z
);

-- Fisher 精确检验 (2x2, 超几何双侧精确 p, lgamma 求组合数避免溢出)
CREATE OR REPLACE MACRO fisher_exact_p(a, b, c, d) AS (
  WITH t AS (SELECT a::BIGINT AS a, b::BIGINT AS b, c::BIGINT AS c, d::BIGINT AS d),
       c0 AS (SELECT t.a+t.b AS r1, t.c+t.d AS r2, t.a+t.c AS k,
                     t.a+t.b+t.c+t.d AS n, t.a AS xobs FROM t),
       bnd AS (SELECT greatest(0, c0.k - c0.r2) AS lo, least(c0.k, c0.r1) AS hi, c0.* FROM c0),
       lg0 AS (SELECT lgamma(bnd.r1+1)+lgamma(bnd.r2+1)-lgamma(bnd.n+1)
                      +lgamma(bnd.k+1)+lgamma(bnd.n-bnd.k+1) AS const FROM bnd),
       obs AS (SELECT bnd.xobs,
               lg0.const - lgamma(bnd.xobs+1) - lgamma(bnd.r1-bnd.xobs+1)
               - lgamma(bnd.k-bnd.xobs+1) - lgamma(bnd.xobs+bnd.r2-bnd.k+1) AS logp_obs
               FROM bnd CROSS JOIN lg0),
       xs AS (SELECT x,
               lg0.const - lgamma(x+1) - lgamma(bnd.r1-x+1)
               - lgamma(bnd.k-x+1) - lgamma(x+bnd.r2-bnd.k+1) AS logp
               FROM bnd CROSS JOIN lg0 CROSS JOIN unnest(range(bnd.lo, bnd.hi+1)) u(x))
  SELECT sum(exp(xs.logp))
  FROM xs CROSS JOIN obs
  WHERE xs.logp <= obs.logp_obs + 1e-9
);

-- ============ H. 时间序列 ============

CREATE OR REPLACE MACRO acf(x, lag) AS (
  WITH z AS (SELECT list_avg(x) AS m, len(x) AS n)
  SELECT list_sum(list_transform(range(1, z.n - lag + 1), i -> (x[i]-z.m)*(x[i+lag]-z.m)))
         / list_sum(list_transform(x, v -> pow(v - z.m, 2)))
  FROM z
);

-- 严格 PACF 需 Durbin-Levinson 递推(要迭代), SQL 宏做不了; 提供 ACF 向量
CREATE OR REPLACE MACRO acf_vector(x, max_lag) AS TABLE
  SELECT k AS lag, acf(x, k) AS r FROM unnest(range(1, max_lag)) u(k)

CREATE OR REPLACE MACRO mae(actual, pred) AS
  list_avg(list_transform(actual, (v, i) -> abs(v - pred[i])));

CREATE OR REPLACE MACRO mse(actual, pred) AS
  list_avg(list_transform(actual, (v, i) -> pow(v - pred[i], 2)));

CREATE OR REPLACE MACRO rmse(actual, pred) AS sqrt(mse(actual, pred));

CREATE OR REPLACE MACRO mape(actual, pred) AS
  list_avg(list_transform(actual, (v, i) -> abs((v - pred[i]) / v))) * 100.0;

-- ============ I. 特殊函数 (statcpp special_functions 里 core/stats_duck 没有的) ============

-- erf(x) = 2*pnorm(x*sqrt(2)) - 1   (数值恒等, 验证到 1e-8 vs scipy)
CREATE OR REPLACE MACRO erf(x) AS 2 * pnorm(x * sqrt(2)) - 1;

CREATE OR REPLACE MACRO erfc(x) AS 1 - erf(x);

-- 正态上尾 (complement CDF; 大 x 时比 1-pnorm 稳)
CREATE OR REPLACE MACRO norm_sf(x) AS 1 - pnorm(x);

-- 注意: betainc 系列 (betainc/betaincinv/gammainc±) 已有分布族覆盖:
--   pbeta(q, a, b) / dbeta / qbeta / pgamma / qgamma 等 15+ 族在 stats_duck。

-- ============ J. 加权统计 + 集中度 (statcpp weighted/concentration 模块) ============
-- 注意: CTE 列名不能与宏参数同名 (会遮蔽 unnest) — 用 vv/ww。

CREATE OR REPLACE MACRO weighted_mean(v, wt) AS TABLE
  WITH u AS (SELECT unnest(v) AS vv, unnest(wt)::DOUBLE AS ww)
  SELECT sum(ww * vv) / sum(ww) AS weighted_mean FROM u;

CREATE OR REPLACE MACRO weighted_stddev(v, wt) AS TABLE
  WITH u AS (SELECT unnest(v) AS vv, unnest(wt)::DOUBLE AS ww),
       m AS (SELECT sum(ww * vv) / sum(ww) AS wm, sum(ww) AS W FROM u),
       d AS (SELECT u.ww * pow(u.vv - m.wm, 2) AS sq, m.W AS W FROM u CROSS JOIN m)
  SELECT sqrt(sum(sq) / (ANY_VALUE(W) - 1)) AS weighted_stddev FROM d;

CREATE OR REPLACE MACRO quantile_weighted(v, wt, p := 0.5) AS TABLE
  WITH u AS (SELECT unnest(v) AS vv, unnest(wt)::DOUBLE AS ww),
       s AS (SELECT vv AS q, sum(ww) OVER (ORDER BY vv) - ww AS cum, sum(ww) OVER () AS W FROM u)
  SELECT min(q) AS quantile FROM s WHERE cum >= p * W;

-- Gini 系数 (0=完全均等, 1=极端集中)
CREATE OR REPLACE MACRO gini(v) AS TABLE
  WITH u AS (SELECT unnest(v) AS vv),
       s AS (SELECT vv, row_number() OVER (ORDER BY vv) AS i, count(*) OVER () AS n FROM u)
  SELECT (2.0 * sum(i * vv)) / (ANY_VALUE(n) * sum(vv)) - (ANY_VALUE(n) + 1.0) / ANY_VALUE(n) AS gini FROM s;

-- HHI 赫芬达尔指数 (基于占比; 1/k=完全竞争, 1=垄断)
CREATE OR REPLACE MACRO hhi(v) AS TABLE
  WITH u AS (SELECT unnest(v) AS vv),
       t AS (SELECT vv, sum(vv) OVER () AS W FROM u)
  SELECT sum(pow(vv / W, 2)) AS hhi FROM t;
