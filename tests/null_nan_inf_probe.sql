-- NULL / NaN / Inf 行为契约的复现脚本（docs/NULL_CONTRACTS.md 的每条断言）
-- 用法: duckdb -unsigned < tests/null_nan_inf_probe.sql   (需要 stats_duck)
-- 实测环境: DuckDB v1.5.6 + stats_duck, 2026-10-13
INSTALL stats_duck FROM community;
LOAD stats_duck;
.read ./stats_fill.sql

-- ===== 1. DuckDB 能存什么 =====
SELECT 'nan_is_null=' || ('nan'::DOUBLE IS NULL);                    -- nan_is_null=false
SELECT 'inf_stored=' || ('inf'::DOUBLE = '9e999'::DOUBLE);           -- inf_stored=true
SELECT 'neg_inf_stored=' || ('-inf'::DOUBLE IS NOT NULL);            -- neg_inf_stored=true

-- ===== 2. LIST 里的 NULL: 宏不跳过, 逐函数继承底层 =====
SELECT 'list_avg_skip=' || list_avg([1.0, NULL, 3.0]);               -- list_avg_skip=2.0 (跳过NULL)
SELECT 'len_counts_null=' || len([1.0, NULL, 3.0]);                  -- len_counts_null=3 (计入NULL)
-- list_stddev_samp 含 NULL -> Out of Range Error (1.5.6); 不在此断言, 见文档
-- array_agg 保留 NULL 元素:
CREATE TEMP TABLE t_null AS SELECT unnest([10.0, NULL, 30.0, 50.0]) AS v;
SELECT 'agg_keeps_null=' || len(array_agg(v)) FROM t_null;           -- agg_keeps_null=4
-- ci_mean 含 NULL -> 底层 stddev 报错 (不在此断言; 文档第 2 节)

-- ===== 3. 参数: NULL 传播 / 越界不校验 =====
SELECT 'pctl_null_p=' || COALESCE(stat_pctl([1.0, 2.0, 3.0], NULL)::VARCHAR, 'NULL');      -- pctl_null_p=NULL
SELECT 'ci_null_se=' || COALESCE(ci_mean([1.0, 2.0, 3.0, 4.0], NULL).lower::VARCHAR, 'NULL'); -- ci_null_se=NULL
SELECT 'pctl_1.5=' || COALESCE(stat_pctl([1.0, 2.0, 3.0], 1.5)::VARCHAR, 'NULL');          -- pctl_1.5=NULL (越界)
SELECT 'ss_moe_neg=' || sample_size_moe_prop(-0.1, 0.05);            -- ss_moe_neg=97.0 (不校验)

-- ===== 4. NaN 进宏: 多数静默错值 =====
SELECT 'iqr_nan=' || stat_iqr([1.0, 2.0, 3.0, 4.0, 5.0, 'nan'::DOUBLE]);   -- iqr_nan=2.5 (真值 2.0, 静默错)
SELECT 'pctl_nan=' || stat_pctl([1.0, 2.0, 'nan'::DOUBLE, 4.0], 0.5);      -- pctl_nan=3.0
SELECT 'trim_nan=' || stat_trimmed_mean([1.0, 2.0, 'nan'::DOUBLE, 4.0, 5.0, 6.0], 0.2); -- trim_nan=4.25
SELECT 'avg_nan=' || list_avg([1.0, 2.0, 'nan'::DOUBLE, 4.0]);             -- avg_nan=nan (传播)
SELECT 'harm_zero=' || stat_harmonic_mean([1.0, 0.0, 4.0]);                -- harm_zero=0.0 (除零被吞)
SELECT 'pctl_inf=' || stat_pctl([1.0, 2.0, 'inf'::DOUBLE, 4.0], 0.5);      -- pctl_inf=3.0
-- stat_cv([1.0,2.0,'inf'::DOUBLE,4.0]) -> Out of Range Error (文档第 4 节)
-- stat_geometric_mean([1.0,-2.0,4.0])  -> ln 域错误 (文档第 4 节)

-- ===== 5. 退化数据 =====
SELECT 'ci_zero=' || ci_mean([5.0, 5.0, 5.0]).se;                     -- ci_zero=0.0 (不报错)
SELECT 'cv_zero=' || stat_cv([5.0, 5.0, 5.0]);                        -- cv_zero=0.0 (静默)
SELECT 'pctl_empty=' || COALESCE(stat_pctl([]::DOUBLE[], 0.5)::VARCHAR, 'NULL');   -- pctl_empty=NULL
SELECT 'ci_empty_point=' || COALESCE(ci_mean([]::DOUBLE[]).point::VARCHAR, 'NULL'); -- ci_empty_point=NULL
SELECT 'mad_empty=' || COALESCE((SELECT mad FROM stat_mad([]::DOUBLE[]))::VARCHAR, 'NULL'); -- mad_empty=NULL
-- n=1 -> list_stddev_samp Out of Range Error (文档第 5 节)

-- ===== 6. 推荐调用模板 (安全姿势) =====
CREATE TEMP TABLE t_mix AS SELECT unnest([1.0, 'nan'::DOUBLE, 3.0, NULL, 5.0]) AS x;
SELECT 'safe_ci=' || round(ci_mean(list_filter(array_agg(x),
    lambda v: v IS NOT NULL AND NOT isnan(v) AND NOT isinf(v))).point, 8) FROM t_mix;
-- safe_ci=3.0  (NaN/NULL 被滤掉, 对 [1,3,5] 取 CI)
