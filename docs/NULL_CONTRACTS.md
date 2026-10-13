# NULL、NaN 与 Inf 的行为契约

`stats_fill.sql` 的 64 个宏是 **SQL 宏层**：每个宏展开为 DuckDB 的
`list_*` / `var_samp` / `stats_duck` 函数调用，**自身没有任何输入校验**。
所以这里的行为不是"统一缺失值策略"，而是**逐宏继承底层 DuckDB 函数的行为**。
与 statcpp（C++ 层，全函数对齐 R 的默认缺失值处理）不同，SQL 宏层会出现
**静默错值**而不是报错——这是最重要的区别，引用本文档前请先读第 6 节。

本文档的每一条都在真 DuckDB 上跑过（`v1.5.6` + `stats_duck`，2026-10-13；
复现脚本：仓库内 `tests/null_nan_inf_probe.sql`）。

---

## 1. DuckDB 能存什么

- **NaN 是值，不是 `NULL`。** `'nan'::DOUBLE IS NULL` 为 `false`，`DOUBLE`
  列可以和 `NULL` 并存 NaN。
- **Inf 可以存。** `'inf'::DOUBLE`、`'-inf'::DOUBLE`、`9e999` 都保留为
  `+Inf` / `-Inf`，按普通值参与计算。

要把 NaN 当缺失值，先转成 `NULL`：

```sql
list_filter(x, lambda v: NOT isnan(v))               -- 过滤 NaN（注意 NULL 也会被过滤）
list_filter(x, lambda v: NOT (isnan(v) OR v IS NULL))
```

---

## 2. 入参 `LIST` 里的 NULL：宏不跳过，行为取决于底层函数

**这是与 statcpp 最关键的区别。** statcpp 的聚合在进入统计前会移除 `NULL`
行；而这里的宏拿到的是 `array_agg(x)` 产出的 LIST——**`array_agg` 保留 NULL
元素**（`len(array_agg(v))` 含 NULL 行也计入长度）。宏把整个 LIST 原样传给
`list_*`，于是：

| 底层函数 | 对 LIST 内 NULL 的行为 | 实测 |
| --- | --- | --- |
| `list_avg` | **跳过**（按非 NULL 元素求均值） | `list_avg([1.0,NULL,3.0]) = 2.0` |
| `len` | **计入** | `len([1.0,NULL,3.0]) = 3` |
| `list_stddev_samp` | **Out of Range Error**（不跳过） | 2 元素 1 NULL 时报错 |
| `median`（`stat_median`） | 跳过 | 与 `list_avg` 同族 |

同一个宏里混用这两类底层函数时，`n`（来自 `len`，含 NULL）与
`avg`/`stddev`（按非 NULL 算）的口径**不一致**。后果：

```sql
-- v = 10, NULL, 30, 50（n 按 4 算，avg 按 3 算，se 直接报错）
SELECT ci_mean(array_agg(v) FROM t);        -- Out of Range Error
-- z_test / f_test_var / kruskal_wallis 的 LIST 里含 NULL 元素时同样触发
-- 底层 stddev/var_samp 报错，不是静默返回错值。
```

**规则：先过滤再聚合。** 所有宏的推荐调用形态：

```sql
SELECT ci_mean(list_filter(array_agg(x), lambda v: v IS NOT NULL)) FROM t;
```

---

## 3. 参数（`alpha`、`p`、`sigma` 等）：`NULL` 参数 → 传播或报错

- **`NULL` 参数按 DuckDB 标量规则传播**：`stat_pctl(l, NULL)` → `NULL`；
  `ci_mean(l, NULL)` 的 `point` 正常但 `se/lower/upper` 为 `NULL`
  （`qnorm(1-NULL/2)` 传播 NULL）。
- **不要传 `NULL` 参数去"使用默认值"。** 默认值只在**省略参数**时生效
  （`alpha := 0.05`）：`ci_mean(array_agg(x))` 用 0.05，
  `ci_mean(array_agg(x), NULL)` 不是。
- **越界数值参数不校验**：`stat_pctl(l, 1.5)` → `NULL`（插值越界），
  `stat_pctl(l, -0.5)` → `NULL`；`sample_size_moe_prop(-0.1, 0.05)` → `97.0`
  （负 `moe` 按绝对值口径算出一个看似合理的数）。

---

## 4. NaN 数据进宏：**多数静默错值，不报错**

statcpp 对 NaN 的策略是按函数类别分流（描述统计给 NULL / 检验丢弃 NaN 行 /
分位数报 "data contains NaN"）。SQL 宏层没有这个分流：

| 宏 | NaN 输入的行为 | 实测 |
| --- | --- | --- |
| `stat_pctl` / `stat_median` / `stat_iqr` | **静默错值**：`list_sort` 把 NaN 排到末尾，插值/排序位置被污染 | `stat_iqr([1,2,3,4,5,nan]) = 2.5`（真值 2.0）；`stat_pctl([1,2,nan,4],0.5) = 3.0` |
| `stat_trimmed_mean` | **静默错值**：NaN 参与排序后被截尾逻辑以正常元素处理 | `stat_trimmed_mean([1,2,nan,4,5,6],0.2) = 4.25` |
| `list_avg` / `list_stddev_samp`（多数 CI、z 检验的底层） | **传播 NaN** | `list_avg([1,2,nan,4]) = nan`（1.5.6 的 `list_stddev_samp` 含 NaN 反而 Out of Range Error；2.0.0-alpha 传播 NaN——版本间行为不同） |
| `stat_cv` | **Out of Range Error** | `stat_cv([1,2,inf,4])` 报错 |
| `stat_geometric_mean` | **硬报错**（`ln` 域错误） | `stat_geometric_mean([1,-2,4])` → "cannot take logarithm of a negative number" |
| `stat_harmonic_mean` | **静默 0.0**：`1/0 = inf`，`sum` 含 inf…实测 `stat_harmonic_mean([1,0,4]) = 0.0` | 除零被吞 |
| `stat_pctl` 含 `inf` | 静默错值（排序把 inf 放末尾） | `stat_pctl([1,2,inf,4],0.5) = 3.0` |

**规则：NaN/Inf 不进宏。** 调用前过滤：

```sql
list_filter(array_agg(x), lambda v: NOT isnan(v) AND NOT isinf(v))
```

---

## 5. 退化数据：错值或硬错误，**没有统一的"报 NULL"约定**

| 输入 | 行为 | 实测 |
| --- | --- | --- |
| 零方差 → `ci_mean` | `se = 0`，`lower = upper = point`（数学上成立但不报错） | `ci_mean([5,5,5]) = point 5.0 / se 0.0 / lower 5.0` |
| 零方差 → `stat_cv` | `0.0`（`0/5`，静默） | `stat_cv([5,5,5]) = 0.0` |
| 空 LIST → 单值宏 | `NULL` | `stat_pctl([],0.5) = NULL`，`stat_trimmed_mean([],0.1) = NULL` |
| 空 LIST → `ci_mean` | `point = NULL`（`list_avg([])` 为 NULL） | 不报错 |
| 空 LIST → TABLE 宏（`stat_mad`） | `NULL` 行 | `stat_mad([])` → `mad = NULL` |
| n=1 → `list_stddev_samp` | **Out of Range Error**（不是 NULL） | 单元素 LIST 的 CI 全部报错 |
| `0` 元素 → 检验 | 底层 `var_samp`/`qt` 报错或 NaN | 需调用方保证 `n>=2` |

与 statcpp 的约定（"退化数据是错误，不是 NULL"，如零方差 t 检验报
`statcpp::t_test: zero variance`）不同：这里**只有部分退化路径报错**，
其余给出数学上成立但语义上误导的值（`se=0` 的 CI、`cv=0`）。

---

## 6. 与 statcpp / DuckDB-stats 的差异（为什么不能照抄它们的契约）

| 维度 | statcpp / DuckDB-stats（C++ 扩展） | 本仓 64 宏（SQL 宏层） |
| --- | --- | --- |
| 入参缺失值 | 聚合在进入统计**前移除 NULL 行**；NaN 按 R 语义分流 | **不处理**：`array_agg` 保留 NULL 进 LIST，行为逐函数继承底层 `list_*` |
| NULL 参数 | 结果 `NULL`（对齐 `lag(v, NULL)`） | 传播或报错（同左，但无文档化约定） |
| NaN 数据 | 按函数类别：NULL / 丢弃 / 报错，三者皆有定义 | **多为静默错值**（排序类），或底层报错（stddev/ln） |
| 退化数据 | 统一报错（"zero variance"） | 零方差 `ci_mean` 给 `se=0` 不报错；n=1 报错 |
| 出错语义 | `Invalid Input Error: statcpp::...`，连接可用 | DuckDB `Out of Range` / `Binder` 错误，或**无错误的错值** |

**一句话：C++ 扩展层把缺失值当契约管，SQL 宏层把缺失值当调用方的责任。**
本文档第 2/3/4 节的"过滤再聚合"形态是唯一安全姿势；不要把 statcpp 的
"NULL 行自动跳过"心智模型带到这里。

---

## 7. 推荐调用模板

```sql
-- 单列
SELECT ci_mean(list_filter(array_agg(x ORDER BY id),
           lambda v: NOT isnan(v) AND NOT isinf(v) AND v IS NOT NULL)) FROM t;

-- 两列（必须同序对齐，过滤后长度仍须一致——两侧用同一 WHERE/过滤器）
SELECT ci_mean_diff_welch(
         list_filter(array_agg(a ORDER BY id),
           lambda v: v IS NOT NULL AND NOT isnan(v)),
         list_filter(array_agg(b ORDER BY id),
           lambda v: v IS NOT NULL AND NOT isnan(v)))
FROM t;
```
