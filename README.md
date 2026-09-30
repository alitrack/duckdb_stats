# duckdb-stats — DuckDB 统计分析层（插件 + SQL 宏）

在 DuckDB 里获得 **R/SciPy 级别的统计分析能力**：两个社区扩展 + 一套 **64 个 SQL 宏**，
覆盖 [statcpp](https://github.com/mitsuruk/statcpp)（386 个 C++ 统计函数）的全部功能模块。

```sql
-- 一个连接里:
INSTALL stats_duck FROM community;  LOAD stats_duck;
-- + 执行 stats_fill.sql
SELECT ci_mean(array_agg(x)) FROM t;                 -- 置信区间
SELECT * FROM oneway_anova((SELECT array_agg(x ORDER BY id) FROM t),
                           (SELECT array_agg(grp ORDER BY id) FROM t));  -- 单因素 ANOVA
SELECT cohens_d_2samp(array_agg(a ORDER BY id), array_agg(b ORDER BY id)) FROM t; -- 效应量
```

## 三层能力（为什么"完成 statcpp 的所有功能"）

| 层 | 提供什么 | 来源 |
|---|---|---|
| **stats_duck**（社区插件, KoliStat） | 15+ 分布族的 d/p/q/r、t 检验族、ANOVA、卡方、非参数检验（Mann-Whitney/Wilcoxon/KS/Shapiro）、Pearson/Spearman/Kendall、`lm` 回归（HC/cluster 稳健 SE）、bootstrap、`adjust_p`（Bonferroni/Holm/Hochberg/FDR/BY）、相关性矩阵、SAS/SPSS/Stata 读写 | 本仓库自动 `INSTALL` |
| **stats_fill.sql**（本仓库, 64 宏） | statcpp 真正"独有"的推断层：置信区间、效应量、功效/样本量、Fisher 精确、Kruskal-Wallis、单/双因素 ANOVA、Tukey/Scheffé/Dunnett 事后、2×2 分类（OR/RR/NNT）、加权统计、Gini/HHI、ACF、erf/erfc | 本仓库 |
| **duckdb-ml**（可选, 你自有的 34 算法 Rust 扩展） | statcpp 里"重"的部分：GLM/Logistic 训练、聚类（kmeans/DBSCAN/层次/FCM/t-SNE）、生存分析（KM/Cox）、ARIMA、ridge/lasso/elastic、PCA/LDA、SVM/XGBoost/RF/MLP/KNN/NB | 按需 `LOAD` |

覆盖对照矩阵见 [`docs/statcpp_vs_duckdb.csv`](docs/statcpp_vs_duckdb.csv)（30 个模块逐项判定）。
**结论**：statcpp 386 函数 = 分布族(→stats_duck) + 检验(→stats_duck) + ML/聚类/生存(→duckdb-ml) +
推断层(→stats_fill.sql)。全部功能有落点；唯一 [APPROX] 是 Tukey/Dunnett 的 p 值
（DuckDB 无 studentized-range CDF，用 z 近似，小样本偏松）。

## 快速开始

### Python（推荐）

```python
# clone 本仓库后:
import duckdb
from duckdb_stats_setup import stats_connect

con = stats_connect()   # 自动: 装 stats_duck → LOAD → 加载 65 宏
con.execute("CREATE TABLE t AS SELECT i::BIGINT id, 100 + 10*(i%3) + random()*8 x, (i%3)::VARCHAR grp FROM range(90) r(i)")
print(con.execute("SELECT ci_mean(array_agg(x ORDER BY id)) FROM t").fetchone())
```

### SQL CLI / 其他语言

```sql
INSTALL stats_duck FROM community;
LOAD stats_duck;
-- 逐条执行 stats_fill.sql (CREATE OR REPLACE MACRO ...)
-- 之后:
SELECT power_2samp(50, 0.5);              -- 功效
SELECT fisher_exact_p(20,10,5,15);        -- Fisher 精确双侧
SELECT * FROM tukey_hsd((SELECT array_agg(x ORDER BY id) FROM t),
                        (SELECT array_agg(grp ORDER BY id) FROM t));
```

完整函数签名与调用约定：**[`docs/FUNCTION_REFERENCE.md`](docs/FUNCTION_REFERENCE.md)**（65 个宏逐一列出）。

## 关键调用约定（踩坑总结）

1. **list 是 1-indexed**：`SELECT l[i] ... FROM unnest(range(1, len(l)+1)) t(i)` 才取全。
2. **多个 `array_agg` 独立聚合、顺序不保证** → 分组函数必须同一 `ORDER BY` 键：
   `array_agg(x ORDER BY id)` / `array_agg(grp ORDER BY id)`。
3. 分组（TABLE）宏的参数必须是**标量子查询**，不能是裸聚合，也不能 `SELECT * FROM m(...) FROM t`（双 FROM 非法）。
4. `stats_duck` 的 `pf/pt/pchisq` 是 **CDF（下尾）**；上尾 p = `1 - pf(...)`。
5. `tukey_hsd`/`dunnett_approx` 的 p 值标 `[APPROX]`（z 近似）。

## 验证

```bash
cd tests
uv venv .venv && uv pip install -p .venv duckdb scipy numpy
.venv/bin/python gen_refs.py      # 生成参照值(固定种子)
.venv/bin/python verify_fill.py   # 64 宏加载 + 58 项断言 vs scipy 1.18 → TOTAL FAILS: 0
```

（Windows: `.venv\Scripts\python ...`；`verify_fill.py` 用临时 DUCKDB_HOME，不污染真实 `~/.duckdb`。）

## 文件

```
stats_fill.sql           64 个统计宏 (A 描述 B 置信区间 C 功效 D 检验 E ANOVA F 效应量 G 分类 H 时序 I 特殊函数 J 加权/集中度)
duckdb_stats_setup.py    stats_connect() 一键启动器
docs/FUNCTION_REFERENCE.md   64 宏签名参考
docs/statcpp_vs_duckdb.csv   30 模块覆盖矩阵
tests/                     scipy 验证 (gen_refs + verify_fill + refs.json)
```

## 为什么是"插件 + 宏"而不是自研 Rust 扩展

- **Lua 路线**（实测否决）：有状态聚合建不出（`must return a function`）、标量慢 ~7×。
- **自研 Rust 包装 statcpp**：只有当产品要求"单一闭源依赖 / 386 函数逐项对齐"才值得；
  且优先包 `special_functions.hpp`（erf 系列已有宏覆盖）与 GLM 推断层。
- 本仓库方案：装一个社区插件 + 一份 SQL 文件，零编译、跨平台、随 `DuckDB >= 1.4`。

## License

MIT
