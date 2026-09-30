# duckdb-stats — Statistical analysis layer for DuckDB (extensions + SQL macros)

**[中文说明 → README_CN.md](README_CN.md)**

Get **R/SciPy-level statistics inside DuckDB**: two plug-ins + **64 SQL macros** that
together cover every functional module of
[statcpp](https://github.com/mitsuruk/statcpp) (386 C++ statistics functions).

```sql
-- one connection:
INSTALL stats_duck FROM community;  LOAD stats_duck;
-- + run stats_fill.sql
SELECT ci_mean(array_agg(x)) FROM t;          -- confidence interval
SELECT * FROM oneway_anova((SELECT array_agg(x ORDER BY id) FROM t),
                           (SELECT array_agg(grp ORDER BY id) FROM t));   -- one-way ANOVA
SELECT cohens_d_2samp(array_agg(a ORDER BY id), array_agg(b ORDER BY id)) FROM t;  -- effect size
```

## Three capability layers (why this "completes statcpp")

| Layer | What it provides | Source |
|---|---|---|
| **stats_duck** (community extension, KoliStat) | d/p/q/r for 15+ distribution families, t-test family, ANOVA, chi-square, non-parametric tests (Mann-Whitney/Wilcoxon/KS/Shapiro), Pearson/Spearman/Kendall, `lm` regression (HC/cluster-robust SEs), bootstrap, `adjust_p` (Bonferroni/Holm/Hochberg/FDR/BY), correlation matrix, SAS/SPSS/Stata readers | auto `INSTALL`ed by this repo |
| **stats_fill.sql** (this repo, 64 macros) | statcpp's genuinely unique inferential layer: CIs, effect sizes, power/sample-size, exact Fisher, Kruskal-Wallis, one/two-way ANOVA, Tukey/Scheffé/Dunnett post-hoc, 2×2 categorical (OR/RR/NNT), weighted stats, Gini/HHI, ACF, erf/erfc | this repo |
| **duckdb-ml** (optional, a 34-algorithm Rust extension) | the "heavy" statcpp modules: GLM/logistic training, clustering (kmeans/DBSCAN/hierarchical/FCM/t-SNE), survival (KM/Cox), ARIMA, ridge/lasso/elastic, PCA/LDA, SVM/XGBoost/RF/MLP/KNN/NB | `LOAD` on demand |

Coverage matrix per statcpp module: [`docs/COVERAGE.md`](docs/COVERAGE.md) (raw CSV: [`statcpp_vs_duckdb.csv`](docs/statcpp_vs_duckdb.csv))
(30 modules, each verified on a real DuckDB build). **Verdict**: statcpp's 386 functions =
distribution families (→ stats_duck) + tests (→ stats_duck) + ML/clustering/survival (→
duckdb-ml) + inferential layer (→ stats_fill.sql). Every function has a home. The only
`[APPROX]` items are Tukey/Dunnett p-values (DuckDB has no studentized-range CDF; z-approx,
anti-conservative for small samples).

## Quick start

### Python (recommended)

```python
# after cloning this repo:
import sys; sys.path.insert(0, '<path-to-repo>')
from duckdb_stats_setup import stats_connect

con = stats_connect()   # auto: INSTALL stats_duck -> LOAD -> load 64 macros
con.execute("CREATE TABLE t AS SELECT i::BIGINT id, 100 + 10*(i%3) + random()*8 x, (i%3)::VARCHAR grp FROM range(90) r(i)")
print(con.execute("SELECT ci_mean(array_agg(x ORDER BY id)) FROM t").fetchone())

# optionally add the ML layer (path to a built duckdb-ml extension):
con = stats_connect(ml_extension=r'/path/to/ml.duckdb_extension')
```

### SQL CLI / other languages

```sql
INSTALL stats_duck FROM community;
LOAD stats_duck;
-- execute stats_fill.sql (CREATE OR REPLACE MACRO ...)
SELECT power_2samp(50, 0.5);              -- power
SELECT fisher_exact_p(20,10,5,15);        -- exact Fisher, two-sided
SELECT * FROM tukey_hsd((SELECT array_agg(x ORDER BY id) FROM t),
                        (SELECT array_agg(grp ORDER BY id) FROM t));
```

Full signatures and calling conventions: **[`docs/FUNCTION_REFERENCE.md`](docs/FUNCTION_REFERENCE.md)**
(all 64 macros listed).

## Key calling conventions (hard-won)

1. **DuckDB lists are 1-indexed**: `SELECT l[i] ... FROM unnest(range(1, len(l)+1)) t(i)` to get all n elements.
2. **Multiple `array_agg` calls are independent aggregates — row order is NOT aligned.**
   Grouped macros must share one `ORDER BY` key: `array_agg(x ORDER BY id)` /
   `array_agg(grp ORDER BY id)`.
3. Table-macro arguments must be **scalar subqueries** — no bare aggregate, no
   `SELECT * FROM m(...) FROM t` (double FROM is illegal).
4. `stats_duck`'s `pf/pt/pchisq` are **CDFs (lower tail)**; upper-tail p = `1 - pf(...)`.
5. `tukey_hsd` / `dunnett_approx` p-values are marked `[APPROX]` (z-approximation).

## Verification

```bash
cd tests
uv venv .venv && uv pip install -p .venv duckdb scipy numpy
.venv/bin/python gen_refs.py      # generate reference values (fixed seeds)
.venv/bin/python verify_fill.py   # loads all 64 macros + 58 assertions vs scipy 1.18 -> TOTAL FAILS: 0
```

(Windows: `.venv\Scripts\python ...`. `verify_fill.py` uses a temporary DUCKDB_HOME and does
not touch your real `~/.duckdb`.)

## Files

```
stats_fill.sql             64 statistics macros (A descriptive, B CIs, C power, D tests,
                           E ANOVA/post-hoc, F effect sizes, G categorical, H time series,
                           I special functions, J weighted + concentration)
duckdb_stats_setup.py      stats_connect() one-line bootstrap
docs/FUNCTION_REFERENCE.md 64 macro signatures (with return fields + examples)
docs/COVERAGE.md           30-module coverage matrix (rendered table)
docs/statcpp_vs_duckdb.csv raw coverage data (machine-readable)
tests/                     scipy verification (gen_refs + verify_fill + refs.json)
```

## Why "plug-ins + macros" instead of a custom Rust extension

- **Lua route** (rejected, measured): stateful aggregates cannot be built
  (`must return a function`); scalars ~7× slower than builtins.
- **Custom Rust wrapper around statcpp**: only worth it when the product demands a single
  closed dependency / exact 386-function parity — then wrap `special_functions.hpp` and the
  GLM inference layer first.
- This repo: one community extension + one SQL file — zero compilation, cross-platform,
  works on DuckDB >= 1.4.

## License

MIT
