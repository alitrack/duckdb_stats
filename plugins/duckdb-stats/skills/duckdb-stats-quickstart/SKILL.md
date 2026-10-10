---
name: duckdb-stats-quickstart
description: >
  Get R/SciPy-level inferential statistics running inside DuckDB with
  duckdb-stats: bootstrap stats_connect() (or the SQL path), the three
  capability layers (stats_duck extension + 64 macros + optional duckdb-ml),
  and the five hard-won calling conventions that make macros work the first
  time. Use before writing any statistical SQL against this stack.
version: 0.1.0
---

# duckdb-stats: quickstart & calling conventions

Layers: **stats_duck** (community extension: d/p/q/r for 15+ distributions,
t-tests, chi-square, non-parametric, correlations, lm, bootstrap, p-adjust) +
**stats_fill.sql** (this repo, 64 macros: the inferential layer) +
**duckdb-ml** optional (GLM/clustering/survival/ARIMA).

## Bootstrap

Python (recommended):

```python
import sys; sys.path.insert(0, '<path-to-repo>')
from duckdb_stats_setup import stats_connect
con = stats_connect()   # auto: INSTALL stats_duck -> LOAD -> load 64 macros
# add ML layer: stats_connect(ml_extension='/path/to/ml.duckdb_extension')
```

SQL CLI:

```sql
INSTALL stats_duck FROM community;
LOAD stats_duck;
-- then execute stats_fill.sql (CREATE OR REPLACE MACRO x64)
SELECT power_2samp(50, 0.5);
SELECT fisher_exact_p(20,10,5,15);
```

## The five calling conventions (break these = wrong results)

1. **DuckDB lists are 1-indexed**: `l[i]` runs 1..len(l).
2. **Multiple `array_agg` calls are independent aggregates — row order is
   NOT aligned.** Grouped macros must share one ORDER BY key:
   `array_agg(x ORDER BY id)` and `array_agg(grp ORDER BY id)`.
3. **Table-macro arguments must be scalar subqueries**:
   `SELECT * FROM oneway_anova((SELECT array_agg(x ORDER BY id) FROM t),
   (SELECT array_agg(grp ORDER BY id) FROM t))` — bare aggregates and
   double-FROM are illegal.
4. **`pf/pt/pchisq` are lower-tail CDFs**; upper-tail p-value = `1 - pf(...)`.
5. **`tukey_hsd` / `dunnett_approx` p-values are `[APPROX]`**
   (z-approximation; anti-conservative for small samples).

## Worked examples

```sql
SELECT ci_mean(array_agg(x)) FROM t;                      -- CI for the mean
SELECT * FROM oneway_anova((SELECT array_agg(x ORDER BY id) FROM t),
                           (SELECT array_agg(grp ORDER BY id) FROM t));
SELECT cohens_d_2samp(array_agg(a ORDER BY id),
                      array_agg(b ORDER BY id)) FROM t;   -- effect size
```

Verification: 58 assertions vs scipy 1.18, 0 fails (`tests/verify_fill.py`).
