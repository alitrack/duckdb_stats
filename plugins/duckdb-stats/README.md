# duckdb-stats — Claude Code Plugin

Two skills that give your AI assistant verified knowledge of
[duckdb_stats](https://github.com/alitrack/duckdb_stats): the stats_connect
bootstrap, the five hard-won calling conventions, and all 64 stats_fill.sql
macros organized by task. R/SciPy-level statistics, scipy-verified (58
assertions, 0 fails).

## Install

```text
/plugin marketplace add alitrack/duckdb_stats
/plugin install duckdb-stats@duckdb-stats
```

## Skills

| Skill | Covers |
|---|---|
| `duckdb-stats-quickstart` | three layers, bootstrap, 5 calling conventions |
| `duckdb-stats-catalog` | all 64 macros by category + picking examples |

## Prerequisite

```sql
INSTALL stats_duck FROM community;   -- extension layer
LOAD stats_duck;
-- then execute stats_fill.sql from the repo (64 macros)
```

MIT, alitrack.
