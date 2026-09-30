"""DuckDB 统计默认环境启动器.

用法 (任意 Python 脚本, 把仓库根加入 sys.path):
    import sys; sys.path.insert(0, '<path-to-repo>')
    from duckdb_stats_setup import stats_connect
    con = stats_connect()

或直接跑本文件做健康检查:  python duckdb_stats_setup.py

功能:
  1. 确保 stats_duck 社区扩展已装入【默认】DuckDB home (~/.duckdb, 不设 DUCKDB_HOME)
  2. LOAD stats_duck
  3. 执行本仓库的 stats_fill.sql (64 个统计宏)
"""
import os
import duckdb

# stats_fill.sql 与本文件同目录 (仓库根)
FILL_SQL = os.path.join(os.path.dirname(os.path.abspath(__file__)), "stats_fill.sql")


def stats_connect(config=None):
    # 不设置 DUCKDB_HOME -> 用真实默认 home ~/.duckdb
    cfg = {'allow_community_extensions': 'true'}
    if config:
        cfg.update(config)
    con = duckdb.connect(config=cfg)
    try:
        installed = con.execute(
            "SELECT installed FROM duckdb_extensions() WHERE extension_name='stats_duck'"
        ).fetchone()
    except Exception:
        installed = None
    if not installed or not installed[0]:
        con.execute("INSTALL stats_duck FROM community")
    con.execute("LOAD stats_duck")
    with open(FILL_SQL, encoding="utf-8") as f:
        sql = f.read()
    import re
    body = re.sub(r"--[^\n]*", "", sql)
    for stmt in re.split(r"(?=CREATE OR REPLACE MACRO)", body):
        if stmt.strip():
            con.execute(stmt)
    return con


if __name__ == "__main__":
    con = stats_connect()
    import numpy as np
    np.random.seed(0)
    x = np.random.normal(100, 12, 60)
    con.execute("CREATE TABLE t(id BIGINT, x DOUBLE, grp VARCHAR)")
    con.executemany("INSERT INTO t VALUES (?,?,?)",
                    [(i, float(v), str(i % 3)) for i, v in enumerate(x)])
    ci = con.execute("SELECT ci_mean(array_agg(x ORDER BY id)) FROM t").fetchone()[0]
    kw = con.execute("""SELECT h_statistic, p_value FROM kruskal_wallis(
        (SELECT array_agg(x ORDER BY id) FROM t),
        (SELECT array_agg(grp ORDER BY id) FROM t))""").fetchone()
    sd = con.execute("SELECT pnorm(1.96)").fetchone()[0]
    print("duckdb", duckdb.__version__)
    print("stats_duck pnorm(1.96) =", round(sd, 4), "(expect 0.9750)")
    print("ci_mean:", {k: round(v, 3) for k, v in ci.items()})
    print("kruskal_wallis: H=%.4f p=%.3e" % kw)
    print("OK - stats_fill 宏已加载:",
          con.execute("SELECT count(*) FROM duckdb_functions() WHERE function_name IN "
                      "('ci_mean','cohens_d_2samp','fisher_exact_p','power_2samp')").fetchone()[0],
          "/ 4 scalar-macro spot-check names (TABLE 宏如 tukey_hsd 不在此表)")
