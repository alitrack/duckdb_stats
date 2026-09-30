"""验证 stats_fill.sql: 58+ 项断言对 scipy/numpy 逐项核对.
用法 (在 tests/ 目录):
    python verify_fill.py
依赖: duckdb, scipy, numpy
"""
import os, sys, json, math, re, duckdb
import numpy as np
from scipy import stats as S

HERE = os.path.dirname(os.path.abspath(__file__))
REPO = os.path.dirname(HERE)
# 用临时 DUCKDB_HOME, 不污染真实 ~/.duckdb; stats_duck 装到这里
os.environ['DUCKDB_HOME'] = os.path.join(HERE, '.duckdb_test_home')
con = duckdb.connect(config={'allow_community_extensions':'true'})
con.execute('INSTALL stats_duck FROM community')
con.execute('LOAD stats_duck')

sql = open(os.path.join(REPO, 'stats_fill.sql'), encoding='utf-8').read()
# split at each 'CREATE OR REPLACE MACRO' boundary
parts = [s for s in re.split(r'(?=^CREATE OR REPLACE MACRO)', sql, flags=re.M) if s.strip().startswith('CREATE OR REPLACE MACRO')]
ok = fail = 0
for s in parts:
    s = s.strip()
    try:
        con.execute(s)
        ok += 1
    except Exception as e:
        fail += 1
        print(f"FAIL: {s.splitlines()[0][:60]}\n    :: {str(e).splitlines()[0][:110]}")
n_macros = len(re.findall(r'^CREATE OR REPLACE MACRO', sql, flags=re.M))
assert ok == n_macros, f"loaded {ok} != {n_macros} declared macros"
print(f"=== macro creation: {ok}/{n_macros} ok, {fail} fail ===\n")

data = json.load(open(os.path.join(HERE, 'refs.json')))
x1 = np.array(data['data']['x1']); x2 = np.array(data['data']['x2'])
pd_ = np.array(data['data']['paired_d']); g = [np.array(c) for c in data['data']['groups']]
refs = data['refs']

con.execute("CREATE TABLE t1 AS SELECT unnest(?::DOUBLE[]) AS x", [x1.tolist()])
con.execute("CREATE TABLE t2 AS SELECT unnest(?::DOUBLE[]) AS x", [x2.tolist()])
con.execute("CREATE TABLE pd AS SELECT unnest(?::DOUBLE[]) AS d", [pd_.tolist()])
# av: row id + value + group (id-aligned for the grouped macros)
vals = [float(v) for k in range(3) for v in g[k]]
aa   = [str(k) for k in range(3) for _ in range(len(g[k]))]
bb   = [str(i % 2) for i in range(len(vals))]
con.execute("CREATE TABLE av(id BIGINT, x DOUBLE, grp VARCHAR)")
con.executemany("INSERT INTO av VALUES (?,?,?)", [(i, v, a) for i,(v,a) in enumerate(zip(vals, aa))])
con.execute("CREATE TABLE av2(id BIGINT, x DOUBLE, a VARCHAR, b VARCHAR)")
con.executemany("INSERT INTO av2 VALUES (?,?,?,?)", [(i, v, a, b) for i,(v,a,b) in enumerate(zip(vals, aa, bb))])

fails = 0
def check(label, got, ref, tol=2e-3):
    global fails
    try: ok = abs(float(got) - float(ref)) <= max(tol, abs(ref)*2e-3)
    except Exception: ok = False
    if not ok: fails += 1
    print(f"{'OK ' if ok else 'BAD'} {label:26s} got={got}  ref={ref}")

def one(sql, params=None): return con.execute(sql, params or []).fetchone()[0]
def trow(sql):
    res = con.execute(sql)
    cols = [d[0] for d in res.description]
    r = res.fetchone()
    return dict(zip(cols, r))
def st(sql):
    # scalar macro returns STRUCT -> python dict
    return con.execute(f"SELECT {sql}").fetchone()[0]

# ---- A descriptive
check("stat_pctl(.5)", one("SELECT stat_pctl(array_agg(x),0.5) FROM t1"), np.median(x1))
check("stat_pctl(.9)", one("SELECT stat_pctl(array_agg(x),0.9) FROM t1"), np.percentile(x1,90))
check("geometric_mean", one("SELECT stat_geometric_mean(array_agg(x)) FROM t1"), refs['geom_mean'])
check("harmonic_mean", one("SELECT stat_harmonic_mean(array_agg(x)) FROM t1"), refs['harm_mean'])
check("stat_cv", one("SELECT stat_cv(array_agg(x)) FROM t1"), refs['cv'])
check("stat_iqr", one("SELECT stat_iqr(array_agg(x)) FROM t1"), refs['iqr'])
check("stat_mad", one("SELECT * FROM stat_mad((SELECT array_agg(x) FROM t1))"), np.median(np.abs(x1-np.median(x1))))
check("trimmed_mean", one("SELECT stat_trimmed_mean(array_agg(x),0.1) FROM t1"), refs['trimmed_mean'])

# ---- B CIs
c = st("ci_mean((SELECT array_agg(x) FROM t1))")
check("ci_mean.lower", c['lower'], refs['ci_mean'][0]); check("ci_mean.upper", c['upper'], refs['ci_mean'][1])
w = st("ci_mean_diff_welch((SELECT array_agg(x) FROM t1),(SELECT array_agg(x) FROM t2))")
dfw = refs['welch_df']; se = math.sqrt(x1.var(ddof=1)/len(x1)+x2.var(ddof=1)/len(x2))
check("welch.df", w['df'], dfw)
check("welch.lower", w['lower'], (x1.mean()-x2.mean())-S.t.ppf(0.975,dfw)*se)
check("welch.upper", w['upper'], (x1.mean()-x2.mean())+S.t.ppf(0.975,dfw)*se)
pp = st("ci_mean_diff_paired((SELECT array_agg(d) FROM pd))")
check("paired.point", pp['point'], pd_.mean())
check("paired.lower", pp['lower'], pd_.mean()-S.t.ppf(0.975,len(pd_)-1)*pd_.std(ddof=1)/math.sqrt(len(pd_)))
wp = st("ci_proportion(12,40)")
check("wilson.lower", wp['lower'], refs['ci_wilson'][0]); check("wilson.upper", wp['upper'], refs['ci_wilson'][1])
wv = st("ci_variance((SELECT array_agg(x) FROM t1))")
n=len(x1); s2=x1.var(ddof=1)
check("varCI.lower", wv['lower'], (n-1)*s2/S.chi2.ppf(0.975,n-1))
check("varCI.upper", wv['upper'], (n-1)*s2/S.chi2.ppf(0.025,n-1))
cp = st("ci_proportion_diff(12,40,30,80)")
check("propdiff.diff", cp['diff'], 12/40-30/80)

# ---- C sample size / power
check("ss_moe_prop(.05)", one("SELECT sample_size_moe_prop(0.05)"), refs['ss_moe_prop_005'], tol=2)
check("ss_moe_mean", one("SELECT sample_size_moe_mean(1.5,4.0)"), refs['ss_moe_mean_1.5'], tol=2)
check("ss_2samp(d=.5)", one("SELECT sample_size_2samp(0.5)"), refs['n_per_group_d0.5'], tol=2)
check("power_2samp(50,.5)", one("SELECT power_2samp(50,0.5)"), refs['power_2samp_n50_d0.5'])
check("power_1samp(100,.5)", one("SELECT power_1samp(100,0.5)"), refs['power_1samp_n100_d0.5'])

# ---- D tests
z = st("z_test((SELECT array_agg(x) FROM t1),100,15.0)")
zt = st("z_test_proportion(12,40,0.5)")
# manual z references (scipy 1.18 removed ztest)
m = x1.mean(); se_z = 15.0/math.sqrt(len(x1))
zref = (m-100)/se_z
check("z_test.stat", z['statistic'], zref)
check("z_test.p", z['p_value'], 2*(1-S.norm.cdf(abs(zref))))
p0 = 0.5
zref2 = (12/40 - p0)/math.sqrt(p0*(1-p0)/40)
check("z_prop.stat", zt['statistic'], zref2)
check("z_prop.p", zt['p_value'], 2*(1-S.norm.cdf(abs(zref2))))
zf = st("f_test_var((SELECT array_agg(x) FROM t1),(SELECT array_agg(x) FROM t2))")
F = x1.var(ddof=1)/x2.var(ddof=1)
pv = 2*min(S.f.cdf(F,len(x1)-1,len(x2)-1), S.f.sf(F,len(x1)-1,len(x2)-1))
check("f_var.stat", zf['statistic'], F); check("f_var.p", zf['p_value'], pv)

# ---- E ANOVA / KW
kw = trow("SELECT * FROM kruskal_wallis((SELECT array_agg(x ORDER BY id) FROM av),(SELECT array_agg(grp ORDER BY id) FROM av))")
rk = S.kruskal(*g)
check("kruskal.H", kw['h_statistic'], rk.statistic); check("kruskal.p", kw['p_value'], rk.pvalue)
av = trow("SELECT * FROM oneway_anova((SELECT array_agg(x ORDER BY id) FROM av),(SELECT array_agg(grp ORDER BY id) FROM av))")
rf = S.f_oneway(*g)
check("anova.F", av['f_stat'], rf.statistic); check("anova.p", av['p_value'], rf.pvalue)
tw = trow("SELECT * FROM twoway_anova((SELECT array_agg(x ORDER BY id) FROM av2),(SELECT array_agg(a ORDER BY id) FROM av2),(SELECT array_agg(b ORDER BY id) FROM av2))")
# reference: two-way ANOVA SS decomposition (fully crossed, pure python)
grand = sum(vals)/len(vals)
_a_levels = sorted(set(aa)); _b_levels = sorted(set(bb))
_cm = {}  # (a,b) -> (mean, count)
for i in range(len(vals)):
    key = (aa[i], bb[i])
    if key not in _cm: _cm[key] = [0.0, 0]
    _cm[key][0] += vals[i]; _cm[key][1] += 1
for key in _cm: _cm[key][0] /= _cm[key][1]
def _mean_of(f, lvl):
    xs = [vals[i] for i in range(len(vals)) if (f=='a' and aa[i]==lvl) or (f=='b' and bb[i]==lvl)]
    return sum(xs)/len(xs)
def _count_of(f, lvl):
    return sum(1 for i in range(len(vals)) if (f=='a' and aa[i]==lvl) or (f=='b' and bb[i]==lvl))
am = {l: _mean_of('a', l) for l in _a_levels}; bm = {l: _mean_of('b', l) for l in _b_levels}
na = {l: _count_of('a', l) for l in _a_levels}; nb = {l: _count_of('b', l) for l in _b_levels}
ss_a = sum((_mean_of('a',l)-grand)**2 * na[l] for l in _a_levels)
ss_b = sum((_mean_of('b',l)-grand)**2 * nb[l] for l in _b_levels)
ss_ab = 0.0
for a_ in _a_levels:
    for b_ in _b_levels:
        if (a_, b_) in _cm:
            ss_ab += (_cm[(a_,b_)][0] - am[a_] - bm[b_] + grand)**2 * _cm[(a_,b_)][1]
sse = 0.0
for i in range(len(vals)):
    sse += (vals[i] - _cm[(aa[i], bb[i])][0])**2
n_tot = len(vals); k_ = len(_a_levels); j_ = len(_b_levels)
check("twoway.ss_a", tw['ss_a'], ss_a); check("twoway.ss_b", tw['ss_b'], ss_b)
check("twoway.ss_ab", tw['ss_ab'], ss_ab); check("twoway.sse", tw['ss_error'], sse)
check("twoway.f_a", tw['f_a'], (ss_a/(k_-1))/(sse/(n_tot-k_*j_)))
check("twoway.f_b", tw['f_b'], (ss_b/(j_-1))/(sse/(n_tot-k_*j_)))
check("twoway.f_ab", tw['f_ab'], (ss_ab/((k_-1)*(j_-1)))/(sse/(n_tot-k_*j_)))
tu = con.execute("SELECT count(*) FROM tukey_hsd((SELECT array_agg(x ORDER BY id) FROM av),(SELECT array_agg(grp ORDER BY id) FROM av))").fetchone()[0]
print("tukey rows:", tu, "(expect 3 pairs)")
sp = one("SELECT scheffe_p(5.0, 15, 15, 3, 42, 100.0)")
print("scheffe p:", sp)
dx = [10.0, 3.0, -2.0]; sech=[1.2,1.2,1.2]
du = con.execute("SELECT count(*) FROM dunnett_approx(?::DOUBLE[],?::DOUBLE[])", [dx, sech]).fetchone()[0]
print("dunnett rows:", du)

# ---- F effect sizes
check("cohens_d_1", one("SELECT cohens_d((SELECT array_agg(x) FROM t1),100)"), refs['cohens_d_1'])
check("cohens_d_2", one("SELECT cohens_d_2samp((SELECT array_agg(x) FROM t1),(SELECT array_agg(x) FROM t2))"), refs['cohens_d_2'])
check("hedges_g_2", one("SELECT hedges_g_2samp((SELECT array_agg(x) FROM t1),(SELECT array_agg(x) FROM t2))"), refs['hedges_g_2'])
check("glass_delta", one("SELECT glass_delta((SELECT array_agg(x) FROM t1),(SELECT array_agg(x) FROM t2))"), refs['glass'])
check("cohens_h", one("SELECT cohens_h(0.3,0.375)"), refs['cohens_h'])
check("t_to_r", one("SELECT t_to_r(1.5,43)"), 1.5/math.sqrt(1.5**2+43))
check("r_to_d", one("SELECT r_to_d(0.5)"), 0.5/math.sqrt(1-0.25))

# ---- G categorical
o = st("odds_ratio(20,10,5,15)")
check("odds_ratio", o['point'], refs['odds_ratio']); check("OR.lower", o['lower'], refs['odds_ratio_ci'][0])
r = st("risk_ratio(20,10,5,15)")
check("risk_ratio", r['point'], refs['risk_ratio'])
rd = st("risk_difference(20,10,5,15)")
check("risk_diff", rd['point'], refs['risk_diff']); check("nnt", rd['nnt'], refs['nnt'])
fp = one("SELECT fisher_exact_p(20,10,5,15)")
rfp = S.fisher_exact([[20,10],[5,15]])[1]
check("fisher_exact_p", fp, rfp, tol=1e-6)

# ---- H time series
x = [3.0,4.0,5.0,6.0,5.0,4.0,3.0,4.0,5.0,6.0]
acf1 = one("SELECT acf(?::DOUBLE[], 1)", [x])
# manual acf lag1
m = np.mean(x); num = sum((x[i]-m)*(x[i+1]-m) for i in range(len(x)-1)); den = sum((v-m)**2 for v in x)
check("acf(1)", acf1, num/den)
act = one("SELECT list_avg(list_transform(?::DOUBLE[], v -> 1.0))", [x])  # sanity
print("\n=== TOTAL FAILS:", fails, "===")
