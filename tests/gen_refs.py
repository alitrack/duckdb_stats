"""Generate scipy/R reference values for stats_fill.sql verification. Standalone (no duckdb needed)."""
import json, math
import numpy as np
from scipy import stats as S
from scipy.stats import distributions as dists

np.random.seed(42)
n1, n2 = 20, 25
x1 = np.random.normal(100, 15, n1)
x2 = np.random.normal(105, 12, n2)
paired_d = np.random.normal(2, 5, 30)
w = np.random.uniform(0.5, 3.0, n1)

def ci_mean(x, alpha=0.05):
    n = len(x); se = x.std(ddof=1)/math.sqrt(n)
    m = x.mean(); t = S.t.ppf(1-alpha/2, n-1)
    return (m - t*se, m + t*se)

ci_welch = S.ttest_ind(x1, x2, equal_var=False)
ci_pooled = S.ttest_ind(x1, x2, equal_var=True)
ci_paired = S.ttest_1samp(paired_d, 0)
def wilson_ci(k, n, alpha=0.05):
    p = k/n; z = S.norm.ppf(1-alpha/2)
    denom = 1 + z*z/n
    center = (p + z*z/(2*n))/denom
    half = z*math.sqrt(p*(1-p)/n + z*z/(4*n*n))/denom
    return (center-half, center+half)
def wald_ci(k, n, alpha=0.05):
    p = k/n; z = S.norm.ppf(1-alpha/2); se = math.sqrt(p*(1-p)/n)
    return (p-z*se, p+z*se)
ci_wilson = wilson_ci(12, 40)
ci_wald = wald_ci(12, 40)
var_data = np.random.normal(0, 4, 40)
s2 = var_data.var(ddof=1)
chi_lo, chi_hi = S.chi2.ppf(0.025, 39), S.chi2.ppf(0.975, 39)
ci_varians = ((39*s2)/chi_hi, (39*s2)/chi_lo)
p1, p2, nA, nB = 12/40, 30/80, 40, 80
z = S.norm.ppf(0.975)
se_pdiff = math.sqrt(p1*(1-p1)/nA + p2*(1-p2)/nB)
ci_pdiff = ((p1-p2)-z*se_pdiff, (p1-p2)+z*se_pdiff)

# sample size for MOE proportion (worst case p=0.5)
def ss_moe_prop(moe, alpha=0.05):
    return (z**2 * 0.25) / (moe**2)
def ss_moe_mean(moe, sd, alpha=0.05):
    t = S.t.ppf(1-alpha/2, 30)  # approx iterative; do exact loop
    n = (S.norm.ppf(1-alpha/2)*sd/moe)**2
    for _ in range(50):
        t = S.t.ppf(1-alpha/2, n-1)
        n2 = (t*sd/moe)**2
        if abs(n2-n) < 1e-9: n = n2; break
        n = n2
    return n

# two-sample power (normal approx): effect d, alpha, power -> n per group
def n_per_group(d, alpha=0.05, power=0.80):
    za = S.norm.ppf(1-alpha/2); zb = S.norm.ppf(power)
    return 2*(za+zb)**2 / d**2
def power_2samp(n, d, alpha=0.05):
    za = S.norm.ppf(1-alpha/2)
    return S.norm.cdf(d*math.sqrt(n/2) - za)
def power_1samp(n, d, alpha=0.05):
    za = S.norm.ppf(1-alpha/2)
    return S.norm.cdf(d*math.sqrt(n) - za)

# one-way ANOVA: 3 groups, means 0/10/20, sd 10, n=15 each
g = [np.random.normal(m, 10, 15) for m in (0, 10, 20)]
allv = np.concatenate(g)
grand = allv.mean()
ss_b = sum(15*(c.mean()-grand)**2 for c in g)
ss_w = sum(((c-c.mean())**2).sum() for c in g)
dfb, dfw = 2, 42
F = (ss_b/dfb)/(ss_w/dfw)
p_anova = S.f.sf(F, dfb, dfw)
eta2 = ss_b/(ss_b+ss_w); omega2 = (ss_b-dfb*(ss_w/ss_w)) / (ss_b+ss_w+ (ss_w/dfw) )
omega2 = (ss_b - dfb*(ss_w/dfw)) / (ss_b + ss_w + (ss_w/dfw))
f_cohen = math.sqrt(eta2/(1-eta2))

# effect sizes
cohens_d_1 = (x1.mean()-100)/x1.std(ddof=1)
s_pooled = math.sqrt(((n1-1)*x1.var(ddof=1)+(n2-1)*x2.var(ddof=1))/(n1+n2-2))
cohens_d_2 = (x1.mean()-x2.mean())/s_pooled
df_tot = n1+n2-2
J = 1 - 3/(4*df_tot-1)
hedges_g_1 = J*(x1.mean()-100)/x1.std(ddof=1)
hedges_g_2 = J*cohens_d_2
glass = (x1.mean()-x2.mean())/x2.std(ddof=1)
tstat = ci_welch.statistic
t_to_r = math.sqrt(tstat**2/(tstat**2+df_tot))
d_to_r = cohens_d_1/math.sqrt(cohens_d_1**2+df_tot)
r_val = np.corrcoef(x1, x2[:n1])[0,1]
r_to_d = r_val/math.sqrt(1-r_val**2)
cohens_h = 2*math.asin(math.sqrt(12/40)) - 2*math.asin(math.sqrt(30/80))
# 2x2 table: a=20,b=10,c=5,d=15
a,b,c,d = 20,10,5,15
odds_ratio = (a*d)/(b*c)
risk_ratio = (a/(a+b))/(c/(c+d))
rr_se = math.sqrt(1/a - 1/(a+b) + 1/c - 1/(c+d))
risk_diff = a/(a+b) - c/(c+d)
rd_se = math.sqrt(a*(a+b)/( (a+b)**2 * a+b* (a+b)**2 ))  # will recompute properly below
rd_se = math.sqrt( a*(a+b)/((a+b)**3) + c*(c+d)/((c+d)**3) )
ldf = math.log(odds_ratio)
lor_se = math.sqrt(1/a+1/b+1/c+1/d)
nnt = 1/risk_diff
# CI for OR via log scale
or_ci = (math.exp(ldf-z*lor_se), math.exp(ldf+z*lor_se))
rr_ci = (risk_ratio*math.exp(-z*rr_se), risk_ratio*math.exp(z*rr_se))
rd_ci = (risk_diff - z*rd_se, risk_diff + z*rd_se)

out = {
 "data": {"x1": x1.tolist(), "x2": x2.tolist(), "paired_d": paired_d.tolist(),
          "w": w.tolist(), "groups": [c.tolist() for c in g]},
 "refs": {
   "ci_mean": ci_mean(x1),
   "ci_mean_se": x1.std(ddof=1)/math.sqrt(n1),
   "ci_mean_welch": (ci_welch.statistic, ci_welch.pvalue),
   "ci_mean_paired": (ci_paired.statistic, ci_paired.pvalue),
   "ci_wilson": ci_wilson,
   "ci_wald": ci_wald,
   "ci_varians": ci_varians,
   "ci_pdiff": ci_pdiff,
   "ss_moe_prop_005": ss_moe_prop(0.05),
   "ss_moe_prop_002": ss_moe_prop(0.02),
   "ss_moe_mean_1.5": ss_moe_mean(1.5, 4.0),
   "n_per_group_d0.5": n_per_group(0.5),
   "n_per_group_d0.8": n_per_group(0.8),
   "power_2samp_n50_d0.5": power_2samp(50, 0.5),
   "power_2samp_n100_d0.3": power_2samp(100, 0.3),
   "power_1samp_n100_d0.5": power_1samp(100, 0.5),
   "anova_F": F, "anova_p": p_anova, "anova_eta2": eta2, "anova_omega2": omega2, "anova_cohens_f": f_cohen,
   "anova_ssb": ss_b, "anova_ssw": ss_w,
   "cohens_d_1": cohens_d_1, "cohens_d_2": cohens_d_2,
   "hedges_g_1": hedges_g_1, "hedges_g_2": hedges_g_2, "glass": glass,
   "t_to_r": t_to_r, "d_to_r": d_to_r, "r_to_d": r_to_d, "cohens_h": cohens_h,
   "odds_ratio": odds_ratio, "odds_ratio_ci": or_ci,
   "risk_ratio": risk_ratio, "risk_ratio_ci": rr_ci,
   "risk_diff": risk_diff, "risk_diff_ci": rd_ci, "nnt": nnt,
   "mean_w": (x1*w).sum()/w.sum(),
   "geom_mean": math.exp(np.log(x1).mean()),
   "harm_mean": n1/(1/x1).sum(),
   "cv": x1.std(ddof=1)/x1.mean(),
   "mad": np.median(np.abs(x1-x1.mean())),
   "iqr": np.percentile(x1,75)-np.percentile(x1,25),
   "trimmed_mean": S.trim_mean(x1, 0.1),
 },
}
# welch df
num = (x1.var(ddof=1)/n1 + x2.var(ddof=1)/n2)**2
den = (x1.var(ddof=1)/n1)**2/(n1-1) + (x2.var(ddof=1)/n2)**2/(n2-1)
out["refs"]["welch_df"] = num/den
import os
json.dump(out, open(os.path.join(os.path.dirname(os.path.abspath(__file__)),"refs.json"),"w"))
print("refs written. sample: ci_mean=", out["refs"]["ci_mean"], "welch_df=", out["refs"]["welch_df"])
print("n_per_group_d0.5=", out["refs"]["n_per_group_d0.5"], "(R pwr.t.test expects 64)")
print("ss_moe_mean_1.5=", out["refs"]["ss_moe_mean_1.5"])
