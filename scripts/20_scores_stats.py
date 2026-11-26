"""Two-stage design: per-cell module scores -> one value per sample x cell type
-> exact sample-level permutation tests with effect sizes and bootstrap CIs.

No pooled-cell differential testing is performed anywhere in this pipeline.
"""
import hashlib
import itertools
import json
import warnings
from pathlib import Path

import numpy as np
import pandas as pd
import scanpy as sc

from common import (ANN, SCORE, STATS, DATASETS, GENE_SETS, QC_THRESHOLDS,
                    STATS_CONFIG)

warnings.filterwarnings("ignore")
sc.settings.verbosity = 0
SEED = STATS_CONFIG["seed"]
T = QC_THRESHOLDS


def bh_fdr(p):
    p = np.asarray(p, dtype=float)
    ok = ~np.isnan(p)
    q = np.full_like(p, np.nan)
    pp = p[ok]
    if len(pp) == 0:
        return q
    order = np.argsort(pp)
    ranked = pp[order] * len(pp) / (np.arange(len(pp)) + 1)
    ranked = np.minimum.accumulate(ranked[::-1])[::-1]
    out = np.empty_like(ranked)
    out[order] = np.clip(ranked, 0, 1)
    q[ok] = out
    return q


def _loc(x):
    """Protocol-specified location statistic: median (Round 2 correction)."""
    return float(np.median(x)) if STATS_CONFIG.get("delta_statistic") == "median" else float(np.mean(x))


def exact_perm(a, b):
    """Two-sided exact permutation p for the difference in location; C(6,3)=20 labelings.

    a = exposed, b = control -> obs = loc(exposed) - loc(control).
    """
    vals = np.concatenate([a, b]).astype(float)
    n = len(a)
    obs = _loc(vals[:n]) - _loc(vals[n:])
    idx = np.arange(len(vals))
    stats = []
    for comb in itertools.combinations(idx, n):
        mask = np.zeros(len(vals), dtype=bool)
        mask[list(comb)] = True
        stats.append(_loc(vals[mask]) - _loc(vals[~mask]))
    stats = np.asarray(stats)
    p = np.mean(np.abs(stats) >= abs(obs) - 1e-12)
    return obs, float(p), len(stats)


def cliffs_delta(a, b):
    gt = sum(1 for x in a for y in b if x > y)
    lt = sum(1 for x in a for y in b if x < y)
    return (gt - lt) / (len(a) * len(b))


def hedges_g(a, b):
    n1, n2 = len(a), len(b)
    s1, s2 = np.var(a, ddof=1), np.var(b, ddof=1)
    sp = np.sqrt(((n1 - 1) * s1 + (n2 - 1) * s2) / (n1 + n2 - 2))
    if sp == 0:
        return np.nan
    d = (np.mean(a) - np.mean(b)) / sp
    J = 1 - 3 / (4 * (n1 + n2) - 9)
    return float(J * d)


def boot_ci(a, b, fn, reps=None, seed=SEED):
    """Percentile CI with an explicit valid-replicate proportion.

    Round 2 (statistics reviewer M5): with n=3 a sizeable fraction of resamples
    can have zero within-group variance (Hedges' g undefined).  Returning a CI
    computed from the surviving subset is misleading, so report the valid
    proportion and mark the CI as NA when < 0.95 of replicates are valid.
    """
    reps = int(reps or STATS_CONFIG["bootstrap_reps"])
    rng = np.random.default_rng(seed)
    out = []
    for _ in range(reps):
        aa = rng.choice(a, len(a), replace=True)
        bb = rng.choice(b, len(b), replace=True)
        out.append(fn(aa, bb))
    arr = np.asarray(out, dtype=float)
    valid = np.isfinite(arr)
    frac = float(valid.mean())
    if not valid.any():
        return [np.nan, np.nan], frac
    lo, hi = float(np.percentile(arr[valid], 2.5)), float(np.percentile(arr[valid], 97.5))
    if frac < 0.95:
        return [np.nan, np.nan], frac
    return [lo, hi], frac


# ---------------------------------------------------------------- scoring
lookup_cache = {}
sample_rows = []
module_names = {
    "ox_core": "oxeiptosis_core",
    "nrf2_antag": "nrf2_antagonistic",
    "ferro_driver": "ferroptosis_driver",
    "ferro_suppressor": "ferroptosis_suppressor",
}
gene_presence = {}

for gse, cfg in DATASETS.items():
    if cfg["tier"] != "INFERENCE":
        continue
    if not (ANN / f"{gse}_annotated.h5ad").exists():
        print("SKIP scoring", gse, "(annotation not ready)")
        continue
    ad = sc.read_h5ad(ANN / f"{gse}_annotated.h5ad")
    lookup = {str(g).lower(): str(g) for g in ad.var_names}
    gene_presence[gse] = {}
    for key, setname in module_names.items():
        wanted = GENE_SETS[setname]
        present = [lookup[g.lower()] for g in wanted if g.lower() in lookup]
        gene_presence[gse][key] = {
            "present": present,
            "missing": [g for g in wanted if g.lower() not in lookup],
        }
        if len(present) > 0:
            sc.tl.score_genes(ad, gene_list=present, score_name=key,
                              ctrl_size=50, n_bins=25, rng=SEED)
        else:
            ad.obs[key] = 0.0
    # sensitivity: plain mean expression of core triad (no control-gene subtraction)
    core = gene_presence[gse]["ox_core"]["present"]
    if core:
        ad.obs["ox_core_meanexpr"] = np.asarray(ad[:, core].X.mean(axis=1)).ravel()
    else:
        ad.obs["ox_core_meanexpr"] = 0.0
    ad.obs["ferro_net"] = ad.obs["ferro_driver"].astype(float) - ad.obs["ferro_suppressor"].astype(float)

    # LowConfidence clusters are excluded from cell-type-specific scores but
    # RETAINED in the per-sample denominator, so fractions remain interpretable
    # as a share of all captured cells (biology reviewer M2; see amendments R2-4).
    _denom = ad.obs["sample"].value_counts().to_dict()
    _obs = ad.obs[ad.obs["celltype"].astype(str) != "LowConfidence"]
    for (ct, smp), grp in _obs.groupby(["celltype", "sample"], observed=True):
        n = len(grp)
        row = {"dataset": gse, "celltype": str(ct), "sample": str(smp),
               "group": str(grp["group"].iloc[0]), "n_cells": n,
               "fraction_of_sample": n / _denom[smp]}
        for key in list(module_names) + ["ferro_net", "ox_core_meanexpr"]:
            vals = grp[key].astype(float).to_numpy()
            if n >= T["min_cells_per_celltype_per_sample"]:
                stable = int(hashlib.md5(f"{gse}|{ct}|{smp}".encode()).hexdigest()[:8], 16)
                rng = np.random.default_rng(SEED + stable)
                if n > T["downsample_per_sample_celltype"]:
                    vals = rng.choice(vals, T["downsample_per_sample_celltype"], replace=False)
                row[key] = float(np.mean(vals))
            else:
                row[key] = np.nan
        sample_rows.append(row)

sdf = pd.DataFrame(sample_rows)
sdf.to_csv(SCORE / "sample_scores.csv", index=False, encoding="utf-8-sig")
(SCORE / "gene_presence.json").write_text(
    json.dumps(gene_presence, indent=2, ensure_ascii=False), encoding="utf-8")
print("WROTE", SCORE / "sample_scores.csv", "rows=", len(sdf))

# ------------------------------------------------------------------ tests
test_rows = []
comp_rows = []
modules = list(module_names) + ["ferro_net", "ox_core_meanexpr"]
for gse, cfg in DATASETS.items():
    if cfg["tier"] != "INFERENCE":
        continue
    sub = sdf[sdf.dataset == gse]
    for ct, g in sub.groupby("celltype"):
        # Composition is a cell-type-level test, computed ONCE per cell type.
        # It must NOT sit inside the module loop: that would duplicate each
        # row once per module and inflate the BH family size.
        ca = g.loc[g.group.isin(["Control"]), "fraction_of_sample"].dropna().to_numpy()
        cb = g.loc[~g.group.isin(["Control"]), "fraction_of_sample"].dropna().to_numpy()
        if len(ca) >= 3 and len(cb) >= 3:
            cobs, cp, cn = exact_perm(cb, ca)
            comp_rows.append({"dataset": gse, "celltype": ct,
                              "fraction_control": float(np.mean(ca)),
                              "fraction_exposed": float(np.mean(cb)),
                              "delta_fraction": cobs,
                              "permutation_p_two_sided": cp,
                              "n_permutations": cn})
        for mod in modules:
            a = g.loc[g.group.isin(["Control"]), mod].dropna().to_numpy()
            b = g.loc[~g.group.isin(["Control"]), mod].dropna().to_numpy()
            if len(a) < 3 or len(b) < 3:
                continue
            obs, p, nperm = exact_perm(b, a)   # b = exposed minus control
            d = cliffs_delta(b, a)
            gval = hedges_g(b, a)
            (ci_d_lo, ci_d_hi), frac_d = boot_ci(b, a, cliffs_delta)
            (ci_g_lo, ci_g_hi), frac_g = boot_ci(b, a, hedges_g)
            # Round 3 (statistics blocking 1): the finite-replicate fraction is
            # almost always >0.95 even when the interval itself is pathological
            # (e.g. [-45, -0.9]) because a single near-zero pooled SD explodes g.
            # The interval is therefore never reported; only the point estimate
            # and the valid fraction are retained, matching the manuscript.
            ci_g_lo, ci_g_hi = np.nan, np.nan
            test_rows.append({
                "dataset": gse, "celltype": ct, "module": mod,
                "n_control": len(a), "n_exposed": len(b),
                "median_control": float(np.median(a)), "median_exposed": float(np.median(b)),
                "mean_control": float(np.mean(a)), "mean_exposed": float(np.mean(b)),
                "delta_exposed_minus_control": obs,
                "permutation_p_two_sided": p, "n_permutations": nperm,
                "cliffs_delta": d, "cliffs_delta_CI95_lo": ci_d_lo, "cliffs_delta_CI95_hi": ci_d_hi,
                "cliffs_delta_boot_valid_fraction": frac_d,
                "hedges_g": gval, "hedges_g_CI95_lo": ci_g_lo, "hedges_g_CI95_hi": ci_g_hi,
                "hedges_g_boot_valid_fraction": frac_g,
                "bootstrap_reps": STATS_CONFIG["bootstrap_reps"],
            })

tdf = pd.DataFrame(test_rows)
# Round 3: family labels now follow the manuscript exactly (A / B1 / B2 / C).
# The pre-analysis plan declared A, B and C; family B was split into B1
# (secondary modules) and B2 (separation endpoint D) during implementation,
# before any test statistic was inspected.  That deviation is amendment A5 in
# 08_audit/amendments.md and must not be described as pre-registration.
t = tdf["module"].to_numpy()
famA = t == "ox_core"
famB = ~famA
tdf["FDR_family"] = np.where(famA, "A_primary_ox_core", "B1_secondary_modules")
tdf["FDR_BH_within_family"] = np.nan
tdf.loc[famA, "FDR_BH_within_family"] = bh_fdr(tdf.loc[famA, "permutation_p_two_sided"].to_numpy())
for ds in tdf.dataset.unique():
    m = (tdf.dataset == ds).to_numpy() & famB
    if m.any():
        tdf.loc[m, "FDR_BH_within_family"] = bh_fdr(tdf.loc[m, "permutation_p_two_sided"].to_numpy())
# Legacy pooled q retained only for comparison with the Round-1 manuscript; never used for claims.
tdf["FDR_BH_all_tests_legacy"] = bh_fdr(tdf["permutation_p_two_sided"].to_numpy())
tdf["FDR_BH_all_tests"] = tdf["FDR_BH_within_family"]
tdf["significant_FDR_005"] = tdf["FDR_BH_within_family"] < STATS_CONFIG["alpha"]
tdf.to_csv(STATS / "sample_level_tests.csv", index=False, encoding="utf-8-sig")
cdf = pd.DataFrame(comp_rows)
if len(cdf):
    cdf["FDR_BH_composition_within_dataset"] = np.nan
    for ds in cdf.dataset.unique():
        m = (cdf.dataset == ds).to_numpy()
        cdf.loc[m, "FDR_BH_composition_within_dataset"] = bh_fdr(cdf.loc[m, "permutation_p_two_sided"].to_numpy())
    cdf["FDR_BH_composition"] = cdf["FDR_BH_composition_within_dataset"]
cdf.to_csv(STATS / "composition_tests.csv", index=False, encoding="utf-8-sig")

meta = {
    "fdr_family": ("A:ox_core(celltype x dataset); B1:secondary modules(per dataset); "
                   "B2:separation D(per dataset); C:LR(per dataset)"),
    "fdr_family_declared_in_plan": "A, B, C (three families); B split into B1/B2 as amendment A5",
    "fdr_family_sizes": {
        "A_primary_ox_core": int(famA.sum()),
        "B1_secondary_by_dataset": {str(d): int(((tdf.dataset == d) & famB).sum()) for d in tdf.dataset.unique()},
    },
    "n_tests": len(tdf),
    "min_possible_two_sided_p": 2 / 20,
    "exact_permutations_per_test": 20,
    "bootstrap_reps": STATS_CONFIG["bootstrap_reps"],
    "delta_statistic": STATS_CONFIG["delta_statistic"],
    "primary_effect_size": "Cliff_delta",
    "seed": SEED,
    "downsample_per_sample_celltype": T["downsample_per_sample_celltype"],
    "max_cells_per_sample": T["max_cells_per_sample"],
}
(STATS / "stats_metadata.json").write_text(json.dumps(meta, indent=2), encoding="utf-8")
print("WROTE", STATS / "sample_level_tests.csv", "tests=", len(tdf),
      "min_p=", meta["min_possible_two_sided_p"])
