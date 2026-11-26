"""Round-3 review responses: annotation-dropout sensitivity and transcript localization.

Three products, all sample-level and none of them a pooled-cell test:

1. `05_stats/sensitivity_lowconf.csv`
   LowConfidence cells are reassigned to their top-scoring panel class
   (`celltype_raw`) instead of being dropped, and the primary ox_core test is
   recomputed with the identical pipeline (per-sample downsample, median
   statistic, exact C(6,3)=20 permutation, Cliff's delta).  This tests whether
   the cell-type-specific primary result is driven by exclusion of low-margin
   clusters, which is a real concern because dropout differed between arms.

2. `08_audit/TableS3_ox_genes_by_celltype.csv`
   Per dataset x cell type detection rate and mean log-normalised expression of
   AIFM1, KEAP1 and PGAM5 among high-confidence cells only.  This is the
   evidence base for the transcript-level "localization" wording in the title.

3. `08_audit/TableS4_lowconf_by_sample.csv`
   Per-sample QC, LowConfidence count and median genes after QC, so that
   group-level differences in annotation dropout are auditable rather than
   summarised as a single pooled percentage.
"""
import hashlib
import itertools
import json
import warnings

import numpy as np
import pandas as pd
import scanpy as sc

from common import (ANN, SCORE, STATS, AUDIT, DATASETS, GENE_SETS,
                    QC_THRESHOLDS, STATS_CONFIG)

warnings.filterwarnings("ignore")
sc.settings.verbosity = 0
SEED = STATS_CONFIG["seed"]
T = QC_THRESHOLDS
CORE_KEY = "oxeiptosis_core"


def _loc(x):
    return float(np.median(x))


def exact_perm(a, b):
    """a = exposed, b = control; median difference; exact over all labelings."""
    vals = np.concatenate([a, b]).astype(float)
    n = len(a)
    obs = _loc(vals[:n]) - _loc(vals[n:])
    stats = []
    for comb in itertools.combinations(np.arange(len(vals)), n):
        mask = np.zeros(len(vals), dtype=bool)
        mask[list(comb)] = True
        stats.append(_loc(vals[mask]) - _loc(vals[~mask]))
    stats = np.asarray(stats)
    return obs, float(np.mean(np.abs(stats) >= abs(obs) - 1e-12))


def cliffs_delta(a, b):
    gt = sum(1 for x in a for y in b if x > y)
    lt = sum(1 for x in a for y in b if x < y)
    return (gt - lt) / (len(a) * len(b))


def bh(p):
    p = np.asarray(p, float)
    ok = ~np.isnan(p)
    q = np.full_like(p, np.nan)
    pp = p[ok]
    if not len(pp):
        return q
    order = np.argsort(pp)
    ranked = pp[order] * len(pp) / (np.arange(len(pp)) + 1)
    ranked = np.minimum.accumulate(ranked[::-1])[::-1]
    out = np.empty_like(ranked)
    out[order] = np.clip(ranked, 0, 1)
    q[ok] = out
    return q


def downsample_mean(vals, key):
    if len(vals) == 0:
        return np.nan
    stable = int(hashlib.md5(key.encode()).hexdigest()[:8], 16)
    rng = np.random.default_rng(SEED + stable)
    if len(vals) > T["downsample_per_sample_celltype"]:
        vals = rng.choice(vals, T["downsample_per_sample_celltype"], replace=False)
    return float(np.mean(vals))


sens_rows, gene_rows, drop_rows = [], [], []

for gse, cfg in DATASETS.items():
    if cfg["tier"] != "INFERENCE":
        continue
    ad = sc.read_h5ad(ANN / f"{gse}_annotated.h5ad")
    lookup = {str(g).lower(): str(g) for g in ad.var_names}
    core = [lookup[g.lower()] for g in GENE_SETS[CORE_KEY] if g.lower() in lookup]
    sc.tl.score_genes(ad, gene_list=core, score_name="ox_core_sens",
                      ctrl_size=50, n_bins=25, rng=SEED)

    obs = ad.obs
    is_low = obs["celltype"].astype(str) == "LowConfidence"
    # fallback label: top-scoring panel class, kept in celltype_raw by 10_qc_annotate.py
    label_sens = np.where(is_low, obs["celltype_raw"].astype(str), obs["celltype"].astype(str))
    obs = obs.assign(celltype_sens=label_sens, ox_core_sens=ad.obs["ox_core_sens"].astype(float))

    # ---- Table S4: per-sample dropout
    for smp, g in obs.groupby("sample", observed=True):
        n_low = int((g["celltype"].astype(str) == "LowConfidence").sum())
        drop_rows.append({
            "dataset": gse, "sample": smp, "group": str(g["group"].iloc[0]),
            "n_cells_qc": int(len(g)), "n_lowconfidence": n_low,
            "pct_lowconfidence": round(100 * n_low / len(g), 2),
            "median_genes_all": float(np.median(g["n_genes_by_counts"])),
            "median_genes_highconf": float(np.median(
                g.loc[g["celltype"].astype(str) != "LowConfidence", "n_genes_by_counts"])),
        })

    # ---- Table S3: transcript detection per cell type (high-confidence only)
    hi = obs[obs["celltype"].astype(str) != "LowConfidence"]
    for ct, g in hi.groupby("celltype", observed=True):
        if len(g) < T["min_cells_per_celltype_per_sample"]:
            continue
        row = {"dataset": gse, "celltype": str(ct), "n_cells": int(len(g))}
        sub = ad[g.index, core]
        mat = np.asarray(sub.X.todense()) if hasattr(sub.X, "todense") else np.asarray(sub.X)
        for j, gene in enumerate(core):
            v = mat[:, j]
            row[f"pct_expressing_{gene}"] = round(100 * float((v > 0).mean()), 2)
            row[f"mean_expr_{gene}"] = round(float(v.mean()), 4)
        gene_rows.append(row)

    # ---- sensitivity test on fallback labels
    prim = pd.read_csv(STATS / "sample_level_tests.csv")
    prim = prim[(prim.dataset == gse) & (prim.module == "ox_core")]
    for ct, g in obs.groupby("celltype_sens", observed=True):
        vals = pd.DataFrame({"sample": g["sample"].astype(str),
                             "group": g["group"].astype(str),
                             "v": g["ox_core_sens"]})
        agg_rows = [{"sample": smp, "group": grp,
                     "v": downsample_mean(s.to_numpy(), f"{gse}|{ct}|{smp}")}
                    for (smp, grp), s in vals.groupby(["sample", "group"])["v"]]
        agg = pd.DataFrame(agg_rows)
        a = agg.loc[agg.group == "Control", "v"].dropna().to_numpy()
        b = agg.loc[agg.group != "Control", "v"].dropna().to_numpy()
        if len(a) < 3 or len(b) < 3:
            continue
        delta, p = exact_perm(b, a)
        prow = prim[prim.celltype == str(ct)]
        delta_primary = float(prow.delta_exposed_minus_control.iloc[0]) if len(prow) else np.nan
        p_primary = float(prow.permutation_p_two_sided.iloc[0]) if len(prow) else np.nan
        sens_rows.append({
            "dataset": gse, "celltype": str(ct),
            "n_control": len(a), "n_exposed": len(b),
            "delta_primary": delta_primary,
            "delta_sensitivity": delta,
            "cliffs_delta_sensitivity": cliffs_delta(b, a),
            "p_sensitivity": p,
            "p_primary": p_primary,
            "sign_agreement": (np.sign(delta) == np.sign(delta_primary))
                              if np.isfinite(delta_primary) else np.nan,
            "celltype_new_in_sensitivity": bool(not len(prow)),
        })

sdf = pd.DataFrame(sens_rows)
sdf["FDR_BH_sensitivity_A"] = bh(sdf["p_sensitivity"].to_numpy())
sdf.to_csv(STATS / "sensitivity_lowconf.csv", index=False, encoding="utf-8-sig")

gdf = pd.DataFrame(gene_rows)
gdf.to_csv(AUDIT / "TableS3_ox_genes_by_celltype.csv", index=False, encoding="utf-8-sig")

ddf = pd.DataFrame(drop_rows)
ddf.to_csv(AUDIT / "TableS4_lowconf_by_sample.csv", index=False, encoding="utf-8-sig")

summary = {
    "n_strata": int(len(sdf)),
    "n_sign_agree": int((sdf.sign_agreement == True).sum()),  # noqa: E712
    "n_sign_disagree": int((sdf.sign_agreement == False).sum()),  # noqa: E712
    "n_new_in_sensitivity": int(sdf.celltype_new_in_sensitivity.sum()),
    "min_p_sensitivity": float(sdf.p_sensitivity.min()) if len(sdf) else None,
    "n_sensitivity_FDR005": int((sdf.FDR_BH_sensitivity_A < STATS_CONFIG["alpha"]).sum()),
    "max_abs_delta_primary": float(sdf.delta_primary.abs().max()) if len(sdf) else None,
    "max_abs_delta_sensitivity": float(sdf.delta_sensitivity.abs().max()) if len(sdf) else None,
    "seed": SEED,
}
(STATS / "sensitivity_summary.json").write_text(
    json.dumps(summary, indent=2), encoding="utf-8")

print("WROTE", STATS / "sensitivity_lowconf.csv", "strata=", len(sdf))
print("WROTE", AUDIT / "TableS3_ox_genes_by_celltype.csv", "rows=", len(gdf))
print("WROTE", AUDIT / "TableS4_lowconf_by_sample.csv", "rows=", len(ddf))
print(json.dumps(summary, indent=2))
