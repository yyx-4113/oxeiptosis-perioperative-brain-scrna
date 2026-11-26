"""Secondary endpoints: oxeiptosis-ferroptosis separation and targeted LR panel.

LR analysis is a per-sample CellChat-style score over a PRE-SPECIFIED 32-pair
panel.  It is explicitly not a full CellChatDB run and is an indirect readout
because AIFM1/KEAP1/PGAM5 are intracellular.
"""
import hashlib
import itertools
import json
import warnings

import numpy as np
import pandas as pd
import scanpy as sc
from scipy.stats import spearmanr

from common import (ANN, STATS, SCORE, DATASETS, LR_PANEL, STATS_CONFIG,
                    QC_THRESHOLDS)

warnings.filterwarnings("ignore")
SEED = STATS_CONFIG["seed"]
T = QC_THRESHOLDS

sdf = pd.read_csv(SCORE / "sample_scores.csv")
tdf = pd.read_csv(STATS / "sample_level_tests.csv")


def _loc(x):
    return float(np.median(x)) if STATS_CONFIG.get("delta_statistic") == "median" else float(np.mean(x))


def exact_perm(a, b):
    vals = np.concatenate([a, b]).astype(float)
    n = len(a)
    obs = _loc(vals[:n]) - _loc(vals[n:])
    idx = np.arange(len(vals))
    stats = []
    for comb in itertools.combinations(idx, n):
        mask = np.zeros(len(vals), dtype=bool)
        mask[list(comb)] = True
        stats.append(_loc(vals[mask]) - _loc(vals[~mask]))
    return obs, float(np.mean(np.abs(np.asarray(stats)) >= abs(obs) - 1e-12)), 20


def boot_ci_two_sample(a, b, reps=None, seed=SEED):
    """Two-sample bootstrap CI for loc(b) - loc(a).

    Round-2 review N5: the previous one-sample implementation resampled only
    the exposed group while treating the control median as fixed, producing
    systematically too-narrow intervals that could exclude zero even when the
    exact permutation p was large.  Both groups are now resampled jointly.
    """
    reps = int(reps or STATS_CONFIG["bootstrap_reps"])
    rng = np.random.default_rng(seed)
    a = np.asarray(a, float); b = np.asarray(b, float)
    draws = np.asarray([_loc(rng.choice(b, len(b), replace=True)) -
                        _loc(rng.choice(a, len(a), replace=True)) for _ in range(reps)])
    return float(np.percentile(draws, 2.5)), float(np.percentile(draws, 97.5))


def bh(p):
    p = np.asarray(p, float)
    ok = ~np.isnan(p)
    q = np.full_like(p, np.nan)
    pp = p[ok]
    if len(pp) == 0:
        return q
    order = np.argsort(pp)
    r = pp[order] * len(pp) / (np.arange(len(pp)) + 1)
    r = np.minimum.accumulate(r[::-1])[::-1]
    out = np.empty_like(r)
    out[order] = np.clip(r, 0, 1)
    q[ok] = out
    return q


# ------------------------------------------------------------- separation
sep_rows = []
for gse, cfg in DATASETS.items():
    if cfg["tier"] != "INFERENCE":
        continue
    for ct, g in sdf[sdf.dataset == gse].groupby("celltype"):
        ctrl = g[g.group == "Control"]
        exp = g[g.group != "Control"]
        if len(ctrl) < 3 or len(exp) < 3:
            continue
        for mod in ("ox_core", "ferro_net", "nrf2_antag"):
            a = ctrl[mod].dropna().to_numpy()
            b = exp[mod].dropna().to_numpy()
            if len(a) < 3 or len(b) < 3:
                continue
            d, p, nperm = exact_perm(b, a)
            sep_rows.append({"dataset": gse, "celltype": ct, "module": mod,
                             "delta": d, "p": p,
                             "control_values": ";".join(f"{x:.4f}" for x in a),
                             "exposed_values": ";".join(f"{x:.4f}" for x in b)})
sep = pd.DataFrame(sep_rows)
piv = sep.pivot_table(index=["dataset", "celltype"], columns="module", values="delta")

# Round 2 (statistics M3 / biology M3): the separation endpoint must carry
# uncertainty and must not compare scores with different control-gene baselines.
# (1) paired per-sample difference after within-dataset module z-scoring:
#     d_i = z(ox_core_i) - z(ferro_net_i), same samples for both modules.
# (2) exact permutation p (20 labelings) and bootstrap 95% CI on the median d_i.
# (3) effect-size-scale sensitivity D = |Cliff's delta(ox)| - |Cliff's delta(ferro)|.
def cliffs_delta(a, b):
    gt = sum(1 for x in a for y in b if x > y)
    lt = sum(1 for x in a for y in b if x < y)
    return (gt - lt) / (len(a) * len(b))


sep_stats = []
for gse, cfg in DATASETS.items():
    if cfg["tier"] != "INFERENCE":
        continue
    sub = sdf[sdf.dataset == gse].copy()
    for mod in ("ox_core", "ferro_net"):
        mu, sd = sub[mod].mean(), sub[mod].std(ddof=0)
        sub["z_" + mod] = (sub[mod] - mu) / sd if sd and np.isfinite(sd) else 0.0
    for ct, g in sub.groupby("celltype"):
        paired = g[["sample", "group", "z_ox_core", "z_ferro_net", "ox_core", "ferro_net"]].dropna()
        ctrl = paired[paired.group == "Control"]
        exp = paired[paired.group != "Control"]
        if len(ctrl) < 3 or len(exp) < 3:
            continue
        paired = paired.assign(d=paired["z_ox_core"] - paired["z_ferro_net"])
        a = paired.loc[paired.group == "Control", "d"].to_numpy()
        b = paired.loc[paired.group != "Control", "d"].to_numpy()
        d_obs, p, nperm = exact_perm(b, a)
        lo, hi = boot_ci_two_sample(a, b)
        c_ox = cliffs_delta(exp["ox_core"].to_numpy(), ctrl["ox_core"].to_numpy())
        c_fe = cliffs_delta(exp["ferro_net"].to_numpy(), ctrl["ferro_net"].to_numpy())
        sep_stats.append({"dataset": gse, "celltype": ct,
                          "D_zscore_exposed_minus_control": d_obs,
                          "permutation_p_two_sided": p, "n_permutations": nperm,
                          "CI95_lo": lo, "CI95_hi": hi,
                          "D_effectsize_sensitivity": abs(c_ox) - abs(c_fe),
                          "n_control": len(a), "n_exposed": len(b)})
sep_df = pd.DataFrame(sep_stats)
if len(sep_df):
    # Family B: BH within each dataset (per protocol §7.3 and amendments R2-1)
    sep_df["FDR_BH_family_B_within_dataset"] = np.nan
    for ds in sep_df.dataset.unique():
        m = (sep_df.dataset == ds).to_numpy()
        sep_df.loc[m, "FDR_BH_family_B_within_dataset"] = bh(sep_df.loc[m, "permutation_p_two_sided"].to_numpy())
    sep_df["significant_FDR_005"] = sep_df["FDR_BH_family_B_within_dataset"] < 0.05
sep_df.to_csv(STATS / "separation_tests.csv", index=False, encoding="utf-8-sig")

piv = piv.reset_index()
if len(sep_df):
    piv = piv.merge(sep_df[["dataset", "celltype", "D_zscore_exposed_minus_control",
                            "permutation_p_two_sided", "CI95_lo", "CI95_hi",
                            "FDR_BH_family_B_within_dataset", "D_effectsize_sensitivity"]],
                    on=["dataset", "celltype"], how="left")
piv["separation_ox_minus_ferro"] = piv["ox_core"] - piv["ferro_net"]
sep.to_csv(STATS / "separation_components.csv", index=False, encoding="utf-8-sig")
piv.to_csv(STATS / "separation_summary.csv", index=False, encoding="utf-8-sig")

corr_rows = []
for gse, cfg in DATASETS.items():
    if cfg["tier"] != "INFERENCE":
        continue
    for ct, g in sdf[sdf.dataset == gse].groupby("celltype"):
        x = g["ox_core"].to_numpy(float)
        y = g["ferro_net"].to_numpy(float)
        ok = ~(np.isnan(x) | np.isnan(y))
        if ok.sum() >= 4:
            r, p = spearmanr(x[ok], y[ok])
            corr_rows.append({"dataset": gse, "celltype": ct, "n_samples": int(ok.sum()),
                              "spearman_rho": float(r), "p": float(p)})
corrdf = pd.DataFrame(corr_rows)
if len(corrdf):
    # Exploratory only, pooled across groups; not reported in the manuscript.
    corrdf["exploratory_unadjusted"] = True
corrdf.to_csv(STATS / "ox_ferro_correlation.csv", index=False, encoding="utf-8-sig")

# ------------------------------------------------------------- LR panel
lr_rows = []
lr_presence = {}
for gse, cfg in DATASETS.items():
    if cfg["tier"] != "INFERENCE" or not (ANN / f"{gse}_annotated.h5ad").exists():
        continue
    ad = sc.read_h5ad(ANN / f"{gse}_annotated.h5ad")
    lookup = {str(g).lower(): str(g) for g in ad.var_names}
    celltypes = sorted(c for c in ad.obs["celltype"].astype(str).unique() if c != "LowConfidence")
    samples = sorted(ad.obs["sample"].astype(str).unique())
    lr_presence[gse] = {}
    # masks per sample x cell type; gene-level means computed lazily on sparse columns
    masks = {}
    for smp in samples:
        for ct in celltypes:
            mask = ((ad.obs["sample"].astype(str) == smp) &
                    (ad.obs["celltype"].astype(str) == ct)).to_numpy()
            if int(mask.sum()) >= QC_THRESHOLDS["min_cells_per_celltype_per_sample"]:
                masks[(smp, ct)] = mask
    gene_idx = {g.lower(): i for i, g in enumerate(ad.var_names)}
    for L, R in LR_PANEL:
        li, ri = gene_idx.get(L.lower()), gene_idx.get(R.lower())
        lr_presence[gse][f"{L}-{R}"] = {"ligand_found": li is not None, "receptor_found": ri is not None}
        if li is None or ri is None:
            continue
        def col_values(idx):
            v = ad.X[:, idx]
            return np.asarray(v.todense() if hasattr(v, "todense") else v).ravel()
        lv = col_values(li)
        rv = col_values(ri)
        lmean = {k: float(lv[m].mean()) for k, m in masks.items()}
        rmean = {k: float(rv[m].mean()) for k, m in masks.items()}
        per_sample = {}
        axes = {}
        for smp in samples:
            best, best_axis = -np.inf, None
            for sender in celltypes:
                for receiver in celltypes:
                    k1, k2 = (smp, sender), (smp, receiver)
                    if k1 not in lmean or k2 not in rmean:
                        continue
                    val = lmean[k1] * rmean[k2]
                    if val > best:
                        best, best_axis = val, f"{sender}->{receiver}"
            per_sample[smp] = best if best_axis else np.nan
            axes[smp] = best_axis
        ctrl = [per_sample[s] for s in samples if DATASETS[gse]["groups"][s] == "Control"]
        exp = [per_sample[s] for s in samples if DATASETS[gse]["groups"][s] != "Control"]
        if len(ctrl) < 3 or len(exp) < 3 or np.isnan(ctrl + exp).any():
            continue
        d, p, nperm = exact_perm(np.array(exp), np.array(ctrl))
        # dominant axis = modal axis across group samples, not a single outlier sample
        from collections import Counter
        dom_ctrl = Counter([axes[s] for s in samples if DATASETS[gse]["groups"][s] == "Control"]).most_common(1)[0][0]
        dom_exp = Counter([axes[s] for s in samples if DATASETS[gse]["groups"][s] != "Control"]).most_common(1)[0][0]
        lr_rows.append({"dataset": gse, "ligand": L, "receptor": R,
                        "pair": f"{L}-{R}", "delta_exposed_minus_control": d,
                        "permutation_p_two_sided": p, "n_permutations": nperm,
                        "dominant_axis_control": dom_ctrl,
                        "dominant_axis_exposed": dom_exp})
lr_cols = ["dataset", "ligand", "receptor", "pair", "delta_exposed_minus_control",
           "permutation_p_two_sided", "n_permutations",
           "dominant_axis_control", "dominant_axis_exposed"]
lrdf = pd.DataFrame(lr_rows, columns=lr_cols)
if len(lrdf):
    # Family C: BH WITHIN dataset (Round 2, statistics reviewer M2)
    lrdf["FDR_BH_LR_family"] = np.nan
    for ds in lrdf.dataset.unique():
        m = (lrdf.dataset == ds).to_numpy()
        lrdf.loc[m, "FDR_BH_LR_family"] = bh(lrdf.loc[m, "permutation_p_two_sided"].to_numpy())
    lrdf["significant_FDR_005"] = lrdf["FDR_BH_LR_family"] < 0.05
lrdf.to_csv(STATS / "lr_panel_tests.csv", index=False, encoding="utf-8-sig")
(STATS / "lr_gene_presence.json").write_text(
    json.dumps(lr_presence, indent=2), encoding="utf-8")
print("WROTE separation/lr: sep rows", len(sep), "lr rows", len(lrdf))
