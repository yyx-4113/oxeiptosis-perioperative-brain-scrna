"""Per-sample QC -> integration -> Leiden clustering -> marker-based annotation.

Outputs (all under 03_annotation / 02_qc):
  <GSE>_annotated.h5ad, <GSE>_cluster_celltype.csv,
  qc_per_sample.csv, qc_summary.json
"""
import hashlib
import json
import warnings
from pathlib import Path

import numpy as np
import pandas as pd
import scanpy as sc

from common import (PROC, QC, ANN, DATASETS, QC_THRESHOLDS, CELLTYPE_MARKERS,
                    GTF_MT_PREFIX, STATS_CONFIG)

warnings.filterwarnings("ignore")
sc.settings.verbosity = 1
SEED = STATS_CONFIG["seed"]
T = QC_THRESHOLDS
qc_rows = []


def load_sample(gse, sample):
    d = PROC / gse / sample
    a = sc.read_10x_mtx(d, var_names="gene_symbols", make_unique=True, cache=False)
    # preserve GEO official symbols; mitochondrial detection is case-insensitive below
    a.obs["sample"] = sample
    a.obs["group"] = DATASETS[gse]["groups"][sample]
    a.obs["dataset"] = gse
    return a


for gse in [g for g, v in DATASETS.items() if v["tier"] == "INFERENCE"]:
    if not all((PROC / gse / s / "matrix.mtx.gz").exists() for s in DATASETS[gse]["groups"]):
        print("=== SKIP", gse, "(matrices not ready)", flush=True)
        continue
    print("=== processing", gse, flush=True)
    adatas = []
    for sample in sorted(DATASETS[gse]["groups"]):
        a = load_sample(gse, sample)
        n_raw = a.n_obs
        a.var["mt"] = pd.Index(a.var_names).str.lower().str.startswith("mt-")
        sc.pp.calculate_qc_metrics(a, qc_vars=("mt",),
                                   percent_top=None, log1p=False, inplace=True)
        before = a.n_obs
        a = a[a.obs["n_genes_by_counts"] >= T["min_genes_per_cell"]].copy()
        a = a[a.obs["total_counts"] >= T["min_counts_per_cell"]].copy()
        a = a[a.obs["pct_counts_mt"] < T["max_mt_pct"]].copy()
        sc.pp.filter_genes(a, min_cells=T["min_cells_per_gene"])
        # computational cap: deterministic random subsample after biological QC
        if a.n_obs > T["max_cells_per_sample"]:
            stable = int(hashlib.md5(f"{gse}|{sample}".encode()).hexdigest()[:8], 16)
            rng = np.random.default_rng(SEED + stable)
            idx = rng.choice(a.n_obs, T["max_cells_per_sample"], replace=False)
            a = a[idx].copy()
            capped = True
        else:
            capped = False
        qc_rows.append({
            "dataset": gse, "sample": sample, "group": a.obs["group"].iloc[0],
            "cells_raw": n_raw, "cells_after_qc": a.n_obs,
            "cells_removed": before - a.n_obs, "genes_after_qc": a.n_vars,
            "median_genes": float(a.obs["n_genes_by_counts"].median()),
            "median_counts": float(a.obs["total_counts"].median()),
            "median_mt_pct": float(a.obs["pct_counts_mt"].median()),
            "computational_cap_applied": capped,
        })
        print(f"  {sample}: {n_raw} -> {a.n_obs} cells, {a.n_vars} genes, capped={capped}", flush=True)
        adatas.append(a)

    # label must not be "sample": concat would overwrite the real GSM sample IDs
    ad = sc.concat(adatas, label="batch_id", index_unique="-", join="inner")
    ad.obs["group"] = ad.obs["group"].astype(str)
    ad.obs["dataset"] = gse
    ad.var_names_make_unique()
    sc.pp.normalize_total(ad, target_sum=1e4)
    sc.pp.log1p(ad)
    ad.raw = ad

    sc.pp.highly_variable_genes(ad, n_top_genes=T["n_hvg"], flavor="seurat", batch_key="sample")
    ad_hvg = ad[:, ad.var["highly_variable"]].copy()
    sc.pp.scale(ad_hvg, max_value=10)
    sc.tl.pca(ad_hvg, n_comps=T["n_pcs"], svd_solver="arpack", rng=SEED)
    ad.obsm["X_pca"] = ad_hvg.obsm["X_pca"]
    sc.pp.neighbors(ad, n_neighbors=15, use_rep="X_pca", rng=SEED)
    try:
        sc.tl.leiden(ad, resolution=T["resolution"], key_added="leiden",
                     flavor="leidenalg", rng=SEED)
    except Exception as exc:
        print("  leiden failed, trying igraph:", exc, flush=True)
        sc.tl.leiden(ad, resolution=T["resolution"], key_added="leiden",
                     flavor="igraph", n_iterations=2, rng=SEED)
    sc.tl.umap(ad, rng=SEED)

    # marker panel scores (locked a priori)
    lookup = {str(g).lower(): str(g) for g in ad.var_names}
    for ct, markers in CELLTYPE_MARKERS.items():
        present = [lookup[m.lower()] for m in markers if m.lower() in lookup]
        missing = [m for m in markers if m.lower() not in lookup]
        if not present:
            ad.obs["score__" + ct] = 0.0
        else:
            sc.tl.score_genes(ad, gene_list=present, score_name="score__" + ct,
                              ctrl_size=50, n_bins=25, rng=SEED)
        print(f"  panel {ct}: present={len(present)} missing={missing}", flush=True)

    # cluster-level assignment; Neuron_pan used only as validation, not a candidate label
    candidates = [ct for ct in CELLTYPE_MARKERS if ct != "Neuron_pan"]
    score_cols = ["score__" + ct for ct in candidates]
    cluster_mean = ad.obs.groupby("leiden", observed=True)[score_cols].mean()
    order = np.argsort(-cluster_mean.values, axis=1)
    top_ct = [candidates[i] for i in order[:, 0]]
    top1 = cluster_mean.values[np.arange(len(cluster_mean)), order[:, 0]]
    top2 = cluster_mean.values[np.arange(len(cluster_mean)), order[:, 1]]
    ct_map = dict(zip(cluster_mean.index.astype(str), top_ct))
    margin_map = dict(zip(cluster_mean.index.astype(str), top1 - top2))
    # Round 2 (biology reviewer M2): clusters with a near-zero top1-top2 margin
    # are arbitrary assignments.  Label them LowConfidence and exclude them from
    # cell-type-specific inference; they remain in UMAP/composition denominators.
    margin_ok = {k: bool(v >= T["min_annotation_margin"]) for k, v in margin_map.items()}
    ad.obs["celltype_raw"] = ad.obs["leiden"].astype(str).map(ct_map).astype(str)
    ad.obs["celltype"] = [raw if margin_ok.get(str(lei), False) else "LowConfidence"
                          for raw, lei in zip(ad.obs["celltype_raw"], ad.obs["leiden"].astype(str))]
    ad.obs["annotation_margin"] = [float(margin_map.get(str(c), np.nan))
                                   for c in ad.obs["leiden"].astype(str)]
    # Neuron_pan validation: fraction of cells in each neuronal candidate cluster
    # that also score positive on the pan-neuronal panel (above its own median).
    pan = ad.obs["score__Neuron_pan"].astype(float)
    pan_high = pan > pan.median()
    validated = {}
    for ct in candidates:
        m = ad.obs["celltype_raw"] == ct
        validated[ct] = float(pan_high[m].mean()) if m.any() else np.nan

    cl = cluster_mean.copy()
    cl["assigned_celltype"] = [ct if margin_ok.get(str(k), False) else "LowConfidence"
                               for ct, k in zip(top_ct, cluster_mean.index.astype(str))]
    cl["raw_top_celltype"] = top_ct
    cl["margin_top1_minus_top2"] = top1 - top2
    cl["n_cells"] = ad.obs.groupby("leiden", observed=True).size().values
    cl["neuron_pan_fraction_positive"] = [validated.get(ct, np.nan) for ct in top_ct]
    cl.to_csv(ANN / f"{gse}_cluster_celltype.csv", encoding="utf-8-sig")
    print("  LowConfidence clusters:",
          int((~np.array(list(margin_ok.values()))).sum()), "/", len(cluster_mean),
          "; cells:", int((ad.obs["celltype"] == "LowConfidence").sum()), flush=True)
    ad.write(ANN / f"{gse}_annotated.h5ad", compression="gzip")
    print("  WROTE", ANN / f"{gse}_annotated.h5ad", "shape", ad.shape, flush=True)

qcdf = pd.DataFrame(qc_rows)
qcdf.to_csv(QC / "qc_per_sample.csv", index=False, encoding="utf-8-sig")
summary = {
    "qc_thresholds": T,
    "dataset_totals": qcdf.groupby("dataset").agg(
        samples=("sample", "nunique"), cells_raw=("cells_raw", "sum"),
        cells_after_qc=("cells_after_qc", "sum")).to_dict("index"),
    "seed": SEED,
}
(QC / "qc_summary.json").write_text(json.dumps(summary, indent=2), encoding="utf-8")
print("WROTE", QC / "qc_per_sample.csv")
