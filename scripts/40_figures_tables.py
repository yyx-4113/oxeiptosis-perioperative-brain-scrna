"""Generate all figures and manuscript-ready tables from result CSVs.

Every number plotted is read from an on-disk result file; nothing is recomputed
ad hoc here except derived display quantities.
"""
import json
import warnings
from pathlib import Path

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import scanpy as sc

from common import (ANN, QC, SCORE, STATS, FIG, AUDIT, DATASETS,
                    GENE_SETS, CELLTYPE_MARKERS)

warnings.filterwarnings("ignore")
plt.rcParams.update({"font.size": 8, "figure.dpi": 150, "savefig.dpi": 300,
                     "axes.spines.top": False, "axes.spines.right": False,
                     "pdf.fonttype": 42})
BLUE, RED, GREY = "#3B6EA5", "#B34242", "#777777"

qcdf = pd.read_csv(QC / "qc_per_sample.csv")
sdf = pd.read_csv(SCORE / "sample_scores.csv")
tdf = pd.read_csv(STATS / "sample_level_tests.csv")
sep = pd.read_csv(STATS / "separation_summary.csv")
lrdf = pd.read_csv(STATS / "lr_panel_tests.csv")
comp = pd.read_csv(STATS / "composition_tests.csv")
gp = json.loads((SCORE / "gene_presence.json").read_text(encoding="utf-8"))

# ---------------------------------------------------------------- Fig S1 QC
fig, axes = plt.subplots(1, 3, figsize=(9, 3))
for ax, col, lab in zip(axes, ["median_genes", "median_counts", "median_mt_pct"],
                        ["Median genes/cell", "Median UMIs/cell", "Median MT %"]):
    x = np.arange(len(qcdf))
    colors = [BLUE if g == "Control" else RED for g in qcdf["group"]]
    ax.bar(x, qcdf[col], color=colors)
    ax.set_xticks(x)
    ax.set_xticklabels([f"{d[-4:]}\n{s}" for d, s in zip(qcdf["dataset"], qcdf["sample"])],
                       rotation=90, fontsize=6)
    ax.set_ylabel(lab)
fig.suptitle("Fig S1. Per-sample QC after filtering", fontsize=10)
fig.tight_layout()
fig.savefig(FIG / "FigS1_QC.png", bbox_inches="tight")
plt.close(fig)

# --------------------------------------------------------- Fig 1 UMAP panels
for gse in [g for g, v in DATASETS.items() if v["tier"] == "INFERENCE"]:
    if not (ANN / f"{gse}_annotated.h5ad").exists():
        continue
    ad = sc.read_h5ad(ANN / f"{gse}_annotated.h5ad")
    fig, ax = plt.subplots(1, 2, figsize=(9, 4))
    sc.pl.umap(ad, color="celltype", ax=ax[0], show=False, legend_loc="on data",
               title=f"{gse}: cell types", palette=None)
    grp_palette = {g: (RED if g != "Control" else BLUE)
                   for g in ad.obs["group"].astype(str).unique()}
    sc.pl.umap(ad, color="group", ax=ax[1], show=False,
               title=f"{gse}: condition", palette=grp_palette)
    fig.suptitle(f"Fig 2. {gse} UMAP (post-QC n={ad.n_obs:,})", fontsize=10)
    fig.tight_layout()
    fig.savefig(FIG / f"Fig2_UMAP_{gse}.png", bbox_inches="tight")
    plt.close(fig)

# ------------------------------------------------- Fig 3 score distributions
for mod, title in [("ox_core", "Oxeiptosis core (AIFM1/KEAP1/PGAM5)"),
                   ("ferro_net", "Ferroptosis net (driver - suppressor)")]:
    fig, axes = plt.subplots(1, 2, figsize=(10, 4))
    for ax, gse in zip(axes, [g for g, v in DATASETS.items() if v["tier"] == "INFERENCE"]):
        sub = sdf[sdf.dataset == gse]
        cts = sorted(sub.celltype.unique())
        pos = np.arange(len(cts))
        for ct_i, ct in enumerate(cts):
            g = sub[(sub.celltype == ct)]
            for grp, col in [("Control", BLUE)]:
                vals = g.loc[g.group == grp, mod].dropna()
                if len(vals):
                    ax.scatter(np.full(len(vals), ct_i), vals, color=col, alpha=.75, s=22,
                               label="Control" if ct_i == 0 else None)
            vals = g.loc[g.group != "Control", mod].dropna()
            if len(vals):
                ax.scatter(np.full(len(vals), ct_i) + .18, vals, color=RED, alpha=.75, s=22,
                           label="Exposed" if ct_i == 0 else None)
            a = g.loc[g.group == "Control", mod].dropna().mean()
            b = g.loc[g.group != "Control", mod].dropna().mean()
            ax.plot([ct_i - .12, ct_i + .3], [a, b], color=GREY, lw=1, alpha=.7)
        ax.set_xticks(pos)
        ax.set_xticklabels(cts, rotation=45, ha="right", fontsize=6)
        ax.set_title(gse)
        ax.set_ylabel("Per-sample module score (mean of equal-depth cells)")
        ax.axhline(0, color="k", lw=.5, alpha=.3)
    ax.legend(fontsize=7)
    fig.suptitle(f"Fig 3. {title}", fontsize=10)
    fig.tight_layout()
    fig.savefig(FIG / f"Fig3_{mod}_sample_scores.png", bbox_inches="tight")
    plt.close(fig)

# ----------------------------------------------------------- Fig 4 forest plot
for mod in ["ox_core", "ferro_net", "nrf2_antag"]:
    sub = tdf[(tdf.module == mod) & tdf.cliffs_delta.notna()].copy()
    if sub.empty:
        continue
    labels = [f"{r.dataset[-4:]} | {r.celltype}" for r in sub.itertuples()]
    y = np.arange(len(sub))[::-1]
    fig, ax = plt.subplots(figsize=(7, max(3, .27 * len(sub))))
    ax.errorbar(sub.cliffs_delta, y,
                xerr=[sub.cliffs_delta - sub.cliffs_delta_CI95_lo,
                      sub.cliffs_delta_CI95_hi - sub.cliffs_delta],
                fmt="o", color=BLUE, ecolor=GREY, capsize=2, ms=4)
    ax.axvline(0, color="k", lw=.6)
    ax.set_yticks(y)
    ax.set_yticklabels(labels, fontsize=6)
    ax.set_xlabel("Cliff's delta (exposed - control), bootstrap 95% CI")
    ax.set_title(f"Fig 4. {mod} effect sizes (exploratory, n=3/group)")
    fig.tight_layout()
    fig.savefig(FIG / f"Fig4_forest_{mod}.png", bbox_inches="tight")
    plt.close(fig)

# ---------------------------------------------------------- Fig 5 separation
fig, ax = plt.subplots(figsize=(5, 4))
for gse, g in sep.groupby("dataset"):
    ax.scatter(g.ox_core, g.ferro_net, label=gse, s=30)
    for r in g.itertuples():
        ax.annotate(r.celltype.split("_")[0][:6], (r.ox_core, r.ferro_net), fontsize=5)
ax.axhline(0, color="k", lw=.5)
ax.axvline(0, color="k", lw=.5)
ax.set_xlabel("Oxeiptosis core delta (exposed - control)")
ax.set_ylabel("Ferroptosis net delta (exposed - control)")
ax.set_title("Fig 5. Separation of two cell-death signatures")
ax.legend(fontsize=7)
fig.tight_layout()
fig.savefig(FIG / "Fig5_separation.png", bbox_inches="tight")
plt.close(fig)

# ------------------------------------------------------------- Fig 6 LR panel
if len(lrdf):
    pivot = lrdf.pivot_table(index="pair", columns="dataset",
                             values="delta_exposed_minus_control")
    fig, ax = plt.subplots(figsize=(5, max(4, .22 * len(pivot))))
    im = ax.imshow(pivot.fillna(0), aspect="auto", cmap="RdBu_r", vmin=-1, vmax=1)
    ax.set_yticks(np.arange(len(pivot)))
    ax.set_yticklabels(pivot.index, fontsize=6)
    ax.set_xticks(np.arange(pivot.shape[1]))
    ax.set_xticklabels(pivot.columns, fontsize=7)
    fig.colorbar(im, ax=ax, label="Delta LR score", shrink=.6)
    ax.set_title("Fig 6. Targeted LR panel (sample-level)")
    fig.tight_layout()
    fig.savefig(FIG / "Fig6_LR_panel.png", bbox_inches="tight")
    plt.close(fig)

# ---------------------------------------------------------------- Fig 1 flow
fig, ax = plt.subplots(figsize=(8, 6))
ax.axis("off")
steps = [
    (0.95, "5 public GEO datasets audited"),
    (0.79, "INFERENCE (n>=3/condition): GSE267933, GSE222430\nLOCALIZATION_ONLY: GSE153598, GSE196239\nBULK_REFERENCE: GSE141242"),
    (0.61, "Per-sample QC -> Leiden clustering -> marker panels\nlow-margin clusters -> LowConfidence (excluded)"),
    (0.44, "LOCKED gene sets: ox_core (AIFM1/KEAP1/PGAM5)\nNRF2 antagonistic | ferroptosis driver/suppressor"),
    (0.27, "One sample-level value per cell type\nExact permutation (20 labelings; floor p=0.1)\nFamilies A (ox_core), B1 (modules), B2 (separation D), C (LR)"),
    (0.09, "All transcript results exploratory\nPROSPECTIVE: p-AIFM1(S116)/total IF + phosphatase\ncontrols + inhibitor rescue matrix"),
]
for y, txt in steps:
    ax.text(0.5, y, txt, ha="center", va="center", fontsize=8,
            bbox=dict(boxstyle="round,pad=0.5", facecolor="#EEF3F8", edgecolor=BLUE))
ax.set_xlim(0, 1); ax.set_ylim(0, 1)
ax.set_title("Study flow and pre-specified inference boundaries", fontsize=10, weight="bold")
fig.savefig(FIG / "Fig1_study_flow.png", bbox_inches="tight")
plt.close(fig)

# --------------------------------------------- Fig 5b: separation with CI
sep_t = pd.read_csv(STATS / "separation_tests.csv")
if len(sep_t):
    fig, ax = plt.subplots(figsize=(8, 4))
    labels = [f"{r.dataset[-4:]}\n{r.celltype[:12]}" for r in sep_t.itertuples()]
    x = np.arange(len(sep_t))
    y = sep_t["D_zscore_exposed_minus_control"].to_numpy(float)
    lo = sep_t["CI95_lo"].to_numpy(float)
    hi = sep_t["CI95_hi"].to_numpy(float)
    yerr = np.vstack([y - lo, hi - y])
    ax.errorbar(x, y, yerr=yerr, fmt="o", capsize=3, color=BLUE)
    ax.axhline(0, color=GREY, lw=0.8, ls="--")
    ax.set_xticks(x); ax.set_xticklabels(labels, rotation=45, ha="right", fontsize=6)
    ax.set_ylabel("Paired z-scored separation D (ox - ferro)")
    ax.set_title("Separation endpoint D: bootstrap 95% CI; exact-permutation p")
    fig.tight_layout()
    fig.savefig(FIG / "Fig5b_separation_CI.png", bbox_inches="tight")
    plt.close(fig)

# ------------------------------------------------------------- Table exports
tab1 = pd.DataFrame([
    {"dataset": k, "tier": v["tier"], "species": v["species"], "tissue": v["tissue"],
     "exposure": v["exposure"], "caveat": v["caveat"]}
    for k, v in DATASETS.items()])
tab1.to_csv(AUDIT / "Table1_dataset_audit.csv", index=False, encoding="utf-8-sig")
qcdf.to_csv(AUDIT / "Table2_sample_QC.csv", index=False, encoding="utf-8-sig")
cl = pd.concat([pd.read_csv(ANN / f"{g}_cluster_celltype.csv").assign(dataset=g)
                for g in [x for x, v in DATASETS.items() if v["tier"] == "INFERENCE"]
                if (ANN / f"{g}_cluster_celltype.csv").exists()])
cl.to_csv(AUDIT / "Table3_annotation.csv", index=False, encoding="utf-8-sig")
tdf.to_csv(AUDIT / "Table4_sample_level_tests.csv", index=False, encoding="utf-8-sig")
sep.to_csv(AUDIT / "Table5_separation.csv", index=False, encoding="utf-8-sig")
lrdf.to_csv(AUDIT / "Table6_LR_tests.csv", index=False, encoding="utf-8-sig")
comp.to_csv(AUDIT / "TableS1_composition.csv", index=False, encoding="utf-8-sig")
pres = pd.DataFrame([
    {"dataset": ds, "module": mod, "present": ";".join(info["present"]),
     "missing": ";".join(info["missing"])}
    for ds, mods in gp.items() for mod, info in mods.items()])
pres.to_csv(AUDIT / "TableS2_gene_presence.csv", index=False, encoding="utf-8-sig")
print("WROTE figures:", sorted(p.name for p in FIG.glob("*.png")))
print("WROTE tables:", sorted(p.name for p in AUDIT.glob("Table*.csv")))
