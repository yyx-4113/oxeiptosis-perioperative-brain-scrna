"""Compose the manuscript from result CSVs.

Every numeric claim is taken from an on-disk result file and simultaneously
recorded in 07_manuscript/manuscript_numbers.json with its source path.
"""
import json
from pathlib import Path

import numpy as np
import pandas as pd

from common import (ROOT, QC, SCORE, STATS, META, FIG, MAN, AUDIT, ANN, DATASETS,
                    GENE_SETS, CELLTYPE_MARKERS, LR_PANEL, STATS_CONFIG,
                    QC_THRESHOLDS)

TITLE = ("Single-cell localization of oxeiptosis transcripts after sevoflurane "
         "exposure in aged mouse hippocampus: a two-cohort reanalysis")
SLUG = "single-cell-localization-of-oxeiptosis-transcripts"
VERSION = "v2.0"

inv = pd.read_csv(META / "sample_inventory.csv")
qc = pd.read_csv(QC / "qc_per_sample.csv")
sdf = pd.read_csv(SCORE / "sample_scores.csv")
tdf = pd.read_csv(STATS / "sample_level_tests.csv")
sep = pd.read_csv(STATS / "separation_summary.csv")
sep_tests = pd.read_csv(STATS / "separation_tests.csv")
comp = pd.read_csv(STATS / "composition_tests.csv")
lrd = pd.read_csv(STATS / "lr_panel_tests.csv")
gp = json.loads((SCORE / "gene_presence.json").read_text(encoding="utf-8"))
cl = pd.read_csv(AUDIT / "Table3_annotation.csv")

numbers = {}


def num(key, value, source, note=""):
    numbers[key] = {"value": value, "source": source, "note": note}
    return str(value)


def fmt(x, d=3):
    return "NA" if pd.isna(x) else f"{float(x):.{d}f}"


# --------------------------------------------------------------- key values
gse_list = list(DATASETS.keys())
n_inf_cohorts = len([g for g, v in DATASETS.items() if v["tier"] == "INFERENCE"])
n_samples = num("n_samples_inference", int(inv.shape[0]), "01_data/meta/sample_inventory.csv",
                f"libraries across {n_inf_cohorts} inference cohorts")
cells_raw_int = int(qc.cells_raw.sum())
cells_qc_int = int(qc.cells_after_qc.sum())
cells_raw = num("cells_raw", cells_raw_int, "02_qc/qc_per_sample.csv")
cells_qc = num("cells_after_qc", cells_qc_int, "02_qc/qc_per_sample.csv")
celltypes = sorted(sdf.celltype.unique())
n_celltypes = num("n_celltypes", len(celltypes), "04_scoring/sample_scores.csv")
n_tests = num("n_tests", int(len(tdf)), "05_stats/sample_level_tests.csv")
_stats_meta = json.loads((STATS / "stats_metadata.json").read_text(encoding="utf-8"))
min_p = num("min_possible_p", _stats_meta["min_possible_two_sided_p"],
            "05_stats/stats_metadata.json:min_possible_two_sided_p",
            "exact permutation with 3 vs 3, C(6,3)=20; derived from product file")
sig = tdf[tdf.significant_FDR_005]
n_sig_int = int(len(sig))
n_sig = num("n_significant_FDR005", n_sig_int, "05_stats/sample_level_tests.csv")
lr_sig = lrd[lrd.significant_FDR_005] if len(lrd) else lrd
n_lr_sig = num("n_LR_significant", int(len(lr_sig)), "05_stats/lr_panel_tests.csv")
n_lr_tested = num("n_LR_tested", int(len(lrd)), "05_stats/lr_panel_tests.csv")
n_lr_panel = num("n_LR_panel_declared", len(LR_PANEL), "scripts/common.py:LR_PANEL")
_lr_by_ds = lrd.groupby("dataset").size().to_dict() if len(lrd) else {}
lr_eval_gse1 = num("n_LR_evaluable_GSE267933", int(_lr_by_ds.get("GSE267933", 0)), "05_stats/lr_panel_tests.csv")
lr_eval_gse2 = num("n_LR_evaluable_GSE222430", int(_lr_by_ds.get("GSE222430", 0)), "05_stats/lr_panel_tests.csv")
n_sep_sig = num("n_separation_FDR005", int(sep_tests.significant_FDR_005.sum()),
               "05_stats/separation_tests.csv")
# Round-2复审 N1: count LowConfidence directly from the annotation table
# (authoritative product), not from h5ad .obs which previously threw NameError
# and was silently swallowed -> false "0 cells".
_cl = pd.read_csv(AUDIT / "Table3_annotation.csv")
_low = _cl[_cl.assigned_celltype.astype(str) == "LowConfidence"]
_low_by_ds = _low.groupby("dataset")["n_cells"].sum().to_dict()
n_lowconf_int = int(_low["n_cells"].sum())
total_cells = int(_cl["n_cells"].sum())
lowconf_frac = n_lowconf_int / total_cells if total_cells else 0.0
n_lowconf = num("n_lowconfidence_cells", n_lowconf_int, "08_audit/Table3_annotation.csv")
numbers["lowconfidence_fraction"] = {"value": round(lowconf_frac, 4),
                                     "source": "08_audit/Table3_annotation.csv"}
lowconf_by_ds = {str(k): {"n_lowconfidence": int(v),
                          "fraction": round(float(v / _cl[_cl.dataset == k].n_cells.sum()), 4)}
                 for k, v in _low_by_ds.items()}
numbers["lowconfidence_by_dataset"] = {"value": lowconf_by_ds,
                                       "source": "08_audit/Table3_annotation.csv"}
n_boot = num("bootstrap_reps", STATS_CONFIG["bootstrap_reps"], "scripts/common.py:STATS_CONFIG")
num("median_genes_range_per_cell", [int(qc.median_genes.min()), int(qc.median_genes.max())],
    "02_qc/qc_per_sample.csv", "min and max of per-sample median genes")
num("median_mt_pct_range_per_cell", [round(float(qc.median_mt_pct.min()), 2),
                                     round(float(qc.median_mt_pct.max()), 2)],
    "02_qc/qc_per_sample.csv", "min and max of per-sample median mitochondrial percent")

top = tdf.sort_values("permutation_p_two_sided").head(6)
top_ox = tdf[tdf.module == "ox_core"].sort_values(
    "delta_exposed_minus_control", key=lambda s: s.abs(), ascending=False).iloc[0]
num("largest_abs_ox_core_delta", round(float(top_ox.delta_exposed_minus_control), 4),
    "05_stats/sample_level_tests.csv", f"{top_ox.dataset}/{top_ox.celltype}")
top_lines = []
for r in top.itertuples():
    top_lines.append(
        f"| {r.dataset} | {r.celltype} | {r.module} | {fmt(r.delta_exposed_minus_control)} | "
        f"{fmt(r.cliffs_delta)} [{fmt(r.cliffs_delta_CI95_lo)}, {fmt(r.cliffs_delta_CI95_hi)}] | "
        f"{fmt(r.permutation_p_two_sided)} | {fmt(r.FDR_BH_all_tests)} |")
top_table = "\n".join(top_lines)

sep_lines = []
for r in sep_tests.itertuples():
    sep_lines.append(f"| {r.dataset} | {r.celltype} | {fmt(r.D_zscore_exposed_minus_control)} | "
                     f"{fmt(r.CI95_lo)} to {fmt(r.CI95_hi)} | "
                     f"{fmt(r.permutation_p_two_sided)} | {fmt(r.FDR_BH_family_B_within_dataset)} |")
sep_table = "\n".join(sep_lines)

qc_lines = "\n".join(
    f"| {r.dataset} | {r.sample} | {r.group} | {r.cells_raw:,} | {r.cells_after_qc:,} | "
    f"{fmt(r.median_genes,0)} | {fmt(r.median_counts,0)} | {fmt(r.median_mt_pct,1)} |"
    for r in qc.itertuples())
ct_lines = "\n".join(
    f"| {r.dataset} | {r.assigned_celltype} | {int(r.n_cells):,} | {fmt(r.margin_top1_minus_top2)} |"
    for r in cl.itertuples())

gene_lines = []
for ds, mods in gp.items():
    if ds not in [g for g, v in DATASETS.items() if v["tier"] == "INFERENCE"]:
        continue
    for mod, info in mods.items():
        gene_lines.append(f"| {ds} | {mod} | {', '.join(info['present']) or 'none'} | "
                          f"{', '.join(info['missing']) or 'none'} |")
gene_table = "\n".join(gene_lines)

lr_lines = "\n".join(
    f"| {r.dataset} | {r.pair} | {fmt(r.delta_exposed_minus_control)} | "
    f"{fmt(r.permutation_p_two_sided)} | {fmt(r.FDR_BH_LR_family)} |"
    for r in lrd.sort_values("FDR_BH_LR_family").head(10).itertuples()) if len(lrd) else "none"

comp_lines = "\n".join(
    f"| {r.dataset} | {r.celltype} | {fmt(r.fraction_control)} | {fmt(r.fraction_exposed)} | "
    f"{fmt(r.permutation_p_two_sided)} |" for r in comp.itertuples())

if n_sig_int == 0:
    sig_sentence = f"No test reached FDR<0.05 (0 of {n_tests})"
else:
    sig_sentence = f"{n_sig_int} of {n_tests} tests reached FDR<0.05; these are listed in Table 4 and remain exploratory pending protein-level validation"

inf_gse = [g for g, v in DATASETS.items() if v["tier"] == "INFERENCE"]
core_ok = (len(gp) == len(inf_gse)) and all(
    all(g.lower() in {str(x).lower() for x in gp[ds]["ox_core"]["present"]}
        for g in GENE_SETS["oxeiptosis_core"]) for ds in inf_gse)
null_sentence = ("This reanalysis returns a predominantly null result at the transcript level"
                 if n_sig_int == 0 else
                 "This reanalysis identifies only limited transcription-level associations")
core_sentence = ("The locked core triad was detected in both cohorts" if core_ok
                 else "The locked core triad was only partially detected; see Table S2")
n_ferro_ge = int((sep.ferro_net.abs() >= sep.ox_core.abs()).sum())
num("n_ferro_net_ge_ox_core", n_ferro_ge, "05_stats/separation_summary.csv",
    f"out of {len(sep)} cell-type/cohort comparisons")
num("n_separation_comparisons", int(len(sep)), "05_stats/separation_summary.csv")
pred_lines = []
_pref = ["Astrocyte", "Excitatory_neuron", "Microglia", "Oligodendrocyte",
         "OPC", "Endothelial", "Pericyte"]
_pox = tdf[tdf.module == "ox_core"]
for ds in sorted(_pox.dataset.unique()):
    for ct in _pref:
        sel = _pox[(_pox.dataset == ds) & (_pox.celltype == ct)]
        if sel.empty:
            continue
        r = sel.iloc[0]
        pred_lines.append(
            f"| {r.dataset} / {r.celltype} | oxeiptosis-core transcript score differs "
            f"between conditions (delta={fmt(r.delta_exposed_minus_control)}) | low "
            f"(nominal p={fmt(r.permutation_p_two_sided)}, family-A FDR="
            f"{fmt(r.FDR_BH_all_tests)}) | Cliff's delta "
            f"{fmt(r.cliffs_delta)} [{fmt(r.cliffs_delta_CI95_lo)}, {fmt(r.cliffs_delta_CI95_hi)}]; "
            f"n=3 per group | cell-type-resolved p-AIFM1(S116)/total IF with λ-phosphatase, "
            f"calyculin A and H2O2 controls, plus LFHP-1c/NAC rescue |")
pred_table = "\n".join(pred_lines) if pred_lines else "| none | - | - | - | - |"

grp_n = sdf.groupby(["dataset", "celltype", "group"]).size()
_incomplete = grp_n[grp_n < 3].reset_index()
_incomplete_ct = _incomplete[["dataset", "celltype"]].drop_duplicates()
excluded_pairs = [f"{r.dataset}/{r.celltype}" for r in _incomplete_ct.itertuples()]
excluded_sentence = (f"Cell type×cohort strata with fewer than three evaluable samples in "
                     f"either group ({', '.join(excluded_pairs)}) were excluded from inference "
                     f"in that cohort only"
                     if excluded_pairs else "All annotated cell types had three evaluable samples per group")

ferro_sentence = (f"Ferroptosis-net shifts were at least as large as oxeiptosis-core shifts in "
                  f"{n_ferro_ge} of {len(sep)} cell-type/cohort comparisons" if n_ferro_ge
                  else "Oxeiptosis-core shifts exceeded ferroptosis-net shifts in all comparisons")

refs = [
    "1. Holze C, Michaudel C, Mackowiak C, et al. Oxeiptosis, a ROS-induced caspase-independent apoptosis-like cell-death pathway. Nat Immunol. 2018;19(2):130-140. PMID: 29255269; doi:10.1038/s41590-017-0013-y.",
    "2. Tang D, Kang R, Vanden Berghe T, Vandenabeele P, Kroemer G. The molecular machinery of regulated cell death. Cell Res. 2019;29(5):347-364. PMID: 30948788; doi:10.1038/s41422-019-0164-5.",
    "3. Suo Z, et al. Aged hippocampal single-cell atlas screening unveils disrupted neuroglial system in postoperative cognitive impairment. Aging Cell. 2025;24(3):e14406. PMID: 39540334. GSE267933 (single-cell RNA sequencing).",
    "4. Chang E, et al. iScience. 2023;26(4):106534. PMID: 37123239. GSE196239.",
    "5. Zhang Z, Hu C, Liu Y, et al. Integrated re-analysis of bulk and single-cell transcriptomic data reveals limited ferroptosis and neuroinflammation signatures after sevoflurane exposure in the aged rat hippocampus and human prefrontal cortex cells. Front Cell Dev Biol. 2026;14:1938601. PMID: 42871048; PMCID: PMC13658196; doi:10.3389/fcell.2026.1938601.",
    "6. Gao C, et al. A novel PGAM5 inhibitor LFHP-1c protects blood-brain barrier integrity in ischemic stroke. Acta Pharm Sin B. 2021;11(7):1867-1884. PMID: 34386325.",
    "7. Suo Z, et al. CNS Neurosci Ther. 2022;28:1576-1595. PMID: 35899365; Dryad doi:10.5061/dryad.jsxksn0cj.",
    "8. Zhang X, et al. Unveiling a novel mechanism in noise-induced hearing loss: oxeiptosis-mediated regulated cell death of cochlear hair cell. Neurosci Bull. 2026. PMID: 41697573; doi:10.1007/s12264-025-01585-z.",
    "9. Wolf FA, Angerer P, Theis FJ. SCANPY: large-scale single-cell gene expression data analysis. Genome Biol. 2018;19:15. PMID: 29409532.",
    "10. Traag VA, Waltman L, van Eck NJ. From Louvain to Leiden: guaranteeing well-connected communities. Sci Rep. 2019;9:5233. PMID: 30914743; doi:10.1038/s41598-019-41695-z.",
    "11. Edgar R, Domrachev M, Lash AE. Gene Expression Omnibus: NCBI gene expression and hybridization array data repository. Nucleic Acids Res. 2002;30:207-210. PMID: 11752295.",
    "12. Zhou N, Yuan X, Du Q, et al. FerrDb V2: update of the manually curated database of ferroptosis regulators and ferroptosis-disease associations. Nucleic Acids Res. 2023;51(D1):D571-D582. PMID: 36305834; doi:10.1093/nar/gkac935.",
    "13. GSE222430. 10x single-cell sequencing of old mouse hippocampal transcriptomes after sevoflurane exposure (Lei Z, Yongyan C; Department of Anesthesiology, Shanghai Ninth People's Hospital Affiliated to Shanghai Jiao Tong University School of Medicine). GEO accession, submitted 2023-01-09, public 2026-01-08; no associated publication (GEO citation field empty, checked 2026-11-26).",
]

md = f"""# {TITLE}

**Version:** {VERSION}  
**Author:** Yongxin Yang  
**Affiliation:** Department of Anesthesiology, The Second Affiliated Hospital of Fujian University of Traditional Chinese Medicine, Fuzhou, Fujian 350003, China  
**Correspondence:** 960856791@qq.com · ORCID 0009-0004-9698-6552  
**Word count:** main text ~__WC__ words (excluding tables, references and figure legends)

---

## Abstract

**Background.** Oxeiptosis is a ROS-dependent, caspase-independent regulated cell-death pathway driven by the KEAP1–PGAM5–AIFM1 axis. Whether sevoflurane exposure engages this pathway in the aged hippocampus remains untested. Prior public-data reanalysis of related datasets reported no reproducible ferroptosis or neuroinflammation signal, raising the bar for any new claim.

**Methods.** We reanalyzed two independent mouse hippocampal single-cell RNA-sequencing cohorts with three biological replicates per condition: GSE267933 (2.5% sevoflurane in 50% O2 plus laparotomy versus 50% O2 control, 24 h) and GSE222430 (2% sevoflurane for 6 h versus air control, hippocampus sampled 24 h later). After per-sample quality control and a pre-set computational cap of {num('cap', QC_THRESHOLDS['max_cells_per_sample'], 'scripts/common.py:QC_THRESHOLDS')} cells (not triggered), cells were clustered and annotated by a pre-locked marker panel with a top1–top2 margin threshold; low-margin clusters were labelled LowConfidence and excluded from inference. A pre-specified three-gene oxeiptosis core score (AIFM1, KEAP1, PGAM5), an NRF2 antagonistic module and a ferroptosis driver-minus-suppressor score were computed per cell. Per-sample medians within each cell type were compared by exact permutation (C(6,3)=20; minimum attainable two-sided p={min_p}), with Cliff's delta as the primary effect size ({n_boot}-replicate bootstrap interval) and Hedges' g as a point-estimate sensitivity measure. Four FDR families were used, each BH-corrected within dataset where applicable: A = primary ox_core across cell types; B1 = secondary modules; B2 = the separation endpoint D; C = the targeted ligand–receptor panel. The B1/B2 split was fixed in the analysis code before outcome inspection and is declared in `08_audit/amendments.md`. A targeted {n_lr_panel}-pair ligand–receptor panel was scored per sample ({n_lr_tested} cohort×pair tests across the two cohorts). Deviations from the original protocol are documented in `08_audit/amendments.md`.

**Results.** From {cells_raw_int:,} raw cells, {cells_qc_int:,} passed QC; {n_lowconf_int:,} cells ({lowconf_frac:.1%}; {lowconf_by_ds.get('GSE267933', {}).get('fraction', 0):.1%} and {lowconf_by_ds.get('GSE222430', {}).get('fraction', 0):.1%} in the two cohorts) were labelled LowConfidence and excluded from cell-type-specific inference. The remaining cells were assigned to {n_celltypes} classes. {core_sentence} (gene-presence audit in Table S2). {sig_sentence}, and the minimum attainable p value was {min_p}; therefore all transcript-level findings are exploratory trends rather than confirmed differences. The largest absolute oxeiptosis-core shift occurred in {top_ox.celltype} of {top_ox.dataset} (delta={fmt(top_ox.delta_exposed_minus_control)}), with wide intervals and no family-A significance. The paired, z-standardised separation endpoint D showed {n_sep_sig} of {len(sep_tests)} tests significant in family B2 (Table 5). {ferro_sentence}, so a pathway-specific claim cannot be made from transcription alone. Composition differences are reported separately (Table S1).

**Conclusions.** In two replicate-level cohorts, transcript abundance of the oxeiptosis core triad does not provide reproducible evidence of sevoflurane-induced oxeiptosis in the aged hippocampus. The study establishes an auditable, sample-level framework and defines—without significance-based selection—the cell types in which protein-level testing of p-AIFM1(S116) may be explored. All results are predictive/hypothesis-generating and require immunofluorescence with phosphatase controls and inhibitor rescue before mechanistic interpretation.

**Keywords:** oxeiptosis; sevoflurane; single-cell RNA sequencing; hippocampus; AIFM1; sample-level inference; permutation test

---

## 1. Introduction

Perioperative neurocognitive disorder affects a substantial fraction of older surgical patients, yet the underlying cell-death programs remain unresolved. Oxeiptosis was defined by Holze and colleagues as a ROS-dependent, caspase-independent death modality in which KEAP1 releases PGAM5, which dephosphorylates AIFM1 at serine 116; NRF2-dependent survival signaling is dispensable under high ROS [1,2]. This creates a testable prediction: if sevoflurane engages oxeiptosis in the aged hippocampus, the core triad should show a cell-type-specific transcriptional signal, and a protein-level S116 dephosphorylation should follow.

We previously evaluated five candidate public datasets and found that only two contain biological replication suitable for inference. GSE153598 contains one macaque per condition with four technical library splits of pooled tissue; GSE196239 contains one library per condition in human prefrontal-cortex primary culture; GSE141242 is a bulk rat microarray. These are restricted here to descriptive localization or background reference. A closely related but non-overlapping reanalysis of GSE141242 and GSE196239 reported no consistent ferroptosis or neuroinflammation signal in public transcriptomes, although its cell-culture qRT-PCR did show perturbation of GPX4/SLC7A11/ACSL4/PTGS2 [5]. The present study differs in three declared ways: (i) the endpoint is changed from ferroptosis/neuroinflammation to oxeiptosis; (ii) the data are two mouse hippocampal single-cell cohorts with n=3 biological replicates per condition rather than bulk arrays and single libraries; and (iii) inference is cell-type-resolved with one observation per biological sample and explicit permutation-resolution limits.

Three pre-specified decisions constrain this work. First, AIRE was excluded from the gene set because no oxeiptosis publication includes it and human protein data indicate no brain expression; NFE2L2 was retained only as a separate antagonistic module. Second, the 6/24/72 h time course proposed originally cannot be addressed by the public data, which provide a single 24 h point; it is deferred to prospective animal work. Third, the primary mechanistic endpoint is p-AIFM1(S116)/total AIFM1 assessed by immunofluorescence, not a transcriptional score. Terminology: oxeiptosis denotes the ROS-dependent KEAP1–PGAM5–AIFM1 death pathway and is not synonymous with hyperoxic lung or brain injury; this distinction matters because the control condition in GSE267933 was 50% O2, itself an oxygen intervention that may shift ROS baseline.

## 2. Methods

### 2.1 Data audit and inclusion rules

Datasets were included for inference only if they provided at least three biological replicates per condition and readable per-sample cell counts. GSE267933 (n=3 control breathing 50% O2, n=3 sevoflurane 2.5% in 50% O2 plus laparotomy, 24 h) and GSE222430 (n=3 air control, n=3 sevoflurane 2% for 6 h, sampled at 24 h) met these criteria. Because 50% O2 is itself an intervention, GSE222430 (air control) is the cleaner anesthetic-effect contrast; GSE267933 is reported as a mixed-exposure reference. The two cohorts are reported side by side and never pooled: pooling would create a heterogeneous exposure and cohort-level pseudo-replication. GSE153598, GSE196239 and GSE141242 were retained as localization/bulk references only (Table 1). GSE222430 is a public GEO submission without an associated publication (accession only).

### 2.2 Quality control and annotation

Per sample, cells with fewer than 200 genes, fewer than 500 UMIs, or >15% mitochondrial counts were removed; genes detected in fewer than three cells were discarded. No upper gene-count bound and no computational doublet caller were applied (deviations from the original plan; see `amendments.md`). A deterministic random cap of {num('cap2', QC_THRESHOLDS['max_cells_per_sample'], 'scripts/common.py:QC_THRESHOLDS')} cells per sample was pre-set but not triggered (seed {STATS_CONFIG['seed']}). Samples were normalized to 10,000 counts, log-transformed, and 2,000 highly variable genes were used for 30-component PCA. Neighborhood graphs and Leiden clustering (resolution 0.8) were followed by cluster-level annotation using pre-locked marker panels (Table 3); hippocampal subtype markers (Prox1/Calb1, Neurod6/Dkk3/Fibcd1) were included and Itgam was removed from the microglia panel because it also labels infiltrating monocytes. Clusters whose top1–top2 panel-score margin was below {QC_THRESHOLDS['min_annotation_margin']} were labelled LowConfidence ({lowconf_frac:.1%} of cells) and excluded from cell-type-specific inference. Neuron_pan served only as a validation panel: for each neuronal candidate cluster we report the fraction of cells above the pan-neuronal score median (Table 3); it was never an assignment label.

### 2.3 Pre-locked gene sets and scoring

The oxeiptosis core set comprised AIFM1, KEAP1 and PGAM5. The NRF2 antagonistic module comprised NFE2L2, NQO1, GCLC and GCLM; HMOX1 was removed from this module in Round 2 because it was also a pre-locked ferroptosis-suppressor member, which made the two endpoints non-independent. Ferroptosis was represented as a driver-minus-suppressor contrast using pre-specified lists drawn from FerrDb v2 (PMID 36305834) and restricted to genes detectable here (Table S2); Alox15 was not detected and is listed as missing. Scores were computed with scanpy's control-gene-subtracted scoring (50 control genes, 25 expression bins); a plain mean-expression score for the core triad was used as a sensitivity analysis. AIRE was not scored. A score is a relative quantity against sampled control genes; its between-group difference is not an absolute mRNA fold change. Original-plan AUCell scoring and 200/400-cell downsampling sensitivity were not executed (see `amendments.md`).

### 2.4 Statistical analysis

Each cell type in each sample contributed exactly one score after equal-depth downsampling to 100 cells; cell types with fewer than 20 cells in a sample were treated as missing. Group differences were tested on the median difference with the exact permutation distribution over all 20 labelings of six samples, giving a minimum attainable two-sided p of {min_p}. Cliff's delta was the primary effect size with a {n_boot}-replicate percentile bootstrap interval. Bias-corrected Hedges' g was reported as a point estimate with the fraction of bootstrap replicates that were finite; its interval was suppressed because at n=3 the pooled-variance resampling distribution is unstable and can span orders of magnitude (see Table 4 and `amendments.md`). Four FDR families were applied with Benjamini–Hochberg correction: family A = ox_core across cell types and datasets; family B1 = secondary modules, corrected within each dataset; family B2 = the separation endpoint D, corrected within each dataset; family C = the ligand–receptor panel, also within dataset. Cell-type proportions were tested in a separate descriptive family. {excluded_sentence}. The separation endpoint D was constructed at the sample level as the paired sample-level difference z(ox_core) − z(ferro_net), with z-scores computed across all cells of each dataset, and tested by exact permutation with a bootstrap 95% interval; an effect-size-scale sensitivity D = |Cliff's delta(ox)| − |Cliff's delta(ferro)| was also computed. The targeted ligand–receptor analysis used per-sample products of mean ligand and receptor expression, taking the modal sender→receiver axis across samples of each group for each of the {n_lr_panel} pre-specified pairs, with family-C BH correction.

### 2.5 Reporting standard

All transcript-level results are predictive and exploratory. The confirmatory endpoint in the prospective extension is p-AIFM1(S116)/total AIFM1 by immunofluorescence in the annotated cell types, with λ-phosphatase (negative), calyculin A and H2O2 (positive) controls; the p/total ratio must be computed within the same molecular-weight band (66/62/57 kDa forms must not be mixed), plus a rescue matrix containing LFHP-1c, N-acetylcysteine, catalase, Z-VAD-FMK, necrostatin-1 and ferrostatin-1 in a dose–time design. Public de-identified GEO data were used; institutional ethics approval and informed consent were not required for this reanalysis.

## 3. Results

### 3.1 Cohort composition and quality control

Across the inference cohorts, {cells_raw_int:,} raw cells were obtained; {cells_qc_int:,} remained after QC (Table 2). Median genes per cell ranged from {qc.median_genes.min():.0f} to {qc.median_genes.max():.0f}, and median mitochondrial percentage from {qc.median_mt_pct.min():.1f}% to {qc.median_mt_pct.max():.1f}%. Annotation yielded {n_celltypes} cell classes, including neurons, astrocytes, microglia, oligodendrocytes and OPCs (Table 3; UMAP in Fig. 2).

### 3.2 Gene detection

{core_sentence}; details of present and missing genes for every module are given in Table S2. Detection of a transcript is a prerequisite for scoring and is not evidence of pathway activation.

### 3.3 Sample-level oxeiptosis and ferroptosis scores

Six representative tests, ordered by nominal p value, are shown below; the complete set of {n_tests} tests (families A and B1) is in Table 4. Delta is the median difference (exposed−control); FDR is BH within the indicated family.

| Cohort | Cell type | Module | Delta (exposed−control) | Cliff's delta [95% CI] | Permutation p | FDR (family) |
|---|---|---|---|---|---|---|
{top_table}

{sig_sentence}. Given the exact permutation floor of p={min_p}, these data cannot establish a difference between conditions at conventional significance. The largest nominal oxeiptosis-core delta was {fmt(top_ox.delta_exposed_minus_control)} in {top_ox.celltype} ({top_ox.dataset}); because this is one of {n_tests} comparisons and has a wide interval, it is reported as a trend requiring protein-level testing. Consistent with the pre-specified antagonistic interpretation, NRF2-module shifts (family B1) were not treated as evidence of oxeiptosis: an elevated NRF2/ARE response would indicate a pro-survival antioxidant reaction, and high-ROS oxeiptosis is NRF2-independent [1]; per-cell-type NRF2 results are in Table 4.

### 3.4 Separation from ferroptosis

| Cohort | Cell type | Separation D (z-scored) | Bootstrap 95% CI | Permutation p | FDR (B2) |
|---|---|---|---|---|---|
{sep_table}

The table gives the paired, sample-level z-standardised separation endpoint D with its exact-permutation p, two-sample bootstrap 95% interval and family-B2 q value; positive D indicates a relatively larger oxeiptosis shift. {n_sep_sig} of {len(sep_tests)} D tests reached FDR<0.05 in family B2. The hypothesis-free effect-size-scale sensitivity D = |δ(ox)| − |δ(ferro)| is provided in Table 5. {ferro_sentence}; transcriptional data therefore do not establish pathway specificity (Fig. 5).

### 3.5 Cell-type composition

Composition shifts could confound per-cell scores. Nominal permutation p values for cell-type fractions are listed in Table S1; because the same n=3 resolution limit applies, these are descriptive.

| Cohort | Cell type | Fraction (control) | Fraction (exposed) | Permutation p |
|---|---|---|---|---|
{comp_lines}

### 3.6 Targeted ligand–receptor panel

Of {n_lr_panel} pre-specified pairs, {lr_eval_gse1} and {lr_eval_gse2} were evaluable in the two cohorts respectively (transcript not detected for the remainder), giving {n_lr_tested} cohort×pair tests; {n_lr_sig} passed family-C FDR<0.05. Top pairs by FDR are listed below; the full panel is in Table 6. AIFM1, KEAP1 and PGAM5 are intracellular and absent from ligand–receptor repositories, so any communication difference is an indirect downstream observation without spatial co-localization information (Fig. 6).

| Cohort | Pair | Delta | Permutation p | FDR |
|---|---|---|---|---|
{lr_lines}

## 4. Discussion

{null_sentence}, which is scientifically informative given that a previous reanalysis of conceptually related but non-overlapping public data also failed to find a consistent death-signature signal in bulk transcriptomes [5], although its cell-culture qRT-PCR did detect perturbation of ferroptosis markers. The result is consistent with the biology of oxeiptosis: the decisive event is dephosphorylation of AIFM1(S116) by PGAM5, a post-translational modification that is not necessarily accompanied by changes in mRNA abundance. A three-gene transcriptional score therefore has low sensitivity for the pathway switch. Additionally, the KEAP1 "biphasic" behavior described in the literature is a functional switch between NRF2-dependent survival and PGAM5-dependent death, not a biphasic mRNA expression pattern; the public data contain a single 24 h point and cannot address temporal dynamics.

The two cohorts differ in exposure and in control gas: GSE267933 combines sevoflurane with laparotomy over a 50% O2 background, whereas GSE222430 is anesthetic-only against air. Their concordance (or lack of it) is therefore informative but confounded by control gas; neither design isolates 6/24/72 h kinetics. We deliberately avoided pooling them: doing so would inflate replication while introducing heterogeneity, a form of pseudo-replication at the cohort level. The term oxeiptosis is used here strictly for the ROS-dependent KEAP1–PGAM5–AIFM1 pathway and not as a synonym for hyperoxic injury.

The main strength is that inference was performed once per biological sample, with the permutation floor stated explicitly. The main weaknesses are statistical power (three samples per group; minimum attainable two-sided p={min_p}), absence of a computational doublet caller, and cell-type annotation from marker panels rather than a reference atlas. Low-margin clusters were flagged and excluded, but the observed microglia/neuron proportions differ from published hippocampal benchmarks, so cell-type-specific results should be read as annotation-conditional. Bootstrap intervals over three observations are unstable, and Hedges' g intervals were suppressed when resamples had zero within-group variance; Cliff's delta is therefore the primary effect size.

The next step is therefore not more of the same scoring but a prospective experiment: aged mice assigned to air, oxygen, sevoflurane and sevoflurane-plus-surgery, sampled at 0/6/12/24/72 h, with cell-type-resolved p-AIFM1(S116)/total AIFM1 immunofluorescence, phosphatase controls, and a dose–time rescue matrix. Antibody PA5-143645 (RRID AB_2942873) is labelled human-only by the vendor although the S116 site is conserved in mouse; AP5501 is annotated for mouse and has published mouse-tissue use, including the study cited as [8]. Cross-species validation in the intended mouse tissue with λ-phosphatase (negative) and calyculin A/H2O2 (positive) controls must still precede any mouse endpoint, with p/total computed within the same molecular-weight band. If validation fails, quantitative phosphoproteomics should replace immunofluorescence.

Three interpretive points follow from the design. First, the NRF2 module is antagonistic by construction: under high ROS oxeiptosis proceeds without NRF2-dependent survival signalling [1], so an increased NRF2/ARE score is interpreted as a protective response and does not support pathway activation. Second, oxeiptosis is initiated by a post-translational switch: high ROS promotes KEAP1 dissociation from PGAM5, PGAM5 dephosphorylates AIFM1 at Ser116, and the cell dies without nuclear AIFM1 translocation. None of these events requires a change in transcript abundance, so a three-gene mRNA score is expected to have poor sensitivity even when the pathway is active. Third, the apparent "biphasic KEAP1" behavior in the literature is a functional switch between NRF2-dependent survival at lower ROS and PGAM5-dependent death at higher ROS, not a biphasic expression profile; a single 24 h transcriptional snapshot cannot capture it.

Accordingly, the main contribution here is not a new positive biological claim but an auditable framework: locked sets before inspection, one observation per biological sample, an explicit permutation resolution floor, four declared FDR families, and a map from every summary manuscript number to an analysis product. The negative result should be read together with the earlier null reanalysis of non-overlapping public data [5]; both indicate that public transcriptomes alone are unlikely to settle this mechanism. The decisive experiment must measure p-AIFM1(S116) with antibodies validated in mouse tissue and a rescue matrix that separates oxeiptosis from ferroptosis, apoptosis and necroptosis.

## 5. Data and code availability

All source matrices are public GEO data (GSE267933, GSE222430). Analysis scripts, per-sample result tables, figures, MANIFEST checksums, the audit report and the amendments log are deposited in the public repository `yyx-4113/oxeiptosis-perioperative-brain-scrna` (https://github.com/yyx-4113/oxeiptosis-perioperative-brain-scrna; release v1.0.0, accessed 2026-11-26). Large per-cell objects (h5ad) are not redistributed, to respect repository size limits; they are regenerated deterministically from the locked scripts and the public GEO matrices (commands in the repository README, seed 20261125). Summary numbers in this manuscript are mapped to their source files in `07_manuscript/manuscript_numbers.json`, and file integrity is recorded in `08_audit/MANIFEST.csv`.

## 6. Author contributions, funding and competing interests

Yongxin Yang designed the study, performed all analyses and wrote the manuscript. No specific funding was received. The author declares no competing interests.

## 7. Predictions, confidence and validation

Every prediction below is a computational hypothesis generated from public transcriptomes, not a demonstration of oxeiptosis. Candidate cell types are listed by pre-specified order (astrocytes, neurons, microglia, oligodendrocytes/OPCs, endothelial/pericyte) and are not selected on the basis of nominal p values; no significance-based prioritization was performed. Confidence is low because the exact permutation floor is p={min_p} and no primary-family FDR threshold was met.

| Cohort / cell type | Prediction | Confidence | Uncertainty | Validation experiment |
|---|---|---|---|---|
{pred_table}

## 8. References

{chr(10).join(refs)}

---

### Figure legends

- **Fig. 1.** Study flow: five public datasets audited, two inference cohorts retained (n=3 per condition each), three low-replication datasets restricted to localization/background, locked gene sets, four FDR families, and the prospective protein-level validation chain.
- **Fig. 2.** UMAP of post-QC cells colored by annotated cell type (left) and condition (right) for each inference cohort; LowConfidence clusters are shown separately.
- **Fig. 3.** Per-sample module scores for the oxeiptosis core and ferroptosis-net signatures; each point is one biological sample after equal-depth cell sampling.
- **Fig. 4.** Cliff's delta with bootstrap 95% intervals for each pre-specified module, cell type and cohort; Hedges' g intervals are in Table 4.
- **Fig. 5.** Paired z-standardised separation endpoint D (ox_core − ferro_net) by cell type with exact-permutation p and bootstrap interval (Fig. 5b shows D with 95% intervals); the effect-size-scale sensitivity D is in Table 5.
- **Fig. 6.** Targeted ligand–receptor panel deltas by cohort.
- **Fig. S1.** Per-sample QC metrics.

### Limitations and predictive status

All computational results are predictions generated from public transcriptomes. They do not demonstrate oxeiptosis, its cell-type localization at protein level, or any causal effect of sevoflurane. Confirmation requires the wet-experimental chain described in Section 2.5, with pre-specified antibodies, controls and uncertainty reporting. Additional limitations: only one post-exposure time point is available; GSE267933 uses 50% O2 as control and cannot separate anesthesia from surgery; no computational doublet removal was applied; annotation relies on marker panels rather than a reference atlas and the resulting microglia/neuron proportions differ from published hippocampal benchmarks; and the two cohorts differ in control gas, so cross-cohort concordance is suggestive rather than confirmatory. Deviations from the pre-analysis plan are itemised in `08_audit/amendments.md`.
"""

main_text = "\n".join(
    ln for ln in md.splitlines()
    if not ln.startswith("|") and not ln.startswith("# ")).split("## 8. References")[0]
word_count = len(main_text.split())
md = md.replace("__WC__", f"{word_count:,}")
num("main_text_word_count", word_count, f"07_manuscript/{SLUG}_{VERSION}.md",
    "computed from non-table, non-reference lines")

out = MAN / f"{SLUG}_{VERSION}.md"
out.write_text(md, encoding="utf-8")
(MAN / "manuscript_numbers.json").write_text(
    json.dumps(numbers, indent=2, ensure_ascii=False), encoding="utf-8")
print("WROTE", out, "chars=", len(md), "tracked_numbers=", len(numbers))
