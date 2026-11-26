# Oxeiptosis localization after sevoflurane: reproducible single-cell reanalysis

This repository contains the complete analysis chain for the manuscript
"Single-cell localization of oxeiptosis transcripts after sevoflurane exposure
in aged mouse hippocampus: a two-cohort reanalysis".

## File map

- `scripts/common.py` — central locked configuration: gene sets, markers, QC and statistics.
- `scripts/00_env_check.py` — execution environment manifest (`00_environment/env.json`).
- `scripts/01b_parallel_download.py` — resumable HTTP-range downloader for GEO supplementary TARs.
- `scripts/01_organize_data.py` — 10x triplet organization and sample inventory with SHA-256.
- `scripts/10_qc_annotate.py` — per-sample QC, integration, Leiden clustering, marker annotation.
- `scripts/20_scores_stats.py` — per-cell module scores, sample-level aggregation, exact permutation tests.
- `scripts/30_lr_separation.py` — oxeiptosis/ferroptosis separation and targeted ligand–receptor panel.
- `scripts/40_figures_tables.py` — all figures and manuscript tables.
- `scripts/50_manifest.py` — checksum manifest of analysis products.
- `scripts/60_manuscript.py` — manuscript assembled from result tables.
- `scripts/90_audit.py` — independent recomputation and consistency checks.
- `08_audit/MANIFEST.csv` — SHA-256 checksums of all analysis products.
- `08_audit/audit_report.json` — PASS/FAIL audit report.

## Reproduction

```bash
python scripts/00_env_check.py
python scripts/01b_parallel_download.py GSE267933
python scripts/01b_parallel_download.py GSE222430
python scripts/01_organize_data.py
python scripts/10_qc_annotate.py
python scripts/20_scores_stats.py
python scripts/30_lr_separation.py
python scripts/40_figures_tables.py
python scripts/60_manuscript.py
python scripts/50_manifest.py
python scripts/90_audit.py
```

## Software

Python 3.13; numpy, scipy, pandas, anndata, scanpy, scikit-learn, h5py,
matplotlib, seaborn, leidenalg, igraph. Exact versions are recorded in
`00_environment/env.json`.

## Caveats

1. All transcript-level results are exploratory predictions and do not
   demonstrate oxeiptosis. The confirmatory readout is p-AIFM1(S116)/total AIFM1
   by immunofluorescence with phosphatase controls.
2. With three biological replicates per group, the exact two-sided permutation
   test has a minimum attainable p of 0.10; absence of significance is not
   evidence of absence of effect.
3. GSE267933 combines sevoflurane with laparotomy; GSE222430 is anesthetic-only.
   Cohorts are reported separately and never pooled.
4. A deterministic computational cap of 20,000 cells per sample is applied after
   biological QC; seed 20261125. Equal-depth sampling of 100 cells per sample
   and cell type is used for scoring.
5. No full CellChatDB run was performed. The ligand–receptor analysis is a
   pre-specified 32-pair, per-sample score and is an indirect downstream readout.
