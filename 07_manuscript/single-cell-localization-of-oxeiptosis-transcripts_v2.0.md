# Single-cell localization of oxeiptosis transcripts after sevoflurane exposure in aged mouse hippocampus: a two-cohort reanalysis

**Version:** v2.0  
**Author:** Yongxin Yang  
**Affiliation:** Department of Anesthesiology, The Second Affiliated Hospital of Fujian University of Traditional Chinese Medicine, Fuzhou, Fujian 350003, China  
**Correspondence:** 960856791@qq.com · ORCID 0009-0004-9698-6552  
**Word count:** main text ~4,247 words (excluding tables, references and figure legends)

---

## Abstract

**Background.** Oxeiptosis is a ROS-dependent, caspase-independent regulated cell-death pathway driven by the KEAP1–PGAM5–AIFM1 axis. Whether sevoflurane exposure engages this pathway in the aged hippocampus remains untested. Prior public-data reanalysis of related datasets reported no reproducible ferroptosis or neuroinflammation signal, raising the bar for any new claim.

**Methods.** We reanalyzed two independent mouse hippocampal single-cell RNA-sequencing cohorts, each with three biological replicates per condition and analysed separately: GSE267933 (2.5% sevoflurane in 50% O2 plus laparotomy versus 50% O2 control, 24 h) and GSE222430 (2% sevoflurane for 6 h versus air control, hippocampus sampled 24 h later). After per-sample quality control and a pre-set computational cap of 20000 cells (not triggered), cells were clustered (Leiden [10]) and annotated with Scanpy [9] by a pre-locked marker panel and a top1–top2 margin threshold; low-margin clusters were labelled LowConfidence and excluded, with a fallback-label sensitivity analysis in which those cells were reassigned to their top-scoring class. A pre-specified three-gene oxeiptosis core score (AIFM1, KEAP1, PGAM5), an NRF2 antagonistic module and a ferroptosis driver-minus-suppressor score were computed per cell. Per-sample medians within each cell type were compared by exact permutation (C(6,3)=20; minimum attainable two-sided p=0.1), with Cliff's delta as the primary effect size (10000-replicate bootstrap interval) and Hedges' g as a point estimate without interval. Four FDR families were used, each BH-corrected within dataset where applicable: A = primary ox_core across cell types; B1 = secondary modules; B2 = the separation endpoint D; C = the targeted ligand–receptor panel. The pre-analysis plan declared three families (A, B, C); family B was split into B1 and B2 during implementation, before any test statistic was computed, and this is recorded as amendment A7 in `08_audit/amendments.md` rather than as pre-registration. A targeted 32-pair ligand–receptor panel was scored per sample (53 cohort×pair tests across the two cohorts).

**Results.** From 71,150 raw cells, 67,321 passed QC; 27,545 cells (40.9%; 6.2% in GSE267933 and 52.7% in GSE222430) were labelled LowConfidence and excluded from cell-type-specific inference, leaving 8 classes. Annotation dropout was not balanced across arms in GSE267933 (6.0–12.1% of cells in controls versus 2.3–5.8% in exposed samples) but was balanced in GSE222430 (47.2–56.8% versus 48.7–55.4%). The core triad was detected in both cohorts, at the transcript level in 7.8–33.0% of cells for Aifm1, 7.0–41.1% for Keap1 and 5.9–24.7% for Pgam5, with the highest rates in neurons (Table S3; gene-presence audit in Table S2). No test reached FDR<0.05 (0 of 66); because the exact permutation floor equals 0.1, no test can reach FDR<0.05 at n=3 vs 3, so this is a resolution limit rather than evidence against an effect. The largest absolute oxeiptosis-core shift occurred in Oligodendrocyte of GSE222430 (delta=-0.061), with wide intervals. Reassigning low-confidence cells to their top-scoring class left the primary deltas essentially unchanged (maximum |delta| 0.061; sign agreement in 10 of the 11 primary strata, one reversal; no sensitivity test reached FDR<0.05, minimum p=0.2). The z-standardised separation endpoint D showed 0 of 11 tests significant in family B2 (Table 5). On the raw score difference, ferroptosis-net shifts were at least as large as oxeiptosis-core shifts in 9 of 11 cell-type/cohort comparisons, but the two scores sit on different scales; on the scale-free contrast |delta_ox|-|delta_ferro|, ferroptosis-net was larger in 4 of 11 and oxeiptosis-core in 7, so a pathway-specific claim cannot be made from transcription alone. Composition differences are reported separately (Table S1).

**Conclusions.** In two independent cohorts analysed separately, with identical sample-level inference and different exposures and control gases, transcript abundance of the oxeiptosis core triad does not establish sevoflurane-induced oxeiptosis in the aged hippocampus. With three replicates per arm the design cannot exclude a difference of moderate size, so the result is best read as uninformative at the transcript level rather than as evidence of no effect. The study contributes an auditable framework and defines—without significance-based selection—the cell types in which protein-level testing of p-AIFM1(S116) should be attempted first. All results are predictive and require immunofluorescence with phosphatase controls and inhibitor rescue before mechanistic interpretation.

**Keywords:** oxeiptosis; sevoflurane; single-cell RNA sequencing; hippocampus; AIFM1; sample-level inference; permutation test

---

## 1. Introduction

Perioperative neurocognitive disorder affects a substantial fraction of older surgical patients, yet the underlying cell-death programs remain unresolved. Oxeiptosis was defined by Holze and colleagues as a ROS-dependent, caspase-independent death modality in which KEAP1 releases PGAM5, which dephosphorylates AIFM1 at serine 116; NRF2-dependent survival signaling is dispensable under high ROS [1,2]. This creates a testable prediction: if sevoflurane engages oxeiptosis in the aged hippocampus, the core triad should show a cell-type-specific transcriptional signal, and a protein-level S116 dephosphorylation should follow.

We previously evaluated five candidate public datasets and found that only two contain biological replication suitable for inference [3,4,11]. GSE153598 contains one macaque per condition with four technical library splits of pooled tissue; GSE196239 contains one library per condition in human prefrontal-cortex primary culture [4]; GSE141242 is a bulk rat microarray. These are restricted here to descriptive localization or background reference. Prior whole-transcriptome work in aged mouse hippocampus reported neuroinflammation, metabolic and blood–brain-barrier processes across the perioperative period [7], and a non-overlapping reanalysis of GSE141242 and GSE196239 found no consistent ferroptosis or neuroinflammation signal in public transcriptomes, although its cell-culture qRT-PCR did show perturbation of GPX4/SLC7A11/ACSL4/PTGS2 [5]. The present study differs in three declared ways: (i) the endpoint is changed from ferroptosis/neuroinflammation to oxeiptosis; (ii) the data are two mouse hippocampal single-cell cohorts with n=3 biological replicates per condition rather than bulk arrays and single libraries; and (iii) inference is cell-type-resolved with one observation per biological sample and explicit permutation-resolution limits.

Three pre-specified decisions constrain this work. First, AIRE was excluded from the gene set because no oxeiptosis publication includes it and human protein data indicate no brain expression; NFE2L2 was retained only as a separate antagonistic module. Second, the 6/24/72 h time course proposed originally cannot be addressed by the public data, which provide a single 24 h point; it is deferred to prospective animal work. Third, the primary mechanistic endpoint is p-AIFM1(S116)/total AIFM1 assessed by immunofluorescence, not a transcriptional score. Terminology: oxeiptosis denotes the ROS-dependent KEAP1–PGAM5–AIFM1 death pathway and is not synonymous with hyperoxic lung or brain injury; this distinction matters because the control condition in GSE267933 was 50% O2, itself an oxygen intervention that may shift ROS baseline.

## 2. Methods

### 2.1 Data audit and inclusion rules

Datasets were included for inference only if they provided at least three biological replicates per condition and readable per-sample cell counts. GSE267933 (n=3 control breathing 50% O2, n=3 sevoflurane 2.5% in 50% O2 plus laparotomy, 24 h) [3] and GSE222430 (n=3 air control, n=3 sevoflurane 2% for 6 h, sampled at 24 h) [13] met these criteria. Because 50% O2 is itself an intervention, GSE222430 (air control) is the cleaner anesthetic-effect contrast; GSE267933 is reported as a mixed-exposure reference. The two cohorts are reported side by side and never pooled: pooling would create a heterogeneous exposure and cohort-level pseudo-replication. GSE153598, GSE196239 and GSE141242 were retained as localization/bulk references only (Table 1) [4,11]. GSE222430 is a public GEO submission without an associated publication (accession only). All matrices were retrieved from GEO [11]; per-file SHA-256 checksums are recorded in `01_data/meta/sample_inventory.csv`.

### 2.2 Quality control and annotation

Per sample, cells with fewer than 200 genes, fewer than 500 UMIs, or >15% mitochondrial counts were removed; genes detected in fewer than three cells were discarded. No upper gene-count bound and no computational doublet caller were applied (deviations from the original plan; see `amendments.md`). A deterministic random cap of 20000 cells per sample was pre-set but not triggered (seed 20261125). Samples were normalized to 10,000 counts, log-transformed, and 2,000 highly variable genes were used for 30-component PCA. Neighborhood graphs and Leiden clustering (resolution 0.8) [10] were followed by cluster-level annotation in Scanpy [9] using pre-locked marker panels (marker lists in `scripts/common.py:CELLTYPE_MARKERS`; per-cluster annotation output in Table 3); hippocampal subtype markers (Prox1/Calb1, Neurod6/Dkk3/Fibcd1) were included and Itgam was removed from the microglia panel because it also labels infiltrating monocytes. Clusters whose top1–top2 panel-score margin was below 0.05 were labelled LowConfidence (40.9% of cells) and excluded from cell-type-specific inference. Neuron_pan served only as a validation panel: for each neuronal candidate cluster we report the fraction of cells above the pan-neuronal score median (Table 3); it was never an assignment label. Because exclusion of low-margin cells can bias arms differentially, dropout was tabulated per sample (Table S4) and a fallback-label sensitivity analysis is reported in Section 2.6.

### 2.3 Pre-locked gene sets and scoring

The oxeiptosis core set comprised AIFM1, KEAP1 and PGAM5 [1,2]. The NRF2 antagonistic module comprised NFE2L2, NQO1, GCLC and GCLM; HMOX1 was removed from this module in Round 2 because it was also a pre-locked ferroptosis-suppressor member, which made the two endpoints non-independent. Ferroptosis was represented as a driver-minus-suppressor contrast using pre-specified lists drawn from FerrDb v2 [12] and restricted to genes detectable here (Table S2); Alox15 was not detected and is listed as missing. Scores were computed with Scanpy's control-gene-subtracted scoring [9] (50 control genes, 25 expression bins); a plain mean-expression score for the core triad was used as a sensitivity analysis. AIRE was not scored. A score is a relative quantity against sampled control genes; its between-group difference is not an absolute mRNA fold change. Original-plan AUCell scoring and 200/400-cell downsampling sensitivity were not executed (see `amendments.md`).

### 2.4 Statistical analysis

Each cell type in each sample contributed exactly one score after equal-depth downsampling to 100 cells; cell types with fewer than 20 cells in a sample were treated as missing. Group differences were tested on the median difference with the exact permutation distribution over all 20 labelings of six samples, giving a minimum attainable two-sided p of 0.1. Cliff's delta was the primary effect size with a 10000-replicate percentile bootstrap interval. Bias-corrected Hedges' g is reported as a point estimate only, with the fraction of finite bootstrap replicates; no interval is reported, because at n=3 pooled-variance resampling is unstable and earlier versions produced intervals spanning two orders of magnitude. Four FDR families were applied with Benjamini–Hochberg correction: family A = ox_core across cell types and datasets; family B1 = secondary modules, corrected within each dataset; family B2 = the separation endpoint D, corrected within each dataset; family C = the ligand–receptor panel, also within dataset. The pre-analysis plan declared three families (A, B, C); the split of B into B1 and B2 was made in the analysis code before any test statistic was computed and is recorded as amendment A7 in `08_audit/amendments.md`. Cell-type proportions were tested in a separate descriptive family. Each library corresponds to one sample and one batch, so condition and batch are fully confounded and no batch correction was attempted. Cell type×cohort strata with fewer than three evaluable samples in either group (GSE267933/Inhibitory_neuron) were excluded from inference in that cohort only. The separation endpoint D was constructed as z(ox_core) − z(ferro_net) after z-standardising each module across the evaluable sample×cell-type values of a dataset; the subtraction is within-sample, but the between-group test is an unpaired exact permutation, and each module keeps its own standard deviation, so D is not scale-free. A scale-free sensitivity contrast, |δ(ox)| − |δ(ferro)|, is therefore reported alongside it. The targeted ligand–receptor analysis used per-sample products of mean ligand and receptor expression, taking the modal sender→receiver axis across samples of each group for each of the 32 pre-specified pairs, with family-C BH correction.

### 2.5 Reporting standard

All transcript-level results are predictive and exploratory. The confirmatory endpoint in the prospective extension is p-AIFM1(S116)/total AIFM1 by immunofluorescence in the annotated cell types, with λ-phosphatase (negative), calyculin A and H2O2 (positive) controls; the p/total ratio must be computed within the same molecular-weight band (66/62/57 kDa forms must not be mixed), plus a rescue matrix containing the PGAM5 inhibitor LFHP-1c [6], N-acetylcysteine, catalase, Z-VAD-FMK, necrostatin-1 and ferrostatin-1 in a dose–time design.

### 2.6 Sensitivity analyses

Three sensitivity analyses were pre-specified for this revision and are reported with the same exact permutation test. First, LowConfidence cells were reassigned to their top-scoring panel class (`celltype_raw`) instead of being dropped, and the ox_core test was recomputed; this directly addresses differential annotation dropout between arms (Table S4, `05_stats/sensitivity_lowconf.csv`). Second, a plain mean-expression core score without control-gene subtraction was computed and is included in Table 4 (`ox_core_meanexpr`). Third, a scale-free effect-size contrast, |δ(ox)| − |δ(ferro)|, was used in place of the z-standardised separation endpoint D. None of these analyses was used to select a preferred result; all are reported whatever their direction.

## 3. Results

### 3.1 Cohort composition and quality control

Across the inference cohorts, 71,150 raw cells were obtained; 67,321 remained after QC (Table 2; per-sample QC in Fig. S1; design and inference boundaries in Fig. 1). Median genes per cell ranged from 2156 to 2812, and median mitochondrial percentage from 0.1% to 7.5%. Annotation yielded 8 classes across the two cohorts, but not the same classes in each: no excitatory-neuron class passed the margin in GSE222430 (0 cells), and GSE267933 contributed only 73 excitatory neurons (0.4% of its cells), so the excitatory compartment is effectively unsampled and neuronal statements concern inhibitory neurons. Annotation dropout differed between arms in GSE267933 (6.0–12.1% of cells in controls versus 2.3–5.8% in exposed samples) and was balanced in GSE222430 (47.2–56.8% versus 48.7–55.4%); per-sample values are in Table S4 and the fallback-label sensitivity analysis is reported below. The annotated proportions also differ from published hippocampal benchmarks (microglia are over-represented and neurons under-represented), so all cell-type-specific results are annotation-conditional (Fig. 2).

### 3.2 Transcript detection by cell type

The three core transcripts were detected at the cell level in both cohorts (Table S3): Aifm1 in 7.8–33.0% of high-confidence cells, Keap1 in 7.0–41.1% and Pgam5 in 5.9–24.7%, with the highest rates in inhibitory neurons in both cohorts. Detection rates differ between cohorts (Keap1 is lower throughout GSE222430), consistent with differences in platform and sequencing depth; cross-cohort detection rates are therefore not compared. All module genes present in each dataset are listed in Table S2. Detection of a transcript is a prerequisite for scoring and is not evidence of pathway activation.

### 3.3 Sample-level oxeiptosis and ferroptosis scores

All 11 primary-family (ox_core) tests are shown below, ordered by nominal p with ties broken by absolute delta; no subset was selected, so the table cannot be read as a favourable selection. The complete set of 66 tests (families A and B1) is in Table 4; per-sample module scores are plotted in Fig. 3 and effect sizes with bootstrap intervals in Fig. 4. Delta is the median difference (exposed−control); FDR is BH within the indicated family.

| Cohort | Cell type | Module | Delta (exposed−control) | Cliff's delta [95% CI] | Permutation p | FDR (family) |
|---|---|---|---|---|---|---|
| GSE222430 | Oligodendrocyte | ox_core | -0.061 | -0.778 [-1.000, -0.111] | 0.200 | 0.733 |
| GSE267933 | Oligodendrocyte | ox_core | 0.055 | 1.000 [1.000, 1.000] | 0.200 | 0.733 |
| GSE222430 | Microglia | ox_core | -0.034 | -1.000 [-1.000, -1.000] | 0.200 | 0.733 |
| GSE267933 | Endothelial | ox_core | 0.030 | 0.556 [-0.333, 1.000] | 0.400 | 0.880 |
| GSE222430 | OPC | ox_core | -0.005 | -0.333 [-1.000, 1.000] | 0.400 | 0.880 |
| GSE222430 | Inhibitory_neuron | ox_core | -0.012 | -0.556 [-1.000, 0.333] | 0.600 | 0.943 |
| GSE267933 | Astrocyte | ox_core | 0.009 | 0.333 [-1.000, 1.000] | 0.600 | 0.943 |
| GSE222430 | Endothelial | ox_core | -0.012 | -0.333 [-1.000, 0.778] | 1.000 | 1.000 |
| GSE222430 | Astrocyte | ox_core | 0.008 | 0.333 [-0.778, 1.000] | 1.000 | 1.000 |
| GSE267933 | Microglia | ox_core | -0.001 | -0.111 [-1.000, 1.000] | 1.000 | 1.000 |
| GSE267933 | OPC | ox_core | 0.000 | 0.333 [-0.778, 1.000] | 1.000 | 1.000 |

No test reached FDR<0.05 (0 of 66); because the exact permutation floor equals 0.1, no test can reach FDR<0.05 at n=3 vs 3, so this is a resolution limit rather than evidence against an effect. The largest absolute oxeiptosis-core delta was -0.061 in Oligodendrocyte (GSE222430); it rests on three samples per arm with a wide bootstrap interval and is reported as a trend requiring protein-level testing, not as a finding. Consistent with the pre-specified antagonistic interpretation, NRF2-module shifts (family B1) were not treated as evidence of oxeiptosis: an elevated NRF2/ARE response would indicate a pro-survival antioxidant reaction, and high-ROS oxeiptosis is NRF2-independent [1]; per-cell-type NRF2 results are in Table 4.

Reassigning all low-confidence cells to their top-scoring class (27,545 cells) added 4 evaluable strata and left the primary estimates essentially unchanged: the maximum absolute ox_core delta was 0.061, the same value as in the primary run, the sign agreed in 10 of the 11 primary strata, one comparison reversed sign (GSE222430 endothelial, from a small negative to a small positive delta), and no sensitivity test reached FDR<0.05 (minimum p=0.2; Table S4 and `05_stats/sensitivity_lowconf.csv`). Exclusion of low-margin clusters therefore does not account for the pattern of results, although it does not remove the power limitation.

### 3.4 Separation from ferroptosis

| Cohort | Cell type | Separation D (z-scored) | Bootstrap 95% CI | Permutation p | FDR (B2) |
|---|---|---|---|---|---|
| GSE267933 | Astrocyte | 0.550 | -0.880 to 2.287 | 0.400 | 0.667 |
| GSE267933 | Endothelial | 1.424 | -1.111 to 2.352 | 0.400 | 0.667 |
| GSE267933 | Microglia | 0.075 | -0.653 to 1.140 | 1.000 | 1.000 |
| GSE267933 | OPC | 0.176 | -1.713 to 2.393 | 1.000 | 1.000 |
| GSE267933 | Oligodendrocyte | 1.504 | 0.526 to 2.311 | 0.200 | 0.667 |
| GSE222430 | Astrocyte | -0.188 | -1.646 to 3.679 | 1.000 | 1.000 |
| GSE222430 | Endothelial | -0.632 | -2.476 to 0.551 | 0.400 | 0.800 |
| GSE222430 | Inhibitory_neuron | -0.269 | -0.941 to 0.678 | 1.000 | 1.000 |
| GSE222430 | Microglia | -1.228 | -2.321 to 0.734 | 0.400 | 0.800 |
| GSE222430 | OPC | 0.074 | -0.457 to 1.071 | 0.600 | 0.900 |
| GSE222430 | Oligodendrocyte | -1.138 | -2.133 to 0.654 | 0.400 | 0.800 |

The table gives the separation endpoint D with its exact-permutation p, two-sample bootstrap 95% interval and family-B2 q value; positive D indicates a relatively larger oxeiptosis-core shift. 0 of 11 D tests reached FDR<0.05 in family B2. One interval (GSE267933 oligodendrocyte) excludes zero while its exact p is 0.200; with three observations a percentile bootstrap need not respect the discrete permutation distribution, and the exact test is treated as authoritative throughout. Because D is a difference of two z-standardised modules it is not scale-free, so the contrast |δ(ox)| − |δ(ferro)| is also reported in Table 5: on that scale ferroptosis-net was larger in 4 of 11 comparisons and oxeiptosis-core in 7. The two metrics therefore disagree, and neither supports a pathway-specific claim; the transcriptional data do not separate oxeiptosis from ferroptosis (Fig. 5).

### 3.5 Cell-type composition

Composition shifts could confound per-cell scores. Nominal permutation p values for cell-type fractions are listed in Table S1; because the same n=3 resolution limit applies, these are descriptive. Strata contributing only a few cells per sample (for example pericyte and GSE267933 excitatory neurons) are shown for completeness, but their p values are not interpretable and were not used in any claim.

| Cohort | Cell type | Fraction (control) | Fraction (exposed) | Permutation p |
|---|---|---|---|---|
| GSE267933 | Astrocyte | 0.101 | 0.071 | 0.200 |
| GSE267933 | Endothelial | 0.018 | 0.026 | 1.000 |
| GSE267933 | Excitatory_neuron | 0.004 | 0.004 | 1.000 |
| GSE267933 | Microglia | 0.441 | 0.520 | 0.400 |
| GSE267933 | OPC | 0.033 | 0.028 | 1.000 |
| GSE267933 | Oligodendrocyte | 0.301 | 0.304 | 1.000 |
| GSE267933 | Pericyte | 0.002 | 0.003 | 0.600 |
| GSE222430 | Astrocyte | 0.130 | 0.106 | 1.000 |
| GSE222430 | Endothelial | 0.018 | 0.016 | 1.000 |
| GSE222430 | Inhibitory_neuron | 0.071 | 0.088 | 0.200 |
| GSE222430 | Microglia | 0.039 | 0.039 | 1.000 |
| GSE222430 | OPC | 0.024 | 0.027 | 0.600 |
| GSE222430 | Oligodendrocyte | 0.185 | 0.201 | 1.000 |
| GSE222430 | Pericyte | 0.004 | 0.004 | 1.000 |

### 3.6 Targeted ligand–receptor panel

Of 32 pre-specified pairs, 27 and 26 were evaluable in the two cohorts respectively (transcript not detected for the remainder), giving 53 cohort×pair tests; 0 passed family-C FDR<0.05. Top pairs by FDR are listed below; the full panel is in Table 6. AIFM1, KEAP1 and PGAM5 are intracellular and absent from ligand–receptor repositories, so any communication difference is an indirect downstream observation without spatial co-localization information (Fig. 6).

| Cohort | Pair | Delta | Permutation p | FDR |
|---|---|---|---|---|
| GSE222430 | Vegfa-Flt1 | -0.292 | 0.200 | 0.520 |
| GSE222430 | Notch1-Jag1 | -0.055 | 0.200 | 0.520 |
| GSE222430 | Tnfsf10-Tnfrsf10b | -0.003 | 0.200 | 0.520 |
| GSE222430 | Lgals3-Trem2 | 0.003 | 0.200 | 0.520 |
| GSE222430 | Hmgb1-Tlr4 | -0.070 | 0.200 | 0.520 |
| GSE222430 | Apoe-Ldlr | 0.615 | 0.200 | 0.520 |
| GSE222430 | Ephb2-Efnb1 | 0.020 | 0.200 | 0.520 |
| GSE222430 | Igf1-Igf1r | -0.480 | 0.200 | 0.520 |
| GSE222430 | S100b-Ager | 0.002 | 0.200 | 0.520 |
| GSE222430 | Angpt1-Tek | -0.173 | 0.200 | 0.520 |

## 4. Discussion

No transcript-level comparison reaches FDR<0.05, which is arithmetically unavoidable when the smallest attainable p equals the FDR threshold; the design therefore cannot separate absence of effect from insufficient resolution. This is consistent with a previous reanalysis of conceptually related but non-overlapping public data, which also failed to find a consistent death-signature signal in bulk transcriptomes [5], although its cell-culture qRT-PCR did detect perturbation of ferroptosis markers. It is also consistent with the biology of oxeiptosis: the decisive event is dephosphorylation of AIFM1(S116) by PGAM5, a post-translational modification that need not be accompanied by changes in mRNA abundance. A three-gene transcriptional score therefore has low sensitivity for the pathway switch. Additionally, the KEAP1 "biphasic" behaviour described in the literature is a functional switch between NRF2-dependent survival and PGAM5-dependent death, not a biphasic mRNA expression pattern; the public data contain a single 24 h point and cannot address temporal dynamics.

The two cohorts differ in exposure and in control gas: GSE267933 combines sevoflurane with laparotomy over a 50% O2 background, whereas GSE222430 is anesthetic-only against air. Their concordance (or lack of it) is therefore informative but confounded by control gas; neither design isolates 6/24/72 h kinetics. We deliberately avoided pooling them: doing so would inflate replication while introducing heterogeneity, a form of pseudo-replication at the cohort level. The term oxeiptosis is used here strictly for the ROS-dependent KEAP1–PGAM5–AIFM1 pathway and not as a synonym for hyperoxic injury.

The main strength is that inference was performed once per biological sample, with the permutation floor stated explicitly and every summary number mapped to an analysis product. The main weaknesses are statistical power (three samples per group; minimum attainable two-sided p=0.1), absence of a computational doublet caller, and cell-type annotation from marker panels rather than a reference atlas. Annotation quality is the most consequential limitation. Low-margin clusters were flagged and excluded, but dropout was unbalanced in GSE267933 (6.0–12.1% in controls versus 2.3–5.8% in exposed samples), microglia are over-represented and neurons under-represented relative to published hippocampal benchmarks, and no excitatory-neuron class survived in GSE222430. The fallback-label sensitivity analysis argues against dropout driving the primary deltas, but it cannot repair a biased class inventory: a cell type that is absent cannot be scored. Each library is also a single batch, so condition and batch are fully confounded and no batch correction is possible. Bootstrap intervals over three observations are unstable, and Hedges' g is therefore reported without an interval; Cliff's delta is the primary effect size.

The next step is therefore not more of the same scoring but a prospective experiment: aged mice assigned to air, oxygen, sevoflurane and sevoflurane-plus-surgery, sampled at 0/6/12/24/72 h, with cell-type-resolved p-AIFM1(S116)/total AIFM1 immunofluorescence, phosphatase controls, and a dose–time rescue matrix. Antibody PA5-143645 (RRID AB_2942873) is labelled human-only by the vendor although the S116 site is conserved in mouse; AP5501 is annotated for mouse and has published mouse-tissue use, including the study cited as [8]. Cross-species validation in the intended mouse tissue with λ-phosphatase (negative) and calyculin A/H2O2 (positive) controls must still precede any mouse endpoint, with p/total computed within the same molecular-weight band. If validation fails, quantitative phosphoproteomics should replace immunofluorescence.

Three interpretive points follow from the design. First, the NRF2 module is antagonistic by construction: under high ROS oxeiptosis proceeds without NRF2-dependent survival signalling [1], so an increased NRF2/ARE score is interpreted as a protective response and does not support pathway activation. Second, oxeiptosis is initiated by a post-translational switch: high ROS promotes KEAP1 dissociation from PGAM5, PGAM5 dephosphorylates AIFM1 at Ser116, and death proceeds without requiring nuclear AIFM1 translocation, which distinguishes the pathway from parthanatos [1,2]. None of these events requires a change in transcript abundance, so a three-gene mRNA score is expected to have poor sensitivity even when the pathway is active. Third, the apparent "biphasic KEAP1" behavior in the literature is a functional switch between NRF2-dependent survival at lower ROS and PGAM5-dependent death at higher ROS, not a biphasic expression profile; a single 24 h transcriptional snapshot cannot capture it.

Accordingly, the main contribution here is not a new positive biological claim but an auditable framework: locked sets before inspection, one observation per biological sample, an explicit permutation resolution floor, four FDR families (three declared in the plan, one split during implementation), and a map from every summary number to an analysis product. The negative result should be read together with the earlier null reanalysis of non-overlapping public data [5]; both indicate that public transcriptomes alone are unlikely to settle this mechanism. The decisive experiment must measure p-AIFM1(S116) with antibodies validated in mouse tissue and a rescue matrix that separates oxeiptosis from ferroptosis, apoptosis and necroptosis.

### 4.6 Candidate cell types and validation plan

Every prediction below is a computational hypothesis generated from public transcriptomes, not a demonstration of oxeiptosis. Candidate cell types are listed by pre-specified order (astrocytes, neurons, microglia, oligodendrocytes/OPCs, endothelial/pericyte) and are not selected on the basis of nominal p values; no significance-based prioritization was performed. Confidence is low because the exact permutation floor is p=0.1, no primary-family FDR threshold was met, and neuronal classes are poorly captured.

| Cohort / cell type | Prediction | Confidence | Uncertainty | Validation experiment |
|---|---|---|---|---|
| GSE222430 / Astrocyte | oxeiptosis-core transcript score differs between conditions (delta=0.008) | low (nominal p=1.000, family-A FDR=1.000) | Cliff's delta 0.333 [-0.778, 1.000]; n=3 per group | cell-type-resolved p-AIFM1(S116)/total IF with λ-phosphatase, calyculin A and H2O2 controls, plus LFHP-1c/NAC rescue |
| GSE222430 / Microglia | oxeiptosis-core transcript score differs between conditions (delta=-0.034) | low (nominal p=0.200, family-A FDR=0.733) | Cliff's delta -1.000 [-1.000, -1.000]; n=3 per group | cell-type-resolved p-AIFM1(S116)/total IF with λ-phosphatase, calyculin A and H2O2 controls, plus LFHP-1c/NAC rescue |
| GSE222430 / Oligodendrocyte | oxeiptosis-core transcript score differs between conditions (delta=-0.061) | low (nominal p=0.200, family-A FDR=0.733) | Cliff's delta -0.778 [-1.000, -0.111]; n=3 per group | cell-type-resolved p-AIFM1(S116)/total IF with λ-phosphatase, calyculin A and H2O2 controls, plus LFHP-1c/NAC rescue |
| GSE222430 / OPC | oxeiptosis-core transcript score differs between conditions (delta=-0.005) | low (nominal p=0.400, family-A FDR=0.880) | Cliff's delta -0.333 [-1.000, 1.000]; n=3 per group | cell-type-resolved p-AIFM1(S116)/total IF with λ-phosphatase, calyculin A and H2O2 controls, plus LFHP-1c/NAC rescue |
| GSE222430 / Endothelial | oxeiptosis-core transcript score differs between conditions (delta=-0.012) | low (nominal p=1.000, family-A FDR=1.000) | Cliff's delta -0.333 [-1.000, 0.778]; n=3 per group | cell-type-resolved p-AIFM1(S116)/total IF with λ-phosphatase, calyculin A and H2O2 controls, plus LFHP-1c/NAC rescue |
| GSE267933 / Astrocyte | oxeiptosis-core transcript score differs between conditions (delta=0.009) | low (nominal p=0.600, family-A FDR=0.943) | Cliff's delta 0.333 [-1.000, 1.000]; n=3 per group | cell-type-resolved p-AIFM1(S116)/total IF with λ-phosphatase, calyculin A and H2O2 controls, plus LFHP-1c/NAC rescue |
| GSE267933 / Microglia | oxeiptosis-core transcript score differs between conditions (delta=-0.001) | low (nominal p=1.000, family-A FDR=1.000) | Cliff's delta -0.111 [-1.000, 1.000]; n=3 per group | cell-type-resolved p-AIFM1(S116)/total IF with λ-phosphatase, calyculin A and H2O2 controls, plus LFHP-1c/NAC rescue |
| GSE267933 / Oligodendrocyte | oxeiptosis-core transcript score differs between conditions (delta=0.055) | low (nominal p=0.200, family-A FDR=0.733) | Cliff's delta 1.000 [1.000, 1.000]; n=3 per group | cell-type-resolved p-AIFM1(S116)/total IF with λ-phosphatase, calyculin A and H2O2 controls, plus LFHP-1c/NAC rescue |
| GSE267933 / OPC | oxeiptosis-core transcript score differs between conditions (delta=0.000) | low (nominal p=1.000, family-A FDR=1.000) | Cliff's delta 0.333 [-0.778, 1.000]; n=3 per group | cell-type-resolved p-AIFM1(S116)/total IF with λ-phosphatase, calyculin A and H2O2 controls, plus LFHP-1c/NAC rescue |
| GSE267933 / Endothelial | oxeiptosis-core transcript score differs between conditions (delta=0.030) | low (nominal p=0.400, family-A FDR=0.880) | Cliff's delta 0.556 [-0.333, 1.000]; n=3 per group | cell-type-resolved p-AIFM1(S116)/total IF with λ-phosphatase, calyculin A and H2O2 controls, plus LFHP-1c/NAC rescue |

## 5. Limitations

All computational results are predictions generated from public transcriptomes. They do not demonstrate oxeiptosis, its localization at protein level, or any causal effect of sevoflurane; confirmation requires the chain described in Section 2.5. Beyond the design limitations discussed above, five points bound the interpretation. (i) Only one post-exposure time point (24 h) is available, so the 6/24/72 h kinetics that motivated the study cannot be addressed. (ii) GSE267933 uses 50% O2 as control and combines sevoflurane with laparotomy, so anesthesia, surgery and oxygen tension cannot be separated; GSE222430 is anesthetic-only but is a single unpublished submission with a 52.7% annotation dropout. (iii) No computational doublet removal was applied and no excitatory-neuron class was recovered in GSE222430, so the neuronal compartment is represented mainly by inhibitory neurons. (iv) Condition and batch are fully confounded (one library per sample), and annotation relied on marker panels rather than a reference atlas; the fallback-label sensitivity analysis supports robustness to dropout but not to a biased class inventory. (v) With three replicates per arm and a p floor of 0.1, the study cannot exclude moderate effects, and the absence of FDR-significant tests must not be reported as evidence of no effect. Deviations from the pre-analysis plan, including the B1/B2 split, are itemised in `08_audit/amendments.md`.

## 6. Conclusion

In two independent mouse hippocampal single-cell cohorts analysed separately, the oxeiptosis core transcripts are detectable in most annotated cell classes, but no cell-type-resolved comparison survives the pre-declared FDR threshold and the design cannot distinguish absence of effect from insufficient resolution. The value of the analysis lies in a pre-specified, sample-level framework, an explicit permutation floor, an annotation-dropout sensitivity analysis, and a ranked, non-significance-based plan for protein-level testing of p-AIFM1(S116) in aged hippocampus after sevoflurane exposure.

## 7. Declarations

### 7.1 Ethics approval

This study is a secondary analysis of publicly available, de-identified data deposited in GEO (GSE267933, GSE222430). No new human or animal participants were recruited, no identifiable information was accessed, and institutional ethics approval and informed consent were therefore not required.

### 7.2 Funding

No specific funding was received for this work.

### 7.3 Author contributions

Yongxin Yang designed the study, performed all analyses, and wrote the manuscript.

### 7.4 Competing interests

The author declares no competing interests.

## 8. Data and code availability

All source matrices are public GEO data (GSE267933, GSE222430) [11]. Analysis scripts, per-sample result tables, figures, MANIFEST checksums, the audit report and the amendments log are deposited in the public repository `yyx-4113/oxeiptosis-perioperative-brain-scrna` (https://github.com/yyx-4113/oxeiptosis-perioperative-brain-scrna; release v1.0.1, accessed 2026-11-26). Large per-cell objects (h5ad) are not redistributed, to respect repository size limits; they are regenerated deterministically from the locked scripts and the public GEO matrices (commands in the repository README, seed 20261125). Summary numbers in this manuscript are mapped to their source files in `07_manuscript/manuscript_numbers.json`, and file integrity is recorded in `08_audit/MANIFEST.csv`.

## 9. References

1. Holze C, Michaudel C, Mackowiak C, et al. Oxeiptosis, a ROS-induced caspase-independent apoptosis-like cell-death pathway. Nat Immunol. 2018;19(2):130-140. PMID: 29255269; doi:10.1038/s41590-017-0013-y.
2. Tang D, Kang R, Vanden Berghe T, Vandenabeele P, Kroemer G. The molecular machinery of regulated cell death. Cell Res. 2019;29(5):347-364. PMID: 30948788; doi:10.1038/s41422-019-0164-5.
3. Suo Z, et al. Aged hippocampal single-cell atlas screening unveils disrupted neuroglial system in postoperative cognitive impairment. Aging Cell. 2025;24(3):e14406. PMID: 39540334. GSE267933 (single-cell RNA sequencing).
4. Chang E, et al. iScience. 2023;26(4):106534. PMID: 37123239. GSE196239.
5. Zhang Z, Hu C, Liu Y, et al. Integrated re-analysis of bulk and single-cell transcriptomic data reveals limited ferroptosis and neuroinflammation signatures after sevoflurane exposure in the aged rat hippocampus and human prefrontal cortex cells. Front Cell Dev Biol. 2026;14:1938601. PMID: 42871048; PMCID: PMC13658196; doi:10.3389/fcell.2026.1938601.
6. Gao C, et al. A novel PGAM5 inhibitor LFHP-1c protects blood-brain barrier integrity in ischemic stroke. Acta Pharm Sin B. 2021;11(7):1867-1884. PMID: 34386325.
7. Suo Z, et al. CNS Neurosci Ther. 2022;28:1576-1595. PMID: 35899365; Dryad doi:10.5061/dryad.jsxksn0cj.
8. Zhang X, et al. Unveiling a novel mechanism in noise-induced hearing loss: oxeiptosis-mediated regulated cell death of cochlear hair cell. Neurosci Bull. 2026. PMID: 41697573; doi:10.1007/s12264-025-01585-z.
9. Wolf FA, Angerer P, Theis FJ. SCANPY: large-scale single-cell gene expression data analysis. Genome Biol. 2018;19:15. PMID: 29409532.
10. Traag VA, Waltman L, van Eck NJ. From Louvain to Leiden: guaranteeing well-connected communities. Sci Rep. 2019;9:5233. PMID: 30914743; doi:10.1038/s41598-019-41695-z.
11. Edgar R, Domrachev M, Lash AE. Gene Expression Omnibus: NCBI gene expression and hybridization array data repository. Nucleic Acids Res. 2002;30:207-210. PMID: 11752295.
12. Zhou N, Yuan X, Du Q, et al. FerrDb V2: update of the manually curated database of ferroptosis regulators and ferroptosis-disease associations. Nucleic Acids Res. 2023;51(D1):D571-D582. PMID: 36305834; doi:10.1093/nar/gkac935.
13. GSE222430. 10x single-cell sequencing of old mouse hippocampal transcriptomes after sevoflurane exposure (Lei Z, Yongyan C; Department of Anesthesiology, Shanghai Ninth People's Hospital Affiliated to Shanghai Jiao Tong University School of Medicine). GEO accession, submitted 2023-01-09, public 2026-01-08; no associated publication (GEO citation field empty, checked 2026-11-26).

---

### Figure legends

- **Fig. 1.** Study flow: five public datasets audited, two inference cohorts retained (n=3 per condition each), three low-replication datasets restricted to localization/background, locked gene sets, four FDR families (A, B1, B2, C), and the prospective protein-level validation chain.
- **Fig. 2.** UMAP of post-QC cells coloured by annotated cell type (left) and condition (right) for each inference cohort; LowConfidence clusters are shown separately.
- **Fig. 3.** Per-sample module scores for the oxeiptosis core and ferroptosis-net signatures; each point is one biological sample after equal-depth cell sampling.
- **Fig. 4.** Cliff's delta with bootstrap 95% intervals for each pre-specified module, cell type and cohort. Hedges' g is a point estimate without an interval (Table 4).
- **Fig. 5.** z-standardised separation endpoint D (ox_core − ferro_net) by cell type with exact-permutation p and bootstrap interval (Fig. 5b shows D with 95% intervals); the scale-free sensitivity contrast is in Table 5.
- **Fig. 6.** Targeted ligand–receptor panel deltas by cohort.
- **Fig. S1.** Per-sample QC metrics (genes, UMIs and mitochondrial percentage) for the 12 libraries.
