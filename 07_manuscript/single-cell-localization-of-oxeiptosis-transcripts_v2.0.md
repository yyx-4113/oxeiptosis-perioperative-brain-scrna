# Single-cell localization of oxeiptosis transcripts after sevoflurane exposure in aged mouse hippocampus: a two-cohort reanalysis

**Version:** v2.0  
**Author:** Yongxin Yang  
**Affiliation:** Department of Anesthesiology, The Second Affiliated Hospital of Fujian University of Traditional Chinese Medicine, Fuzhou, Fujian 350003, China  
**Correspondence:** 960856791@qq.com · ORCID 0009-0004-9698-6552  
**Word count:** main text ~2,978 words (excluding tables, references and figure legends)

---

## Abstract

**Background.** Oxeiptosis is a ROS-dependent, caspase-independent regulated cell-death pathway driven by the KEAP1–PGAM5–AIFM1 axis. Whether sevoflurane exposure engages this pathway in the aged hippocampus remains untested. Prior public-data reanalysis of related datasets reported no reproducible ferroptosis or neuroinflammation signal, raising the bar for any new claim.

**Methods.** We reanalyzed two independent mouse hippocampal single-cell RNA-sequencing cohorts with three biological replicates per condition: GSE267933 (2.5% sevoflurane in 50% O2 plus laparotomy versus 50% O2 control, 24 h) and GSE222430 (2% sevoflurane for 6 h versus air control, hippocampus sampled 24 h later). After per-sample quality control and a pre-set computational cap of 20000 cells (not triggered), cells were clustered and annotated by a pre-locked marker panel with a top1–top2 margin threshold; low-margin clusters were labelled LowConfidence and excluded from inference. A pre-specified three-gene oxeiptosis core score (AIFM1, KEAP1, PGAM5), an NRF2 antagonistic module and a ferroptosis driver-minus-suppressor score were computed per cell. Per-sample medians within each cell type were compared by exact permutation (C(6,3)=20; minimum attainable two-sided p=0.1), with Cliff's delta as the primary effect size (10000-replicate bootstrap interval) and Hedges' g as a point-estimate sensitivity measure. Four FDR families were used, each BH-corrected within dataset where applicable: A = primary ox_core across cell types; B1 = secondary modules; B2 = the separation endpoint D; C = the targeted ligand–receptor panel. The B1/B2 split was fixed in the analysis code before outcome inspection and is declared in `08_audit/amendments.md`. A targeted 32-pair ligand–receptor panel was scored per sample (53 cohort×pair tests across the two cohorts). Deviations from the original protocol are documented in `08_audit/amendments.md`.

**Results.** From 71,150 raw cells, 67,321 passed QC; 27,545 cells (40.9%; 6.2% and 52.7% in the two cohorts) were labelled LowConfidence and excluded from cell-type-specific inference. The remaining cells were assigned to 8 classes. The locked core triad was detected in both cohorts (gene-presence audit in Table S2). No test reached FDR<0.05 (0 of 66), and the minimum attainable p value was 0.1; therefore all transcript-level findings are exploratory trends rather than confirmed differences. The largest absolute oxeiptosis-core shift occurred in Oligodendrocyte of GSE222430 (delta=-0.061), with wide intervals and no family-A significance. The paired, z-standardised separation endpoint D showed 0 of 11 tests significant in family B2 (Table 5). Ferroptosis-net shifts were at least as large as oxeiptosis-core shifts in 9 of 11 cell-type/cohort comparisons, so a pathway-specific claim cannot be made from transcription alone. Composition differences are reported separately (Table S1).

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

Per sample, cells with fewer than 200 genes, fewer than 500 UMIs, or >15% mitochondrial counts were removed; genes detected in fewer than three cells were discarded. No upper gene-count bound and no computational doublet caller were applied (deviations from the original plan; see `amendments.md`). A deterministic random cap of 20000 cells per sample was pre-set but not triggered (seed 20261125). Samples were normalized to 10,000 counts, log-transformed, and 2,000 highly variable genes were used for 30-component PCA. Neighborhood graphs and Leiden clustering (resolution 0.8) were followed by cluster-level annotation using pre-locked marker panels (Table 3); hippocampal subtype markers (Prox1/Calb1, Neurod6/Dkk3/Fibcd1) were included and Itgam was removed from the microglia panel because it also labels infiltrating monocytes. Clusters whose top1–top2 panel-score margin was below 0.05 were labelled LowConfidence (40.9% of cells) and excluded from cell-type-specific inference. Neuron_pan served only as a validation panel: for each neuronal candidate cluster we report the fraction of cells above the pan-neuronal score median (Table 3); it was never an assignment label.

### 2.3 Pre-locked gene sets and scoring

The oxeiptosis core set comprised AIFM1, KEAP1 and PGAM5. The NRF2 antagonistic module comprised NFE2L2, NQO1, GCLC and GCLM; HMOX1 was removed from this module in Round 2 because it was also a pre-locked ferroptosis-suppressor member, which made the two endpoints non-independent. Ferroptosis was represented as a driver-minus-suppressor contrast using pre-specified lists drawn from FerrDb v2 (PMID 36305834) and restricted to genes detectable here (Table S2); Alox15 was not detected and is listed as missing. Scores were computed with scanpy's control-gene-subtracted scoring (50 control genes, 25 expression bins); a plain mean-expression score for the core triad was used as a sensitivity analysis. AIRE was not scored. A score is a relative quantity against sampled control genes; its between-group difference is not an absolute mRNA fold change. Original-plan AUCell scoring and 200/400-cell downsampling sensitivity were not executed (see `amendments.md`).

### 2.4 Statistical analysis

Each cell type in each sample contributed exactly one score after equal-depth downsampling to 100 cells; cell types with fewer than 20 cells in a sample were treated as missing. Group differences were tested on the median difference with the exact permutation distribution over all 20 labelings of six samples, giving a minimum attainable two-sided p of 0.1. Cliff's delta was the primary effect size with a 10000-replicate percentile bootstrap interval. Bias-corrected Hedges' g was reported as a point estimate with the fraction of bootstrap replicates that were finite; its interval was suppressed because at n=3 the pooled-variance resampling distribution is unstable and can span orders of magnitude (see Table 4 and `amendments.md`). Four FDR families were applied with Benjamini–Hochberg correction: family A = ox_core across cell types and datasets; family B1 = secondary modules, corrected within each dataset; family B2 = the separation endpoint D, corrected within each dataset; family C = the ligand–receptor panel, also within dataset. Cell-type proportions were tested in a separate descriptive family. Cell type×cohort strata with fewer than three evaluable samples in either group (GSE267933/Inhibitory_neuron) were excluded from inference in that cohort only. The separation endpoint D was constructed at the sample level as the paired sample-level difference z(ox_core) − z(ferro_net), with z-scores computed across all cells of each dataset, and tested by exact permutation with a bootstrap 95% interval; an effect-size-scale sensitivity D = |Cliff's delta(ox)| − |Cliff's delta(ferro)| was also computed. The targeted ligand–receptor analysis used per-sample products of mean ligand and receptor expression, taking the modal sender→receiver axis across samples of each group for each of the 32 pre-specified pairs, with family-C BH correction.

### 2.5 Reporting standard

All transcript-level results are predictive and exploratory. The confirmatory endpoint in the prospective extension is p-AIFM1(S116)/total AIFM1 by immunofluorescence in the annotated cell types, with λ-phosphatase (negative), calyculin A and H2O2 (positive) controls; the p/total ratio must be computed within the same molecular-weight band (66/62/57 kDa forms must not be mixed), plus a rescue matrix containing LFHP-1c, N-acetylcysteine, catalase, Z-VAD-FMK, necrostatin-1 and ferrostatin-1 in a dose–time design. Public de-identified GEO data were used; institutional ethics approval and informed consent were not required for this reanalysis.

## 3. Results

### 3.1 Cohort composition and quality control

Across the inference cohorts, 71,150 raw cells were obtained; 67,321 remained after QC (Table 2). Median genes per cell ranged from 2156 to 2812, and median mitochondrial percentage from 0.1% to 7.5%. Annotation yielded 8 cell classes, including neurons, astrocytes, microglia, oligodendrocytes and OPCs (Table 3; UMAP in Fig. 2).

### 3.2 Gene detection

The locked core triad was detected in both cohorts; details of present and missing genes for every module are given in Table S2. Detection of a transcript is a prerequisite for scoring and is not evidence of pathway activation.

### 3.3 Sample-level oxeiptosis and ferroptosis scores

Six representative tests, ordered by nominal p value, are shown below; the complete set of 66 tests (families A and B1) is in Table 4. Delta is the median difference (exposed−control); FDR is BH within the indicated family.

| Cohort | Cell type | Module | Delta (exposed−control) | Cliff's delta [95% CI] | Permutation p | FDR (family) |
|---|---|---|---|---|---|---|
| GSE267933 | Astrocyte | ferro_suppressor | 0.082 | 1.000 [1.000, 1.000] | 0.200 | 0.938 |
| GSE267933 | Astrocyte | ferro_net | -0.073 | -0.778 [-1.000, -0.111] | 0.200 | 0.938 |
| GSE222430 | Astrocyte | nrf2_antag | 0.034 | 1.000 [1.000, 1.000] | 0.200 | 0.600 |
| GSE267933 | Oligodendrocyte | ox_core_meanexpr | 0.056 | 1.000 [1.000, 1.000] | 0.200 | 0.938 |
| GSE267933 | Oligodendrocyte | ox_core | 0.055 | 1.000 [1.000, 1.000] | 0.200 | 0.733 |
| GSE222430 | Inhibitory_neuron | ferro_driver | 0.034 | 1.000 [1.000, 1.000] | 0.200 | 0.600 |

No test reached FDR<0.05 (0 of 66). Given the exact permutation floor of p=0.1, these data cannot establish a difference between conditions at conventional significance. The largest nominal oxeiptosis-core delta was -0.061 in Oligodendrocyte (GSE222430); because this is one of 66 comparisons and has a wide interval, it is reported as a trend requiring protein-level testing. Consistent with the pre-specified antagonistic interpretation, NRF2-module shifts (family B1) were not treated as evidence of oxeiptosis: an elevated NRF2/ARE response would indicate a pro-survival antioxidant reaction, and high-ROS oxeiptosis is NRF2-independent [1]; per-cell-type NRF2 results are in Table 4.

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

The table gives the paired, sample-level z-standardised separation endpoint D with its exact-permutation p, two-sample bootstrap 95% interval and family-B2 q value; positive D indicates a relatively larger oxeiptosis shift. 0 of 11 D tests reached FDR<0.05 in family B2. The hypothesis-free effect-size-scale sensitivity D = |δ(ox)| − |δ(ferro)| is provided in Table 5. Ferroptosis-net shifts were at least as large as oxeiptosis-core shifts in 9 of 11 cell-type/cohort comparisons; transcriptional data therefore do not establish pathway specificity (Fig. 5).

### 3.5 Cell-type composition

Composition shifts could confound per-cell scores. Nominal permutation p values for cell-type fractions are listed in Table S1; because the same n=3 resolution limit applies, these are descriptive.

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

This reanalysis returns a predominantly null result at the transcript level, which is scientifically informative given that a previous reanalysis of conceptually related but non-overlapping public data also failed to find a consistent death-signature signal in bulk transcriptomes [5], although its cell-culture qRT-PCR did detect perturbation of ferroptosis markers. The result is consistent with the biology of oxeiptosis: the decisive event is dephosphorylation of AIFM1(S116) by PGAM5, a post-translational modification that is not necessarily accompanied by changes in mRNA abundance. A three-gene transcriptional score therefore has low sensitivity for the pathway switch. Additionally, the KEAP1 "biphasic" behavior described in the literature is a functional switch between NRF2-dependent survival and PGAM5-dependent death, not a biphasic mRNA expression pattern; the public data contain a single 24 h point and cannot address temporal dynamics.

The two cohorts differ in exposure and in control gas: GSE267933 combines sevoflurane with laparotomy over a 50% O2 background, whereas GSE222430 is anesthetic-only against air. Their concordance (or lack of it) is therefore informative but confounded by control gas; neither design isolates 6/24/72 h kinetics. We deliberately avoided pooling them: doing so would inflate replication while introducing heterogeneity, a form of pseudo-replication at the cohort level. The term oxeiptosis is used here strictly for the ROS-dependent KEAP1–PGAM5–AIFM1 pathway and not as a synonym for hyperoxic injury.

The main strength is that inference was performed once per biological sample, with the permutation floor stated explicitly. The main weaknesses are statistical power (three samples per group; minimum attainable two-sided p=0.1), absence of a computational doublet caller, and cell-type annotation from marker panels rather than a reference atlas. Low-margin clusters were flagged and excluded, but the observed microglia/neuron proportions differ from published hippocampal benchmarks, so cell-type-specific results should be read as annotation-conditional. Bootstrap intervals over three observations are unstable, and Hedges' g intervals were suppressed when resamples had zero within-group variance; Cliff's delta is therefore the primary effect size.

The next step is therefore not more of the same scoring but a prospective experiment: aged mice assigned to air, oxygen, sevoflurane and sevoflurane-plus-surgery, sampled at 0/6/12/24/72 h, with cell-type-resolved p-AIFM1(S116)/total AIFM1 immunofluorescence, phosphatase controls, and a dose–time rescue matrix. Antibody PA5-143645 (RRID AB_2942873) is labelled human-only by the vendor although the S116 site is conserved in mouse; AP5501 is annotated for mouse and has published mouse-tissue use, including the study cited as [8]. Cross-species validation in the intended mouse tissue with λ-phosphatase (negative) and calyculin A/H2O2 (positive) controls must still precede any mouse endpoint, with p/total computed within the same molecular-weight band. If validation fails, quantitative phosphoproteomics should replace immunofluorescence.

Three interpretive points follow from the design. First, the NRF2 module is antagonistic by construction: under high ROS oxeiptosis proceeds without NRF2-dependent survival signalling [1], so an increased NRF2/ARE score is interpreted as a protective response and does not support pathway activation. Second, oxeiptosis is initiated by a post-translational switch: high ROS promotes KEAP1 dissociation from PGAM5, PGAM5 dephosphorylates AIFM1 at Ser116, and the cell dies without nuclear AIFM1 translocation. None of these events requires a change in transcript abundance, so a three-gene mRNA score is expected to have poor sensitivity even when the pathway is active. Third, the apparent "biphasic KEAP1" behavior in the literature is a functional switch between NRF2-dependent survival at lower ROS and PGAM5-dependent death at higher ROS, not a biphasic expression profile; a single 24 h transcriptional snapshot cannot capture it.

Accordingly, the main contribution here is not a new positive biological claim but an auditable framework: locked sets before inspection, one observation per biological sample, an explicit permutation resolution floor, four declared FDR families, and a map from every summary manuscript number to an analysis product. The negative result should be read together with the earlier null reanalysis of non-overlapping public data [5]; both indicate that public transcriptomes alone are unlikely to settle this mechanism. The decisive experiment must measure p-AIFM1(S116) with antibodies validated in mouse tissue and a rescue matrix that separates oxeiptosis from ferroptosis, apoptosis and necroptosis.

## 5. Data and code availability

All source matrices are public GEO data (GSE267933, GSE222430). Analysis scripts, per-sample result tables, figures, MANIFEST checksums, the audit report and the amendments log are deposited in the public repository `yyx-4113/oxeiptosis-perioperative-brain-scrna` (https://github.com/yyx-4113/oxeiptosis-perioperative-brain-scrna; release v1.0.0, accessed 2026-11-26). Large per-cell objects (h5ad) are not redistributed, to respect repository size limits; they are regenerated deterministically from the locked scripts and the public GEO matrices (commands in the repository README, seed 20261125). Summary numbers in this manuscript are mapped to their source files in `07_manuscript/manuscript_numbers.json`, and file integrity is recorded in `08_audit/MANIFEST.csv`.

## 6. Author contributions, funding and competing interests

Yongxin Yang designed the study, performed all analyses and wrote the manuscript. No specific funding was received. The author declares no competing interests.

## 7. Predictions, confidence and validation

Every prediction below is a computational hypothesis generated from public transcriptomes, not a demonstration of oxeiptosis. Candidate cell types are listed by pre-specified order (astrocytes, neurons, microglia, oligodendrocytes/OPCs, endothelial/pericyte) and are not selected on the basis of nominal p values; no significance-based prioritization was performed. Confidence is low because the exact permutation floor is p=0.1 and no primary-family FDR threshold was met.

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

## 8. References

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

- **Fig. 1.** Study flow: five public datasets audited, two inference cohorts retained (n=3 per condition each), three low-replication datasets restricted to localization/background, locked gene sets, four FDR families, and the prospective protein-level validation chain.
- **Fig. 2.** UMAP of post-QC cells colored by annotated cell type (left) and condition (right) for each inference cohort; LowConfidence clusters are shown separately.
- **Fig. 3.** Per-sample module scores for the oxeiptosis core and ferroptosis-net signatures; each point is one biological sample after equal-depth cell sampling.
- **Fig. 4.** Cliff's delta with bootstrap 95% intervals for each pre-specified module, cell type and cohort; Hedges' g intervals are in Table 4.
- **Fig. 5.** Paired z-standardised separation endpoint D (ox_core − ferro_net) by cell type with exact-permutation p and bootstrap interval (Fig. 5b shows D with 95% intervals); the effect-size-scale sensitivity D is in Table 5.
- **Fig. 6.** Targeted ligand–receptor panel deltas by cohort.
- **Fig. S1.** Per-sample QC metrics.

### Limitations and predictive status

All computational results are predictions generated from public transcriptomes. They do not demonstrate oxeiptosis, its cell-type localization at protein level, or any causal effect of sevoflurane. Confirmation requires the wet-experimental chain described in Section 2.5, with pre-specified antibodies, controls and uncertainty reporting. Additional limitations: only one post-exposure time point is available; GSE267933 uses 50% O2 as control and cannot separate anesthesia from surgery; no computational doublet removal was applied; annotation relies on marker panels rather than a reference atlas and the resulting microglia/neuron proportions differ from published hippocampal benchmarks; and the two cohorts differ in control gas, so cross-cohort concordance is suggestive rather than confirmatory. Deviations from the pre-analysis plan are itemised in `08_audit/amendments.md`.
