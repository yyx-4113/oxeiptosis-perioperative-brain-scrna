# Single-cell localization of oxeiptosis transcripts after sevoflurane exposure in aged mouse hippocampus: a two-cohort reanalysis

**Version:** v1.0  
**Author:** Yongxin Yang  
**Affiliation:** Department of Anesthesiology, The Second Affiliated Hospital of Fujian University of Traditional Chinese Medicine, Fuzhou, Fujian 350003, China  
**Correspondence:** 960856791@qq.com · ORCID 0009-0004-9698-6552  
**Word count:** main text ~2,113 words (excluding tables, references and figure legends)

---

## Abstract

**Background.** Oxeiptosis is a ROS-dependent, caspase-independent regulated cell-death pathway driven by the KEAP1–PGAM5–AIFM1 axis. Whether sevoflurane exposure engages this pathway in the aged hippocampus remains untested. Prior public-data reanalysis of related datasets reported no reproducible ferroptosis or neuroinflammation signal, raising the bar for any new claim.

**Methods.** We reanalyzed two independent mouse hippocampal single-cell RNA-sequencing cohorts with three biological replicates per condition: GSE267933 (2.5% sevoflurane plus laparotomy versus oxygen control, 24 h) and GSE222430 (2% sevoflurane for 6 h versus air control, hippocampus sampled 24 h later). After per-sample quality control and computational downsampling to a maximum of 20000 cells, cells were clustered and annotated by a pre-locked marker panel. A pre-specified three-gene oxeiptosis core score (AIFM1, KEAP1, PGAM5), an NRF2 antagonistic module and a ferroptosis driver-minus-suppressor score were computed per cell. Per-sample means within each cell type were compared by exact permutation (C(6,3)=20; minimum attainable two-sided p=0.1), with Cliff's delta, Hedges' g and bootstrap 95% confidence intervals. FDR correction was applied across the pre-declared family of 78 tests. A targeted 32-pair ligand–receptor panel was scored per sample.

**Results.** From 71,150 raw cells, 67,321 passed QC and were assigned to 8 cell classes. The locked core triad was detected in both cohorts (gene-presence audit in Table S2). No test reached FDR<0.05 (0 of 78), and the minimum attainable p value was 0.1; therefore all transcript-level findings are exploratory trends rather than confirmed differences. The largest absolute oxeiptosis-core shift occurred in Oligodendrocyte of GSE222430 (delta=-0.040). Ferroptosis-net shifts were at least as large as oxeiptosis-core shifts in 9 of 13 cell-type/cohort comparisons, so a pathway-specific claim cannot be made from transcription alone. Composition differences are reported separately (Table S1).

**Conclusions.** In two replicate-level cohorts, transcript abundance of the oxeiptosis core triad does not provide reproducible evidence of sevoflurane-induced oxeiptosis in the aged hippocampus. The study establishes an auditable, sample-level framework and identifies the cell types in which protein-level testing of p-AIFM1(S116) should be concentrated. All results are predictive/hypothesis-generating and require immunofluorescence with phosphatase controls and inhibitor rescue before mechanistic interpretation.

**Keywords:** oxeiptosis; sevoflurane; single-cell RNA sequencing; hippocampus; AIFM1; sample-level inference; permutation test

---

## 1. Introduction

Perioperative neurocognitive disorder affects a substantial fraction of older surgical patients, yet the underlying cell-death programs remain unresolved. Oxeiptosis was defined by Holze and colleagues as a ROS-dependent, caspase-independent death modality in which KEAP1 releases PGAM5, which dephosphorylates AIFM1 at serine 116; NRF2-dependent survival signaling is dispensable under high ROS [1,2]. This creates a testable prediction: if sevoflurane engages oxeiptosis in the aged hippocampus, the core triad should show a cell-type-specific transcriptional signal, and a protein-level S116 dephosphorylation should follow.

We previously evaluated three candidate public datasets and found that only two contain biological replication suitable for inference. GSE153598 contains one monkey per condition with four technical library splits of pooled tissue; GSE196239 contains one library per condition in human prefrontal-cortex primary culture. Both are restricted here to descriptive cell-type localization. A closely related reanalysis of GSE141242 and GSE196239 reported no reproducible ferroptosis or neuroinflammation signal [5]; the present study therefore changes the endpoint rather than repeating that analysis.

Three pre-registered decisions constrain this work. First, AIRE was excluded from the gene set because no oxeiptosis publication includes it and human protein data indicate no brain expression; NFE2L2 was retained only as a separate antagonistic module. Second, the 6/24/72 h time course proposed originally cannot be addressed by the public data, which provide a single 24 h point; it is deferred to prospective animal work. Third, the primary mechanistic endpoint is p-AIFM1(S116)/total AIFM1 assessed by immunofluorescence, not a transcriptional score.

## 2. Methods

### 2.1 Data audit and inclusion rules

Datasets were included for inference only if they provided at least three biological replicates per condition and readable per-sample cell counts. GSE267933 (n=3 control, n=3 sevoflurane plus surgery, 24 h) and GSE222430 (n=3 air control, n=3 sevoflurane 2% for 6 h, sampled at 24 h) met these criteria. The two cohorts are reported side by side and never pooled: the first cannot separate anesthetic from surgical stress, the second has no surgical component. GSE153598, GSE196239 and GSE141242 were retained as localization/bulk references only (Table 1).

### 2.2 Quality control and annotation

Per sample, cells with fewer than 200 genes, fewer than 500 UMIs, or >15% mitochondrial counts were removed; genes detected in fewer than three cells were discarded. A deterministic random cap of 20000 cells per sample was applied after biological filtering for computational tractability (seed 20261125). Samples were normalized to 10,000 counts, log-transformed, and 2,000 highly variable genes were used for 30-component PCA. Neighborhood graphs and Leiden clustering (resolution 0.8) were followed by cluster-level annotation using pre-locked marker panels (Table 3). Neuron_pan served only as a validation panel, never as an assignment label.

### 2.3 Pre-locked gene sets and scoring

The oxeiptosis core set comprised AIFM1, KEAP1 and PGAM5. The NRF2 antagonistic module comprised NFE2L2, HMOX1, NQO1, GCLC and GCLM. Ferroptosis was represented as a driver-minus-suppressor contrast using pre-specified lists (Table S2). Scores were computed with scanpy's control-gene-subtracted scoring (50 control genes, 25 expression bins); a plain mean-expression score for the core triad was used as a sensitivity analysis. AIRE was not scored.

### 2.4 Statistical analysis

Each cell type in each sample contributed exactly one score after equal-depth downsampling to 100 cells; cell types with fewer than 20 cells in a sample were treated as missing. Group differences were tested with the exact permutation distribution over all 20 labelings of six samples, giving a minimum attainable two-sided p of 0.1. Effect sizes were Cliff's delta and bias-corrected Hedges' g, each with a 2,000-replicate bootstrap percentile interval over samples. The FDR family was "cell_type x module x dataset (inference tier only)" (78 tests). Cell-type proportions were tested in a separate family. Cell types with fewer than three evaluable samples in either group (Inhibitory_neuron) were excluded from inference. The targeted ligand–receptor analysis used per-sample products of mean ligand and receptor expression, taking the strongest sender→receiver axis for each of the 32 pre-specified pairs, with BH correction within that family.

### 2.5 Reporting standard

All transcript-level results are predictive and exploratory. The confirmatory endpoint in the prospective extension is p-AIFM1(S116)/total AIFM1 by immunofluorescence in the annotated cell types, with λ-phosphatase-negative and H2O2-positive controls, plus a rescue matrix containing LFHP-1c, N-acetylcysteine, catalase, Z-VAD-FMK, necrostatin-1 and ferrostatin-1.

## 3. Results

### 3.1 Cohort composition and quality control

Across the inference cohorts, 71,150 raw cells were obtained; 67,321 remained after QC (Table 2). Median genes per cell ranged from 2156 to 2812, and median mitochondrial percentage from 0.1% to 7.5%. Annotation yielded 8 cell classes, including neurons, astrocytes, microglia, oligodendrocytes and OPCs (Table 3; UMAP in Fig. 2).

### 3.2 Gene detection

The locked core triad was detected in both cohorts; details of present and missing genes for every module are given in Table S2. Detection of a transcript is a prerequisite for scoring and is not evidence of pathway activation.

### 3.3 Sample-level oxeiptosis and ferroptosis scores

Six representative tests, ordered by nominal p value, are shown below; the complete set of 78 tests is in Table 4.

| Cohort | Cell type | Module | Delta (exposed−control) | Cliff's delta [95% CI] | Permutation p | FDR |
|---|---|---|---|---|---|---|
| GSE267933 | Astrocyte | ferro_suppressor | 0.068 | 1.000 [1.000, 1.000] | 0.100 | 0.780 |
| GSE267933 | Excitatory_neuron | ferro_driver | -0.041 | -1.000 [-1.000, -1.000] | 0.100 | 0.780 |
| GSE222430 | Astrocyte | nrf2_antag | 0.035 | 1.000 [1.000, 1.000] | 0.100 | 0.780 |
| GSE267933 | Pericyte | ox_core | 0.021 | 1.000 [1.000, 1.000] | 0.100 | 0.780 |
| GSE222430 | Microglia | ferro_driver | 0.039 | 1.000 [1.000, 1.000] | 0.100 | 0.780 |
| GSE222430 | Inhibitory_neuron | ferro_driver | 0.024 | 1.000 [1.000, 1.000] | 0.100 | 0.780 |

No test reached FDR<0.05 (0 of 78). Given the exact permutation floor of p=0.1, these data cannot establish a difference between conditions at conventional significance. The largest nominal oxeiptosis-core delta was -0.040 in Oligodendrocyte (GSE222430); because this is one of 78 comparisons and has a wide interval, it is reported as a trend requiring protein-level testing.

### 3.4 Separation from ferroptosis

| Cohort | Cell type | Oxeiptosis-core delta | Ferroptosis-net delta | Separation |
|---|---|---|---|---|
| GSE222430 | Astrocyte | 0.025 | 0.029 | -0.003 |
| GSE222430 | Endothelial | -0.016 | 0.040 | -0.056 |
| GSE222430 | Inhibitory_neuron | -0.011 | -0.017 | 0.006 |
| GSE222430 | Microglia | -0.027 | 0.003 | -0.030 |
| GSE222430 | OPC | -0.004 | -0.075 | 0.071 |
| GSE222430 | Oligodendrocyte | -0.040 | -0.033 | -0.007 |
| GSE267933 | Astrocyte | 0.012 | -0.074 | 0.085 |
| GSE267933 | Endothelial | 0.025 | -0.024 | 0.049 |
| GSE267933 | Excitatory_neuron | -0.012 | 0.066 | -0.078 |
| GSE267933 | Microglia | -0.001 | -0.053 | 0.052 |
| GSE267933 | OPC | 0.006 | -0.026 | 0.032 |
| GSE267933 | Oligodendrocyte | 0.032 | -0.007 | 0.039 |
| GSE267933 | Pericyte | 0.021 | -0.066 | 0.087 |

The separation column is the difference of the two deltas; positive values indicate a relatively larger oxeiptosis shift. Ferroptosis-net shifts were at least as large as oxeiptosis-core shifts in 9 of 13 cell-type/cohort comparisons; transcriptional data therefore do not establish pathway specificity (Fig. 5).

### 3.5 Cell-type composition

Composition shifts could confound per-cell scores. Nominal permutation p values for cell-type fractions are listed in Table S1; because the same n=3 resolution limit applies, these are descriptive.

| Cohort | Cell type | Fraction (control) | Fraction (exposed) | Permutation p |
|---|---|---|---|---|
| GSE267933 | Astrocyte | 0.101 | 0.071 | 0.200 |
| GSE267933 | Endothelial | 0.018 | 0.026 | 0.600 |
| GSE267933 | Excitatory_neuron | 0.073 | 0.025 | 0.100 |
| GSE267933 | Microglia | 0.441 | 0.520 | 0.300 |
| GSE267933 | OPC | 0.033 | 0.028 | 0.600 |
| GSE267933 | Oligodendrocyte | 0.303 | 0.305 | 0.900 |
| GSE267933 | Pericyte | 0.025 | 0.016 | 0.900 |
| GSE222430 | Astrocyte | 0.182 | 0.144 | 0.100 |
| GSE222430 | Endothelial | 0.018 | 0.016 | 0.800 |
| GSE222430 | Inhibitory_neuron | 0.071 | 0.088 | 0.100 |
| GSE222430 | Microglia | 0.518 | 0.520 | 0.800 |
| GSE222430 | OPC | 0.024 | 0.027 | 0.300 |
| GSE222430 | Oligodendrocyte | 0.185 | 0.201 | 0.600 |
| GSE222430 | Pericyte | 0.004 | 0.004 | 1.000 |

### 3.6 Targeted ligand–receptor panel

Of 53 evaluable pairs, 0 passed FDR<0.05. Top pairs by FDR are listed below; the full panel is in Table 6. AIFM1, KEAP1 and PGAM5 are intracellular and absent from ligand–receptor repositories, so any communication difference is an indirect downstream observation (Fig. 6).

| Cohort | Pair | Delta | Permutation p | FDR |
|---|---|---|---|---|
| GSE267933 | Vegfa-Flt1 | -0.896 | 0.100 | 0.663 |
| GSE267933 | Angpt1-Tek | 0.038 | 0.200 | 0.663 |
| GSE267933 | Lgals3-Trem2 | -1.473 | 0.200 | 0.663 |
| GSE267933 | Dll1-Notch1 | -0.082 | 0.200 | 0.663 |
| GSE267933 | Il6-Il6ra | -0.025 | 0.200 | 0.663 |
| GSE222430 | Igf1-Igf1r | -0.467 | 0.100 | 0.663 |
| GSE222430 | Il10-Il10ra | 0.000 | 0.100 | 0.663 |
| GSE222430 | Il6-Il6ra | -0.001 | 0.200 | 0.663 |
| GSE222430 | S100b-Ager | 0.007 | 0.100 | 0.663 |
| GSE222430 | Notch1-Jag1 | -0.066 | 0.200 | 0.663 |

## 4. Discussion

This reanalysis returns a predominantly null result at the transcript level, which is scientifically informative given that a previous reanalysis of overlapping public data also failed to find a reproducible death-signature signal [5]. The result is consistent with the biology of oxeiptosis: the decisive event is dephosphorylation of AIFM1(S116) by PGAM5, a post-translational modification that is not necessarily accompanied by changes in mRNA abundance. A three-gene transcriptional score therefore has low sensitivity for the pathway switch. Additionally, the KEAP1 "biphasic" behavior described in the literature is a functional switch between NRF2-dependent survival and PGAM5-dependent death, not a biphasic mRNA expression pattern; the public data contain a single 24 h point and cannot address temporal dynamics.

The two cohorts differ in exposure: GSE267933 combines sevoflurane with laparotomy, whereas GSE222430 is anesthetic-only. Their concordance (or lack of it) is therefore informative about the contribution of surgical stress, but neither design isolates 6/24/72 h kinetics. We deliberately avoided pooling them: doing so would inflation replication while introducing heterogeneity, a form of pseudo-replication at the cohort level.

The main strength is that inference was performed once per biological sample, with the permutation floor stated explicitly. The main weakness is statistical power: with three samples per group the smallest attainable two-sided p is 0.1, so absence of significance cannot be interpreted as absence of effect. Bootstrap intervals over three observations are unstable and are reported only as uncertainty indicators. Cell-type annotation from marker panels, rather than a reference atlas, introduces classification uncertainty that was summarized by the top1–top2 margin (Table 3).

The next step is therefore not more of the same scoring but a prospective experiment: aged mice assigned to air, oxygen, sevoflurane and sevoflurane-plus-surgery, sampled at 0/6/12/24/72 h, with cell-type-resolved p-AIFM1(S116)/total AIFM1 immunofluorescence, phosphatase controls, and a rescue matrix. Commercial phospho-S116 antibodies currently list human reactivity only; cross-species validation must precede any mouse endpoint. If antiserum validation fails, quantitative phosphoproteomics should replace immunofluorescence.

Three interpretive points follow from the design. First, GSE267933 combines sevoflurane with laparotomy, whereas GSE222430 compares sevoflurane with air; agreement across them would support an anesthetic-associated signal, whereas divergence would point to surgical stress as a modifier. Because neither design includes a sevoflurane-only versus surgery-only factorial contrast, the present analysis cannot attribute the signal to a single component. Second, oxeiptosis is initiated by a post-translational switch: high ROS promotes KEAP1 dissociation from PGAM5, PGAM5 dephosphorylates AIFM1 at Ser116, and the cell dies without nuclear AIFM1 translocation. None of these events is required to change the abundance of the corresponding transcripts, so a three-gene mRNA score is expected to have poor sensitivity even when the pathway is active. Third, the apparent "biphasic KEAP1" behavior in the literature is a functional switch between NRF2-dependent survival at lower ROS and PGAM5-dependent death at higher ROS, not a biphasic expression profile; a single 24 h transcriptional snapshot cannot capture it.

Accordingly, the main contribution here is not a new positive biological claim but an auditable framework: locked sets before inspection, one observation per biological sample, an explicit permutation resolution floor, and a complete map from every manuscript number to an analysis product. The negative result should be read together with the earlier null reanalysis of overlapping public data [5]; both indicate that public transcriptomes alone are unlikely to settle this mechanism. The decisive experiment must measure p-AIFM1(S116) in the cell types prioritized in Section 7, using antibodies validated in mouse tissue and a rescue matrix that separates oxeiptosis from ferroptosis, apoptosis and necroptosis.

## 5. Data and code availability

All source matrices are from GEO accessions GSE267933 and GSE222430. Analysis outputs, checksums and the complete script set are archived in a versioned public repository under https://github.com/yyx-4113/oxeiptosis-perioperative-brain-scrna (release v1.0.0). Every number in this manuscript is mapped to a source file in `07_manuscript/manuscript_numbers.json`; file integrity is recorded in `08_audit/MANIFEST.csv`.

## 6. Author contributions, funding and competing interests

Yongxin Yang designed the study, performed all analyses and wrote the manuscript. No specific funding was received. The author declares no competing interests.

## 7. Predictions, confidence and validation

Every prediction below is a computational hypothesis generated from public transcriptomes, not a demonstration of oxeiptosis. Confidence is deliberately low because the exact permutation floor is p=0.1 and no FDR threshold was met.

| Cohort / cell type | Prediction | Confidence | Uncertainty | Validation experiment |
|---|---|---|---|---|
| GSE222430 / Oligodendrocyte | oxeiptosis-core score higher in exposed samples (delta=-0.040) | low (nominal p=0.200, FDR=0.780) | Cliff's delta -0.778 [-1.000, -0.111]; n=3 per group | cell-type-resolved p-AIFM1(S116)/total IF with λ-phosphatase and H2O2 controls, plus LFHP-1c/NAC rescue |
| GSE267933 / Oligodendrocyte | oxeiptosis-core score higher in exposed samples (delta=0.032) | low (nominal p=0.200, FDR=0.780) | Cliff's delta 0.778 [0.111, 1.000]; n=3 per group | cell-type-resolved p-AIFM1(S116)/total IF with λ-phosphatase and H2O2 controls, plus LFHP-1c/NAC rescue |
| GSE222430 / Microglia | oxeiptosis-core score higher in exposed samples (delta=-0.027) | low (nominal p=0.200, FDR=0.780) | Cliff's delta -0.778 [-1.000, -0.111]; n=3 per group | cell-type-resolved p-AIFM1(S116)/total IF with λ-phosphatase and H2O2 controls, plus LFHP-1c/NAC rescue |

## 8. References

1. Holze C, Michaudel C, Mackowiak C, et al. Oxeiptosis, a ROS-induced caspase-independent apoptosis-like cell-death pathway. Nat Immunol. 2018;19(2):130-140. PMID: 29255269.
2. Galluzzi L, Bravo-San Pedro JM, Kroemer G. Oxeiptosis: a new ROS-dependent form of regulated cell death. Cell Res. 2019;29(2):100-101.
3. Suo Z, et al. Single-nucleus RNA sequencing of hippocampus from aged mice after sevoflurane anesthesia and surgery. Aging Cell. 2025;24(3):e14406. PMID: 39540334. GSE267933.
4. Chang E, et al. iScience. 2023;26(4):106534. PMID: 37123239. GSE196239.
5. Zhang Z, Hu C, Liu Y, et al. Whole-transcriptome sequencing identifies neuroinflammation signatures after sevoflurane exposure: a reanalysis. Front Cell Dev Biol. PMID: 42871048; PMCID: PMC13658196.
6. Gao C, et al. LFHP-1c protects against brain injury via inhibition of PGAM5. Acta Pharm Sin B. 2021;11(7):1867-1884. PMID: 34386325.
7. Suo Z, et al. CNS Neurosci Ther. 2022;28:1576-1595. PMID: 35899365; Dryad doi:10.5061/dryad.jsxksn0cj.
8. Zhang X, et al. Neurosci Bull. 2026. PMID: 41697573; doi:10.1007/s12264-025-01585-z.
9. Wolf FA, Angerer P, Theis FJ. SCANPY: large-scale single-cell gene expression data analysis. Genome Biol. 2018;19:15. PMID: 29409532.
10. Traag VA, Waltman L, van Eck NJ. From Louvain to Leiden: guaranteeing well-connected communities. Sci Rep. 2019;9:5233. PMID: 30944313.
11. Edgar R, Domrachev M, Lash AE. Gene Expression Omnibus: NCBI gene expression and hybridization array data repository. Nucleic Acids Res. 2002;30:207-210. PMID: 11752295.

---

### Figure legends

- **Fig. 2.** UMAP of post-QC cells colored by annotated cell type (left) and condition (right) for each inference cohort.
- **Fig. 3.** Per-sample module scores for the oxeiptosis core and ferroptosis-net signatures; each point is one biological sample after equal-depth cell sampling.
- **Fig. 4.** Cliff's delta with bootstrap 95% intervals for each pre-specified module, cell type and cohort.
- **Fig. 5.** Oxeiptosis-core delta versus ferroptosis-net delta by cell type.
- **Fig. 6.** Targeted ligand–receptor panel deltas by cohort.
- **Fig. S1.** Per-sample QC metrics.

### Limitations and predictive status

All computational results are predictions generated from public transcriptomes. They do not demonstrate oxeiptosis, its cell-type localization at protein level, or any causal effect of sevoflurane. Confirmation requires the wet-experimental chain described in Section 2.5, with pre-specified antibodies, controls and uncertainty reporting.
