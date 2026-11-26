# Pre-analysis plan amendments and deviations log

This file records every deviation between the locked protocol
(`output/oxeiptosis-scrna-protocol-20261125/stage3/oxeiptosis-localization-perioperative-brain-injury_v1.0.docx`)
and the executed analysis. Entries are written at the time of revision.

## A. Protocol-specified items not executed as written (environment-driven)

| # | Protocol | Executed | Reason | Consequence / mitigation |
|---|---|---|---|---|
| A1 | AUCell primary scoring; ssGSEA and mean-expression Z-score as sensitivity | scanpy `score_genes` (control-gene subtraction, ctrl_size=50, n_bins=25); plain mean expression as sensitivity | AUCell requires additional dependency and ranking infrastructure; `score_genes` is a deterministic control-gene-subtracted score with fixed seed | Disclosed in Methods; results are explicitly exploratory; sensitivity score retained |
| A2 | Downsampling sensitivity at 100/200/400 cells per sample×cell type | 100 cells only | n=3 per group; additional depths do not change the permutation resolution (floor p=0.1) | Disclosed; labelled as partial sensitivity |
| A3 | QC: genes 200–6000, MT%<20, DoubletFinder/Scrublet | genes ≥200 (no upper bound), MT%<15, no doublet caller | Doublet callers unavailable in environment; upper gene bound risk assessed as low in these 10x libraries | Disclosed in Methods and Limitations; low-confidence annotation category added as partial mitigation |
| A4 | CellChat pathway-level ligand–receptor testing | targeted 32-pair ligand×receptor product panel (27–26 evaluable) | Core triad proteins are intracellular and absent from CellChatDB; targeted panel is transparent and sample-level | Disclosed in Methods; results labelled indirect |
| A5 | Bootstrap B≥10000 | 10000 (revised in Round 2; originally 2000) | Round-1 peer review identified B=2000 as a protocol deviation | Corrected before resubmission; see entry R2 |
| A6 | Delta defined as median difference | mean difference (revised in Round 2 to median difference) | Round-1 review identified deviation | Corrected before resubmission; see entry R2 |
| A7 | Three FDR families (A, B, C) in protocol §7.3 | Four families: A, B1 (secondary modules), B2 (separation endpoint D), C | B was split during implementation because the separation endpoint was added after the protocol was locked; the split was written into the executed code before any test statistic was computed | Declared in Methods and Discussion as an amendment, **not** as pre-registration; see section B |

## B. Statistical analysis plan — family structure (corrected in Round 2)

The locked protocol §7.3 specified three FDR families. The Round-1 implementation
used one pooled 42-test family plus separate composition and LR families. This
was identified as a Major concern by the independent statistical reviewer and is
corrected here **before** re-submission:

- **Family A (primary)**: `ox_core` only, tested across cell types × datasets.
- **Family B (secondary)**: split into two sub-families, each BH-corrected
  **within each dataset** (this split is amendment **A7**):
  - **B1** = `nrf2_antag`, `ferro_driver`, `ferro_suppressor`, `ferro_net`, `ox_core_meanexpr`;
  - **B2** = the separation endpoint D alone.
  The Round-2 reviewer noted that pooling B1 and B2 into one BH family would change
  individual q values (though no q crosses 0.05 either way).  Because keeping them
  separate is the more conservative option for the separation endpoint (smaller family
  size gives smaller q), and because this structure was already fixed in the executed
  code before round-2 inspection, it is retained and explicitly declared here rather
  than silently merged.  Manuscript text uses "families B1 and B2" consistently.
- **Family C (communication)**: targeted LR panel, BH-corrected **within each dataset**.
- Composition proportions are reported as exploratory descriptive tests in a
  separate explicitly labelled family (not part of the primary inferential claim).

Two different facts must not be merged. (i) The **Round-1** correction, which
replaced one pooled 42-test family with families A/B/C, was made **after** seeing
Round-1 results and is therefore a post hoc refinement. (ii) The **A7 split of B
into B1/B2** was written into the executed code before any test statistic was
computed, but after the protocol was locked; it is an amendment, not a
pre-registered analysis. The manuscript states exactly this and never describes
the four-family structure as pre-registered.

## C. Round-2 revisions triggered by peer review

| # | Reviewer concern | Action | Timestamp |
|---|---|---|---|
| R2-1 | Statistics M2: FDR family mismatch | Re-implement three families as described in B; all q values regenerated from result CSVs | 2026-11-25 |
| R2-2 | Statistics M3 / Biology M3: separation D had no p or CI | Compute per-sample paired d_i = ox_core_i − ferro_net_i (standardised within sample to remove scale artifact), exact permutation p and bootstrap 95% CI; report effect-size-scale D as sensitivity | 2026-11-25 |
| R2-3 | Statistics M5: B=2000 and pathological Hedges' g CI | B=10000; Hedges' g reported only as point estimate with bootstrap-valid-replicate proportion; Cliff's delta is primary effect size; percentile CI reported only when ≥95% of replicates are valid | 2026-11-25 |
| R2-4 | Biology M2: low-margin annotation | Add margin threshold; clusters below threshold labelled `LowConfidence` and excluded from cell-type-specific inference; add hippocampal subtype markers; implement Neuron_pan validation | 2026-11-25 |
| R2-5 | Biology M5: NRF2 direction not interpreted | Add explicit reverse-interpretation rule in Results/Discussion; remove Hmox1 overlap between NRF2 and ferroptosis-suppressor modules | 2026-11-25 |
| R2-6 | Biology M6 / M7 / M8; Editor M6 | Correct 50% O2 control discussion, terminology, antibody fact, three reference titles/PMIDs, add ethics statement | 2026-11-25 |
| R2-7 | Editor M2/M4 | Honest data-availability statement (no URL for a non-existent repository); complete number-to-product map | 2026-11-25 |

## D. Items deliberately not changed (with justification)

- **No pooling of two cohorts**: explicitly prohibited by the locked protocol because
  GSE267933 contains surgical stress and GSE222430 does not; pooling would create
  heterogeneous exposure and pseudo-replication.
- **AIRE remains excluded**: no oxeiptosis literature supports inclusion; Human Protein
  Atlas reports no brain expression.
- **No cross-species inference**: GSE153598 (macaque) and GSE196239 (human iPSC-derived
  culture) remain descriptive only.
- **No claim of pathway activation from transcription**: p-AIFM1(S116) immunofluorescence
  remains the proposed confirmatory endpoint.

## E. Round-3 revisions (2026-11-26)

Three independent reviewers (statistics, biology, editor) reviewed the Round-2
manuscript and all returned "major revision". Actions:

| # | Reviewer concern | Action | Timestamp |
|---|---|---|---|
| R3-1 | Hedges' g intervals were reported in Table 4 while the text said they were suppressed | Interval suppressed in the code for every row (`hedges_g_CI95_lo/hi` set to NA); text and Fig. 4 legend now state point estimate only; audit check added | 2026-11-26 |
| R3-2 | Four FDR families presented as pre-registered although the plan declared three | Amendment A7 added; Methods, Abstract and Discussion now state that A/B/C were declared and B was split before any statistic was computed; `stats_metadata.json` records both facts | 2026-11-26 |
| R3-3 | "0 of 66 FDR<0.05" framed as a negative finding although it is arithmetically forced by the p floor | Results, Discussion and Conclusion rewritten: the design cannot separate absence of effect from insufficient resolution | 2026-11-26 |
| R3-4 | Differential annotation dropout between arms in GSE267933 not disclosed and untested | Per-sample dropout table added (Table S4); fallback-label sensitivity analysis added (`scripts/25_sensitivity.py`, `05_stats/sensitivity_lowconf.csv`) | 2026-11-26 |
| R3-5 | "9 of 11" compared raw score magnitudes across two non-comparable scales | Both counts now reported with their definitions: raw 9/11, scale-free 4/11; neither supports a pathway-specific claim | 2026-11-26 |
| R3-6 | Separation D described as computed across all cells and as a paired test | Methods corrected to the executed code: z across sample×cell-type values, within-sample subtraction, unpaired permutation; scale-free contrast reported | 2026-11-26 |
| R3-7 | §3.3 table selected on nominal p with undeclared tie-breaking and omitted the largest shift | Ties broken by absolute delta; top-6 now contains the largest ox_core shift | 2026-11-26 |
| R3-8 | Title claims "localization" without per-cell-type transcript evidence; excitatory neurons essentially absent | Table S3 added (per cell type detection and mean expression of Aifm1/Keap1/Pgam5); per-cohort class composition and the missing excitatory class now stated | 2026-11-26 |
| R3-9 | Nine of 13 references never cited in text; Limitations placed after References; no Conclusion or Ethics heading; Figs. 1/3/4/S1 uncited | Citations inserted; sections reordered (Discussion → Limitations → Conclusion → Declarations → Data availability → References); figure citations added | 2026-11-26 |
| R3-10 | Data-availability URL returned 404 at review time | Repository created and pushed, tag v1.0.0, then v1.0.1 after these revisions; URL and release verified anonymous HTTP 200 | 2026-11-26 |
| R3-11 | Manifest included the deposit directory, so local and repository manifests diverged | `repro_repo`, `.git` and `__pycache__` excluded from manifest generation; repository manifest regenerated over its own file set | 2026-11-26 |
