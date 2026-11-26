"""Central configuration for the oxeiptosis scRNA-seq reanalysis.

All gene sets, markers, QC thresholds and statistical constants are locked here
before any result is inspected.  No result-derived adjustment is permitted
without recording an amendment in 08_audit/amendments.md.
"""
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
RAW = ROOT / "01_data" / "raw"
PROC = ROOT / "01_data" / "processed"
META = ROOT / "01_data" / "meta"
QC = ROOT / "02_qc"
ANN = ROOT / "03_annotation"
SCORE = ROOT / "04_scoring"
STATS = ROOT / "05_stats"
FIG = ROOT / "06_figures"
MAN = ROOT / "07_manuscript"
AUDIT = ROOT / "08_audit"
for _d in (RAW, PROC, META, QC, ANN, SCORE, STATS, FIG, MAN, AUDIT):
    _d.mkdir(parents=True, exist_ok=True)

# ---------------------------------------------------------------- datasets
# tier: INFERENCE = n>=3 biological replicates, eligible for sample-level tests
#       LOCALIZATION_ONLY = no biological replication, descriptive only
DATASETS = {
    "GSE267933": {
        "tier": "INFERENCE",
        "species": "mouse",
        "tissue": "hippocampus",
        "group_col": "group",
        "groups": {"GSM8281758_Hippocampus_C1": "Control",
                   "GSM8281759_Hippocampus_C2": "Control",
                   "GSM8281760_Hippocampus_C3": "Control",
                   "GSM8281761_Hippocampus_S1": "Sevo+Surgery",
                   "GSM8281762_Hippocampus_S2": "Sevo+Surgery",
                   "GSM8281763_Hippocampus_S3": "Sevo+Surgery"},
        "exposure": "2.5% sevoflurane + 50% O2, 30 min + laparotomy; 24 h",
        "caveat": "drug effect not separable from surgical stress",
    },
    "GSE222430": {
        "tier": "INFERENCE",
        "species": "mouse",
        "tissue": "hippocampus",
        "group_col": "group",
        "groups": {"GSM6923198_mCtrl-1": "Control", "GSM6923199_mCtrl-2": "Control",
                   "GSM6923200_mCtrl-3": "Control", "GSM6923201_mSEV-1": "Sevoflurane",
                   "GSM6923202_mSEV-2": "Sevoflurane", "GSM6923203_mSEV-3": "Sevoflurane"},
        "exposure": "2% sevoflurane 6 h vs air control; hippocampus harvested 24 h later",
        "caveat": "24 h single time point",
    },
    "GSE153598": {"tier": "LOCALIZATION_ONLY", "species": "rhesus macaque", "tissue": "hippocampus",
                  "groups": {}, "exposure": "P7/P14/P28 4 h sevoflurane; 4 libraries per pooled tissue",
                  "caveat": "n=1 animal per condition; libraries are technical splits"},
    "GSE196239": {"tier": "LOCALIZATION_ONLY", "species": "human", "tissue": "prefrontal cortex primary culture",
                  "groups": {}, "exposure": "3.3% sevoflurane 6 h / 2 h recovery; one library per condition",
                  "caveat": "one library per condition; descriptive only"},
    "GSE141242": {"tier": "BULK_REFERENCE", "species": "rat", "tissue": "hippocampus",
                  "groups": {}, "exposure": "2.5% sevoflurane 4 h vs 100% O2; n=3 microarray",
                  "caveat": "bulk array, no cell-type resolution"},
}

# ------------------------------------------------------------- locked sets
# Core triad: Holze et al. Nat Immunol 2018 (PMID 29255269) and Galluzzi 2019.
# AIRE deliberately excluded: no oxeiptosis literature support, not detected in brain (HPA).
GENE_SETS = {
    "oxeiptosis_core": ["Aifm1", "Keap1", "Pgam5"],
    # Hmox1 REMOVED from the NRF2 module in Round 2: it was also a member of
    # the ferroptosis-suppressor module, making the two endpoints non-independent.
    # Retained in the ferroptosis-suppressor module (its canonical death-modifier role).
    # See 08_audit/amendments.md entry R2-5.
    "nrf2_antagonistic": ["Nfe2l2", "Nqo1", "Gclc", "Gclm"],
    "ferroptosis_driver": ["Acsl4", "Lpcat3", "Pebp1", "Ncoa4",
                           "Tfrc", "Alox15", "Sat1", "Elavl1"],
    "ferroptosis_suppressor": ["Gpx4", "Slc7a11", "Fth1", "Ftl1", "Aifm2",
                               "Hmox1", "Cisd1", "Prdx6"],
}
# Provenance for the ferroptosis lists: FerrDb v2 (Zhou N et al. 2023, Nucleic
# Acids Res 51:D571-D582, PMID 36305834) driver/suppressor annotations, restricted to genes
# detectable in these datasets.  Elavl1/Pebp1/Sat1 are context-dependent and are
# retained only as pre-locked members; see Limitations.

# Cell identity marker panels (mouse symbols).  Scored per cell with
# scanpy.tl.score_genes; cluster assignment = highest mean panel score.
CELLTYPE_MARKERS = {
    # Hippocampal subtype markers added in Round 2 (DG granule: Prox1/Calb1;
    # CA pyramidal: Neurod6/Dkk3/Fibcd1).  Itgam removed from Microglia
    # (also labels infiltrating monocytes).  Tbr1 retained but is not a
    # hippocampus-specific marker; see Limitations.
    "Excitatory_neuron": ["Slc17a7", "Slc17a6", "Neurod2", "Prox1", "Calb1",
                          "Neurod6", "Dkk3", "Fibcd1"],
    "Inhibitory_neuron": ["Gad1", "Gad2", "Slc6a1"],
    "Neuron_pan": ["Rbfox3", "Snap25", "Syt1", "Nefl"],
    "Astrocyte": ["Aqp4", "Gfap", "Slc1a2", "Aldh1l1", "Slc1a3"],
    "Microglia": ["Cx3cr1", "P2ry12", "Tmem119", "C1qa", "Inpp5d"],
    "Oligodendrocyte": ["Mbp", "Plp1", "Mog", "Mal"],
    "OPC": ["Pdgfra", "Cspg4", "Sox10"],
    "Endothelial": ["Pecam1", "Cldn5", "Flt1"],
    "Pericyte": ["Pdgfrb", "Rgs5", "Acta2"],
}

# Targeted ligand-receptor panel (mouse nomenclature, CellChatDB-style pairs).
# Explicitly a targeted panel, NOT a full CellChatDB run.
LR_PANEL = [
    ("Cx3cl1", "Cx3cr1"), ("Csf1", "Csf1r"), ("Il1b", "Il1r1"),
    ("Tnf", "Tnfrsf1a"), ("Tgfb1", "Tgfbr1"), ("Pdgfb", "Pdgfrb"),
    ("Vegfa", "Flt1"), ("Angpt1", "Tek"), ("Sema3a", "Nrp1"),
    ("Ephb2", "Efnb1"), ("Apoe", "Ldlr"), ("Apoe", "Lrp8"),
    ("C3", "C3ar1"), ("Ccl2", "Ccr2"), ("Cxcl10", "Cxcr3"),
    ("Mif", "Cd74"), ("Hmgb1", "Tlr4"), ("S100b", "Ager"),
    ("Notch1", "Jag1"), ("Dll1", "Notch1"), ("Faslg", "Tnfrsf6"),
    ("Tnfsf10", "Tnfrsf10b"), ("Lgals3", "Trem2"), ("Bdnf", "Ntrk2"),
    ("Ntf3", "Ntrk3"), ("Gdnf", "Ret"), ("Fgf2", "Fgfr1"),
    ("Igf1", "Igf1r"), ("Cxcl12", "Cxcr4"), ("Il6", "Il6ra"),
    ("Il10", "Il10ra"), ("Ccl3", "Ccr1"),
]

# --------------------------------------------------------------------- QC
QC_THRESHOLDS = {
    "min_genes_per_cell": 200,
    "min_counts_per_cell": 500,
    "max_mt_pct": 15.0,
    "min_cells_per_gene": 3,
    "n_hvg": 2000,
    "n_pcs": 30,
    "resolution": 0.8,
    "min_cells_per_celltype_per_sample": 20,   # below this the sample is NA for that cell type
    "downsample_per_sample_celltype": 100,     # equal-depth sensitivity analysis
    "max_cells_per_sample": 20000,             # computational cap, random seed fixed; disclosed in manuscript
    "min_annotation_margin": 0.05,             # clusters below this are LowConfidence and excluded from inference
}

# ------------------------------------------------------------------ stats
STATS_CONFIG = {
    "permutation": "exact",       # C(6,3)=20 labelings, two-sided min p = 0.10
    "alpha": 0.05,
    # Three families per locked protocol §7.3 (corrected in Round 2; see amendments.md):
    # A = ox_core x celltype x dataset (primary);
    # B = secondary modules + separation D, per dataset;
    # C = targeted LR panel, per dataset.
    "fdr_family": "A:ox_core(celltype x dataset); B1:secondary modules(per dataset); B2:separation D(per dataset); C:LR(per dataset)",
    "bootstrap_reps": 10000,
    "seed": 20261125,
    "effect": "Cliff_delta_primary_Hedges_g_sensitivity",
    "delta_statistic": "median",  # Round 2: protocol-specified median difference
}

GTF_MT_PREFIX = ("mt-", "MT-")
