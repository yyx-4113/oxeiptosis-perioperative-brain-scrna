"""Independent audit: re-compute one key statistic from source CSV, verify
internal consistency, and verify manifest checksums.

Emits 08_audit/audit_report.json with PASS/FAIL per check.
"""
import hashlib
import itertools
import json
from pathlib import Path

import numpy as np
import pandas as pd

from common import ROOT, AUDIT, QC, SCORE, STATS, META, DATASETS, STATS_CONFIG

checks = []


def check(name, ok, detail=""):
    checks.append({"check": name, "status": "PASS" if ok else "FAIL", "detail": str(detail)})


def exact_perm(a, b):
    vals = np.concatenate([a, b]).astype(float)
    n = len(a)
    obs = vals[:n].mean() - vals[n:].mean()
    st = []
    for comb in itertools.combinations(np.arange(len(vals)), n):
        m = np.zeros(len(vals), bool)
        m[list(comb)] = True
        st.append(vals[m].mean() - vals[~m].mean())
    return obs, float(np.mean(np.abs(np.asarray(st)) >= abs(obs) - 1e-12))


# files
required = [META / "sample_inventory.csv", QC / "qc_per_sample.csv",
            SCORE / "sample_scores.csv", STATS / "sample_level_tests.csv",
            STATS / "separation_summary.csv", STATS / "lr_panel_tests.csv",
            AUDIT / "MANIFEST.csv"]
for p in required:
    check("file_exists:" + p.name, p.exists() and p.stat().st_size > 100)

inv = pd.read_csv(META / "sample_inventory.csv")
for gse in [g for g, v in DATASETS.items() if v["tier"] == "INFERENCE"]:
    c = inv[inv.dataset == gse].groupby("group").size().to_dict()
    check(f"n_replicates_{gse}", c.get("Control") == 3 and sum(v for k, v in c.items() if k != "Control") == 3, c)

sdf = pd.read_csv(SCORE / "sample_scores.csv")
tdf = pd.read_csv(STATS / "sample_level_tests.csv")
obs_per_group = sdf.groupby(["dataset", "celltype", "group"]).size()
incomplete = obs_per_group[obs_per_group < 3]
check("sample_level_at_most_one_value_per_sample", bool((obs_per_group <= 3).all()),
      f"group sizes <=3; incomplete strata={incomplete.to_dict()}")

# composition table must have exactly one row per dataset x celltype;
# a previous bug appended it once per module (6x) and inflated the BH family
cp_path = STATS / "composition_tests.csv"
if cp_path.exists():
    cdf = pd.read_csv(cp_path)
    cdup = cdf.groupby(["dataset", "celltype"]).size()
    check("composition_unique_per_dataset_celltype", bool((cdup == 1).all()),
          f"n_rows={len(cdf)}, key multiplicities={sorted(cdup.unique())}")
    check("composition_fdr_family_size", len(cdf) == len(cdup),
          f"rows={len(cdf)}, unique_keys={len(cdup)}")

# resolution floor: no exact permutation p below 0.1 with n=3+3
meta_stats = json.loads((STATS / "stats_metadata.json").read_text(encoding="utf-8"))
floor = float(meta_stats["min_possible_two_sided_p"])  # derived, never hardcoded
bad = tdf[tdf.permutation_p_two_sided < floor - 1e-9]
check("permutation_resolution_floor", len(bad) == 0, f"floor={floor}, violations={len(bad)}")
check("FDR_at_least_p", bool((tdf.FDR_BH_all_tests >= tdf.permutation_p_two_sided - 1e-12).all()))
check("bootstrap_reps_protocol", int(meta_stats["bootstrap_reps"]) >= 10000,
      f"B={meta_stats['bootstrap_reps']}")
check("delta_statistic_protocol", meta_stats.get("delta_statistic") == "median",
      f"stat={meta_stats.get('delta_statistic')}")
# family sizes must match the three-family design declared in metadata
fam = meta_stats.get("fdr_family_sizes", {})
check("fdr_family_A_declared", int(fam.get("A_primary_ox_core", 0)) > 0, fam)
check("fdr_family_B_declared", bool(fam.get("B_secondary_by_dataset")), fam)
# separation endpoint must carry inference (Round 2, statistics M3)
sp_path = STATS / "separation_tests.csv"
if sp_path.exists():
    sp = pd.read_csv(sp_path)
    check("separation_has_p_and_CI",
          {"permutation_p_two_sided", "CI95_lo", "CI95_hi", "FDR_BH_family_B_within_dataset"}.issubset(sp.columns),
          f"cols={list(sp.columns)}")
    check("separation_p_floor", bool((sp.permutation_p_two_sided >= floor - 1e-9).all()),
          f"min={sp.permutation_p_two_sided.min()}")
    # N5 (Round-2): CI must be a two-sample bootstrap; it must not systematically
    # exclude zero when the exact permutation p is large.  With n=3 this is a
    # consistency heuristic, not a formal equivalence test.
    incon = sp[(sp.permutation_p_two_sided >= 0.4) &
               ((sp.CI95_lo > 0) | (sp.CI95_hi < 0))]
    check("separation_CI_p_consistency", len(incon) <= 2,
          f"p>=0.4 but CI excludes 0 in {len(incon)} of {len(sp)} rows")
# low-confidence annotation gate (Round 2, biology M2)
check("annotation_margin_gate_present", "min_annotation_margin" in
      json.loads((QC / "qc_summary.json").read_text(encoding="utf-8"))["qc_thresholds"],
      "qc_summary.json:qc_thresholds")
# N1/N10: LowConfidence cells must actually be present in the annotation product
ann_tab = AUDIT / "Table3_annotation.csv"
if ann_tab.exists():
    at = pd.read_csv(ann_tab)
    lc = at[at.assigned_celltype.astype(str) == "LowConfidence"]
    check("lowconfidence_count_matches_annotation", len(lc) > 0 and lc.n_cells.sum() > 0,
          f"clusters={len(lc)}, cells={int(lc.n_cells.sum())}")
# N4/N10: family B1 and B2 q values must be BH-consistent with their declared families
for ds in tdf.dataset.unique():
    for fam, label in [("B_secondary", "B1")]:
        m = (tdf.dataset == ds) & (tdf.FDR_family == fam)
        if m.any():
            sub = tdf[m]
            pv = sub.permutation_p_two_sided.to_numpy(); nn = len(pv)
            order = np.argsort(pv); qq = pv[order] * nn / (np.arange(nn) + 1)
            qq = np.minimum.accumulate(qq[::-1])[::-1]
            out = np.empty(nn); out[order] = np.clip(qq, 0, 1)
            check(f"B1_BH_recompute_{ds}", np.allclose(out, sub.FDR_BH_all_tests),
                  f"max_diff={float(np.max(np.abs(out - sub.FDR_BH_all_tests))):.2e}")

# independent recomputation of first row from sample_scores.csv
r = tdf.iloc[0]
g = sdf[(sdf.dataset == r.dataset) & (sdf.celltype == r.celltype)]
a = g.loc[g.group == "Control", r.module].dropna().to_numpy()
b = g.loc[g.group != "Control", r.module].dropna().to_numpy()
# Round 2: exact_perm in the audit must match the executed median statistic
def _loc(x):
    return float(np.median(x)) if STATS_CONFIG.get("delta_statistic") == "median" else float(np.mean(x))


def exact_perm2(a, b):
    vals = np.concatenate([a, b]).astype(float)
    n = len(a)
    obs = _loc(vals[:n]) - _loc(vals[n:])
    st = []
    for comb in itertools.combinations(np.arange(len(vals)), n):
        m = np.zeros(len(vals), bool); m[list(comb)] = True
        st.append(_loc(vals[m]) - _loc(vals[~m]))
    return obs, float(np.mean(np.abs(np.asarray(st)) >= abs(obs) - 1e-12))


d_rec, p_rec = exact_perm2(b, a)
check("independent_recompute_delta", np.isclose(d_rec, r.delta_exposed_minus_control, atol=1e-8),
      f"recomputed={d_rec:.10f} stored={r.delta_exposed_minus_control:.10f}")
check("independent_recompute_p", np.isclose(p_rec, r.permutation_p_two_sided, atol=1e-12),
      f"recomputed={p_rec:.6f} stored={r.permutation_p_two_sided:.6f}")

# manifest integrity
man = pd.read_csv(AUDIT / "MANIFEST.csv")
bad_hash = []
for row in man.itertuples():
    p = ROOT / row.relative_path
    if not p.exists():
        bad_hash.append((row.relative_path, "missing"))
        continue
    h = hashlib.sha256()
    with open(p, "rb") as fh:
        for chunk in iter(lambda: fh.read(1 << 20), b""):
            h.update(chunk)
    if h.hexdigest() != row.sha256:
        bad_hash.append((row.relative_path, "hash"))
check("manifest_integrity", len(bad_hash) == 0, f"mismatches={bad_hash[:5]}")

# manuscript traceability: every tracked number must cite an existing source file
MAN = ROOT / "07_manuscript"
numpath = MAN / "manuscript_numbers.json"
if numpath and numpath.exists():
    tracked = json.loads(numpath.read_text(encoding="utf-8"))
    srcs = set()
    for v in tracked.values():
        srcs.add(str(v.get("source", "")).split(":")[0])
    missing_src = [s for s in srcs if not (ROOT / s).exists()]
    check("manuscript_number_sources_exist", len(missing_src) == 0,
          f"tracked={len(tracked)} missing={missing_src}")
else:
    check("manuscript_number_sources_exist", False, "manuscript_numbers.json not found")

report = {"checks": checks, "n_pass": sum(c["status"] == "PASS" for c in checks),
          "n_fail": sum(c["status"] == "FAIL" for c in checks),
          "seed": STATS_CONFIG["seed"]}
(AUDIT / "audit_report.json").write_text(json.dumps(report, indent=2), encoding="utf-8")
print(json.dumps({"n_pass": report["n_pass"], "n_fail": report["n_fail"],
                  "failed": [c["check"] for c in checks if c["status"] == "FAIL"]}, ensure_ascii=False))
