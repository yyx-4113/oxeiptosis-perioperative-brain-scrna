"""Organize downloaded 10x triplets into scanpy-readable per-sample folders.

Writes 01_data/meta/sample_inventory.csv with sha256 of each source archive/matrix.
No biological inference in this script.
"""
import csv
import gzip
import hashlib
import re
import shutil
from pathlib import Path

import numpy as np
from scipy.io import mmread

from common import RAW, PROC, META, DATASETS

INFERENCE = [d for d, v in DATASETS.items() if v["tier"] == "INFERENCE"]
rows = []

for gse in INFERENCE:
    raw_dir = RAW / gse
    files = sorted(raw_dir.glob("*.gz"))
    if not files:
        print(gse, "SKIP: no extracted files yet")
        continue
    suffixes = {
        "barcodes": ("_barcodes.tsv.gz", ".barcodes.tsv.gz"),
        "features": ("_features.tsv.gz", ".features.tsv.gz"),
        "matrix": ("_matrix.mtx.gz", ".matrix.mtx.gz"),
    }
    parsed = {}
    for f in files:
        for kind, sufs in suffixes.items():
            for suf in sufs:
                if f.name.endswith(suf):
                    parsed.setdefault(f.name[:-len(suf)], {})[kind] = f
                    break
    samples = sorted(parsed)
    print(gse, "detected prefixes:", samples)
    for s in samples:
        src = {kind: parsed[s].get(kind) for kind in ("barcodes", "features", "matrix")}
        if any(v is None for v in src.values()):
            print("  SKIP incomplete triplet:", s, src)
            continue
        dst_dir = PROC / gse / s
        dst_dir.mkdir(parents=True, exist_ok=True)
        targets = {
            "barcodes": dst_dir / "barcodes.tsv.gz",
            "features": dst_dir / "features.tsv.gz",
            "matrix": dst_dir / "matrix.mtx.gz",
        }
        for kind, dst in targets.items():
            if not dst.exists():
                shutil.copy2(src[kind], dst)
        with gzip.open(targets["barcodes"], "rt") as fh:
            n_bar = sum(1 for _ in fh)
        with gzip.open(targets["features"], "rt") as fh:
            n_feat = sum(1 for _ in fh)
        with gzip.open(targets["matrix"], "rb") as fh:
            mat = mmread(fh).tocsr()
        h = hashlib.sha256()
        with open(targets["matrix"], "rb") as fh:
            for chunk in iter(lambda: fh.read(1 << 20), b""):
                h.update(chunk)
        grp = DATASETS[gse]["groups"].get(s, "UNMAPPED")
        rows.append({
            "dataset": gse, "sample": s, "group": grp, "tier": DATASETS[gse]["tier"],
            "n_barcodes": n_bar, "n_features": n_feat,
            "matrix_shape_0": mat.shape[0], "matrix_shape_1": mat.shape[1],
            "nnz": int(mat.nnz),
            "matrix_sha256": h.hexdigest(),
            "source_matrix": src["matrix"].name,
            "exposure": DATASETS[gse]["exposure"],
        })
        print(f"  {s:12s} {grp:14s} cells={n_bar:6d} genes={n_feat:6d} "
              f"shape={mat.shape} nnz={mat.nnz}")

if not rows:
    raise SystemExit("No samples organized; run the download/extract steps first.")
out = META / "sample_inventory.csv"
with open(out, "w", newline="", encoding="utf-8-sig") as fh:
    w = csv.DictWriter(fh, fieldnames=list(rows[0].keys()))
    w.writeheader()
    w.writerows(rows)
print("WROTE", out, "rows=", len(rows))
