"""Create a checksum manifest for every analysis product (excluding raw TARs)."""
import csv
import hashlib
from pathlib import Path

from common import ROOT, AUDIT

# Round 3 (editor non-blocking 2): the deposit directory and VCS metadata must
# not enter the manifest, otherwise local and repository manifests diverge.
SKIP_DIRS = {"raw", "repro_repo", ".git", "__pycache__"}
# Round 2 (editor minor 8): include scripts and repository metadata so the
# "complete script set under checksum" statement is accurate.  h5ad objects
# are large and reproducible from locked scripts + source matrices; they are
# deliberately excluded from redistribution (documented in the manuscript).
EXT = {".csv", ".json", ".md", ".txt", ".png", ".py", ".cff", ".yml", ".yaml"}
SKIP_NAMES = {"MANIFEST.csv", "audit_report.json"}
SKIP_EXACT = {"03_annotation/GSE267933_annotated.h5ad",
              "03_annotation/GSE222430_annotated.h5ad"}

rows = []
for p in sorted(ROOT.rglob("*")):
    if not p.is_file() or p.suffix not in EXT:
        continue
    if any(part in SKIP_DIRS for part in p.relative_to(ROOT).parts):
        continue
    if p.name in SKIP_NAMES:
        continue
    if str(p.relative_to(ROOT)).replace("\\", "/") in SKIP_EXACT:
        continue
    h = hashlib.sha256()
    with open(p, "rb") as fh:
        for chunk in iter(lambda: fh.read(1 << 20), b""):
            h.update(chunk)
    rows.append({"relative_path": str(p.relative_to(ROOT)).replace("\\", "/"),
                 "bytes": p.stat().st_size, "sha256": h.hexdigest()})

out = AUDIT / "MANIFEST.csv"
with open(out, "w", newline="", encoding="utf-8-sig") as fh:
    w = csv.DictWriter(fh, fieldnames=["relative_path", "bytes", "sha256"])
    w.writeheader()
    w.writerows(rows)
print("WROTE", out, "files=", len(rows))
