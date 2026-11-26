"""Record the execution environment.  Run first; writes 00_environment/env.json."""
import importlib.metadata
import json
import platform
import shutil
import subprocess
import sys
from datetime import datetime, timezone

from common import ROOT, DATASETS

PKGS = ["numpy", "scipy", "pandas", "anndata", "scanpy", "sklearn",
        "statsmodels", "matplotlib", "seaborn", "h5py", "requests",
        "leidenalg", "igraph"]

out = {
    "timestamp_utc": datetime.now(timezone.utc).isoformat(),
    "python": sys.version.split()[0],
    "platform": platform.platform(),
    "executable": sys.executable,
    "R_in_PATH": shutil.which("R") is not None,
    "packages": {},
}
for name in PKGS:
    try:
        mod = __import__(name)
        try:
            out["packages"][name] = importlib.metadata.version(name)
        except Exception:
            out["packages"][name] = getattr(mod, "__version__", "unknown")
    except Exception as exc:
        out["packages"][name] = f"UNAVAILABLE ({type(exc).__name__})"

out["datasets_declared"] = {
    k: {"tier": v["tier"], "species": v["species"], "exposure": v["exposure"],
        "caveat": v["caveat"]}
    for k, v in DATASETS.items()
}
p = ROOT / "00_environment" / "env.json"
p.parent.mkdir(parents=True, exist_ok=True)
p.write_text(json.dumps(out, indent=2, ensure_ascii=False), encoding="utf-8")
print(json.dumps({"written": str(p), "python": out["python"],
                  "missing": [k for k, v in out["packages"].items()
                              if str(v).startswith("UNAVAILABLE")]}, ensure_ascii=False))
