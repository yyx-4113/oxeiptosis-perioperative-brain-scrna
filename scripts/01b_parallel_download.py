"""Parallel HTTP range downloader for large GEO supplementary TAR files.

Usage: python 01b_parallel_download.py GSE267933
       python 01b_parallel_download.py GSE222430
Writes <raw>/<GSE>_RAW.tar.part then atomically renames to .tar.
"""
import os
import random
import sys
import time
import urllib.error
import urllib.request
from concurrent.futures import ThreadPoolExecutor, as_completed
from pathlib import Path

from common import RAW

CHUNK = 2 * 1024 * 1024
WORKERS = 8


def fetch(url, start, end, dest, idx):
    req = urllib.request.Request(url, headers={"Range": f"bytes={start}-{end}",
                                                "User-Agent": "Mozilla/5.0"})
    for attempt in range(10):
        try:
            with urllib.request.urlopen(req, timeout=180) as resp:
                if resp.status != 206:
                    raise IOError(f"unexpected status {resp.status}")
                data = resp.read()
            if len(data) != end - start + 1:
                raise IOError(f"short read {len(data)} vs {end-start+1}")
            with open(dest, "r+b") as fh:
                fh.seek(start)
                fh.write(data)
            return idx, len(data)
        except Exception as exc:
            if attempt == 9:
                raise
            time.sleep(min(30, 2 ** attempt) + random.random() * 2)
    return idx, 0


def main(gse):
    pre = gse[:6] + "nnn"
    url = f"https://ftp.ncbi.nlm.nih.gov/geo/series/{pre}/{gse}/suppl/{gse}_RAW.tar"
    with urllib.request.urlopen(urllib.request.Request(url, method="HEAD"), timeout=60) as r:
        size = int(r.headers["Content-Length"])
        assert r.headers.get("Accept-Ranges") == "bytes", "server does not support ranges"
    out = RAW / gse / f"{gse}_RAW.tar"
    part = RAW / gse / f"{gse}_RAW.tar.part"
    out.parent.mkdir(parents=True, exist_ok=True)
    with open(part, "wb") as fh:
        fh.truncate(size)
    ranges = [(s, min(s + CHUNK - 1, size - 1)) for s in range(0, size, CHUNK)]
    print(f"{gse}: {size} bytes, {len(ranges)} chunks, {WORKERS} workers", flush=True)
    t0 = time.time()
    done = 0
    with ThreadPoolExecutor(max_workers=WORKERS) as ex:
        futs = [ex.submit(fetch, url, s, e, part, i) for i, (s, e) in enumerate(ranges)]
        for fut in as_completed(futs):
            idx, n = fut.result()
            done += n
            if (idx + 1) % 10 == 0:
                print(f"  {idx+1}/{len(ranges)} chunks, {done/1e6:.1f} MB, "
                      f"{time.time()-t0:.0f}s", flush=True)
    part.replace(out)
    print(f"DONE {out} {out.stat().st_size} bytes in {time.time()-t0:.0f}s", flush=True)


if __name__ == "__main__":
    main(sys.argv[1])
