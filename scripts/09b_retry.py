#!/usr/bin/env python3
"""Stage 09b-retry: refetch CADD for cache entries with null scores,
excluding deletions (API systematically returns [] for them — documented).

Appends corrected entries to the main cache (load_cache: later lines win).
Usage: python3 scripts/09b_retry.py [--workers 8]
"""

from __future__ import annotations

import argparse
import gzip
import json
import time
import urllib.request
from concurrent.futures import ThreadPoolExecutor, as_completed
from datetime import date
from pathlib import Path

REPO = Path(__file__).resolve().parents[1]
API = "https://cadd.gs.washington.edu/api/v1.0/GRCh38-v1.7/{chrom}:{pos}_{ref}_{alt}"


def fetch_one(chrom: str, pos: int, ref: str, alt: str):
    url = API.format(chrom=chrom.replace("chr", ""), pos=pos, ref=ref, alt=alt)
    req = urllib.request.Request(url, headers={"User-Agent": "IsoformVUS/0.1.0 (research)"})
    try:
        with urllib.request.urlopen(req, timeout=30) as resp:
            payload = json.load(resp)
    except Exception:
        return None
    if not payload:
        return None
    rec = payload[0] if isinstance(payload, list) else payload
    try:
        return {"phred": float(rec["PHRED"]), "raw": float(rec["RawScore"])}
    except (KeyError, ValueError, TypeError):
        return None


def main() -> int:
    ap = argparse.ArgumentParser()
    ap.add_argument("--workers", type=int, default=8)
    ap.add_argument("--cache", default=str(REPO / "data/interim/cadd_cache.jsonl"))
    args = ap.parse_args()

    cache, order = {}, []
    with open(args.cache) as fh:
        for line in fh:
            try:
                r = json.loads(line)
            except json.JSONDecodeError:
                continue
            if r["variant_key"] not in cache:
                order.append(r["variant_key"])
            cache[r["variant_key"]] = r

    # variant coords for retry candidates (null phred, not deletion)
    coords, vartype = {}, {}
    with gzip.open(REPO / "data/processed/variant_labels.tsv.gz", "rt") as fh:
        fh.readline()
        for line in fh:
            c = line.rstrip("\n").split("\t")
            coords[c[0]] = (c[1], int(c[2]), c[3], c[4])
            vartype[c[0]] = c[5]

    todo = [k for k, v in cache.items()
            if v.get("phred") is None and vartype.get(k) != "deletion" and k in coords]
    print(f"[09b-retry] {len(todo)} null non-deletion variants to refetch", flush=True)

    results = {}
    with ThreadPoolExecutor(max_workers=args.workers) as ex:
        futs = {ex.submit(fetch_one, *coords[k]): k for k in todo}
        for i, fut in enumerate(as_completed(futs)):
            k = futs[fut]
            try:
                results[k] = fut()
            except Exception:
                results[k] = None
            if (i + 1) % 200 == 0:
                print(f"[09b-retry] {i+1}/{len(todo)}", flush=True)

    n_ok = 0
    with open(args.cache, "a") as cf:
        for k in todo:
            rec = results[k]
            entry = {"variant_key": k, "query_date": date.today().isoformat(),
                     "cadd_version": "v1.7", "build": "GRCh38",
                     "retry": True}
            if rec:
                entry.update(rec); n_ok += 1
            else:
                entry.update({"phred": None, "raw": None})
            cf.write(json.dumps(entry) + "\n")
    print(f"[09b-retry] recovered {n_ok}/{len(todo)}", flush=True)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
