#!/usr/bin/env python3
"""Stage 09b: CADD v1.7 (GRCh38) enrichment via the public REST API.

Fetches PHRED + RawScore per variant with a resume-safe JSONL cache
(data/interim/cadd_cache.jsonl). Polite rate limiting; skips cached keys.

Default scope is validation+test variants only (the benchmark sets), because
60k sequential API calls would take ~12h. Train-set CADD can be backfilled
later; the cache makes reruns incremental.

Writes data/processed/variant_cadd.tsv.gz:
    variant_key, cadd_phred, cadd_raw

Provenance: CADD v1.7, GRCh38, REST endpoint GRCh38-v1.7, query dates in cache.
"""

from __future__ import annotations

import argparse
import gzip
import json
import sys
import time
import urllib.request
from datetime import date
from pathlib import Path

REPO = Path(__file__).resolve().parents[1]

API = "https://cadd.gs.washington.edu/api/v1.0/GRCh38-v1.7/{chrom}:{pos}_{ref}_{alt}"


def load_cache(path: Path) -> dict:
    cache = {}
    if path.exists():
        with open(path) as fh:
            for line in fh:
                try:
                    r = json.loads(line)
                    cache[r["variant_key"]] = r
                except (json.JSONDecodeError, KeyError):
                    continue
    return cache


def fetch_one(chrom: str, pos: int, ref: str, alt: str) -> dict | None:
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
    ap.add_argument("--labels", default=str(REPO / "data/processed/variant_labels.tsv.gz"))
    ap.add_argument("--splits", default=str(REPO / "results/splits.json"))
    ap.add_argument("--scope", default="val_test", choices=["val_test", "all"],
                    help="val_test = validation+test variants only (default)")
    ap.add_argument("--cache", default=str(REPO / "data/interim/cadd_cache.jsonl"))
    ap.add_argument("--out", default=str(REPO / "data/processed/variant_cadd.tsv.gz"))
    ap.add_argument("--sleep", type=float, default=0.4, help="seconds between calls")
    ap.add_argument("--limit", type=int, default=0)
    args = ap.parse_args()

    splits = {a["gene"]: a["split"]
              for a in json.loads(Path(args.splits).read_text())["assignments"]}

    targets = []
    with gzip.open(args.labels, "rt") as fh:
        fh.readline()
        for line in fh:
            c = line.rstrip("\n").split("\t")
            vkey, chrom, pos, ref, alt, gene = c[0], c[1], int(c[2]), c[3], c[4], c[6]
            sp = splits.get(gene, "")
            if args.scope == "val_test" and sp not in ("validation", "test"):
                continue
            targets.append((vkey, chrom, pos, ref, alt))
            if args.limit and len(targets) >= args.limit:
                break
    print(f"[09b] {len(targets)} variants to score (scope={args.scope})")

    cache_path = Path(args.cache)
    cache_path.parent.mkdir(parents=True, exist_ok=True)
    cache = load_cache(cache_path)
    todo = [t for t in targets if t[0] not in cache]
    print(f"[09b] {len(cache)} cached, {len(todo)} to fetch")

    n_ok, n_fail = 0, 0
    with open(cache_path, "a") as cf:
        for i, (vkey, chrom, pos, ref, alt) in enumerate(todo):
            rec = fetch_one(chrom, pos, ref, alt)
            entry = {"variant_key": vkey, "query_date": date.today().isoformat(),
                     "cadd_version": "v1.7", "build": "GRCh38"}
            if rec:
                entry.update(rec); n_ok += 1
            else:
                entry.update({"phred": None, "raw": None}); n_fail += 1
            cf.write(json.dumps(entry) + "\n")
            if (i + 1) % 200 == 0:
                cf.flush()
                print(f"[09b] {i+1}/{len(todo)} fetched ({n_ok} ok, {n_fail} failed)")
            time.sleep(args.sleep)

    # Write the merged output table.
    cache = load_cache(cache_path)
    with gzip.open(args.out, "wt") as fh:
        fh.write("variant_key\tcadd_phred\tcadd_raw\n")
        for vkey, *_ in targets:
            e = cache.get(vkey, {})
            ph, ra = e.get("phred"), e.get("raw")
            fh.write(f"{vkey}\t{'' if ph is None else ph}\t{'' if ra is None else ra}\n")
    print(f"[09b] fetched {n_ok} ok / {n_fail} failed; wrote {args.out}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
