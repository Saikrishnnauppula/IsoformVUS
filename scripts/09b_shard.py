#!/usr/bin/env python3
"""Stage 09b-shard: parallel CADD v1.7 fetcher.

Run N copies with different --shard-idx (0..N-1). Each handles a disjoint
subset of variants (todo[i] where i % n_shards == shard_idx) and appends to
its own shard cache file, so workers never contend on one file.

After all shards finish, run once with --merge-only to combine the main
cache + all shards into data/processed/variant_cadd.tsv.gz.

Usage:
    python3 scripts/09b_shard.py --n-shards 4 --shard-idx 0 --sleep 0.35
    python3 scripts/09b_shard.py --n-shards 4 --merge-only
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
sys.path.insert(0, str(REPO / "scripts"))

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


def compute_targets(labels: Path, splits: Path, scope: str):
    sp = {a["gene"]: a["split"]
          for a in json.loads(splits.read_text())["assignments"]}
    targets = []
    with gzip.open(labels, "rt") as fh:
        fh.readline()
        for line in fh:
            c = line.rstrip("\n").split("\t")
            vkey, chrom, pos, ref, alt, gene = c[0], c[1], int(c[2]), c[3], c[4], c[6]
            if scope == "val_test" and sp.get(gene, "") not in ("validation", "test"):
                continue
            targets.append((vkey, chrom, pos, ref, alt))
    return targets


def main() -> int:
    ap = argparse.ArgumentParser()
    ap.add_argument("--n-shards", type=int, default=4)
    ap.add_argument("--shard-idx", type=int, default=0)
    ap.add_argument("--sleep", type=float, default=0.35)
    ap.add_argument("--scope", default="val_test")
    ap.add_argument("--merge-only", action="store_true")
    ap.add_argument("--cache", default=str(REPO / "data/interim/cadd_cache.jsonl"))
    ap.add_argument("--out", default=str(REPO / "data/processed/variant_cadd.tsv.gz"))
    args = ap.parse_args()

    if args.merge_only:
        merged: dict = load_cache(Path(args.cache))
        for i in range(args.n_shards):
            sp = Path(args.cache).parent / f"cadd_shard_{i}.jsonl"
            for k, v in load_cache(sp).items():
                merged.setdefault(k, v)
        with gzip.open(args.out, "wt") as fh:
            fh.write("variant_key\tcadd_phred\tcadd_raw\n")
            for vkey, e in merged.items():
                ph, ra = e.get("phred"), e.get("raw")
                fh.write(f"{vkey}\t{'' if ph is None else ph}\t{'' if ra is None else ra}\n")
        n_scored = sum(1 for e in merged.values() if e.get("phred") is not None)
        print(f"[09b-shard] merged {len(merged)} entries ({n_scored} scored) -> {args.out}")
        return 0

    targets = compute_targets(
        REPO / "data/processed/variant_labels.tsv.gz",
        REPO / "results/splits.json", args.scope)
    cache = load_cache(Path(args.cache))
    todo = [t for t in targets if t[0] not in cache]
    mine = [t for i, t in enumerate(todo) if i % args.n_shards == args.shard_idx]
    print(f"[09b-shard {args.shard_idx}/{args.n_shards}] {len(mine)} variants "
          f"(total todo {len(todo)})", flush=True)

    shard_path = Path(args.cache).parent / f"cadd_shard_{args.shard_idx}.jsonl"
    have = load_cache(shard_path)
    mine = [t for t in mine if t[0] not in have]
    print(f"[09b-shard {args.shard_idx}] {len(have)} already in shard, "
          f"{len(mine)} to fetch", flush=True)

    n_ok, n_fail = 0, 0
    with open(shard_path, "a") as cf:
        for i, (vkey, chrom, pos, ref, alt) in enumerate(mine):
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
                print(f"[09b-shard {args.shard_idx}] {i+1}/{len(mine)} "
                      f"({n_ok} ok, {n_fail} failed)", flush=True)
            time.sleep(args.sleep)
    print(f"[09b-shard {args.shard_idx}] DONE: {n_ok} ok / {n_fail} failed", flush=True)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
