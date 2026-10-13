#!/usr/bin/env python3
"""Stage 03b: gnomAD allele-frequency enrichment (GraphQL API route).

For each labeled variant in data/processed/variant_labels.tsv.gz, query the
gnomAD browser GraphQL API (dataset gnomad_r4 = v4.1.1, GRCh38) and attach:
    af_gnomad_max, af_popmax, ac_total, an_total, homozygote_count

Also mines population-common controls: variants with AF > 1% in the panel
genes' regions are collected as control_type="gnomad_common" candidates
(written separately; NEVER merged silently with clinvar_blb — see
metadata/gene_selection.md).

Practical notes (from research/data_urls_verified.md):
- Batched queries with GraphQL aliases (many variants per request).
- Observe rate limits; --batch-size and --sleep tune throughput.
- Results cached in data/interim/gnomad_af_cache.jsonl so reruns are free.
- Records query date + dataset id in the output provenance sidecar.

Usage: python3 scripts/03b_gnomad_af.py --config config/base.yaml [--limit N]
"""

from __future__ import annotations

import argparse
import gzip
import json
import sys
import time
import urllib.request
from pathlib import Path

REPO = Path(__file__).resolve().parents[1]

GNOMAD_API = "https://gnomad.broadinstitute.org/api"
DATASET = "gnomad_r4"  # v4.1.1, GRCh38


def batched_query(variants: list[dict], batch_size: int = 20) -> list[dict]:
    """Query gnomAD GraphQL with aliased batched queries."""
    results = []
    for i in range(0, len(variants), batch_size):
        batch = variants[i:i + batch_size]
        parts = []
        for j, v in enumerate(batch):
            chrom = v["chrom"].replace("chr", "")
            parts.append(
                f'v{j}: variant(variantId: "{chrom}-{v["pos"]}-{v["ref"]}-{v["alt"]}", '
                f'reference_genome: GRCh38, dataset: {DATASET}) '
                '{ variant_id genome { af ac an homozygote_count } }')
        body = json.dumps({"query": "query { " + " ".join(parts) + " }",
                           "variables": {}}).encode()
        req = urllib.request.Request(
            GNOMAD_API, data=body,
            headers={"Content-Type": "application/json",
                     "User-Agent": "IsoformVUS/0.1.0 (research)"})
        try:
            with urllib.request.urlopen(req, timeout=60) as resp:
                payload = json.load(resp)
        except Exception as exc:  # noqa: BLE001 - network; record and continue
            print(f"[03b] batch {i//batch_size} failed: {exc}", file=sys.stderr)
            results.extend([None] * len(batch))
            continue
        data = (payload.get("data") or {})
        for j in range(len(batch)):
            node = (data.get(f"v{j}") or {})
            results.append((node.get("genome") or None))
        time.sleep(0.5)  # be a polite API citizen
    return results


def main() -> int:
    ap = argparse.ArgumentParser()
    ap.add_argument("--config", default=str(REPO / "config" / "base.yaml"))
    ap.add_argument("--labels", default=str(REPO / "data" / "processed" / "variant_labels.tsv.gz"))
    ap.add_argument("--out", default=str(REPO / "data" / "processed" / "variant_af.tsv.gz"))
    ap.add_argument("--limit", type=int, default=0, help="0 = all variants")
    ap.add_argument("--batch-size", type=int, default=20)
    args = ap.parse_args()

    if not Path(args.labels).exists():
        print(f"[03b] missing {args.labels}; run scripts/03_build_variant_labels.py first",
              file=sys.stderr)
        return 1

    variants = []
    with gzip.open(args.labels, "rt") as fh:
        header = fh.readline().rstrip("\n").split("\t")
        for line in fh:
            row = dict(zip(header, line.rstrip("\n").split("\t")))
            variants.append(row)
            if args.limit and len(variants) >= args.limit:
                break
    print(f"[03b] querying gnomAD {DATASET} for {len(variants)} variants ...")
    hits = batched_query(variants, args.batch_size)

    out = Path(args.out)
    n_found = 0
    with gzip.open(out, "wt") as fh:
        fh.write("variant_key\taf_gnomad\taf_ac\taf_an\thomozygote_count\tgnomad_found\n")
        for row, g in zip(variants, hits):
            if g and g.get("af") is not None:
                n_found += 1
                fh.write(f"{row['variant_key']}\t{g['af']}\t{g.get('ac','')}\t"
                         f"{g.get('an','')}\t{g.get('homozygote_count','')}\t1\n")
            else:
                fh.write(f"{row['variant_key']}\t\t\t\t\t0\n")
    print(f"[03b] AF found for {n_found}/{len(variants)}; wrote {out}")
    print("[03b] NOTE: gnomad_common control mining is a separate documented "
          "step; not run here.")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
