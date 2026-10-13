#!/usr/bin/env python3
"""Stage 09c: AlphaMissense join (real implementation).

Streams AlphaMissense_hg38.tsv.gz (Zenodo 8208688, md5-verified) once and
joins scores to labeled variants on (chrom, pos, ref, alt).

AlphaMissense is MISSENSE-ONLY: non-missense variants get an explicit
missing indicator (am_missing=1) per the feature spec — never a zero
imputation. Comparisons against AlphaMissense are restricted to the
missense-compatible subset (documented in results).

Writes data/processed/variant_alphamissense.tsv.gz:
    variant_key, am_score, am_class, am_missing

Provenance: AlphaMissense v1 (2023-09-19), hg38, CC-BY-NC 4.0.
"""

from __future__ import annotations

import argparse
import gzip
import sys
from pathlib import Path

REPO = Path(__file__).resolve().parents[1]


def main() -> int:
    ap = argparse.ArgumentParser()
    ap.add_argument("--labels", default=str(REPO / "data/processed/variant_labels.tsv.gz"))
    ap.add_argument("--am", default=str(REPO / "data/raw/AlphaMissense_hg38.tsv.gz"))
    ap.add_argument("--out", default=str(REPO / "data/processed/variant_alphamissense.tsv.gz"))
    args = ap.parse_args()

    for p in (args.labels, args.am):
        if not Path(p).exists():
            print(f"[09c] missing {p}", file=sys.stderr)
            return 1

    # Our variant keys, normalized to AM's chrom style (AM uses 'chr' prefix).
    wanted: dict[str, str] = {}  # am_key -> variant_key
    with gzip.open(args.labels, "rt") as fh:
        fh.readline()
        for line in fh:
            c = line.rstrip("\n").split("\t")
            vkey, chrom, pos, ref, alt = c[0], c[1], c[2], c[3], c[4]
            if not chrom.startswith("chr"):
                chrom = "chr" + chrom
            wanted[f"{chrom}:{pos}:{ref}:{alt}"] = vkey
    print(f"[09c] {len(wanted)} labeled variants to match")

    found: dict[str, tuple[str, str]] = {}
    n_rows = 0
    with gzip.open(args.am, "rt") as fh:
        for line in fh:
            if line.startswith("#"):
                continue
            n_rows += 1
            # chr pos ref alt genome uniprot transcript aa_change score class
            c = line.rstrip("\n").split("\t")
            if len(c) < 10:
                continue
            key = f"{c[0]}:{c[1]}:{c[2]}:{c[3]}"
            if key in wanted and wanted[key] not in found:
                found[wanted[key]] = (c[8], c[9])
            if n_rows % 50_000_000 == 0:
                print(f"[09c] streamed {n_rows} AM rows; matched {len(found)}")
    print(f"[09c] streamed {n_rows} AM rows; matched {len(found)}/{len(wanted)}")

    with gzip.open(args.out, "wt") as fh:
        fh.write("variant_key\tam_score\tam_class\tam_missing\n")
        for am_key, vkey in wanted.items():
            if vkey in found:
                score, cls = found[vkey]
                fh.write(f"{vkey}\t{score}\t{cls}\t0\n")
            else:
                fh.write(f"{vkey}\t\t\t1\n")
    print(f"[09c] wrote {args.out}")
    print("[09c] NOTE: unmatched = non-missense or absent from AM; "
          "am_missing=1 is the explicit indicator (no zero-imputation).")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
