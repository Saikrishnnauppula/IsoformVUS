#!/usr/bin/env python3
"""Merge chunked SpliceAI scored VCFs into the final enrichment TSV.

Finalization step of the isoformvus-baseline-watch: collects variant lines
with SpliceAI= annotations from every data/interim/spliceai*.scored.vcf,
dedupes by variant ID (keeping the max delta across all annotated copies),
and writes data/processed/variant_spliceai.tsv.gz with columns
    variant_key, spliceai_max_delta, spliceai_type

Variants that SpliceAI could not score (outside annotated transcript
windows -- SpliceAI v1.3 returns no delta scores for these) have no
SpliceAI= annotation in any chunk and are excluded; the baseline is
evaluated on its scored subset by scripts/14_add_baselines.py.

Parse logic mirrors scripts/09e_spliceai_run.py.
"""

from __future__ import annotations

import gzip
import glob
from pathlib import Path

REPO = Path(__file__).resolve().parents[1]
LABELS = ["AG", "AL", "DG", "DL"]


def parse_best(info: str) -> tuple[float, str] | None:
    if "SpliceAI=" not in info:
        return None
    pred = info.split("SpliceAI=")[1].split(";")[0]
    best, stype = 0.0, "none"
    for entry in pred.split(","):
        f = entry.split("|")
        if len(f) >= 6:
            try:
                ds = [float(f[2]), float(f[3]), float(f[4]), float(f[5])]
            except ValueError:
                continue
            m = max(ds)
            if m > best:
                best = m
                stype = LABELS[ds.index(m)]
    return best, stype


def main() -> int:
    scored: dict[str, tuple[float, str]] = {}
    files = sorted(glob.glob(str(REPO / "data/interim/spliceai*.scored.vcf")))
    print(f"[merge] {len(files)} scored VCFs")
    n_lines = n_ann = 0
    for path in files:
        with open(path) as fh:
            for line in fh:
                if line.startswith("#"):
                    continue
                n_lines += 1
                c = line.rstrip("\n").split("\t")
                if len(c) < 8:
                    continue
                vid, info = c[2], c[7]
                parsed = parse_best(info)
                if parsed is None:
                    continue
                n_ann += 1
                if vid not in scored or parsed[0] > scored[vid][0]:
                    scored[vid] = parsed
    out = REPO / "data/processed/variant_spliceai.tsv.gz"
    with gzip.open(out, "wt") as fh:
        fh.write("variant_key\tspliceai_max_delta\tspliceai_type\n")
        for vid in sorted(scored):
            best, stype = scored[vid]
            fh.write(f"{vid}\t{best:.4f}\t{stype}\n")
    print(f"[merge] variant lines={n_lines} annotated={n_ann} "
          f"unique scored={len(scored)} -> {out}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
