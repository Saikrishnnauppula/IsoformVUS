#!/usr/bin/env python3
"""Build a VCF of test-set variants for SpliceAI scoring.

SpliceAI input: standard VCF with CHROM/POS/ID/REF/ALT. Only SNVs and
short indels (<=50bp) are passed through; complex/symbolic alleles are
logged and skipped (SpliceAI limitation, disclosed).
"""

from __future__ import annotations

import gzip
import json
import sys
from pathlib import Path

REPO = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(REPO / "src"))
from isoformvus.utils.variants import parse_variant_string  # noqa: E402


def main() -> int:
    proc = REPO / "data" / "processed"
    splits = {a["gene"]: a["split"]
              for a in json.loads((REPO / "results/splits.json").read_text())["assignments"]}
    out = REPO / "data/interim/spliceai_test.vcf"
    skipped = 0
    written = 0
    with gzip.open(proc / "variant_labels.tsv.gz", "rt") as fh, \
            open(out, "w") as outfh:
        header = fh.readline()
        outfh.write("##fileformat=VCFv4.2\n")
        outfh.write("##reference=GRCh38\n")
        for i in list(range(1, 23)) + ["X", "Y", "M"]:
            outfh.write(f"##contig=<ID=chr{i}>\n")
        outfh.write("#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\n")
        for line in fh:
            c = line.rstrip("\n").split("\t")
            key, gene = c[0], c[6]
            if splits.get(gene) != "test":
                continue
            try:
                v = parse_variant_string(key)
            except Exception:  # noqa: BLE001
                skipped += 1
                continue
            if len(v.ref) > 50 or len(v.alt) > 50 or "<" in v.alt:
                skipped += 1
                continue
            outfh.write(f"{v.chrom}\t{v.pos}\t{key}\t{v.ref}\t{v.alt}\t.\t.\t.\n")
            written += 1
    print(f"[spliceai] wrote {out}: {written} variants, {skipped} skipped")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
