#!/usr/bin/env python3
"""Run SpliceAI on the test-set VCF and parse max delta scores.

Requires: GRCh38 FASTA (data/raw/GRCh38.primary_assembly.genome.fa.gz,
uncompressed to .fa), GENCODE v46 GTF annotation, spliceai venv at
~/workspace/.spliceai_env.

Output: data/processed/variant_spliceai.tsv.gz with columns
    variant_key, spliceai_max_delta, spliceai_type
spliceai_max_delta = max over (DS-AG, DS-AL, DS-DG, DS-DL).
"""

from __future__ import annotations

import gzip
import subprocess
import sys
from pathlib import Path

REPO = Path(__file__).resolve().parents[1]
VENV_PY = Path.home() / "workspace/.spliceai_env/bin/python"


def main() -> int:
    raw = REPO / "data/raw"
    fa_gz = raw / "GRCh38.primary_assembly.genome.fa.gz"
    fa = raw / "GRCh38.primary_assembly.genome.fa"
    gtf = raw / "gencode.v46.primary_assembly.annotation.gtf.gz"
    vcf = REPO / "data/interim/spliceai_test.vcf"
    out_vcf = REPO / "data/interim/spliceai_test.scored.vcf"

    if not fa.exists():
        print(f"[spliceai] decompressing {fa_gz.name} (~3GB uncompressed)...")
        with gzip.open(fa_gz, "rb") as fin, open(fa, "wb") as fout:
            while chunk := fin.read(1 << 20):
                fout.write(chunk)

    # SpliceAI ships its own GRCh38 canonical annotation ("grch38"); the
    # GENCODE GTF is not in the UCSC table format SpliceAI requires.
    cmd = [str(VENV_PY), "-m", "spliceai", "-I", str(vcf), "-O", str(out_vcf),
           "-R", str(fa), "-A", "grch38", "-D", "500"]
    print("[spliceai] running:", " ".join(cmd))
    r = subprocess.run(cmd, capture_output=True, text=True, timeout=86400)
    print(r.stdout[-2000:])
    if r.returncode != 0:
        print(r.stderr[-3000:], file=sys.stderr)
        return 1

    # Parse SpliceAI VCF INFO: SpliceAI=ALLELE|SYMBOL|DS_AG|DS_AL|DS_DG|DS_DL|...
    rows = []
    with open(out_vcf) as fh:
        for line in fh:
            if line.startswith("#"):
                continue
            c = line.rstrip("\n").split("\t")
            vid = c[2]
            info = c[7]
            best = 0.0
            stype = "none"
            if "SpliceAI=" in info:
                pred = info.split("SpliceAI=")[1].split(";")[0]
                for entry in pred.split(","):
                    f = entry.split("|")
                    if len(f) >= 6:
                        try:
                            ds = [float(f[2]), float(f[3]),
                                  float(f[4]), float(f[5])]
                        except ValueError:
                            continue
                        m = max(ds)
                        if m > best:
                            best = m
                            stype = ["AG", "AL", "DG", "DL"][ds.index(m)]
            rows.append((vid, best, stype))
    out = REPO / "data/processed/variant_spliceai.tsv.gz"
    with gzip.open(out, "wt") as fh:
        fh.write("variant_key\tspliceai_max_delta\tspliceai_type\n")
        for vid, best, stype in rows:
            fh.write(f"{vid}\t{best:.4f}\t{stype}\n")
    print(f"[spliceai] wrote {out} ({len(rows)} variants)")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
