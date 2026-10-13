#!/usr/bin/env python3
"""Stage 06: isoform-specific variant effects (real implementation).

For each labeled variant: load its gene's GENCODE v46 transcripts (once per
gene), run src/isoformvus/data/transcripts.py classification per transcript,
and summarize per gene: all / some / isoform_specific / none.

Also records, per variant:
    n_transcripts, n_transcripts_overlapped, n_isoform_specific_hits,
    has_junction_specific_hit, has_coding_specific_hit,
    mane_select_effect (effect on the MANE Select transcript, if annotated)

Output: data/processed/variant_isoform_effects.tsv.gz

Usage: python3 scripts/06_variant_isoform_effects.py --config config/base.yaml
"""

from __future__ import annotations

import argparse
import gzip
import sys
from collections import defaultdict
from pathlib import Path

REPO = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(REPO / "src"))

from isoformvus.data.gencode import load_gene_transcripts  # noqa: E402
from isoformvus.data.transcripts import (IsoformEffect,  # noqa: E402
                                         classify_transcript_effects,
                                         summarize_gene_effect)
from isoformvus.utils.variants import parse_variant_string  # noqa: E402


def load_mane_select(path: Path) -> dict[str, str]:
    """gene symbol -> MANE Select Ensembl transcript id (versioned)."""
    mane = {}
    with gzip.open(path, "rt") as fh:
        for line in fh:
            if line.startswith("#"):
                continue
            cols = line.rstrip("\n").split("\t")
            if len(cols) < 10:
                continue
            symbol, ens_nuc, status = cols[3], cols[7], cols[9]
            if status == "MANE Select":
                mane[symbol] = ens_nuc
    return mane


def main() -> int:
    ap = argparse.ArgumentParser()
    ap.add_argument("--config", default=str(REPO / "config" / "base.yaml"))
    ap.add_argument("--labels", default=str(REPO / "data" / "processed" / "variant_labels.tsv.gz"))
    ap.add_argument("--out", default=str(REPO / "data" / "processed" / "variant_isoform_effects.tsv.gz"))
    args = ap.parse_args()

    raw = REPO / "data" / "raw"
    gtf = raw / "gencode.v46.primary_assembly.annotation.gtf.gz"
    mane_path = raw / "MANE.GRCh38.v1.5.summary.txt.gz"
    for p in (Path(args.labels), gtf):
        if not p.exists():
            print(f"[06] missing {p}", file=sys.stderr)
            return 1

    mane = load_mane_select(mane_path) if mane_path.exists() else {}
    print(f"[06] MANE Select transcripts for {len(mane)} genes")

    # Group variants by gene so each gene's transcripts load once.
    by_gene: dict[str, list[tuple[str, str]]] = defaultdict(list)
    with gzip.open(args.labels, "rt") as fh:
        fh.readline()
        for line in fh:
            cols = line.rstrip("\n").split("\t")
            by_gene[cols[6]].append((cols[0], cols[0]))  # (variant_key, variant_key)

    print(f"[06] loading GENCODE transcripts for {len(by_gene)} genes ...")
    tx_by_gene = load_gene_transcripts(gtf, gene_symbols=set(by_gene))

    out = Path(args.out)
    n_var, n_iso_specific = 0, 0
    effect_counts: dict[str, int] = defaultdict(int)
    with gzip.open(out, "wt") as outfh:
        outfh.write("\t".join(["variant_key", "gene_symbol", "n_transcripts",
                                "n_overlapped", "gene_effect_summary",
                                "n_isoform_specific_hits",
                                "has_junction_specific_hit",
                                "has_coding_specific_hit",
                                "mane_select_effect"]) + "\n")
        for gene in sorted(by_gene):
            txs = tx_by_gene.get(gene, [])
            mane_tx = mane.get(gene, "")
            for vkey, _ in by_gene[gene]:
                v = parse_variant_string(vkey)
                effs = classify_transcript_effects(v, txs)
                summary = summarize_gene_effect(effs)
                effect_counts[summary] += 1
                overlapped = sum(1 for e in effs if e.effect != IsoformEffect.NONE)
                specific = {IsoformEffect.EXON_SPECIFIC,
                            IsoformEffect.CODING_SPECIFIC,
                            IsoformEffect.JUNCTION_SPECIFIC}
                n_specific = sum(1 for e in effs if e.effect in specific)
                has_junc = any(e.effect == IsoformEffect.JUNCTION_SPECIFIC for e in effs)
                has_coding = any(e.effect == IsoformEffect.CODING_SPECIFIC for e in effs)
                mane_eff = next((e.effect.value for e in effs
                                 if e.transcript_id == mane_tx), "")
                if summary == "isoform_specific":
                    n_iso_specific += 1
                outfh.write("\t".join([vkey, gene, str(len(txs)),
                                        str(overlapped), summary,
                                        str(n_specific), str(int(has_junc)),
                                        str(int(has_coding)), mane_eff]) + "\n")
                n_var += 1
    print(f"[06] {n_var} variants classified; gene_effect distribution: "
          f"{dict(effect_counts)}")
    print(f"[06] isoform-specific variants: {n_iso_specific} "
          f"({100*n_iso_specific/max(n_var,1):.1f}%)")
    print(f"[06] wrote {out}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
