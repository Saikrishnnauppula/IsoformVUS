#!/usr/bin/env python3
"""Stage 03: build labeled variant set from the ClinVar GRCh38 VCF.

Labels (per docs/leakage_prevention.md — VUSs are NEVER a training class):
    1 = Pathogenic / Likely pathogenic, review status >= 1 star, no conflicts
    0 = Benign / Likely benign, review status >= 1 star, no conflicts
        (control_type = "clinvar_blb")

gnomAD population-common controls (control_type = "gnomad_common") are added
by scripts/03b_gnomad_af.py, which also fills af_gnomad_max.

Review-status mapping (ClinVar VCF CLNREVSTAT):
    "criteria_provided,_single_submitter"            -> 1 star
    "criteria_provided,_multiple_submitters,_no_conflicts" -> 2 stars
    "reviewed_by_expert_panel"                      -> 3 stars
    "practice_guideline"                            -> 4 stars
Anything else (no_assertion, conflicting) is excluded from labels.

Output: data/processed/variant_labels.tsv.gz
    variant_key, chrom, pos, ref, alt, vartype, gene_symbol,
    label, control_type, review_stars, clnsig, condition, rsid

Usage: python3 scripts/03_build_variant_labels.py --config config/base.yaml
       [--panel metadata/gene_panel.tsv] [--vcf data/raw/clinvar_20261004.vcf.gz]
"""

from __future__ import annotations

import argparse
import gzip
import sys
from pathlib import Path

REPO = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(REPO / "src"))

from isoformvus.utils.variants import vcf_to_variant  # noqa: E402

STARS = {
    "criteria_provided,_single_submitter": 1,
    "criteria_provided,_multiple_submitters,_no_conflicts": 2,
    "reviewed_by_expert_panel": 3,
    "practice_guideline": 4,
}

PATHOGENIC_TOKENS = ("Pathogenic", "Likely_pathogenic")
BENIGN_TOKENS = ("Benign", "Likely_benign")


def parse_info(info: str) -> dict:
    d: dict[str, str] = {}
    for kv in info.split(";"):
        if "=" in kv:
            k, v = kv.split("=", 1)
            d[k] = v
    return d


def classify(clnsig: str) -> str | None:
    """Return 'pathogenic' / 'benign' / None (excluded)."""
    s = clnsig
    if "Conflicting" in s or "Uncertain" in s or not s:
        return None
    is_p = any(t in s for t in PATHOGENIC_TOKENS)
    is_b = any(t in s for t in BENIGN_TOKENS)
    if is_p and not is_b:
        return "pathogenic"
    if is_b and not is_p:
        return "benign"
    return None


def load_panel_genes(path: Path) -> set[str]:
    genes = set()
    with open(path) as fh:
        header = fh.readline()
        for line in fh:
            genes.add(line.split("\t")[0])
    return genes


def main() -> int:
    ap = argparse.ArgumentParser()
    ap.add_argument("--config", default=str(REPO / "config" / "base.yaml"))
    ap.add_argument("--panel", default=str(REPO / "metadata" / "gene_panel.tsv"))
    ap.add_argument("--vcf", default=str(REPO / "data" / "raw" / "clinvar_20261004.vcf.gz"))
    ap.add_argument("--out", default=str(REPO / "data" / "processed" / "variant_labels.tsv.gz"))
    ap.add_argument("--max-per-gene-per-class", type=int, default=300,
                    help="cap variants per (gene, class); 0 = no cap. "
                         "Prevents TTN-scale genes dominating the benchmark.")
    ap.add_argument("--seed", type=int, default=20261006)
    args = ap.parse_args()

    panel_genes = load_panel_genes(Path(args.panel)) if Path(args.panel).exists() else None
    if panel_genes is None:
        print("[03] no gene panel; labeling all genes (panel filter skipped)",
              file=sys.stderr)
    else:
        print(f"[03] panel genes: {len(panel_genes)}")

    out = Path(args.out)
    out.parent.mkdir(parents=True, exist_ok=True)

    import random
    rng = random.Random(args.seed)

    # Buffer rows per (gene, label) so the per-gene cap can be applied with
    # a seeded subsample instead of truncation order.
    buffered: dict[tuple[str, int], list[str]] = {}
    n_in, n_skip_conflict, n_skip_nostar = 0, 0, 0
    with gzip.open(args.vcf, "rt") as fh:
        for line in fh:
            if line.startswith("#"):
                continue
            n_in += 1
            cols = line.rstrip("\n").split("\t")
            chrom, pos, _id, ref, alts = cols[0], int(cols[1]), cols[2], cols[3], cols[4]
            info = parse_info(cols[7])
            bucket = classify(info.get("CLNSIG", ""))
            if bucket is None:
                n_skip_conflict += 1
                continue
            stars = STARS.get(info.get("CLNREVSTAT", ""), 0)
            if stars < 1:
                n_skip_nostar += 1
                continue
            geneinfo = info.get("GENEINFO", "")
            gene = geneinfo.split(":")[0] if geneinfo else ""
            if panel_genes is not None and gene not in panel_genes:
                continue
            label = 1 if bucket == "pathogenic" else 0
            for alt in alts.split(","):
                if alt == "<NON_REF>" or "*" in alt:
                    continue
                try:
                    v = vcf_to_variant(chrom, pos, ref, alt)
                except ValueError:
                    continue
                row = "\t".join([
                    v.key, v.chrom, str(v.pos), v.ref, v.alt, v.variant_type,
                    gene, str(label), "clinvar_blb", str(stars),
                    info.get("CLNSIG", ""), info.get("CLNDN", ""),
                    info.get("RS", "")])
                buffered.setdefault((gene, label), []).append(row)

    n_out, n_capped = 0, 0
    with gzip.open(out, "wt") as outfh:
        outfh.write("\t".join(["variant_key", "chrom", "pos", "ref", "alt",
                                "vartype", "gene_symbol", "label",
                                "control_type", "review_stars", "clnsig",
                                "condition", "rsid"]) + "\n")
        for (gene, label), rows in sorted(buffered.items()):
            if args.max_per_gene_per_class and len(rows) > args.max_per_gene_per_class:
                rows = rng.sample(rows, args.max_per_gene_per_class)
                n_capped += 1
            for row in rows:
                outfh.write(row + "\n")
                n_out += 1
    print(f"[03] scanned {n_in} VCF records; wrote {n_out} labeled variants "
          f"({n_capped} gene-class groups capped at {args.max_per_gene_per_class}); "
          f"skipped conflicting/uncertain={n_skip_conflict}, no-stars={n_skip_nostar}")
    print(f"[03] wrote {out}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
