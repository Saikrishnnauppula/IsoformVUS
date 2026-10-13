#!/usr/bin/env python3
"""Stage 02: build the ~100-gene Mendelian panel (real implementation).

Inputs (data/raw/):
    variant_summary.txt.gz          ClinVar monthly summary
    phenotype_to_genes.txt          HPO gene -> phenotype annotations
    gencode.v46.primary_assembly.annotation.gtf.gz

Algorithm (per metadata/gene_selection.md):
    1. Stream variant_summary; keep Assembly == GRCh38 rows.
    2. Per gene: count P/LP (review = criteria provided, no conflicts) and
       B/LB variants. Exclude conflicting classifications.
    3. Keep genes with >= min_pathogenic (default 10) P/LP and >= 10 B/LB.
    4. Keep genes with >= 2 protein-coding GENCODE v46 transcripts having
       distinct exon structures (multi-isoform requirement).
    5. Keep genes with >= 1 HPO annotation (phenotype link for the
       phenotype -> cell-type chain).
    6. Assign a coarse disease category from HPO annotations; cap any single
       category at 40% of the panel (diversity guardrail).
    7. Prefer genes with coding-alternative exons; flag known
       isoform-selective genes (TRNT1 + literature list).

Output: metadata/gene_panel.tsv (see PANEL_COLUMNS).

Usage: python3 scripts/02_build_gene_panel.py --config config/base.yaml
"""

from __future__ import annotations

import argparse
import gzip
import sys
from collections import Counter, defaultdict
from pathlib import Path

REPO = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(REPO / "src"))

from isoformvus.data.gencode import load_gene_transcripts  # noqa: E402
from isoformvus.data.hpo import categorize_gene, parse_obo  # noqa: E402

PATHOGENIC = {"Pathogenic", "Likely pathogenic", "Pathogenic/Likely pathogenic"}
BENIGN = {"Benign", "Likely benign", "Benign/Likely benign"}

# Coarse disease categories are derived from the HPO ontology itself
# (src/isoformvus/data/hpo.py): each gene's HPO terms are walked up to the
# direct children of HP:0000118 (Phenotypic abnormality); the most frequent
# child names the category. No hand-written term->category hints.

KNOWN_ISOFORM_SELECTIVE = {
    # gene -> (source, note); literature-backed isoform-selective mechanisms
    "TRNT1": ("Ly et al. 2025", "alternative start codons; SIFD"),
    "DNM1": ("AlphaGenome Atlas 2026 preprint", "brain-specific exon 10a"),
}


def has_criteria(review_status: str) -> bool:
    return "criteria provided" in review_status.lower()


def is_conflicting(clinsig: str) -> bool:
    return "conflict" in clinsig.lower()


def stream_variant_summary(path: Path):
    """Yield (gene_symbol, bucket) for GRCh38 rows with usable labels."""
    with gzip.open(path, "rt") as fh:
        header = None
        for line in fh:
            if line.startswith("#"):
                if line.startswith("#AlleleID"):
                    header = line.lstrip("#").rstrip("\n").split("\t")
                continue
            if header is None:
                continue
            cols = line.rstrip("\n").split("\t")
            if len(cols) < len(header):
                continue
            row = dict(zip(header, cols))
            if row.get("Assembly") != "GRCh38":
                continue
            gene = row.get("GeneSymbol", "").strip()
            if not gene or gene == "-":
                continue
            clinsig = row.get("ClinicalSignificance", "")
            if is_conflicting(clinsig) or not has_criteria(row.get("ReviewStatus", "")):
                continue
            if clinsig in PATHOGENIC:
                yield gene, "pathogenic"
            elif clinsig in BENIGN:
                yield gene, "benign"


def load_hpo_gene_annotations(path: Path) -> dict[str, set[str]]:
    """gene symbol -> set of HPO ids (from phenotype_to_genes.txt)."""
    gene_hpo: dict[str, set[str]] = defaultdict(set)
    with open(path) as fh:
        for line in fh:
            if line.startswith("#") or not line.strip():
                continue
            parts = line.rstrip("\n").split("\t")
            if len(parts) < 4:
                continue
            # Verified 2026-10-06 (v2026-02-16 release):
            # hpo_id, hpo_name, ncbi_gene_id, gene_symbol, disease_id
            if parts[0] == "hpo_id":
                continue  # header row
            gene_hpo[parts[3]].add(parts[0])
    return dict(gene_hpo)


def main() -> int:
    ap = argparse.ArgumentParser()
    ap.add_argument("--config", default=str(REPO / "config" / "base.yaml"))
    ap.add_argument("--min-pathogenic", type=int, default=10)
    ap.add_argument("--min-benign", type=int, default=10)
    ap.add_argument("--target-size", type=int, default=100)
    ap.add_argument("--max-category-frac", type=float, default=0.40)
    ap.add_argument("--out", default=str(REPO / "metadata" / "gene_panel.tsv"))
    args = ap.parse_args()

    raw = REPO / "data" / "raw"
    vs_path = raw / "variant_summary.txt.gz"
    hpo_path = raw / "phenotype_to_genes.txt"
    gtf_path = raw / "gencode.v46.primary_assembly.annotation.gtf.gz"
    for p in (vs_path, gtf_path):
        if not p.exists():
            print(f"[02] missing {p}; run scripts/01_download_sources.py first",
                  file=sys.stderr)
            return 1
    hpo_available = hpo_path.exists()
    if not hpo_available:
        print("[02] WARNING: phenotype_to_genes.txt missing; HPO criterion "
              "will be marked pending and genes kept provisionally.")

    print("[02] streaming variant_summary ...")
    counts: dict[str, Counter] = defaultdict(Counter)
    for gene, bucket in stream_variant_summary(vs_path):
        counts[gene][bucket] += 1
    print(f"[02] {len(counts)} genes with GRCh38 labeled variants")

    candidates = [g for g, c in counts.items()
                  if c["pathogenic"] >= args.min_pathogenic
                  and c["benign"] >= args.min_benign]
    print(f"[02] {len(candidates)} genes pass count thresholds "
          f"(P>={args.min_pathogenic}, B>={args.min_benign})")

    gene_hpo = load_hpo_gene_annotations(hpo_path) if hpo_available else {}
    hpo_names, hpo_parents = {}, {}
    obo_path = raw / "hp.obo"
    if hpo_available and obo_path.exists():
        print("[02] parsing hp.obo for ontology-based categorization ...")
        hpo_names, hpo_parents = parse_obo(obo_path)
        print(f"[02] hp.obo: {len(hpo_names)} terms")
    if hpo_available:
        print(f"[02] HPO annotations for {len(gene_hpo)} genes")
        candidates = [g for g in candidates if g in gene_hpo]
        print(f"[02] {len(candidates)} candidates have HPO annotations")

    print("[02] loading GENCODE transcripts for candidates ...")
    tx_by_gene = load_gene_transcripts(gtf_path, gene_symbols=set(candidates))

    rows = []
    for gene in sorted(candidates):
        txs = tx_by_gene.get(gene, [])
        pc = [t for t in txs if t.biotype == "protein_coding"]
        structs = {tuple((e.start, e.end) for e in t.exons) for t in pc}
        if len(structs) < 2:
            continue  # multi-isoform requirement
        # coding-alternative exons: distinct CDS exon structures
        cds_structs = set()
        for t in pc:
            if t.cds_start and t.cds_end:
                cds_exons = tuple((max(e.start, t.cds_start), min(e.end, t.cds_end))
                                  for e in t.exons
                                  if e.end >= t.cds_start and e.start <= t.cds_end)
                cds_structs.add(cds_exons)
        hpos = gene_hpo.get(gene, set())
        category = (categorize_gene(hpos, hpo_parents, hpo_names)
                    if hpo_parents else "unmapped")
        rows.append({
            "gene_symbol": gene,
            "disease_category": category,
            "n_pathogenic_clinvar": counts[gene]["pathogenic"],
            "n_benign_clinvar": counts[gene]["benign"],
            "n_protein_coding_transcripts": len(pc),
            "n_distinct_exon_structures": len(structs),
            "n_distinct_cds_structures": len(cds_structs),
            "has_coding_alternative_exon": len(cds_structs) >= 2,
            "n_hpo_terms": len(hpos),
            "known_isoform_selective": gene in KNOWN_ISOFORM_SELECTIVE,
            "known_isoform_note": KNOWN_ISOFORM_SELECTIVE.get(gene, ("", ""))[1],
        })

    print(f"[02] {len(rows)} genes pass the multi-isoform requirement")

    # Diversity cap: sort by pathogenic count desc, greedily take while no
    # category exceeds max_category_frac of target size.
    rows.sort(key=lambda r: r["n_pathogenic_clinvar"], reverse=True)
    cap = int(args.target_size * args.max_category_frac)
    cat_counts: Counter = Counter()
    panel = []
    overflow = []
    for r in rows:
        if len(panel) >= args.target_size:
            overflow.append(r)
            continue
        if cat_counts[r["disease_category"]] >= cap:
            overflow.append(r)
            continue
        panel.append(r)
        cat_counts[r["disease_category"]] += 1

    # Always include known isoform-selective genes if they passed filters.
    for r in overflow:
        if r["known_isoform_selective"] and r not in panel:
            panel.append(r)
    print(f"[02] panel size {len(panel)} "
          f"(target {args.target_size}, overflow {len(overflow)})")
    print(f"[02] category distribution: {dict(cat_counts)}")

    out = Path(args.out)
    out.parent.mkdir(parents=True, exist_ok=True)
    cols = ["gene_symbol", "disease_category", "n_pathogenic_clinvar",
            "n_benign_clinvar", "n_protein_coding_transcripts",
            "n_distinct_exon_structures", "n_distinct_cds_structures",
            "has_coding_alternative_exon", "n_hpo_terms",
            "known_isoform_selective", "known_isoform_note"]
    with open(out, "w") as fh:
        fh.write("\t".join(cols) + "\n")
        for r in panel:
            fh.write("\t".join(str(r[c]) for c in cols) + "\n")
    print(f"[02] wrote {out}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
