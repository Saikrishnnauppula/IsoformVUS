#!/usr/bin/env python3
"""Build the annotation-based MANE Select exon atlas (workstream 1c).

Panel finding C-2 / questions Q3, Q4, Q5, Q10: script 18's 4x2x2 tensor is
built on exons *reconstructed from observed variant positions* (clusters of
coding variants split at >90 bp gaps; s_pos/e_pos are the first/last
observed variant positions), not on annotated MANE/GENCODE coordinates.
This script builds the annotated counterpart from primary sources:

  - GENCODE v46 comprehensive annotation GTF (GRCh38), EBI release 46:
    https://ftp.ebi.ac.uk/pub/databases/gencode/Gencode_human/release_46/gencode.v46.annotation.gtf.gz
  - MANE Select v1.5 summary (NCBI RefSeq MANE release 1.5):
    https://ftp.ncbi.nlm.nih.gov/refseq/MANE/MANE_human/release_1.5/MANE.GRCh38.v1.5.summary.txt.gz

Both downloads live in data/raw/ (gitignored). If the comprehensive GTF is
absent but the primary-assembly GTF shipped with the repo is present, the
script falls back to it and says so on stdout (MANE Select transcripts sit
on the primary assembly, so the extracted exons are identical; script 21
verifies that claim exon-by-exon when both files are present).

For each of the 102 panel genes (metadata/gene_panel.tsv) the MANE Select
transcript is taken from the MANE summary (gene symbol -> Ensembl
transcript), and every exon of that transcript is extracted from the GTF:
exon number, chrom, start, end (GENCODE GTF 1-based inclusive), strand,
CDS overlap (from CDS features), and the GTF CDS phase.

Internal-coding definitions (both shipped, so the comparison with
script 18 is explicit):
  - is_internal_coding          : CDS exon that is neither the first nor
                                  the last CDS exon of the transcript
                                  (transcript orientation). This is the
                                  primary "internal coding exon"
                                  definition used by script 21.
  - is_internal_coding_s18rule  : CDS exon that is not the last CDS exon.
                                  Script 18's only exon-level exclusion is
                                  its variant-level in_last_exon == 0
                                  filter; it has NO first-exon exclusion,
                                  so this column mirrors what script 18
                                  actually keeps.

Writes: data/processed/mane_exon_atlas.tsv.gz (one row per exon).

Usage: PYTHONNOUSERSITE=1 .venv/bin/python scripts/20_build_exon_atlas.py
"""

from __future__ import annotations

import gzip
import re
import sys
from pathlib import Path

import pandas as pd

REPO = Path(__file__).resolve().parents[1]

GENCODE_URL = (
    "https://ftp.ebi.ac.uk/pub/databases/gencode/Gencode_human/"
    "release_46/gencode.v46.annotation.gtf.gz"
)
GENCODE_FALLBACK_URL = (
    "https://ftp.ebi.ac.uk/pub/databases/gencode/Gencode_human/"
    "release_46/gencode.v46.primary_assembly.annotation.gtf.gz"
)
MANE_URL = (
    "https://ftp.ncbi.nlm.nih.gov/refseq/MANE/MANE_human/"
    "release_1.5/MANE.GRCh38.v1.5.summary.txt.gz"
)

ATTR_RE = re.compile(r'(\S+) "([^"]*)";')


def parse_attrs(field: str) -> dict:
    return dict(ATTR_RE.findall(field))


def base_id(versioned: str) -> str:
    return str(versioned).split(".")[0]


def load_mane_select(mane_path: Path) -> pd.DataFrame:
    mane = pd.read_csv(mane_path, sep="\t", comment="#", dtype=str)
    # The file's header line starts with '#NCBI_GeneID'; pandas comment='#'
    # would eat it, so re-read properly if the columns look wrong.
    if "symbol" not in mane.columns:
        with gzip.open(mane_path, "rt") as fh:
            header = fh.readline().lstrip("#").rstrip("\n").split("\t")
        mane = pd.read_csv(
            mane_path, sep="\t", header=None, names=header, dtype=str, skiprows=1
        )
    sel = mane[mane["MANE_status"] == "MANE Select"].copy()
    sel["ensembl_tx_base"] = sel["Ensembl_nuc"].map(base_id)
    return sel


def extract_transcript_features(gtf_path: Path, wanted_tx: set[str]) -> pd.DataFrame:
    """Stream the GTF once, keeping exon/CDS rows for wanted transcripts."""
    rows = []
    opener = gzip.open if str(gtf_path).endswith(".gz") else open
    with opener(gtf_path, "rt") as fh:
        for line in fh:
            if line.startswith("#"):
                continue
            parts = line.rstrip("\n").split("\t")
            if len(parts) < 9 or parts[2] not in ("exon", "CDS"):
                continue
            attrs = parse_attrs(parts[8])
            tx = base_id(attrs.get("transcript_id", ""))
            if tx not in wanted_tx:
                continue
            rows.append(
                {
                    "chrom": parts[0],
                    "feature": parts[2],
                    "start": int(parts[3]),
                    "end": int(parts[4]),
                    "strand": parts[6],
                    "phase": parts[7],
                    "gene_id": attrs.get("gene_id", ""),
                    "gene_name": attrs.get("gene_name", ""),
                    "transcript_id": attrs.get("transcript_id", ""),
                    "exon_number": attrs.get("exon_number", ""),
                    "exon_id": attrs.get("exon_id", ""),
                }
            )
    return pd.DataFrame(rows)


def build_atlas(
    gtf_path: Path, mane_path: Path, panel_path: Path
) -> tuple[pd.DataFrame, dict]:
    panel = pd.read_csv(panel_path, sep="\t")
    panel_genes = panel["gene_symbol"].tolist()
    mane = load_mane_select(mane_path)

    mane_by_symbol = {}
    multi_select = []
    for sym, grp in mane.groupby("symbol"):
        if len(grp) > 1:
            multi_select.append((sym, len(grp)))
        mane_by_symbol[sym] = grp.iloc[0]

    missing_in_mane = [g for g in panel_genes if g not in mane_by_symbol]
    wanted_tx = {
        mane_by_symbol[g]["ensembl_tx_base"] for g in panel_genes if g in mane_by_symbol
    }
    feats = extract_transcript_features(gtf_path, wanted_tx)

    atlas_rows = []
    genes_missing_in_gtf = []
    for gene in panel_genes:
        if gene not in mane_by_symbol:
            genes_missing_in_gtf.append(gene)
            continue
        mrow = mane_by_symbol[gene]
        tx_base = mrow["ensembl_tx_base"]
        sub = feats[feats["transcript_id"].map(base_id) == tx_base]
        exons = sub[sub["feature"] == "exon"].copy()
        if len(exons) == 0:
            genes_missing_in_gtf.append(gene)
            continue
        cds = sub[sub["feature"] == "CDS"]
        # GENCODE numbers exons in transcript order, but a few exon rows
        # carry no exon_number attribute; fall back to genomic order in
        # transcript orientation for those.
        exons["exon_number"] = pd.to_numeric(exons["exon_number"], errors="coerce")
        if exons["exon_number"].isna().any():
            strand = str(exons["strand"].iloc[0])
            genomic_rank = exons["start"].rank(
                method="first", ascending=(strand == "+")
            )
            exons["exon_number"] = exons["exon_number"].fillna(genomic_rank)
        exons["exon_number"] = exons["exon_number"].astype(int)
        exons = exons.sort_values("exon_number")
        # CDS rank in transcript order is assigned after CDS overlap is known.
        cds_flags = []
        for _, ex in exons.iterrows():
            ov = cds[(cds["start"] <= ex["end"]) & (cds["end"] >= ex["start"])]
            if len(ov) == 0:
                cds_flags.append((None, None, 0, "."))
            else:
                c_start = int(max(ov["start"].min(), ex["start"]))
                c_end = int(min(ov["end"].max(), ex["end"]))
                # Total CDS bases inside this exon (CDS features of one
                # transcript do not overlap each other).
                c_len = int(
                    sum(
                        min(r["end"], ex["end"]) - max(r["start"], ex["start"]) + 1
                        for _, r in ov.iterrows()
                    )
                )
                # GTF CDS phase of the overlapping CDS feature (for a
                # minus-strand transcript this is the phase of the CDS
                # segment as annotated; kept verbatim).
                phase = str(ov.iloc[0]["phase"])
                cds_flags.append((c_start, c_end, c_len, phase))
        exons = exons.reset_index(drop=True)
        cds_exon_positions = [i for i, f in enumerate(cds_flags) if f[2] > 0]
        first_cds_pos = cds_exon_positions[0] if cds_exon_positions else -1
        last_cds_pos = cds_exon_positions[-1] if cds_exon_positions else -1
        rank_of_pos = {pos: r + 1 for r, pos in enumerate(cds_exon_positions)}
        for i, (_, ex) in enumerate(exons.iterrows()):
            c_start, c_end, c_len, phase = cds_flags[i]
            is_cds = c_len > 0
            is_first = i == first_cds_pos
            is_last = i == last_cds_pos
            atlas_rows.append(
                {
                    "gene_symbol": gene,
                    "ensembl_gene_id": mrow["Ensembl_Gene"],
                    "ensembl_transcript_id": tx_base,
                    "ensembl_transcript_versioned": mrow["Ensembl_nuc"],
                    "refseq_transcript": mrow["RefSeq_nuc"],
                    "mane_chr": mrow["GRCh38_chr"],
                    "mane_strand": mrow["chr_strand"],
                    "chrom": ex["chrom"],
                    "strand": ex["strand"],
                    "exon_number": int(ex["exon_number"]),
                    "exon_id": ex["exon_id"],
                    "start": int(ex["start"]),
                    "end": int(ex["end"]),
                    "exon_length": int(ex["end"] - ex["start"] + 1),
                    "cds_start": c_start,
                    "cds_end": c_end,
                    "cds_length": int(c_len),
                    "cds_phase_gtf": phase,
                    "is_cds_exon": int(is_cds),
                    "cds_exon_rank": int(rank_of_pos.get(i, 0)),
                    "is_first_cds_exon": int(is_cds and is_first),
                    "is_last_cds_exon": int(is_cds and is_last),
                    "is_internal_coding": int(is_cds and not is_first and not is_last),
                    "is_internal_coding_s18rule": int(is_cds and not is_last),
                }
            )
    atlas = pd.DataFrame(atlas_rows)
    info = {
        "panel_genes": len(panel_genes),
        "genes_missing_in_mane_summary": missing_in_mane,
        "genes_missing_in_gtf": genes_missing_in_gtf,
        "genes_with_multiple_mane_select_rows": multi_select,
        "genes_covered": int(atlas["gene_symbol"].nunique()) if len(atlas) else 0,
    }
    return atlas, info


def main() -> int:
    raw = REPO / "data/raw"
    gtf_path = raw / "gencode.v46.annotation.gtf.gz"
    gtf_kind = "comprehensive"
    if not gtf_path.exists():
        gtf_path = raw / "gencode.v46.primary_assembly.annotation.gtf.gz"
        gtf_kind = "primary_assembly (fallback)"
    mane_path = raw / "MANE.GRCh38.v1.5.summary.txt.gz"
    for p in (gtf_path, mane_path):
        if not p.exists():
            print(f"MISSING download: {p}", file=sys.stderr)
            return 2
    print(f"[20] GTF: {gtf_path.name} ({gtf_kind})")
    print(f"[20] MANE summary: {mane_path.name}")
    atlas, info = build_atlas(gtf_path, mane_path, REPO / "metadata/gene_panel.tsv")
    out = REPO / "data/processed/mane_exon_atlas.tsv.gz"
    atlas.to_csv(out, sep="\t", index=False, compression="gzip")
    print(f"[20] wrote {out} ({len(atlas):,} exon rows)")
    print(f"[20] genes covered: {info['genes_covered']}/{info['panel_genes']}")
    print(f"[20] missing in MANE summary: {info['genes_missing_in_mane_summary']}")
    print(f"[20] missing in GTF: {info['genes_missing_in_gtf']}")
    print(
        f"[20] multi MANE Select rows: {info['genes_with_multiple_mane_select_rows']}"
    )
    print(f"[20] total exons: {len(atlas):,}")
    print(f"[20] CDS exons: {int(atlas['is_cds_exon'].sum()):,}")
    print(
        f"[20] internal coding (exclude first+last CDS): "
        f"{int(atlas['is_internal_coding'].sum()):,}"
    )
    print(
        f"[20] internal coding (script-18 rule, exclude last CDS only): "
        f"{int(atlas['is_internal_coding_s18rule'].sum()):,}"
    )
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
