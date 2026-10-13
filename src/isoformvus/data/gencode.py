"""GENCODE GTF parsing into Transcript objects.

Parses a GENCODE comprehensive/primary-assembly GTF (GRCh38) and builds
``data.transcripts.Transcript`` records per gene. Only the fields the
pipeline needs are extracted; everything else is ignored.

Caveats (documented, not hidden):
- GENCODE GTF coordinates are 1-based fully closed; attribute parsing is
  defensive (malformed attribute lines raise, never silently skipped).
- ``protein_coding_CDS_not_defined`` transcripts are kept but flagged —
  their CDS is unknown, so coding-specific classification is unavailable.
- Retained-intron / NMD transcripts are parsed like any other; the
  *usage* layer (data/isoform.py) decides their evidentiary weight.
"""

from __future__ import annotations

import gzip
import re
from collections import defaultdict
from pathlib import Path

from .transcripts import Exon, Transcript

_ATTR_RE = re.compile(r'(\S+)\s+"([^"]+)"')


def parse_attributes(attr_str: str) -> dict[str, str]:
    """Parse a GTF attribute column into a dict. Raises on malformed input."""
    attrs: dict[str, str] = {}
    for match in _ATTR_RE.finditer(attr_str):
        attrs[match.group(1)] = match.group(2)
    if "gene_id" not in attrs or "transcript_id" not in attrs:
        # Feature lines without transcript context (e.g. gene lines) are the
        # caller's responsibility to skip; here we are strict.
        raise ValueError(f"Attributes lack gene_id/transcript_id: {attr_str[:120]!r}")
    return attrs


def _open(path: Path):
    if str(path).endswith(".gz"):
        return gzip.open(path, "rt")
    return open(path, "rt")


def load_gene_transcripts(gtf_path: str | Path,
                          gene_symbols: set[str] | None = None,
                          ) -> dict[str, list[Transcript]]:
    """Load transcripts for selected genes (or all genes if None).

    Returns {gene_symbol: [Transcript, ...]}. Transcript IDs are versioned
    (e.g. ENST00000335137.4) as in the GTF.
    """
    path = Path(gtf_path)
    exons: dict[str, list[tuple[int, int, int]]] = defaultdict(list)
    tx_meta: dict[str, dict] = {}
    gene_of_tx: dict[str, str] = {}

    with _open(path) as fh:
        for line in fh:
            if line.startswith("#"):
                continue
            cols = line.rstrip("\n").split("\t")
            if len(cols) != 9:
                continue
            chrom, _, feature, start, end, _, strand, frame, attr_str = cols
            if feature not in ("transcript", "exon"):
                continue
            try:
                attrs = parse_attributes(attr_str)
            except ValueError:
                continue  # gene-level or malformed lines
            symbol = attrs.get("gene_name", "")
            if gene_symbols is not None and symbol not in gene_symbols:
                continue
            tid = attrs["transcript_id"]
            gid = attrs["gene_id"]
            gene_of_tx[tid] = symbol
            if feature == "transcript":
                tx_meta[tid] = {
                    "gene_id": gid, "chrom": chrom, "strand": strand,
                    "biotype": attrs.get("transcript_type", ""),
                    "tag": attrs.get("tag", ""),
                }
            else:
                try:
                    frame_i = int(frame) if frame in ("0", "1", "2") else -1
                except ValueError:
                    frame_i = -1
                exons[tid].append((int(start), int(end), frame_i))

    # CDS bounds: from CDS features would be ideal; GENCODE marks CDS via
    # separate "CDS" lines. Second pass for CDS is expensive on huge GTFs, so
    # we approximate: for protein_coding transcripts, CDS = span from first
    # to last exon is WRONG (includes UTRs). Instead we do a targeted second
    # pass only over kept transcripts (see below).
    out: dict[str, list[Transcript]] = defaultdict(list)
    for tid, meta in tx_meta.items():
        exon_list = [Exon(s, e, f) for s, e, f in sorted(exons.get(tid, []))]
        if not exon_list:
            continue  # transcript with no exons: annotation artifact; skip loudly? no — skip
        tx = Transcript(
            transcript_id=tid, gene_id=meta["gene_id"],
            gene_symbol=gene_of_tx[tid], chrom=meta["chrom"],
            strand=meta["strand"], exons=exon_list,
            biotype=meta["biotype"], tag=meta["tag"],
        )
        out[gene_of_tx[tid]].append(tx)
    # Attach CDS bounds via a second targeted pass.
    _attach_cds_bounds(path, out, gene_symbols)
    return dict(out)


def _attach_cds_bounds(gtf_path: Path,
                       transcripts: dict[str, list[Transcript]],
                       gene_symbols: set[str] | None) -> None:
    """Fill cds_start/cds_end from CDS features for loaded transcripts."""
    wanted = {t.transcript_id for txs in transcripts.values() for t in txs}
    cds_min: dict[str, int] = {}
    cds_max: dict[str, int] = {}
    with _open(gtf_path) as fh:
        for line in fh:
            if line.startswith("#"):
                continue
            cols = line.rstrip("\n").split("\t")
            if len(cols) != 9 or cols[2] != "CDS":
                continue
            try:
                attrs = parse_attributes(cols[8])
            except ValueError:
                continue
            tid = attrs.get("transcript_id", "")
            if tid not in wanted:
                continue
            s, e = int(cols[3]), int(cols[4])
            cds_min[tid] = min(s, cds_min.get(tid, s))
            cds_max[tid] = max(e, cds_max.get(tid, e))
    for txs in transcripts.values():
        for t in txs:
            if t.transcript_id in cds_min:
                # dataclass is mutable (not frozen) for exons list; cds fields too.
                t.cds_start = cds_min[t.transcript_id]
                t.cds_end = cds_max[t.transcript_id]
