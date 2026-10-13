"""Transcript structures and isoform-specific variant-effect classification.

Given a variant and the set of transcripts of its gene, classify *per transcript*
whether the variant overlaps the transcript at all, and if so whether the overlap
is shared across isoforms or isoform-specific (exon / coding region / splice
junction).

This module is purely structural: it answers "where does this variant land in
each annotated transcript?" It does NOT predict functional impact — that is the
job of the feature layer (SpliceAI, consequence annotation, etc.), which consumes
these structural labels.
"""

from __future__ import annotations

from dataclasses import dataclass, field
from enum import Enum

from ..utils.variants import Variant, overlaps


class IsoformEffect(str, Enum):
    """Per-transcript structural effect of a variant."""

    NONE = "none"                        # variant does not overlap this transcript
    SHARED = "shared"                    # overlaps an exon/region present in all transcripts
    SOME = "some"                        # overlaps the transcript but not all isoforms
    EXON_SPECIFIC = "exon_specific"      # overlaps an exon unique to a subset of isoforms
    CODING_SPECIFIC = "coding_specific"  # overlaps CDS unique to a subset of isoforms
    JUNCTION_SPECIFIC = "junction_specific"  # within splice-junction window of a subset-specific junction
    INTRONIC = "intronic"                # inside transcript bounds but not in any exon


@dataclass(frozen=True)
class Exon:
    start: int  # 1-based inclusive
    end: int    # 1-based inclusive
    frame: int = -1  # CDS frame offset, -1 if non-coding / unknown

    def __post_init__(self) -> None:
        if not (1 <= self.start <= self.end):
            raise ValueError(f"Invalid exon interval [{self.start}, {self.end}]")


@dataclass
class Transcript:
    """One annotated transcript (e.g. a GENCODE transcript of a gene)."""

    transcript_id: str          # versioned, e.g. "ENST00000335137.4"
    gene_id: str
    gene_symbol: str
    chrom: str
    strand: str                 # "+" or "-"
    exons: list[Exon] = field(default_factory=list)
    cds_start: int | None = None  # 1-based inclusive; None if non-coding
    cds_end: int | None = None
    biotype: str = "protein_coding"
    is_mane_select: bool = False
    tag: str = ""               # e.g. "dominant" / "supported" / "low_usage" / "uncertain"

    def __post_init__(self) -> None:
        if self.strand not in ("+", "-"):
            raise ValueError("strand must be '+' or '-'")
        self.exons = sorted(self.exons, key=lambda e: e.start)

    @property
    def tx_start(self) -> int:
        return self.exons[0].start

    @property
    def tx_end(self) -> int:
        return self.exons[-1].end

    def exon_index_at(self, var: Variant) -> int | None:
        """0-based exon index overlapped by the variant, or None."""
        if var.chrom != self.chrom:
            return None
        for i, exon in enumerate(self.exons):
            if overlaps(var, self.chrom, exon.start, exon.end):
                return i
        return None

    def in_cds(self, var: Variant) -> bool:
        if self.cds_start is None or self.cds_end is None:
            return False
        return overlaps(var, self.chrom, self.cds_start, self.cds_end)


# Splice-junction window: variants within this many bp of an exon boundary are
# flagged as junction-proximal (SpliceAI itself uses wider context; this window
# is only for *structural* classification, not effect prediction).
JUNCTION_WINDOW_BP = 10


@dataclass(frozen=True)
class TranscriptEffect:
    transcript_id: str
    effect: IsoformEffect
    exon_index: int | None = None
    junction_proximal: bool = False
    in_cds: bool = False


def _exon_signature(exon: Exon) -> tuple[int, int]:
    return (exon.start, exon.end)


def classify_transcript_effects(var: Variant, transcripts: list[Transcript]) -> list[TranscriptEffect]:
    """Classify the structural effect of ``var`` on each transcript.

    Algorithm
    ---------
    1. For each transcript, find overlapped exons / intronic / outside.
    2. An exon is "shared" if its (start, end) signature appears in *every*
       transcript of the gene; otherwise it is isoform-specific.
    3. A variant is "shared" if it overlaps only shared exons in every
       transcript it touches; "exon_specific"/"coding_specific" if it touches
       an exon (or CDS) absent from at least one isoform.
    4. Junction proximity is evaluated per transcript against that transcript's
       own exon boundaries.

    Returns one :class:`TranscriptEffect` per input transcript, in input order.
    """
    if not transcripts:
        return []
    # Exon signatures present in every transcript -> shared exons.
    sig_counts: dict[tuple[int, int], int] = {}
    for tx in transcripts:
        for exon in tx.exons:
            sig = _exon_signature(exon)
            sig_counts[sig] = sig_counts.get(sig, 0) + 1
    n_tx = len(transcripts)

    results: list[TranscriptEffect] = []
    for tx in transcripts:
        if var.chrom != tx.chrom or not (tx.tx_start <= var.end and var.pos <= tx.tx_end):
            results.append(TranscriptEffect(tx.transcript_id, IsoformEffect.NONE))
            continue
        exon_idx = tx.exon_index_at(var)
        if exon_idx is None:
            results.append(TranscriptEffect(tx.transcript_id, IsoformEffect.INTRONIC,
                                           in_cds=tx.in_cds(var)))
            continue
        exon = tx.exons[exon_idx]
        shared = sig_counts[_exon_signature(exon)] == n_tx
        # Junction proximity: within window of either boundary of this exon.
        junction_proximal = (
            abs(var.pos - exon.start) <= JUNCTION_WINDOW_BP
            or abs(var.end - exon.end) <= JUNCTION_WINDOW_BP
            or abs(var.pos - exon.end) <= JUNCTION_WINDOW_BP
            or abs(var.end - exon.start) <= JUNCTION_WINDOW_BP
        )
        in_cds = tx.in_cds(var)
        if shared:
            effect = IsoformEffect.SHARED
        else:
            # Is the *coding* overlap also isoform-specific? Conservative check:
            # coding-specific if the variant is in CDS here and this exon is not shared.
            effect = IsoformEffect.CODING_SPECIFIC if in_cds else IsoformEffect.EXON_SPECIFIC
            if junction_proximal:
                # Junction specificity takes precedence for splicing interpretation,
                # but we keep the exon label recoverable via the flags.
                effect = IsoformEffect.JUNCTION_SPECIFIC
        results.append(TranscriptEffect(
            transcript_id=tx.transcript_id, effect=effect, exon_index=exon_idx,
            junction_proximal=junction_proximal, in_cds=in_cds))
    return results


def summarize_gene_effect(effects: list[TranscriptEffect]) -> str:
    """Gene-level summary: 'all', 'some', 'none', or 'isoform_specific'.

    - ``all``: variant affects every annotated transcript (SHARED everywhere it lands).
    - ``isoform_specific``: at least one transcript carries an isoform-specific
      effect (exon/coding/junction) — the counterfactual-relevant case.
    - ``some``: touches a strict subset of transcripts without isoform-specific exons.
    - ``none``: no transcript overlap.
    """
    touched = [e for e in effects if e.effect != IsoformEffect.NONE]
    if not touched:
        return "none"
    specific = {IsoformEffect.EXON_SPECIFIC, IsoformEffect.CODING_SPECIFIC,
                IsoformEffect.JUNCTION_SPECIFIC}
    if any(e.effect in specific for e in touched):
        return "isoform_specific"
    if all(e.effect == IsoformEffect.SHARED for e in touched) and len(touched) == len(effects):
        return "all"
    return "some"
