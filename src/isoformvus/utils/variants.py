"""Variant normalization and coordinate utilities.

Genome build: GRCh38 throughout. All coordinates are 1-based, fully-closed
intervals (VCF/ClinVar convention) unless a function documents otherwise.

Design notes
------------
* Normalization here is *syntactic* (parsing, contig canonicalization, allele
  trimming). True left-alignment / parsimonious normalization requires a
  reference FASTA; that lives in ``scripts/`` where the reference is available,
  not in this import-safe module.
* No function here invents data: missing reference context raises or returns a
  documented sentinel, never a guess.
"""

from __future__ import annotations

from dataclasses import dataclass
from typing import Iterable, Optional


# Canonical contig names: "1" -> "chr1". GRCh38 primary assembly.
_AUTOSOMES = [str(i) for i in range(1, 23)]
_CONTIGS = _AUTOSOMES + ["X", "Y", "M"]
_ALIAS = {"MT": "M", "chrM": "M", "chrMT": "M"}


def canonical_contig(contig: str) -> str:
    """Return ``chr``-prefixed canonical contig name for GRCh38.

    Raises:
        ValueError: if the contig is not a canonical GRCh38 chromosome.
    """
    c = contig.strip()
    c = _ALIAS.get(c, c)
    if c.startswith("chr"):
        c = c[3:]
    c = _ALIAS.get(c, c)
    if c not in _CONTIGS:
        raise ValueError(f"Non-canonical contig for GRCh38: {contig!r}")
    return f"chr{c}"


@dataclass(frozen=True)
class Variant:
    """A normalized genomic variant (1-based, fully closed, GRCh38)."""

    chrom: str   # canonical, e.g. "chr1"
    pos: int     # 1-based position of first reference base
    ref: str     # reference allele (upper-case)
    alt: str     # alternate allele (upper-case)

    def __post_init__(self) -> None:
        if self.pos < 1:
            raise ValueError("pos must be >= 1")
        if not self.ref or not self.alt:
            raise ValueError("ref and alt must be non-empty")
        if any(b not in "ACGTN" for b in self.ref + self.alt):
            raise ValueError("alleles must be IUPAC DNA (ACGTN)")

    @property
    def end(self) -> int:
        """1-based inclusive end of the reference allele."""
        return self.pos + len(self.ref) - 1

    @property
    def variant_type(self) -> str:
        if len(self.ref) == 1 and len(self.alt) == 1:
            return "SNV"
        if len(self.ref) > len(self.alt):
            return "deletion"
        if len(self.alt) > len(self.ref):
            return "insertion"
        return "complex"

    @property
    def key(self) -> str:
        return f"{self.chrom}:{self.pos}:{self.ref}>{self.alt}"


def trim_alleles(ref: str, alt: str) -> tuple[str, str, int]:
    """Trim shared suffix/prefix bases shared by ref and alt.

    Returns ``(ref_trimmed, alt_trimmed, left_shift)`` where ``left_shift`` is
    the number of 5' bases removed (caller must add it to pos). Shared-prefix
    trimming is right-to-left minimal; callers needing full left-alignment
    must use a reference FASTA (see module docstring).
    """
    ref, alt = ref.upper(), alt.upper()
    shift = 0
    # Trim common suffix first (keeps 5' anchor stable).
    while len(ref) > 1 and len(alt) > 1 and ref[-1] == alt[-1]:
        ref, alt = ref[:-1], alt[:-1]
    # Trim common prefix.
    while len(ref) > 1 and len(alt) > 1 and ref[0] == alt[0]:
        ref, alt = ref[1:], alt[1:]
        shift += 1
    return ref, alt, shift


def parse_variant_string(s: str) -> Variant:
    """Parse ``chr:pos:ref>alt`` (e.g. ``chr7:117199645:G>A``)."""
    try:
        chrom, pos, alleles = s.split(":")
        ref, alt = alleles.split(">")
    except ValueError as exc:
        raise ValueError(f"Cannot parse variant string {s!r}; expected chr:pos:ref>alt") from exc
    ref_t, alt_t, shift = trim_alleles(ref, alt)
    return Variant(chrom=canonical_contig(chrom), pos=int(pos) + shift, ref=ref_t, alt=alt_t)


def vcf_to_variant(chrom: str, pos: int, ref: str, alt: str) -> Variant:
    """Build a normalized :class:`Variant` from VCF-style fields (1-based POS)."""
    ref_t, alt_t, shift = trim_alleles(ref, alt)
    return Variant(chrom=canonical_contig(chrom), pos=int(pos) + shift, ref=ref_t, alt=alt_t)


def to_zero_based_half_open(var: Variant) -> tuple[str, int, int]:
    """Convert to 0-based half-open ``(chrom, start, end)`` (BED convention)."""
    return var.chrom, var.pos - 1, var.end


def from_zero_based_half_open(chrom: str, start: int, end: int,
                              ref: str, alt: str) -> Variant:
    """Build a Variant from 0-based half-open coordinates.

    Note: ``ref`` must be the actual reference sequence over ``[start, end)``;
    this function cannot verify that without a FASTA — pass ``"N"``-padded
    alleles only if you document the missing reference explicitly upstream.
    """
    if not (0 <= start < end):
        raise ValueError(f"Invalid half-open interval [{start}, {end})")
    return Variant(chrom=canonical_contig(chrom), pos=start + 1, ref=ref.upper(), alt=alt.upper())


def overlaps(var: Variant, chrom: str, start_1b: int, end_1b: int) -> bool:
    """Whether the variant's reference span overlaps a 1-based closed interval."""
    return canonical_contig(chrom) == var.chrom and var.pos <= end_1b and var.end >= start_1b


def batch_deduplicate(variants: Iterable[Variant]) -> list[Variant]:
    """Deduplicate preserving first-seen order."""
    seen: set[str] = set()
    out: list[Variant] = []
    for v in variants:
        if v.key not in seen:
            seen.add(v.key)
            out.append(v)
    return out
