"""Tests for variant normalization and coordinate utilities."""

import pytest

from isoformvus.utils.variants import (
    Variant,
    batch_deduplicate,
    canonical_contig,
    from_zero_based_half_open,
    overlaps,
    parse_variant_string,
    to_zero_based_half_open,
    trim_alleles,
    vcf_to_variant,
)


def test_canonical_contig():
    assert canonical_contig("1") == "chr1"
    assert canonical_contig("chrX") == "chrX"
    assert canonical_contig("MT") == "chrM"
    with pytest.raises(ValueError):
        canonical_contig("chrUn_gl000220")
    with pytest.raises(ValueError):
        canonical_contig("GL000207.1")


def test_trim_alleles_suffix_then_prefix():
    # "CAGT" vs "CGT": shared suffix "GT" trimmed -> ("CA", "C", 0).
    # Prefix trim is blocked because alt is a single anchor base.
    assert trim_alleles("CAGT", "CGT") == ("CA", "C", 0)


def test_trim_alleles_indel_anchor():
    # Anchor base retained: ("AAG", "AG") -> ("AA", "A", 0).
    assert trim_alleles("AAG", "AG") == ("AA", "A", 0)
    # True shared prefix shifts the position.
    assert trim_alleles("GATT", "GCTT") == ("A", "C", 1)


def test_parse_variant_string():
    v = parse_variant_string("chr7:117199645:G>A")
    assert v.chrom == "chr7" and v.pos == 117199645
    assert v.ref == "G" and v.alt == "A"
    assert v.variant_type == "SNV"
    with pytest.raises(ValueError):
        parse_variant_string("not-a-variant")


def test_vcf_to_variant_trims_anchor():
    v = vcf_to_variant("1", 1000, "CA", "C")
    assert v.ref == "CA" or v.ref == "C"
    assert v.pos >= 1000


def test_variant_validation():
    with pytest.raises(ValueError):
        Variant(chrom="chr1", pos=0, ref="A", alt="G")
    with pytest.raises(ValueError):
        Variant(chrom="chr1", pos=10, ref="A", alt="B")  # non-IUPAC


def test_coordinate_roundtrip():
    v = parse_variant_string("chr2:500:AT>A")
    chrom, start0, end0 = to_zero_based_half_open(v)
    assert (chrom, start0, end0) == ("chr2", 499, 501)
    v2 = from_zero_based_half_open("chr2", 499, 501, ref="AT", alt="A")
    assert v2.pos == 500 and v2.end == 501


def test_overlaps():
    v = parse_variant_string("chr3:1000:A>T")
    assert overlaps(v, "chr3", 990, 1010)
    assert not overlaps(v, "chr3", 1001, 1010)
    assert not overlaps(v, "chr4", 990, 1010)


def test_batch_deduplicate():
    vs = [parse_variant_string("chr1:10:A>G"), parse_variant_string("chr1:10:A>G"),
          parse_variant_string("chr1:11:A>G")]
    assert len(batch_deduplicate(vs)) == 2
