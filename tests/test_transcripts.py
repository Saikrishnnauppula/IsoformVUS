"""Tests for isoform-specific variant-effect classification."""

from isoformvus.data.transcripts import (
    Exon,
    IsoformEffect,
    Transcript,
    classify_transcript_effects,
    summarize_gene_effect,
)
from isoformvus.utils.variants import parse_variant_string


def _gene_two_isoforms():
    # Isoform 1: exons [100-200], [300-400], [500-600]
    # Isoform 2: exons [100-200], [500-600]  (skips the middle exon)
    tx1 = Transcript(
        transcript_id="ENST000001.1", gene_id="ENSG000001", gene_symbol="GENEX",
        chrom="chr1", strand="+",
        exons=[Exon(100, 200), Exon(300, 400), Exon(500, 600)],
        cds_start=100, cds_end=600,
    )
    tx2 = Transcript(
        transcript_id="ENST000002.1", gene_id="ENSG000001", gene_symbol="GENEX",
        chrom="chr1", strand="+",
        exons=[Exon(100, 200), Exon(500, 600)],
        cds_start=100, cds_end=600,
    )
    return [tx1, tx2]


def test_shared_exon_variant():
    txs = _gene_two_isoforms()
    v = parse_variant_string("chr1:150:A>G")  # in shared exon [100-200]
    effs = classify_transcript_effects(v, txs)
    assert all(e.effect == IsoformEffect.SHARED for e in effs)
    assert summarize_gene_effect(effs) == "all"


def test_isoform_specific_exon_variant():
    txs = _gene_two_isoforms()
    v = parse_variant_string("chr1:350:A>G")  # in skipped exon [300-400]
    effs = classify_transcript_effects(v, txs)
    by_tx = {e.transcript_id: e for e in effs}
    assert by_tx["ENST000001.1"].effect == IsoformEffect.CODING_SPECIFIC
    # tx2 skips the exon: position 350 lies within its bounds but in an intron.
    assert by_tx["ENST000002.1"].effect == IsoformEffect.INTRONIC
    assert summarize_gene_effect(effs) == "isoform_specific"


def test_intronic_variant():
    txs = _gene_two_isoforms()
    v = parse_variant_string("chr1:250:A>G")  # intron between exon 1 and 2
    effs = classify_transcript_effects(v, txs)
    by_tx = {e.transcript_id: e for e in effs}
    assert by_tx["ENST000001.1"].effect == IsoformEffect.INTRONIC
    assert by_tx["ENST000002.1"].effect == IsoformEffect.INTRONIC


def test_junction_proximal_flagged():
    txs = _gene_two_isoforms()
    v = parse_variant_string("chr1:305:A>G")  # 5bp inside skipped exon -> junction proximal
    effs = classify_transcript_effects(v, txs)
    by_tx = {e.transcript_id: e for e in effs}
    assert by_tx["ENST000001.1"].junction_proximal is True
    assert by_tx["ENST000001.1"].effect == IsoformEffect.JUNCTION_SPECIFIC


def test_outside_gene():
    txs = _gene_two_isoforms()
    v = parse_variant_string("chr1:900:A>G")
    effs = classify_transcript_effects(v, txs)
    assert all(e.effect == IsoformEffect.NONE for e in effs)
    assert summarize_gene_effect(effs) == "none"


def test_empty_transcript_list():
    v = parse_variant_string("chr1:150:A>G")
    assert classify_transcript_effects(v, []) == []
