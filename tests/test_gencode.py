"""Tests for GENCODE GTF parsing."""

import gzip
import tempfile
from pathlib import Path

from isoformvus.data.gencode import load_gene_transcripts, parse_attributes
from isoformvus.data.transcripts import IsoformEffect, classify_transcript_effects
from isoformvus.utils.variants import parse_variant_string

TINY_GTF = """\
chr1\tHAVANA\ttranscript\t100\t600\t.\t+\t.\tgene_id "ENSG000001"; transcript_id "ENST000001.1"; gene_name "GENEX"; transcript_type "protein_coding";
chr1\tHAVANA\texon\t100\t200\t.\t+\t.\tgene_id "ENSG000001"; transcript_id "ENST000001.1"; gene_name "GENEX"; transcript_type "protein_coding"; exon_number "1";
chr1\tHAVANA\texon\t300\t400\t.\t+\t.\tgene_id "ENSG000001"; transcript_id "ENST000001.1"; gene_name "GENEX"; transcript_type "protein_coding"; exon_number "2";
chr1\tHAVANA\texon\t500\t600\t.\t+\t.\tgene_id "ENSG000001"; transcript_id "ENST000001.1"; gene_name "GENEX"; transcript_type "protein_coding"; exon_number "3";
chr1\tHAVANA\tCDS\t120\t390\t.\t+\t0\tgene_id "ENSG000001"; transcript_id "ENST000001.1"; gene_name "GENEX"; transcript_type "protein_coding";
chr1\tHAVANA\ttranscript\t100\t600\t.\t+\t.\tgene_id "ENSG000001"; transcript_id "ENST000002.1"; gene_name "GENEX"; transcript_type "protein_coding";
chr1\tHAVANA\texon\t100\t200\t.\t+\t.\tgene_id "ENSG000001"; transcript_id "ENST000002.1"; gene_name "GENEX"; transcript_type "protein_coding"; exon_number "1";
chr1\tHAVANA\texon\t500\t600\t.\t+\t.\tgene_id "ENSG000001"; transcript_id "ENST000002.1"; gene_name "GENEX"; transcript_type "protein_coding"; exon_number "2";
chr1\tHAVANA\tCDS\t120\t190\t.\t+\t0\tgene_id "ENSG000001"; transcript_id "ENST000002.1"; gene_name "GENEX"; transcript_type "protein_coding";
chr1\tHAVANA\ttranscript\t700\t900\t.\t+\t.\tgene_id "ENSG000002"; transcript_id "ENST000003.1"; gene_name "OTHER"; transcript_type "protein_coding";
chr1\tHAVANA\texon\t700\t900\t.\t+\t.\tgene_id "ENSG000002"; transcript_id "ENST000003.1"; gene_name "OTHER"; transcript_type "protein_coding"; exon_number "1";
"""


def _write_gtf(tmp_path: Path, gz: bool = True) -> Path:
    p = tmp_path / ("tiny.gtf.gz" if gz else "tiny.gtf")
    if gz:
        with gzip.open(p, "wt") as fh:
            fh.write(TINY_GTF)
    else:
        p.write_text(TINY_GTF)
    return p


def test_parse_attributes():
    a = parse_attributes('gene_id "ENSG1"; transcript_id "ENST1.1"; gene_name "GENEX";')
    assert a == {"gene_id": "ENSG1", "transcript_id": "ENST1.1", "gene_name": "GENEX"}


def test_load_gene_transcripts_filtered(tmp_path):
    p = _write_gtf(tmp_path)
    out = load_gene_transcripts(p, gene_symbols={"GENEX"})
    assert set(out) == {"GENEX"}
    assert len(out["GENEX"]) == 2
    tx1 = next(t for t in out["GENEX"] if t.transcript_id == "ENST000001.1")
    assert len(tx1.exons) == 3
    assert tx1.cds_start == 120 and tx1.cds_end == 390
    assert tx1.biotype == "protein_coding"


def test_end_to_end_effect_classification(tmp_path):
    p = _write_gtf(tmp_path, gz=False)
    txs = load_gene_transcripts(p, gene_symbols={"GENEX"})["GENEX"]
    v = parse_variant_string("chr1:350:A>G")  # skipped exon in tx2
    effs = {e.transcript_id: e for e in classify_transcript_effects(v, txs)}
    assert effs["ENST000001.1"].effect == IsoformEffect.CODING_SPECIFIC
    assert effs["ENST000001.1"].in_cds is True
    assert effs["ENST000002.1"].effect == IsoformEffect.INTRONIC
