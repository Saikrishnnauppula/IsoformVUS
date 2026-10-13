"""Tests for leakage-safe splitting."""

import json

import pytest

from isoformvus.data.splits import (
    assign_splits,
    splits_to_json,
    verify_no_gene_overlap,
    verify_variant_split_integrity,
)


GENES = [f"GENE{i:03d}" for i in range(1, 51)]
CATS = {g: ("neuro" if i % 2 == 0 else "cardiac") for i, g in enumerate(GENES)}


def test_no_gene_in_two_splits():
    a = assign_splits(GENES, CATS)
    verify_no_gene_overlap(a)  # raises on violation
    splits = {x.split for x in a}
    assert splits == {"train", "validation", "test"}


def test_deterministic():
    a1 = assign_splits(GENES, CATS, seed=42)
    a2 = assign_splits(GENES, CATS, seed=42)
    assert [(x.gene_symbol, x.split) for x in a1] == [(x.gene_symbol, x.split) for x in a2]
    a3 = assign_splits(GENES, CATS, seed=43)
    assert [(x.gene_symbol, x.split) for x in a3] != [(x.gene_symbol, x.split) for x in a1]


def test_stratification_balanced():
    a = assign_splits(GENES, CATS)
    for split in ("train", "validation", "test"):
        cats = [x.disease_category for x in a if x.split == split]
        assert "neuro" in cats and "cardiac" in cats


def test_duplicate_genes_rejected():
    with pytest.raises(ValueError):
        assign_splits(["G1", "G1"])


def test_variant_integrity_check():
    a = assign_splits(GENES, CATS)
    gene_split = {x.gene_symbol: x.split for x in a}
    variant_gene = {f"var{i}": GENES[i % len(GENES)] for i in range(200)}
    assert verify_variant_split_integrity(variant_gene, gene_split) == []
    bad = dict(variant_gene)
    bad["varX"] = "NOT_A_GENE"
    assert len(verify_variant_split_integrity(bad, gene_split)) == 1


def test_json_roundtrip():
    a = assign_splits(GENES, CATS)
    payload = json.loads(splits_to_json(a))
    assert len(payload["assignments"]) == len(GENES)
