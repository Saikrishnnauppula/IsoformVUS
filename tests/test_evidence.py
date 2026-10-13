"""Tests for phenotype->cell mapping and isoform usage categorization."""

import pytest

from isoformvus.data.isoform import (
    IsoformUsage,
    UsageCategory,
    categorize_usage,
    check_usage_sums,
)
from isoformvus.data.phenotype_cell import (
    CellTypeRelevance,
    HpoTerm,
    combine_evidence,
    normalize_score,
    rank_cell_types,
)


def test_hpo_term_validation():
    HpoTerm(hpo_id="HP:0001250", name="Seizure")
    with pytest.raises(ValueError):
        HpoTerm(hpo_id="0001250", name="Seizure")


def test_combine_evidence_weighted():
    s = combine_evidence(phenotype_tissue=1.0, gene_expression=0.0)
    assert s == pytest.approx(0.6)  # w1=0.6, w2=0.4 renormalized
    assert normalize_score(2.5) == 1.0
    assert normalize_score(-1.0) == 0.0
    with pytest.raises(ValueError):
        combine_evidence(0.5, 0.5, w1=0.0, w2=0.0)


def test_rank_cell_types():
    cs = [
        CellTypeRelevance("CL:1", "a", 0.2),
        CellTypeRelevance("CL:2", "b", 0.9),
        CellTypeRelevance("CL:3", "c", 0.5),
    ]
    ranked = rank_cell_types(cs)
    assert [c.cell_type_name for c in ranked] == ["b", "c", "a"]
    assert len(rank_cell_types(cs, top_k=2)) == 2


def test_categorize_usage():
    assert categorize_usage(0.8, 0.1, 2, False) == UsageCategory.DOMINANT
    assert categorize_usage(0.8, 0.1, 1, True) == UsageCategory.DOMINANT  # long-read backs it
    assert categorize_usage(0.3, 0.1, 1, False) == UsageCategory.SUPPORTED
    assert categorize_usage(0.02, 0.01, 3, False) == UsageCategory.LOW_USAGE
    assert categorize_usage(0.02, 0.01, 1, False) == UsageCategory.UNCERTAIN
    assert categorize_usage(0.5, 0.2, 0, False) == UsageCategory.UNCERTAIN  # no datasets


def test_isoform_usage_validation():
    with pytest.raises(ValueError):
        IsoformUsage("G", "T", "CL:1", "c", usage=1.5)


def test_check_usage_sums():
    us = [
        IsoformUsage("GENEX", "T1", "CL:1", "neuron", 0.7, n_datasets=2),
        IsoformUsage("GENEX", "T2", "CL:1", "neuron", 0.3, n_datasets=2),
        IsoformUsage("GENEY", "T1", "CL:1", "neuron", 0.4, n_datasets=1),
    ]
    warnings = check_usage_sums(us)
    assert len(warnings) == 1 and "GENEY" in warnings[0]
