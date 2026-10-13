"""Phenotype -> cell-type relevance mapping.

Pipeline: HPO terms -> disease/phenotype relevance -> tissues -> cell types,
producing a *cell-relevance score* per (phenotype-set, cell type) pair.

The score is built from measurable evidence, not hard-coded disease->tissue
rules:

    cell_relevance(hpo_terms, cell_type) =
        w1 * phenotype_tissue_match   (HPO/phenotype -> tissue via curated
                                       phenotype-tissue annotations + GTEx)
      + w2 * gene_expression_support (panel genes' expression in that cell type
                                       from single-cell data)
      + w3 * literature_cooccurrence (placeholder: 0 until a text-mining source
                                       is wired; never silently imputed)

Weights (w1, w2) are set in config and *ablated* (Model 3 vs Model 2 tests
whether this layer adds signal at all).

Evidence tiers for every mapping edge are recorded: 'curated' (HPO/Cell
Ontology annotations), 'expression' (measured RNA), 'inferred' (ontology
propagation, e.g. HPO ancestor terms). Inferred edges are down-weighted and
flagged.
"""

from __future__ import annotations

from dataclasses import dataclass, field


@dataclass(frozen=True)
class HpoTerm:
    """A Human Phenotype Ontology term."""

    hpo_id: str   # e.g. "HP:0001250"
    name: str     # e.g. "Seizure"

    def __post_init__(self) -> None:
        if not self.hpo_id.startswith("HP:"):
            raise ValueError(f"Not an HPO id: {self.hpo_id!r}")


@dataclass
class CellTypeRelevance:
    """Relevance of one cell type to a phenotype set."""

    cell_ontology_id: str   # e.g. "CL:0000540" (neuron)
    cell_type_name: str
    score: float            # 0..1
    evidence: dict = field(default_factory=dict)
    # evidence maps source -> contribution, e.g.
    # {"phenotype_tissue": 0.6, "gene_expression": 0.3, "tier": "curated"}


def normalize_score(raw: float) -> float:
    """Clamp a raw relevance score into [0, 1]."""
    return max(0.0, min(1.0, raw))


def combine_evidence(phenotype_tissue: float,
                     gene_expression: float,
                     w1: float = 0.6,
                     w2: float = 0.4) -> float:
    """Weighted combination of the two measured evidence channels.

    Literature co-occurrence (w3) is 0 until wired; the weights are
    renormalized over the active channels so a missing channel cannot
    silently dilute the score.
    """
    total_w = w1 + w2
    if total_w <= 0:
        raise ValueError("At least one evidence weight must be positive")
    return normalize_score((w1 * phenotype_tissue + w2 * gene_expression) / total_w)


def rank_cell_types(candidates: list[CellTypeRelevance],
                    top_k: int | None = None) -> list[CellTypeRelevance]:
    """Rank cell types by relevance score, descending."""
    ranked = sorted(candidates, key=lambda c: c.score, reverse=True)
    return ranked[:top_k] if top_k is not None else ranked
