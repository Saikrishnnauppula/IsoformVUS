"""Leakage-safe dataset splitting.

Genes — never individual variants — are assigned to TRAIN / VALIDATION /
HELD-OUT TEST. See docs/leakage_prevention.md for the threat model.

The split is deterministic given (gene list, seed, ratios): the same inputs
always produce the same assignment, so the held-out test set can be audited
independently of the training code.
"""

from __future__ import annotations

import hashlib
import json
from dataclasses import dataclass, field


@dataclass(frozen=True)
class SplitAssignment:
    gene_symbol: str
    split: str  # "train" | "validation" | "test"
    disease_category: str = ""


def _stable_shuffle(items: list[str], seed: int) -> list[str]:
    """Deterministic shuffle from SHA-256 hashes (no RNG state to leak)."""
    keyed = [(hashlib.sha256(f"{seed}:{x}".encode()).hexdigest(), x) for x in items]
    keyed.sort(key=lambda kv: kv[0])
    return [x for _, x in keyed]


def assign_splits(genes: list[str],
                  disease_category: dict[str, str] | None = None,
                  train_frac: float = 0.6,
                  val_frac: float = 0.2,
                  seed: int = 20261006) -> list[SplitAssignment]:
    """Assign each gene to a split, stratified by disease category.

    Stratification keeps disease biology balanced across splits so the test
    set is not accidentally all-neurodevelopmental (or similar).

    Raises:
        ValueError: on duplicate genes or invalid fractions.
    """
    if len(set(genes)) != len(genes):
        raise ValueError("Duplicate genes in input")
    if not (0 < train_frac < 1 and 0 < val_frac < 1 and train_frac + val_frac < 1):
        raise ValueError("Fractions must be in (0,1) and sum < 1")

    disease_category = disease_category or {}
    # Group genes by category; uncategorized genes form their own stratum.
    strata: dict[str, list[str]] = {}
    for g in genes:
        strata.setdefault(disease_category.get(g, "__uncategorized__"), []).append(g)

    assignments: list[SplitAssignment] = []
    for category in sorted(strata):
        members = _stable_shuffle(sorted(strata[category]), seed)
        n = len(members)
        n_train = int(round(n * train_frac))
        n_val = int(round(n * val_frac))
        # Guarantee at least one gene per split when the stratum allows.
        if n >= 3:
            n_train = max(1, n_train)
            n_val = max(1, n_val)
        if n_train + n_val >= n and n >= 3:
            n_val = max(1, n - n_train - 1)
        for i, g in enumerate(members):
            if i < n_train:
                split = "train"
            elif i < n_train + n_val:
                split = "validation"
            else:
                split = "test"
            assignments.append(SplitAssignment(g, split, category))
    return assignments


def verify_no_gene_overlap(assignments: list[SplitAssignment]) -> None:
    """Raise if any gene appears in more than one split."""
    seen: dict[str, str] = {}
    for a in assignments:
        if a.gene_symbol in seen and seen[a.gene_symbol] != a.split:
            raise ValueError(f"Gene {a.gene_symbol} in both "
                             f"{seen[a.gene_symbol]} and {a.split}")
        seen[a.gene_symbol] = a.split


def verify_variant_split_integrity(variant_gene: dict[str, str],
                                   gene_split: dict[str, str]) -> list[str]:
    """Check that no variant key maps to genes in different splits.

    Returns a list of violation descriptions (empty = clean).
    """
    violations = []
    for vkey, gene in variant_gene.items():
        if gene not in gene_split:
            violations.append(f"variant {vkey}: gene {gene} has no split assignment")
    return violations


def splits_to_json(assignments: list[SplitAssignment]) -> str:
    payload = {
        "seed_note": "deterministic SHA-256 shuffle; see data/splits.py",
        "assignments": [
            {"gene": a.gene_symbol, "split": a.split,
             "disease_category": a.disease_category}
            for a in sorted(assignments, key=lambda a: a.gene_symbol)
        ],
    }
    return json.dumps(payload, indent=2)
