#!/usr/bin/env python3
"""Stage 04: phenotype -> cell-type relevance mapping.

For each gene's HPO phenotype set (from phenotype_to_genes.txt + hp.obo),
compute cell-relevance scores across Cell Ontology cell types using:
    - phenotype->tissue evidence (HPO annotations + GTEx)
    - gene expression in cell types (single-cell atlases)
Writes data/processed/cell_relevance.tsv.gz
(see src/isoformvus/data/phenotype_cell.py for the scoring model).
"""

from __future__ import annotations

import argparse
from pathlib import Path

REPO = Path(__file__).resolve().parents[1]


def main() -> int:
    ap = argparse.ArgumentParser()
    ap.add_argument("--config", default=str(REPO / "config" / "base.yaml"))
    args = ap.parse_args()
    print("[04] phenotype->cell mapping — NOT YET IMPLEMENTED "
          "(awaits HPO + single-cell inputs)")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
