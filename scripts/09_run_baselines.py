#!/usr/bin/env python3
"""Stage 09: baselines on the held-out test set.

Baselines (test-set only scoring; no training on test):
    - CADD PHRED (v1.7, GRCh38)
    - SpliceAI max delta (v1.3.1 run)
    - AlphaMissense class/score (missense variants only; other classes flagged
      incompatible, not silently scored)
    - consequence + AF rule baseline (e.g. rare + protein-altering)
    - Exomiser-style phenotype rank (where runnable)
    - IMPPROVE-class cell-type model == our M3 (gene-level expression features)

Writes results/baselines.json with per-variant baseline scores + the
incompatible-class flags.
"""

from __future__ import annotations

import argparse
from pathlib import Path

REPO = Path(__file__).resolve().parents[1]


def main() -> int:
    ap = argparse.ArgumentParser()
    ap.add_argument("--config", default=str(REPO / "config" / "base.yaml"))
    args = ap.parse_args()
    print("[09] baselines — NOT YET IMPLEMENTED (awaits feature matrix + splits)")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
