#!/usr/bin/env python3
"""Stage 12: ablation study.

Compares M1..M5 on the held-out test set (metrics from stage 11 artifacts):
does each added biological layer (phenotype, cell type, isoform, isoform
effect) actually improve performance? Reports deltas with bootstrap CIs and
paired tests; writes results/ablation.md + figures/05_ablation.png.

A null result for any layer is reported as a null result.
"""

from __future__ import annotations

import argparse
from pathlib import Path

REPO = Path(__file__).resolve().parents[1]


def main() -> int:
    ap = argparse.ArgumentParser()
    ap.add_argument("--config", default=str(REPO / "config" / "base.yaml"))
    args = ap.parse_args()
    print("[12] ablation — NOT YET IMPLEMENTED (awaits benchmark)")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
