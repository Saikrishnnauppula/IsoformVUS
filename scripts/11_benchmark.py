#!/usr/bin/env python3
"""Stage 11: held-out benchmark (FIRST and ONLY test-set access).

Scores M1..M5 + baselines on the locked test set; computes the full metric
panel (src/isoformvus/evaluation/metrics.py) with bootstrap CIs; computes
the headline 'recovered pathogenic variants missed by baselines' per k.
Writes results/benchmark.md (human-readable) + results/benchmark.json.

Protocol: this script refuses to run if results/splits.json is missing or if
any model artifact has a modification time newer than the split file without
an accompanying protocol-deviation log entry.
"""

from __future__ import annotations

import argparse
from pathlib import Path

REPO = Path(__file__).resolve().parents[1]


def main() -> int:
    ap = argparse.ArgumentParser()
    ap.add_argument("--config", default=str(REPO / "config" / "base.yaml"))
    args = ap.parse_args()
    print("[11] held-out benchmark — NOT YET IMPLEMENTED (awaits trained models)")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
