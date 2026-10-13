#!/usr/bin/env python3
"""Stage 10: train tiered IsoformVUS models (M1..M5).

Uses src/isoformvus/models/tiered.py::fit_tiered_resolved with the
FEATURE_TIERS registry from stage 07. Trains on train split, selects the
learner (logreg vs lightgbm) on validation AUROC, and writes
results/tiered_models.pkl + results/validation_metrics.json.

The TEST split is not loaded in this stage (enforced by only reading the
train/validation row masks).
"""

from __future__ import annotations

import argparse
from pathlib import Path

REPO = Path(__file__).resolve().parents[1]


def main() -> int:
    ap = argparse.ArgumentParser()
    ap.add_argument("--config", default=str(REPO / "config" / "base.yaml"))
    args = ap.parse_args()
    print("[10] tiered training — NOT YET IMPLEMENTED (awaits feature matrix)")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
