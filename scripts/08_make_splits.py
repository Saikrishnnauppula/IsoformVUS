#!/usr/bin/env python3
"""Stage 08: leakage-safe splits (real implementation).

Assigns panel genes to train/validation/test via
src/isoformvus/data/splits.py (gene-level, disease-stratified, deterministic).
Writes results/splits.json ONCE; refuses to overwrite an existing file unless
--force is given (overwriting is logged as a protocol deviation).

The held-out test set must not be touched again until stage 11.
"""

from __future__ import annotations

import argparse
import sys
from pathlib import Path

REPO = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(REPO / "src"))

from isoformvus.data.splits import (assign_splits, splits_to_json,  # noqa: E402
                                    verify_no_gene_overlap)


def main() -> int:
    ap = argparse.ArgumentParser()
    ap.add_argument("--config", default=str(REPO / "config" / "base.yaml"))
    ap.add_argument("--train-frac", type=float, default=0.6)
    ap.add_argument("--val-frac", type=float, default=0.2)
    ap.add_argument("--seed", type=int, default=20261006)
    ap.add_argument("--force", action="store_true")
    args = ap.parse_args()

    out = REPO / "results" / "splits.json"
    if out.exists() and not args.force:
        print(f"[08] {out} exists; refusing to overwrite (leakage protocol).")
        return 0

    panel = REPO / "metadata" / "gene_panel.tsv"
    genes, cats = [], {}
    with open(panel) as fh:
        fh.readline()
        for line in fh:
            parts = line.rstrip("\n").split("\t")
            genes.append(parts[0])
            cats[parts[0]] = parts[1]

    assignments = assign_splits(genes, cats, train_frac=args.train_frac,
                                val_frac=args.val_frac, seed=args.seed)
    verify_no_gene_overlap(assignments)

    from collections import Counter
    dist = Counter(a.split for a in assignments)
    print(f"[08] split distribution: {dict(dist)} (seed {args.seed})")
    for split in ("train", "validation", "test"):
        sc = Counter(a.disease_category for a in assignments if a.split == split)
        print(f"[08]   {split}: {dict(sc)}")

    out.write_text(splits_to_json(assignments))
    print(f"[08] wrote {out} — TEST SET LOCKED until stage 11")
    if args.force:
        print("[08] PROTOCOL DEVIATION LOGGED: splits regenerated with --force")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
