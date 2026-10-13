#!/usr/bin/env python3
"""Prototype ablation: M4 without MANE-derived features (ascertainment check).

MANE Select transcripts are chosen partly for clinical relevance, so
'mane_affected' may carry ascertainment signal. This script re-runs the
M1 vs M4 comparison excluding all MANE-derived features.
"""

from __future__ import annotations

import json
import sys
from pathlib import Path

import pandas as pd
from sklearn.linear_model import LogisticRegression
from sklearn.metrics import average_precision_score, roc_auc_score

REPO = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(REPO / "scripts"))


def main() -> int:
    labels = pd.read_csv(REPO / "data/processed/variant_labels.tsv.gz", sep="\t",
                         usecols=["variant_key", "gene_symbol", "label",
                                  "vartype", "review_stars"])
    effects = pd.read_csv(REPO / "data/processed/variant_isoform_effects.tsv.gz", sep="\t",
                          usecols=["variant_key", "gene_effect_summary",
                                   "n_isoform_specific_hits",
                                   "has_junction_specific_hit",
                                   "has_coding_specific_hit"])
    df = labels.merge(effects, on="variant_key", how="inner")
    splits = {a["gene"]: a["split"]
              for a in json.loads((REPO / "results/splits.json").read_text())["assignments"]}
    df["split"] = df["gene_symbol"].map(splits)

    X = pd.get_dummies(df[["vartype", "gene_effect_summary"]].astype(str),
                       columns=["vartype", "gene_effect_summary"])
    for c in ["review_stars", "n_isoform_specific_hits",
              "has_junction_specific_hit", "has_coding_specific_hit"]:
        X[c] = df[c]
    y = df["label"].astype(int)
    train = df["split"] == "train"
    val = df["split"] == "validation"

    m1 = [c for c in X.columns if c == "review_stars" or c.startswith("vartype_")]
    m4_nomane = m1 + [c for c in X.columns
                      if c.startswith("gene_effect_summary_") or c in
                      ("n_isoform_specific_hits", "has_junction_specific_hit",
                       "has_coding_specific_hit")]
    print(f"train={train.sum()} val={val.sum()}")
    for name, cols in [("M1 variant-only", m1), ("M4-noMANE +isoform", m4_nomane)]:
        clf = LogisticRegression(C=1.0, max_iter=5000, class_weight="balanced",
                                 random_state=20261006)
        clf.fit(X.loc[train, cols], y[train])
        p = clf.predict_proba(X.loc[val, cols])[:, 1]
        print(f"{name}: n_cols={len(cols)} "
              f"val AUROC={roc_auc_score(y[val], p):.4f} "
              f"AUPRC={average_precision_score(y[val], p):.4f}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
