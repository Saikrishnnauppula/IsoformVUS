#!/usr/bin/env python3
"""Prototype experiment: does isoform structure add signal over variant-only features?

M1 (variant-only): vartype one-hot + review_stars
M4 (variant+isoform): M1 + gene_effect_summary one-hot + n_isoform_specific_hits
    + has_junction_specific_hit + has_coding_specific_hit + mane_select_affected

Trains logistic regression on TRAIN genes, evaluates on VALIDATION genes
(test set untouched). This is a prototype — full tiered models with
CADD/SpliceAI/AlphaMissense come in stages 07/10.
"""

from __future__ import annotations

import gzip
import json
import sys
from pathlib import Path

import numpy as np
import pandas as pd
from sklearn.linear_model import LogisticRegression
from sklearn.metrics import average_precision_score, roc_auc_score

REPO = Path(__file__).resolve().parents[1]


def load_labels(path: Path) -> pd.DataFrame:
    return pd.read_csv(path, sep="\t",
                       usecols=["variant_key", "gene_symbol", "label",
                                "vartype", "review_stars"])


def load_effects(path: Path) -> pd.DataFrame:
    return pd.read_csv(path, sep="\t",
                       usecols=["variant_key", "gene_effect_summary",
                                "n_isoform_specific_hits",
                                "has_junction_specific_hit",
                                "has_coding_specific_hit",
                                "mane_select_effect"])


def load_splits(path: Path) -> dict[str, str]:
    payload = json.loads(path.read_text())
    return {a["gene"]: a["split"] for a in payload["assignments"]}


def main() -> int:
    labels = load_labels(REPO / "data/processed/variant_labels.tsv.gz")
    effects = load_effects(REPO / "data/processed/variant_isoform_effects.tsv.gz")
    df = labels.merge(effects, on="variant_key", how="inner")
    splits = load_splits(REPO / "results/splits.json")
    df["split"] = df["gene_symbol"].map(splits)
    assert df["split"].notna().all(), "variants without split assignment!"

    df["mane_affected"] = (df["mane_select_effect"] != "none").astype(int)
    df["mane_affected"] = df["mane_affected"].fillna(0)

    m1_cols = ["review_stars"] + [f"vartype_{t}" for t in
                                  ["SNV", "deletion", "insertion", "complex"]]
    m4_extra = ["n_isoform_specific_hits", "has_junction_specific_hit",
                "has_coding_specific_hit", "mane_affected"] + \
               [f"geff_{s}" for s in ["all", "some", "isoform_specific", "none"]]

    X = pd.get_dummies(df[["vartype", "gene_effect_summary"]].astype(str),
                       columns=["vartype", "gene_effect_summary"])
    for c in m1_cols + m4_extra:
        if c not in X.columns and c in df.columns:
            X[c] = df[c]
    # Align dummy names
    X = X.reindex(columns=[c for c in X.columns], fill_value=0)
    y = df["label"].astype(int)

    def cols_for(prefixes):
        return [c for c in X.columns
                if any(c == p or c.startswith(p + "_") or c.startswith(p)
                       for p in prefixes)]

    m1 = [c for c in X.columns if c == "review_stars" or c.startswith("vartype_")]
    m4 = m1 + [c for c in X.columns if c.startswith("geff_") or c in
               ("n_isoform_specific_hits", "has_junction_specific_hit",
                "has_coding_specific_hit", "mane_affected")]

    train = df["split"] == "train"
    val = df["split"] == "validation"
    print(f"train={train.sum()} val={val.sum()} "
          f"(test={(df['split']=='test').sum()} untouched)")

    for name, cols in [("M1 variant-only", m1), ("M4 +isoform", m4)]:
        clf = LogisticRegression(C=1.0, max_iter=5000, class_weight="balanced",
                                 random_state=20261006)
        clf.fit(X.loc[train, cols], y[train])
        p = clf.predict_proba(X.loc[val, cols])[:, 1]
        print(f"{name}: cols={len(cols)} "
              f"val AUROC={roc_auc_score(y[val], p):.4f} "
              f"AUPRC={average_precision_score(y[val], p):.4f}")

    # Feature coefficients for the M4 model (interpretability)
    clf = LogisticRegression(C=1.0, max_iter=5000, class_weight="balanced",
                             random_state=20261006).fit(X.loc[train, m4], y[train])
    coef = pd.Series(clf.coef_[0], index=m4).sort_values(key=abs, ascending=False)
    print("\nM4 top coefficients (log-odds):")
    print(coef.head(10).to_string())
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
