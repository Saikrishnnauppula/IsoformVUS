#!/usr/bin/env python3
"""Post-benchmark analysis: paired statistical tests + recovered-variant extraction.

1. Paired bootstrap CI for AUROC differences (M4-M1, M5-M4, M3-M1).
2. Extract the top-k pathogenic variants M5 ranks highly that M1 misses,
   with full annotation chains for case-study selection.
3. Isoform-selective subgroup evaluation (pre-registered subgroup).
"""

from __future__ import annotations

import gzip
import json
import sys
from pathlib import Path

import numpy as np
import pandas as pd
from sklearn.metrics import roc_auc_score

REPO = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(REPO / "src"))
from isoformvus.evaluation.metrics import recovered_not_in_baseline  # noqa: E402


def paired_auroc_diff_ci(y, s1, s2, n_boot=2000, seed=20261006):
    rng = np.random.default_rng(seed)
    diffs = []
    n = len(y)
    for _ in range(n_boot):
        idx = rng.integers(0, n, n)
        if len(np.unique(y[idx])) < 2:
            continue
        diffs.append(roc_auc_score(y[idx], s1[idx]) - roc_auc_score(y[idx], s2[idx]))
    lo, hi = np.percentile(diffs, [2.5, 97.5])
    return float(np.mean(diffs)), float(lo), float(hi)


def main() -> int:
    proc = REPO / "data" / "processed"
    df = pd.read_csv(proc / "feature_matrix.tsv.gz", sep="\t")
    tiers = json.loads((proc / "feature_tiers.json").read_text())["tiers"]
    splits = {a["gene"]: a["split"]
              for a in json.loads((REPO / "results/splits.json").read_text())["assignments"]}
    df["split"] = df["gene_symbol"].map(splits)
    test = df["split"] == "test"
    y = df.loc[test, "label"].astype(int).values

    # Re-fit M1/M4/M5 on train, score test (same protocol as benchmark).
    from sklearn.linear_model import LogisticRegression
    train = df["split"] == "train"
    probas = {}
    for m in ["M1", "M3", "M4", "M5"]:
        cols = tiers[m]
        clf = LogisticRegression(C=1.0, max_iter=5000, class_weight="balanced",
                                 random_state=20261006)
        clf.fit(df.loc[train, cols].values, df.loc[train, "label"].values)
        probas[m] = clf.predict_proba(df.loc[test, cols].values)[:, 1]

    print("=== Paired AUROC differences (test, bootstrap 95% CI) ===")
    # Fit M2 as well (needed for the M2-M1 comparison).
    cols2 = tiers["M2"]
    clf2 = LogisticRegression(C=1.0, max_iter=5000, class_weight="balanced",
                              random_state=20261006)
    clf2.fit(df.loc[train, cols2].values, df.loc[train, "label"].values)
    probas["M2"] = clf2.predict_proba(df.loc[test, cols2].values)[:, 1]
    for a, b in [("M4", "M1"), ("M5", "M4"), ("M3", "M1"), ("M2", "M1")]:
        mean, lo, hi = paired_auroc_diff_ci(y, probas[a], probas[b])
        sig = "SIGNIFICANT" if lo > 0 or hi < 0 else "n.s."
        print(f"  {a}-{b}: {mean:+.4f} [{lo:+.4f},{hi:+.4f}] {sig}")

    # --- Isoform-selective subgroup (pre-registered) ---
    eff = pd.read_csv(proc / "variant_isoform_effects.tsv.gz", sep="\t",
                      usecols=["variant_key", "gene_effect_summary"])
    dft = df[test].merge(eff, on="variant_key", how="left")
    sub = dft["gene_effect_summary"] == "isoform_specific"
    print(f"\n=== Isoform-selective subgroup (n={sub.sum()}/{test.sum()}) ===")
    ys = dft.loc[sub, "label"].astype(int).values
    # Map back to probas indices
    idx = np.where(test.values)[0][sub.values]
    for m in ["M1", "M4", "M5"]:
        pm = probas[m][sub.values]
        print(f"  {m}: AUROC={roc_auc_score(ys, pm):.4f}")

    # --- Recovered variants for case studies ---
    test_df = df[test].reset_index(drop=True)
    r = recovered_not_in_baseline(y, probas["M5"], probas["M1"], k=50)
    rec_idx = r["recovered_indices"]
    print(f"\n=== {len(rec_idx)} variants in M5 top-50 missed by M1 top-50 ===")
    # Annotate with labels + effects
    lab = pd.read_csv(proc / "variant_labels.tsv.gz", sep="\t",
                      usecols=["variant_key", "gene_symbol", "vartype",
                               "clnsig", "condition"])
    rec = test_df.iloc[rec_idx][["variant_key", "gene_symbol"]].merge(
        lab, on=["variant_key", "gene_symbol"], how="left")
    rec = rec.merge(eff, on="variant_key", how="left")
    rec["m5_rank"] = range(1, len(rec) + 1)
    rec.to_csv(REPO / "results/recovered_variants.tsv", sep="\t", index=False)
    print(rec[["m5_rank", "gene_symbol", "variant_key", "vartype",
               "gene_effect_summary"]].head(15).to_string(index=False))
    print(f"\n[wrote results/recovered_variants.tsv]")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
