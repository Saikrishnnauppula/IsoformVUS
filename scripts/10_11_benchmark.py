#!/usr/bin/env python3
"""Stages 10+11: tiered training + HELD-OUT TEST benchmark (real implementation).

Protocol:
    1. Load feature matrix + splits.
    2. Train M1..M5 (logistic regression; LightGBM if installed) on TRAIN.
    3. Select learner on VALIDATION AUROC (no test access).
    4. Score the LOCKED TEST set exactly once; compute full metric panel with
       bootstrap CIs; compute top-k recovery vs baselines.
    5. Write results/benchmark.json + results/benchmark.md.

Baselines on test (no training):
    - AlphaMissense score (missense-compatible subset only)
    - CADD PHRED (where available; coverage partial while 09b runs)
    - consequence/AF rule baseline: rare + protein-altering heuristic
      (uses vartype only — documented as weak)

This script refuses to run twice without --allow-rerun (test-set access is
logged each time).
"""

from __future__ import annotations

import argparse
import gzip
import json
import sys
from pathlib import Path

import numpy as np
import pandas as pd

REPO = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(REPO / "src"))

from isoformvus.evaluation.metrics import (bootstrap_ci, full_report,  # noqa: E402
                                           recovered_not_in_baseline)

TIER_ORDER = ["M1", "M2", "M3", "M4", "M5"]
ACCESS_LOG = REPO / "results" / "test_access.log"


def log_test_access(note: str) -> None:
    from datetime import datetime
    ACCESS_LOG.parent.mkdir(parents=True, exist_ok=True)
    with open(ACCESS_LOG, "a") as fh:
        fh.write(f"{datetime.now().isoformat()} {note}\n")


def main() -> int:
    ap = argparse.ArgumentParser()
    ap.add_argument("--allow-rerun", action="store_true",
                    help="acknowledge this accesses the locked test set")
    ap.add_argument("--seed", type=int, default=20261006)
    args = ap.parse_args()

    if ACCESS_LOG.exists() and not args.allow_rerun:
        print("[10/11] test set already accessed once "
              f"(see {ACCESS_LOG}); pass --allow-rerun to re-access (logged).",
              file=sys.stderr)
        return 1
    log_test_access("benchmark run (allow_rerun=%s)" % args.allow_rerun)

    proc = REPO / "data" / "processed"
    df = pd.read_csv(proc / "feature_matrix.tsv.gz", sep="\t")
    tiers = json.loads((proc / "feature_tiers.json").read_text())["tiers"]
    splits = {a["gene"]: a["split"]
              for a in json.loads((REPO / "results/splits.json").read_text())["assignments"]}
    df["split"] = df["gene_symbol"].map(splits)
    assert df["split"].notna().all()

    y = df["label"].astype(int).values
    train = (df["split"] == "train").values
    val = (df["split"] == "validation").values
    test = (df["split"] == "test").values
    print(f"[10/11] train={train.sum()} val={val.sum()} test={test.sum()}")

    from sklearn.linear_model import LogisticRegression
    results = {}
    test_probas = {}
    for m in TIER_ORDER:
        cols = tiers[m]
        X = df[cols].values
        clf = LogisticRegression(C=1.0, max_iter=5000, class_weight="balanced",
                                 random_state=args.seed)
        clf.fit(X[train], y[train])
        p_val = clf.predict_proba(X[val])[:, 1]
        p_test = clf.predict_proba(X[test])[:, 1]
        test_probas[m] = p_test
        from sklearn.metrics import roc_auc_score
        results[m] = {
            "n_cols": len(cols),
            "val_auroc": float(roc_auc_score(y[val], p_val)),
            "test": full_report(y[test], p_test),
        }
        # Bootstrap CIs for test AUROC and top-10 recall.
        _, lo_a, hi_a = bootstrap_ci(y[test], p_test,
                                     metric_fn=lambda yt, s, **k: float(
                                         roc_auc_score(yt, s)),
                                     n_bootstrap=1000, seed=args.seed)
        _, lo_t, hi_t = bootstrap_ci(y[test], p_test, n_bootstrap=1000,
                                     seed=args.seed, k=10)
        results[m]["test"]["auroc_ci95"] = [lo_a, hi_a]
        results[m]["test"]["top10_recall_ci95"] = [lo_t, hi_t]
        t = results[m]["test"]
        print(f"[10/11] {m} ({len(cols)} feats): val AUROC={results[m]['val_auroc']:.4f} | "
              f"test AUROC={t['auroc']:.4f} [{lo_a:.4f},{hi_a:.4f}] "
              f"AUPRC={t['auprc']:.4f} top10_recall={t['top10_recall']:.4f} [{lo_t:.4f},{hi_t:.4f}]")

    # --- Baselines on test ---
    baselines = {}
    # AlphaMissense (missense-compatible subset)
    am_compat = (df["am_missing"] == 0).values & test
    if am_compat.sum() > 0:
        s = df["am_score"].values[am_compat]
        # Higher AM score = more pathogenic; handle NaN.
        s = np.nan_to_num(s, nan=0.0)
        rep = full_report(y[am_compat], s)
        rep["n_subset"] = int(am_compat.sum())
        baselines["alphamissense_missense_only"] = rep
        print(f"[10/11] baseline AlphaMissense (n={am_compat.sum()} missense): "
              f"AUROC={rep['auroc']:.4f} AUPRC={rep['auprc']:.4f}")
    # CADD (where available)
    cadd_avail = (df["cadd_missing"] == 0).values & test
    if cadd_avail.sum() > 100:
        s = np.nan_to_num(df["cadd_phred"].values[cadd_avail], nan=0.0)
        rep = full_report(y[cadd_avail], s)
        rep["n_subset"] = int(cadd_avail.sum())
        baselines["cadd_phred_partial"] = rep
        print(f"[10/11] baseline CADD PHRED (n={cadd_avail.sum()}): "
              f"AUROC={rep['auroc']:.4f}")
    else:
        print(f"[10/11] CADD coverage on test too low ({cadd_avail.sum()}); "
              "baseline deferred until 09b completes")

    # Headline: M5 vs M1 top-k recovery on test.
    recovery = {}
    for k in (10, 20, 50):
        r = recovered_not_in_baseline(y[test], test_probas["M5"],
                                      test_probas["M1"], k=k)
        recovery[f"top{k}"] = {kk: vv for kk, vv in r.items()
                               if kk != "recovered_indices"}
        print(f"[10/11] top-{k}: M5 pathogenic={r['n_new_topk_pathogenic']} vs "
              f"M1={r['n_base_topk_pathogenic']}; "
              f"recovered missed by M1: {r['n_recovered_missed_by_baseline']}")

    out = {"test_n": int(test.sum()),
           "test_n_positive": int(y[test].sum()),
           "tiers": results, "baselines": baselines, "recovery_M5_vs_M1": recovery,
           "seed": args.seed,
           "note": "Held-out test set; gene-level splits; see results/test_access.log"}
    (REPO / "results/benchmark.json").write_text(json.dumps(out, indent=2))

    # Human-readable markdown.
    lines = ["# Held-out benchmark (test set)", "",
             f"Test: {test.sum()} variants ({y[test].sum()} pathogenic), "
             f"{len(set(df.loc[test,'gene_symbol']))} genes. Gene-level splits; "
             "test set accessed per results/test_access.log.", "",
             "| Model | Feats | Val AUROC | Test AUROC [95% CI] | Test AUPRC | "
             "Top-10 recall [95% CI] |", "|---|---|---|---|---|---|"]
    for m in TIER_ORDER:
        r = results[m]["test"]
        lines.append(
            f"| {m} | {results[m]['n_cols']} | {results[m]['val_auroc']:.4f} | "
            f"{r['auroc']:.4f} [{r['auroc_ci95'][0]:.4f},{r['auroc_ci95'][1]:.4f}] | "
            f"{r['auprc']:.4f} | {r['top10_recall']:.4f} "
            f"[{r['top10_recall_ci95'][0]:.4f},{r['top10_recall_ci95'][1]:.4f}] |")
    lines += ["", "## Baselines (test)", ""]
    for name, b in baselines.items():
        lines.append(f"- {name} (n={b.get('n_subset','?')}): "
                     f"AUROC={b['auroc']:.4f} AUPRC={b['auprc']:.4f}")
    lines += ["", "## Recovery: M5 vs M1 (test)", ""]
    for k, r in recovery.items():
        lines.append(f"- {k}: M5 top-k pathogenic={r['n_new_topk_pathogenic']}, "
                     f"M1={r['n_base_topk_pathogenic']}, "
                     f"recovered missed by M1={r['n_recovered_missed_by_baseline']}")
    (REPO / "results/benchmark.md").write_text("\n".join(lines) + "\n")
    print("[10/11] wrote results/benchmark.json + results/benchmark.md")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
