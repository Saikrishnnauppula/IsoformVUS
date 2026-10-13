"""Evaluation metrics for variant prioritization.

Reports AUROC, AUPRC, precision, recall, F1, calibration, and top-k recall
(k = 1, 5, 10, 20). Top-k recall is the headline clinical-relevance metric:
of all pathogenic variants, how many does the model surface in its top-k
ranked candidates — and crucially, how many of those did conventional
baselines miss?

Conventions
-----------
* ``y_true``: 1 = pathogenic/likely-pathogenic, 0 = benign/likely-benign or
  population control (control type recorded separately; never mixed silently).
* ``scores``: higher = more likely pathogenic.
* Top-k recall is computed per *test set* ranking (global), not per gene.
* Confidence intervals: bootstrap over the test set (seeded); DeLong for
  AUROC comparison where implemented.
"""

from __future__ import annotations

import numpy as np
from sklearn.metrics import (auc, average_precision_score, brier_score_loss,
                             f1_score, precision_recall_curve, precision_score,
                             recall_score, roc_auc_score, roc_curve)


def top_k_recall(y_true: np.ndarray, scores: np.ndarray, k: int) -> float:
    """Fraction of positives ranked within the top-k scores."""
    y_true = np.asarray(y_true)
    scores = np.asarray(scores)
    n_pos = int(y_true.sum())
    if n_pos == 0:
        return float("nan")
    order = np.argsort(-scores, kind="stable")
    topk = order[:k]
    return float(y_true[topk].sum() / n_pos)


def recovered_not_in_baseline(y_true: np.ndarray,
                              scores_new: np.ndarray,
                              scores_base: np.ndarray,
                              k: int) -> dict:
    """Headline metric: pathogenic variants in new model's top-k that the
    baseline's top-k missed.

    Returns counts and the variant indices, so every claimed "recovery" is
    auditable back to a real variant.
    """
    y_true = np.asarray(y_true)
    order_new = np.argsort(-np.asarray(scores_new), kind="stable")[:k]
    order_base = np.asarray(np.argsort(-np.asarray(scores_base), kind="stable")[:k])
    new_hits = {i for i in order_new if y_true[i] == 1}
    base_hits = {i for i in order_base if y_true[i] == 1}
    recovered = sorted(new_hits - base_hits)
    return {
        "k": k,
        "n_new_topk_pathogenic": len(new_hits),
        "n_base_topk_pathogenic": len(base_hits),
        "n_recovered_missed_by_baseline": len(recovered),
        "recovered_indices": recovered,
    }


def full_report(y_true: np.ndarray, scores: np.ndarray,
                ks: tuple[int, ...] = (1, 5, 10, 20)) -> dict:
    """Compute the full metric panel."""
    y_true = np.asarray(y_true, dtype=int)
    scores = np.asarray(scores, dtype=float)
    preds = (scores >= 0.5).astype(int)
    out = {
        "n": len(y_true),
        "n_positive": int(y_true.sum()),
        "auroc": float(roc_auc_score(y_true, scores)),
        "auprc": float(average_precision_score(y_true, scores)),
        "precision": float(precision_score(y_true, preds, zero_division=0)),
        "recall": float(recall_score(y_true, preds, zero_division=0)),
        "f1": float(f1_score(y_true, preds, zero_division=0)),
        "brier": float(brier_score_loss(y_true, scores)),
    }
    for k in ks:
        out[f"top{k}_recall"] = top_k_recall(y_true, scores, k)
    return out


def bootstrap_ci(y_true: np.ndarray, scores: np.ndarray,
                 metric_fn=top_k_recall, n_bootstrap: int = 1000,
                 alpha: float = 0.05, seed: int = 20261006, **kwargs) -> tuple[float, float, float]:
    """Seeded bootstrap CI for a metric. Returns (point, lo, hi)."""
    rng = np.random.default_rng(seed)
    y_true = np.asarray(y_true)
    scores = np.asarray(scores)
    n = len(y_true)
    point = float(metric_fn(y_true, scores, **kwargs))
    boots = []
    for _ in range(n_bootstrap):
        idx = rng.integers(0, n, n)
        # Skip resamples with no positives (metric undefined there).
        if y_true[idx].sum() == 0:
            continue
        boots.append(metric_fn(y_true[idx], scores[idx], **kwargs))
    lo, hi = np.percentile(boots, [100 * alpha / 2, 100 * (1 - alpha / 2)])
    return point, float(lo), float(hi)
