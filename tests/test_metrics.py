"""Tests for evaluation metrics."""

import numpy as np
import pytest

from isoformvus.evaluation.metrics import (
    bootstrap_ci,
    full_report,
    recovered_not_in_baseline,
    top_k_recall,
)


def _toy():
    # 4 positives out of 10; scores rank 3 of 4 positives in top-4.
    y = np.array([1, 1, 0, 1, 0, 0, 1, 0, 0, 0])
    s = np.array([0.9, 0.8, 0.75, 0.7, 0.6, 0.5, 0.4, 0.3, 0.2, 0.1])
    return y, s


def test_top_k_recall():
    y, s = _toy()
    assert top_k_recall(y, s, 4) == pytest.approx(0.75)  # 3 of 4 positives
    assert top_k_recall(y, s, 10) == pytest.approx(1.0)
    assert np.isnan(top_k_recall(np.zeros(5), np.ones(5), 3))


def test_recovered_not_in_baseline():
    y, s_new = _toy()
    # Baseline ranks positives poorly (reversed-ish).
    s_base = np.array([0.1, 0.2, 0.9, 0.3, 0.8, 0.7, 0.4, 0.6, 0.5, 0.35])
    r = recovered_not_in_baseline(y, s_new, s_base, k=4)
    assert r["n_recovered_missed_by_baseline"] >= 1
    assert all(y[i] == 1 for i in r["recovered_indices"])


def test_full_report_keys():
    y, s = _toy()
    rep = full_report(y, s)
    for key in ("auroc", "auprc", "precision", "recall", "f1", "brier",
                "top1_recall", "top5_recall", "top10_recall", "top20_recall"):
        assert key in rep


def test_bootstrap_ci_seeded():
    y, s = _toy()
    p1, lo1, hi1 = bootstrap_ci(y, s, n_bootstrap=200, seed=7, k=5)
    p2, lo2, hi2 = bootstrap_ci(y, s, n_bootstrap=200, seed=7, k=5)
    assert (p1, lo1, hi1) == (p2, lo2, hi2)
    assert lo1 <= p1 <= hi1
