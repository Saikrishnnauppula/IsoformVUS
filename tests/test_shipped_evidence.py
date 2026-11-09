"""Data-grounded tests over the shipped evidence files.

The other test modules exercise the package on toy fixtures. These tests
read the shipped artifacts (feature matrix, splits, benchmark, statistics,
test-access log) and check the claims the manuscript and docs make about
them, so the suite proves something about the science, not only the code.

Choices, stated:
- CADD counts and split sizes are recomputed from the shipped files on
  every run (gene->split mapping is read from results/splits.json, never
  hardcoded); the expected literals are the values verified from those
  files on 2026-11-09.
- The benchmark check is a subsample refit, not a full refit: a full refit
  under the locked protocol (LogisticRegression, C=1.0,
  class_weight="balanced", max_iter=5000, 35,921 train rows) is too slow
  for pytest. The subsample refit (5,000 train rows, max_iter=1000) checks
  that benchmark.json is reproducible in kind: each refit AUPRC lands
  within 0.02 of the shipped value, and the tier ordering matches the
  shipped benchmark (on the 2026-11-09 no-CADD refit, M4 above M1).
- The statistics checks target the 2026-11-09 refit structure written by
  scripts/19_refit_statistics.py: live diffs in `auprc_diffs` and
  `tie_corrected_refit_no_cadd`; the dead pre-rebuild +0.0509 gain may
  appear only under `auprc_diffs_superseded_pre_rebuild`, labelled
  superseded.
- The test-access log is the original record: 16 timestamped entries
  (1 on 2026-10-06, 8 on 2026-10-07, 4 on 2026-11-08, 3 on 2026-11-09),
  matching the counts in paper/manuscript.md Section 2.4. The 2026-11-09
  entries are a script-18 verification rerun, the pinned-stack benchmark
  refit rerun, and a baseline-enrichment run. It was untracked until
  2026-11-09 because .gitignore carried a bare *.log rule; the negation
  !results/test_access.log now keeps it in the repository.
"""

import json
import re
import warnings
from pathlib import Path

import pandas as pd
import pytest

REPO = Path(__file__).resolve().parents[1]
PROCESSED = REPO / "data" / "processed"
RESULTS = REPO / "results"


@pytest.fixture(scope="module")
def matrix():
    return pd.read_csv(PROCESSED / "feature_matrix.tsv.gz", sep="\t")


@pytest.fixture(scope="module")
def splits():
    return json.loads((RESULTS / "splits.json").read_text())


@pytest.fixture(scope="module")
def gene_split(splits):
    return {a["gene"]: a["split"] for a in splits["assignments"]}


@pytest.fixture(scope="module")
def matrix_by_split(matrix, gene_split):
    df = matrix.copy()
    df["split"] = df["gene_symbol"].map(gene_split)
    return df


def test_feature_matrix_shape_and_no_nans(matrix):
    assert matrix.shape == (60521, 29)
    assert int(matrix.isna().sum().sum()) == 0


def test_splits_have_zero_gene_overlap(splits, matrix_by_split):
    by_split = {}
    for a in splits["assignments"]:
        by_split.setdefault(a["split"], set()).add(a["gene"])
    assert {k: len(v) for k, v in by_split.items()} == {
        "train": 61,
        "validation": 17,
        "test": 24,
    }
    names = sorted(by_split)
    for i, a in enumerate(names):
        for b in names[i + 1 :]:
            assert by_split[a].isdisjoint(by_split[b]), f"{a} overlaps {b}"
    # The matrix's own gene sets agree with splits.json.
    for split, genes in by_split.items():
        assert (
            set(matrix_by_split.loc[matrix_by_split["split"] == split, "gene_symbol"])
            == genes
        )
    counts = matrix_by_split.groupby("split").size().to_dict()
    assert counts == {"train": 35921, "validation": 10200, "test": 14400}


def test_cadd_scored_counts_by_split(matrix_by_split):
    scored = matrix_by_split[matrix_by_split["cadd_missing"] == 0]
    per_split = scored.groupby("split").size().to_dict()
    assert per_split.get("train", 0) == 0
    assert per_split["validation"] == 7774
    assert per_split["test"] == 11240
    assert len(scored) == 19014


def test_train_rows_all_carry_cadd_sentinel(matrix_by_split):
    train = matrix_by_split[matrix_by_split["split"] == "train"]
    assert len(train) == 35921
    assert (train["cadd_phred"] == -1.0).all()
    assert (train["cadd_missing"] == 1).all()


def test_fitted_tiers_do_not_use_cadd():
    tiers = json.loads((PROCESSED / "feature_tiers.json").read_text())["tiers"]
    for tier, cols in tiers.items():
        assert "cadd_phred" not in cols, f"{tier} still fits cadd_phred"
        assert "cadd_missing" not in cols, f"{tier} still fits cadd_missing"


def test_benchmark_json_consistent_with_statistics():
    bench = json.loads((RESULTS / "benchmark.json").read_text())
    stats = json.loads((RESULTS / "statistics.json").read_text())
    assert bench["test_n"] == 14400
    assert bench["test_n_positive"] == 7200
    step = stats["tie_corrected_refit_no_cadd"]["step_auprc"]
    auprcs = {}
    for tier in ["M1", "M2", "M3", "M4", "M5"]:
        auprcs[tier] = bench["tiers"][tier]["test"]["auprc"]
        assert auprcs[tier] == pytest.approx(step[tier], abs=0.001)
    # On the 2026-11-09 no-CADD refit, M4 has the highest step AUPRC.
    assert max(auprcs, key=auprcs.get) == "M4"
    assert auprcs["M4"] > auprcs["M1"]


def test_benchmark_tiers_match_subsample_refit(matrix_by_split):
    from sklearn.linear_model import LogisticRegression
    from sklearn.metrics import average_precision_score

    tiers = json.loads((PROCESSED / "feature_tiers.json").read_text())["tiers"]
    bench = json.loads((RESULTS / "benchmark.json").read_text())
    train = matrix_by_split[matrix_by_split["split"] == "train"].sample(
        n=5000, random_state=20261006
    )
    test = matrix_by_split[matrix_by_split["split"] == "test"]
    refit = {}
    for tier in ["M1", "M4"]:
        cols = tiers[tier]
        clf = LogisticRegression(C=1.0, max_iter=1000, class_weight="balanced")
        with warnings.catch_warnings():
            warnings.simplefilter("ignore")
            clf.fit(train[cols], train["label"])
        proba = clf.predict_proba(test[cols])[:, 1]
        refit[tier] = average_precision_score(test["label"], proba)
        shipped = bench["tiers"][tier]["test"]["auprc"]
        assert refit[tier] == pytest.approx(shipped, abs=0.02), (
            f"{tier}: subsample refit AUPRC {refit[tier]:.4f} is more than "
            f"0.02 from shipped benchmark.json value {shipped:.4f}"
        )
    # The refit ordering matches the shipped benchmark: M4 above M1.
    assert refit["M4"] > refit["M1"]


def test_statistics_has_no_live_significant_0509():
    stats = json.loads((RESULTS / "statistics.json").read_text())
    live = stats["auprc_diffs"]["M4_minus_M1"]
    refit = stats["tie_corrected_refit_no_cadd"]["diffs"]["M4_minus_M1"]
    # The live top-level diff is the 2026-11-09 refit value, not the dead one.
    assert live["mean"] == pytest.approx(refit["mean"], abs=0.0005)
    # The dead pre-rebuild +0.0509 survives only under the superseded key,
    # explicitly labelled as such.
    sup = stats["auprc_diffs_superseded_pre_rebuild"]
    assert "SUPERSEDED" in sup["note"]
    assert "DEAD" in sup["note"]
    old = sup["values"]["M4_minus_M1"]
    assert old["mean"] == pytest.approx(0.0509, abs=0.0005)


def test_statistics_gene_clustered_m2_verdict_present():
    stats = json.loads((RESULTS / "statistics.json").read_text())
    clustered = stats["gene_clustered_cis"]
    assert clustered["n_clusters"] == 24
    m2 = clustered["diffs"]["M2_minus_M1"]
    assert len(m2["cluster_ci95"]) == 2
    assert "M2" in clustered["m2_minus_m1_verdict"]


def test_statistics_significance_flags_match_rule():
    # Defect fixed 2026-11-09: the permutation p-value and the CI are both
    # tie-aware AUPRC, and significant requires perm_p < alpha AND the CI
    # to exclude zero. Before the fix M2-M1 and M3-M1 were flagged
    # significant:true on a step-AUPRC permutation p while their
    # tie-aware CIs included zero.
    stats = json.loads((RESULTS / "statistics.json").read_text())
    diffs = stats["auprc_diffs"]
    for name, row in diffs.items():
        expected = row["perm_p"] < row["alpha"] and (
            row["ci_lo"] > 0 or row["ci_hi"] < 0
        )
        assert row["significant"] == expected, name
    assert diffs["M2_minus_M1"]["significant"] is False
    assert diffs["M3_minus_M1"]["significant"] is False
    assert diffs["M4_minus_M1"]["significant"] is True
    assert diffs["M5_minus_M1"]["significant"] is True
    assert diffs["M5_minus_M4"]["significant"] is False


def test_access_log_matches_manuscript_record():
    # The committed record is 16 entries. The log is append-only by
    # design: any local rerun of a test-scoring script (the benchmark,
    # 12b, 14, 18, or a full verify-mode run) appends real entries, so
    # the committed prefix is what this test pins; extra lines must be
    # later, timestamped entries, never edits to the record.
    lines = [
        line
        for line in (RESULTS / "test_access.log").read_text().splitlines()
        if line.strip()
    ]
    assert len(lines) >= 16
    committed = lines[:16]
    dates = [re.match(r"\d{4}-\d{2}-\d{2}", line).group(0) for line in committed]
    assert dates.count("2026-10-06") == 1
    assert dates.count("2026-10-07") == 8
    assert dates.count("2026-11-08") == 4
    assert dates.count("2026-11-09") == 3
    joined = "\n".join(committed)
    assert joined.count("benchmark run") == 4
    assert joined.count("baseline enrichment") == 10
    assert joined.count("script-18") == 2
    for extra in lines[16:]:
        assert re.match(r"\d{4}-\d{2}-\d{2}T", extra), extra
