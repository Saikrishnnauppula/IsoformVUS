#!/usr/bin/env python3
"""Statistical strengthening for IsoformVUS (post-hoc, test-set-only).

1. Bootstrap 95% CIs for AUPRC *differences* (paired, test set):
   M4-M1, M5-M1, M2-M1, M3-M1, M5-M4. The headline claim is the
   M1->M4 AUPRC gain (0.802 -> 0.852); AUROC CIs already exist in
   scripts/12b_paired_analysis.py.
2. Permutation p-values for the AUPRC gains (paired label-score shuffle).
3. M2/M3 clarification: paired AUROC *and* AUPRC CIs in one table.
4. CADD complementarity with a PROPER held-out protocol:
   - Threshold T = median CADD PHRED among pathogenic VALIDATION variants
     (defined without touching the test set).
   - Test subgroup = test variants with CADD PHRED < T ("CADD-weak").
   - Report M5 discrimination *within* the subgroup (AUROC/AUPRC with
     bootstrap CIs), mean M5 score for pathogenic variants in/out of the
     subgroup, and top-decile enrichment.

All model fits use the locked benchmark protocol (logistic regression,
C=1.0, class_weight=balanced, seed=20261006, train split only).

Significance rule: 5 paired tier comparisons are tested (M4-M1, M5-M1,
M2-M1, M3-M1, M5-M4); a comparison is marked significant only if its
permutation p-value clears the Bonferroni-corrected per-comparison
alpha = 0.05/5 = 0.01. E.g. M5-M4 p=0.017 is nominally < 0.05 but
non-significant under this rule (and its CI is entirely negative).

Writes results/statistics.json.
"""

from __future__ import annotations

import json
import sys
from pathlib import Path

import numpy as np
import pandas as pd
from sklearn.linear_model import LogisticRegression
from sklearn.metrics import average_precision_score, roc_auc_score

REPO = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(REPO / "src"))

SEED = 20261006
N_BOOT = 5000
N_PERM = 2000
# Significance rule: 5 paired tier comparisons -> Bonferroni-corrected
# per-comparison alpha. A result is "significant" only if its permutation
# p-value clears this threshold AND the effect points in the tested direction.
N_COMPARISONS = 5
BONF_ALPHA = 0.05 / N_COMPARISONS  # = 0.01


def fit_predict(df, tiers, split_col):
    train = df[split_col] == "train"
    out = {}
    for m, cols in tiers.items():
        clf = LogisticRegression(C=1.0, max_iter=5000, class_weight="balanced",
                                 random_state=SEED)
        clf.fit(df.loc[train, cols].values, df.loc[train, "label"].values)
        out[m] = clf.predict_proba(df[cols].values)[:, 1]
    return out


def paired_diff_ci(y, s1, s2, metric, n_boot=N_BOOT, seed=SEED):
    rng = np.random.default_rng(seed)
    fn = average_precision_score if metric == "auprc" else roc_auc_score
    diffs = np.empty(n_boot)
    n = len(y)
    for b in range(n_boot):
        idx = rng.integers(0, n, n)
        if len(np.unique(y[idx])) < 2:
            diffs[b] = np.nan
            continue
        diffs[b] = fn(y[idx], s1[idx]) - fn(y[idx], s2[idx])
    diffs = diffs[~np.isnan(diffs)]
    return float(np.mean(diffs)), float(np.percentile(diffs, 2.5)), float(np.percentile(diffs, 97.5))


def paired_perm_pvalue(y, s1, s2, metric, n_perm=N_PERM, seed=SEED + 1):
    """Two-sided permutation p-value for paired score difference.

    Under H0 the two models are exchangeable: randomly swap s1/s2 per
    variant and recompute the metric difference.
    """
    rng = np.random.default_rng(seed)
    fn = average_precision_score if metric == "auprc" else roc_auc_score
    obs = fn(y, s1) - fn(y, s2)
    ge = 0
    n = len(y)
    for _ in range(n_perm):
        swap = rng.random(n) < 0.5
        a = np.where(swap, s2, s1)
        b = np.where(swap, s1, s2)
        if abs(fn(y, a) - fn(y, b)) >= abs(obs):
            ge += 1
    return float((ge + 1) / (n_perm + 1)), float(obs)


def boot_ci_metric(y, s, metric, n_boot=N_BOOT, seed=SEED + 2):
    rng = np.random.default_rng(seed)
    fn = average_precision_score if metric == "auprc" else roc_auc_score
    vals = np.empty(n_boot)
    n = len(y)
    for b in range(n_boot):
        idx = rng.integers(0, n, n)
        if len(np.unique(y[idx])) < 2:
            vals[b] = np.nan
            continue
        vals[b] = fn(y[idx], s[idx])
    vals = vals[~np.isnan(vals)]
    return float(np.mean(vals)), float(np.percentile(vals, 2.5)), float(np.percentile(vals, 97.5))


def main() -> int:
    proc = REPO / "data" / "processed"
    df = pd.read_csv(proc / "feature_matrix.tsv.gz", sep="\t")
    tiers = json.loads((proc / "feature_tiers.json").read_text())["tiers"]
    splits = {a["gene"]: a["split"]
              for a in json.loads((REPO / "results/splits.json").read_text())["assignments"]}
    df["split"] = df["gene_symbol"].map(splits)

    probas = fit_predict(df, tiers, "split")
    for m in probas:
        df[f"p_{m}"] = probas[m]

    test = df["split"] == "test"
    val = df["split"] == "validation"
    yt = df.loc[test, "label"].astype(int).values

    stats = {"seed": SEED, "n_boot": N_BOOT, "n_perm": N_PERM,
             "test_n": int(test.sum())}

    # --- 1 & 2: AUPRC difference CIs + permutation p-values ---
    print("=== Paired AUPRC differences (test, bootstrap 95% CI, permutation p) ===")
    stats["auprc_diffs"] = {}
    for a, b in [("M4", "M1"), ("M5", "M1"), ("M2", "M1"), ("M3", "M1"), ("M5", "M4")]:
        mean, lo, hi = paired_diff_ci(yt, df.loc[test, f"p_{a}"].values,
                                      df.loc[test, f"p_{b}"].values, "auprc")
        pval, obs = paired_perm_pvalue(yt, df.loc[test, f"p_{a}"].values,
                                       df.loc[test, f"p_{b}"].values, "auprc")
        sig = pval < BONF_ALPHA and mean > 0
        print(f"  {a}-{b}: {mean:+.4f} [{lo:+.4f},{hi:+.4f}] perm p={pval:.4f} "
              f"{'SIGNIFICANT' if sig else 'n.s.'} (Bonferroni α={BONF_ALPHA:.3f})")
        stats["auprc_diffs"][f"{a}_minus_{b}"] = {
            "mean": mean, "ci_lo": lo, "ci_hi": hi,
            "perm_p": pval, "significant": bool(sig),
            "alpha": BONF_ALPHA}

    # --- 3: M2/M3 clarification table (AUROC + AUPRC) ---
    print("\n=== M2/M3 clarification (paired diffs vs M1) ===")
    stats["m2_m3_clarification"] = {}
    for m in ["M2", "M3"]:
        row = {}
        for metric in ["auroc", "auprc"]:
            mean, lo, hi = paired_diff_ci(yt, df.loc[test, f"p_{m}"].values,
                                          df.loc[test, "p_M1"].values, metric)
            pval, _ = paired_perm_pvalue(yt, df.loc[test, f"p_{m}"].values,
                                         df.loc[test, "p_M1"].values, metric)
            sig = pval < BONF_ALPHA and (mean > 0 if metric == "auprc" else lo > 0 or hi < 0)
            row[metric] = {"mean": mean, "ci_lo": lo, "ci_hi": hi,
                           "perm_p": pval, "significant": bool(sig),
                           "alpha": BONF_ALPHA}
            print(f"  {m}-M1 {metric}: {mean:+.4f} [{lo:+.4f},{hi:+.4f}] "
                  f"p={pval:.4f} {'SIGNIFICANT' if row[metric]['significant'] else 'n.s.'}")
        stats["m2_m3_clarification"][m] = row

    # --- 4: CADD complementarity, held-out protocol ---
    print("\n=== CADD complementarity (threshold from VALIDATION, applied to TEST) ===")
    cadd = pd.read_csv(proc / "variant_cadd.tsv.gz", sep="\t")
    cadd = cadd[cadd["cadd_phred"].notna()].rename(columns={"cadd_phred": "cadd"})
    dfc = df.merge(cadd[["variant_key", "cadd"]], on="variant_key", how="left")

    # Threshold from validation pathogenic variants only (no test peeking).
    valv = dfc[val & dfc["cadd"].notna()]
    T = float(valv.loc[valv["label"] == 1, "cadd"].median())
    print(f"  T (median CADD among pathogenic validation variants): {T:.2f}")

    testc = dfc[test & dfc["cadd"].notna()].copy()
    weak = testc["cadd"] < T
    y_w = testc.loc[weak, "label"].astype(int).values
    s_w = testc.loc[weak, "p_M5"].values
    y_all = testc["label"].astype(int).values

    auroc_w, alo, ahi = boot_ci_metric(y_w, s_w, "auroc")
    auprc_w, plo, phi = boot_ci_metric(y_w, s_w, "auprc")
    print(f"  CADD-weak test subgroup: n={weak.sum()} "
          f"({(testc.loc[weak,'label']==1).sum()} pathogenic)")
    print(f"    M5 AUROC within subgroup: {auroc_w:.4f} [{alo:.4f},{ahi:.4f}]")
    print(f"    M5 AUPRC within subgroup: {auprc_w:.4f} [{plo:.4f},{phi:.4f}]")

    # Mean M5 score: pathogenic in weak subgroup vs pathogenic overall (test, CADD-scored).
    m_path_weak = float(testc.loc[weak & (testc["label"] == 1), "p_M5"].mean())
    m_path_all = float(testc.loc[testc["label"] == 1, "p_M5"].mean())
    m_ben_weak = float(testc.loc[weak & (testc["label"] == 0), "p_M5"].mean())
    print(f"    mean M5 | pathogenic, CADD-weak: {m_path_weak:.3f}")
    print(f"    mean M5 | pathogenic, all CADD-scored test: {m_path_all:.3f}")
    print(f"    mean M5 | benign, CADD-weak: {m_ben_weak:.3f}")

    # Top-decile enrichment within the weak subgroup.
    k = max(1, int(0.1 * weak.sum()))
    order = np.argsort(-s_w)[:k]
    prec = float(y_w[order].mean())
    base_rate = float(y_w.mean())
    print(f"    top-decile M5 precision in weak subgroup: {prec:.3f} "
          f"(base rate {base_rate:.3f}, enrichment {prec/base_rate:.2f}x)")

    stats["cadd_complementarity"] = {
        "threshold_source": "validation pathogenic median",
        "threshold_T": T,
        "test_cadd_scored_n": int(len(testc)),
        "weak_n": int(weak.sum()),
        "weak_pathogenic_n": int((testc.loc[weak, "label"] == 1).sum()),
        "m5_auroc_weak": auroc_w, "m5_auroc_weak_ci": [alo, ahi],
        "m5_auprc_weak": auprc_w, "m5_auprc_weak_ci": [plo, phi],
        "mean_m5_pathogenic_weak": m_path_weak,
        "mean_m5_pathogenic_all_cadd_test": m_path_all,
        "mean_m5_benign_weak": m_ben_weak,
        "top_decile_precision": prec,
        "top_decile_base_rate": base_rate,
        "top_decile_enrichment_ratio": round(prec / base_rate, 2),
    }
    stats["significance_rule"] = {
        "method": "Bonferroni correction across the 5 paired tier comparisons",
        "n_comparisons": N_COMPARISONS,
        "per_comparison_alpha": BONF_ALPHA,
        "note": ("A tier comparison is marked significant only if its permutation "
                 "p-value clears the Bonferroni-corrected per-comparison alpha "
                 f"({BONF_ALPHA}) AND the effect points in the tested (positive) "
                 "direction. E.g. M5-M4 p=0.017 is nominally < 0.05 but does not "
                 "clear alpha=0.01, hence non-significant."),
    }

    out = REPO / "results" / "statistics.json"
    out.write_text(json.dumps(stats, indent=2))
    print(f"\n[wrote {out}]")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
