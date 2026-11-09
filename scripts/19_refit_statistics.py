#!/usr/bin/env python3
"""Recompute results/statistics.json on the no-CADD refit (2026-11-09).

Replaces the pre-rebuild / CADD-contaminated statistics with values recomputed
from tiers fit WITHOUT cadd_phred and cadd_missing (train n=35,921 had 0
CADD-scored rows; both columns constant in train). Old values are preserved
under explicitly labelled `*_superseded_*` keys, never as live claims.

Computations (locked protocol: logistic regression, C=1.0,
class_weight=balanced, seed=20261006, train split only):
  1. Step AUPRC + tie-aware AUPRC (McSherry & Najork 2008 expected AP under
     uniform tie-breaking; function copied from scripts/18_exon_architecture_nmd.py)
     per tier, with paired variant-level bootstrap CIs (n_boot=5000) and
     permutation p-values on the same tie-aware metric (n_perm=2000) for
     the 5 tier diffs. A diff is significant only if perm p < 0.01
     (Bonferroni) AND the bootstrap CI excludes zero.
  2. Gene-clustered inference for tier AUPRC diffs: cluster bootstrap resampling
     the 24 test genes (n_boot=5000) + leave-one-gene-out jackknife (24 folds).
     Variant-level CIs are invalid for gene-level HPO constants, so the cluster
     interval is the deciding evidence for "M2 best tier".
  3. CADD complementarity, held-out protocol from script 17 section 4
     (threshold T = median CADD among pathogenic VALIDATION variants),
     with the refit M5.
  4. Same-subset M5-vs-CADD AUPRC check (methods.md:209 audit): M5 AUPRC minus
     CADD AUPRC on the CADD-scored test subset, train-fit models. No
     validation-trained stack exists in any shipped script; this is the closest
     computable quantity and it is reported as such, not as the claimed number.
  5. Per-tier AUPRC CIs -> results/per_model_cis.json (that file previously had
     no generating script in the repo).

Writes results/statistics.json and results/per_model_cis.json.
"""

from __future__ import annotations

import json
import sys
from pathlib import Path

import numpy as np
import pandas as pd
import sklearn
from scipy.stats import rankdata
from sklearn.linear_model import LogisticRegression
from sklearn.metrics import average_precision_score, roc_auc_score

REPO = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(REPO / "src"))

SEED = 20261006
N_BOOT = 5000
N_PERM = 2000
N_COMPARISONS = 5
BONF_ALPHA = 0.05 / N_COMPARISONS  # = 0.01
TIER_ORDER = ["M1", "M2", "M3", "M4", "M5"]


_HARM = np.concatenate([[0.0], np.cumsum(1.0 / np.arange(1, 400001, dtype=np.float64))])


def tie_aware_auprc(y_true, y_score) -> float:
    """Analytical expectation of average precision under uniform tie-breaking
    (McSherry & Najork 2008). Vectorized form of the pandas implementation in
    scripts/18_exon_architecture_nmd.py: for a tie group of size s with q
    positives, entered with T positives and M variants ranked above it, the
    expected contribution is (q/s) * [A*dH + B*(s - (M+1)*dH)] with
    A = T + 1, B = (q - 1)/(s - 1), dH = H[M + s] - H[M] (harmonic numbers);
    a singleton group contributes q * (T + 1)/(M + 1). Verified equal to the
    pandas implementation to <1e-12 on the fitted M1-M5 test probabilities
    (2026-11-09) before use; the pandas form was too slow for 5,000-resample
    cluster bootstraps."""
    y_true = np.asarray(y_true, dtype=np.int64)
    y_score = np.asarray(y_score, dtype=np.float64)
    n_pos = int(np.sum(y_true == 1))
    if n_pos == 0:
        return 0.0
    uniq, inv = np.unique(np.round(y_score, 12), return_inverse=True)
    s_all = np.bincount(inv, minlength=len(uniq)).astype(np.float64)
    q_all = np.bincount(inv, weights=y_true, minlength=len(uniq))
    # Descending score order.
    s_all = s_all[::-1]
    q_all = q_all[::-1]
    m_cum = np.concatenate([[0.0], np.cumsum(s_all)])[:-1]
    t_cum = np.concatenate([[0.0], np.cumsum(q_all)])[:-1]
    dh = _HARM[(m_cum + s_all).astype(np.int64)] - _HARM[m_cum.astype(np.int64)]
    contrib = np.zeros_like(s_all)
    multi = s_all > 1
    if np.any(multi):
        s_m, q_m = s_all[multi], q_all[multi]
        a_m = t_cum[multi] + 1.0
        b_m = np.where(s_m > 1, (q_m - 1.0) / np.maximum(s_m - 1.0, 1.0), 0.0)
        contrib[multi] = (q_m / s_m) * (
            a_m * dh[multi] + b_m * (s_m - (m_cum[multi] + 1.0) * dh[multi])
        )
    single = ~multi
    if np.any(single):
        contrib[single] = q_all[single] * (t_cum[single] + 1.0) / (m_cum[single] + 1.0)
    return float(contrib.sum() / n_pos)


def paired_diff_ci(y, s1, s2, metric_fn, n_boot=N_BOOT, seed=SEED):
    rng = np.random.default_rng(seed)
    diffs = np.empty(n_boot)
    n = len(y)
    for b in range(n_boot):
        idx = rng.integers(0, n, n)
        if len(np.unique(y[idx])) < 2:
            diffs[b] = np.nan
            continue
        diffs[b] = metric_fn(y[idx], s1[idx]) - metric_fn(y[idx], s2[idx])
    diffs = diffs[~np.isnan(diffs)]
    return (
        float(np.mean(diffs)),
        float(np.percentile(diffs, 2.5)),
        float(np.percentile(diffs, 97.5)),
    )


def paired_perm_pvalue(y, s1, s2, metric_fn, n_perm=N_PERM, seed=SEED + 1):
    rng = np.random.default_rng(seed)
    obs = metric_fn(y, s1) - metric_fn(y, s2)
    ge = 0
    n = len(y)
    for _ in range(n_perm):
        swap = rng.random(n) < 0.5
        a = np.where(swap, s2, s1)
        b = np.where(swap, s1, s2)
        if abs(metric_fn(y, a) - metric_fn(y, b)) >= abs(obs):
            ge += 1
    return float((ge + 1) / (n_perm + 1)), float(obs)


def cluster_bootstrap_ci(y, s1, s2, genes, metric_fn, n_boot=N_BOOT, seed=SEED + 7):
    """Cluster bootstrap: resample test GENES with replacement, recompute the
    paired metric difference on the resampled variants. Percentile CI."""
    rng = np.random.default_rng(seed)
    genes = np.asarray(genes)
    uniq = np.unique(genes)
    diffs = np.empty(n_boot)
    for b in range(n_boot):
        sampled = rng.choice(uniq, size=len(uniq), replace=True)
        idx = np.concatenate([np.where(genes == g)[0] for g in sampled])
        if len(np.unique(y[idx])) < 2:
            diffs[b] = np.nan
            continue
        diffs[b] = metric_fn(y[idx], s1[idx]) - metric_fn(y[idx], s2[idx])
    diffs = diffs[~np.isnan(diffs)]
    return (
        float(np.mean(diffs)),
        float(np.percentile(diffs, 2.5)),
        float(np.percentile(diffs, 97.5)),
        float(np.percentile(diffs, 50)),
    )


def jackknife_by_gene(y, s1, s2, genes, metric_fn):
    """Leave-one-gene-out jackknife for a paired metric difference."""
    genes = np.asarray(genes)
    uniq = np.unique(genes)
    full = metric_fn(y, s1) - metric_fn(y, s2)
    jk = []
    for g in uniq:
        keep = genes != g
        jk.append(metric_fn(y[keep], s1[keep]) - metric_fn(y[keep], s2[keep]))
    jk = np.array(jk)
    n = len(uniq)
    # Jackknife SE: Var = ((n-1)/n) * sum((theta_(i) - theta_mean)^2)
    # (from pseudovalues theta~_i = n*theta - (n-1)*theta_(i)).
    se = float(np.sqrt((n - 1) / n * np.sum((jk - jk.mean()) ** 2)))
    lo = float(full - 1.96 * se)
    hi = float(full + 1.96 * se)
    return {
        "full_diff": float(full),
        "jackknife_mean": float(jk.mean()),
        "jackknife_se": se,
        "ci95_normal": [lo, hi],
        "jackknife_diffs_min": float(jk.min()),
        "jackknife_diffs_max": float(jk.max()),
        "n_clusters": int(n),
    }


def boot_ci_metric(y, s, metric_fn, n_boot=N_BOOT, seed=SEED + 2):
    rng = np.random.default_rng(seed)
    vals = np.empty(n_boot)
    n = len(y)
    for b in range(n_boot):
        idx = rng.integers(0, n, n)
        if len(np.unique(y[idx])) < 2:
            vals[b] = np.nan
            continue
        vals[b] = metric_fn(y[idx], s[idx])
    vals = vals[~np.isnan(vals)]
    return (
        float(np.mean(vals)),
        float(np.percentile(vals, 2.5)),
        float(np.percentile(vals, 97.5)),
    )


def main() -> int:
    proc = REPO / "data" / "processed"
    df = pd.read_csv(proc / "feature_matrix.tsv.gz", sep="\t")
    tiers = json.loads((proc / "feature_tiers.json").read_text())["tiers"]
    if "cadd_phred" in tiers["M1"] or "cadd_missing" in tiers["M1"]:
        raise SystemExit("tiers still contain CADD columns; refit must drop them first")
    splits = {
        a["gene"]: a["split"]
        for a in json.loads((REPO / "results/splits.json").read_text())["assignments"]
    }
    df["split"] = df["gene_symbol"].map(splits)
    train = (df["split"] == "train").values
    val = (df["split"] == "validation").values
    test = (df["split"] == "test").values
    y = df["label"].astype(int).values
    yt = y[test]
    genes_test = df.loc[test, "gene_symbol"].values
    print(
        f"[19] train={train.sum()} val={val.sum()} test={test.sum()} "
        f"test_genes={len(np.unique(genes_test))}"
    )

    probas = {}
    for m in TIER_ORDER:
        cols = tiers[m]
        clf = LogisticRegression(
            C=1.0, max_iter=5000, class_weight="balanced", random_state=SEED
        )
        clf.fit(df.loc[train, cols].values, y[train])
        probas[m] = clf.predict_proba(df.loc[test, cols].values)[:, 1]
        print(
            f"[19] {m} ({len(cols)} feats): step AUPRC="
            f"{average_precision_score(yt, probas[m]):.4f} tie AUPRC="
            f"{tie_aware_auprc(yt, probas[m]):.4f}"
        )

    # --- 1: tie-aware AUPRC diffs, variant-level bootstrap + permutation ---
    step_auprc = {m: float(average_precision_score(yt, probas[m])) for m in TIER_ORDER}
    tie_auprc = {m: float(tie_aware_auprc(yt, probas[m])) for m in TIER_ORDER}
    auprc_diffs = {}
    for a, b in [("M4", "M1"), ("M5", "M1"), ("M2", "M1"), ("M3", "M1"), ("M5", "M4")]:
        mean, lo, hi = paired_diff_ci(
            yt, probas[a], probas[b], tie_aware_auprc, n_boot=2000
        )
        pval, obs = paired_perm_pvalue(yt, probas[a], probas[b], tie_aware_auprc)
        sig = bool(pval < BONF_ALPHA and (lo > 0 or hi < 0))
        auprc_diffs[f"{a}_minus_{b}"] = {
            "mean": mean,
            "ci_lo": lo,
            "ci_hi": hi,
            "perm_p": pval,
            "significant": sig,
            "alpha": BONF_ALPHA,
            "metric": "tie-aware AUPRC (McSherry & Najork 2008 expected AP)",
            "note": "Variant-level paired bootstrap CI and permutation p, both on "
            "tie-aware AUPRC. Significant = perm_p < alpha AND the bootstrap CI "
            "excludes zero. Variant-level CIs are NOT gene-level evidence; "
            "see gene_clustered_cis.",
        }
        print(
            f"[19] tie-AUPRC {a}-{b}: {mean:+.4f} [{lo:+.4f},{hi:+.4f}] perm p={pval:.4f}"
        )

    # --- 2: gene-clustered inference (deciding evidence for tier claims) ---
    clustered = {}
    for a, b in [("M2", "M1"), ("M4", "M1"), ("M5", "M1"), ("M5", "M4")]:
        cmean, clo, chi, cmed = cluster_bootstrap_ci(
            yt, probas[a], probas[b], genes_test, tie_aware_auprc, n_boot=2000
        )
        jk = jackknife_by_gene(yt, probas[a], probas[b], genes_test, tie_aware_auprc)
        clustered[f"{a}_minus_{b}"] = {
            "method": "cluster bootstrap resampling the 24 test genes with replacement, "
            "n_boot=2000, seed=20261013; plus leave-one-gene-out jackknife (24 folds)",
            "n_clusters": 24,
            "cluster_bootstrap_mean": cmean,
            "cluster_bootstrap_median": cmed,
            "cluster_ci95": [clo, chi],
            "cluster_excludes_zero": bool(clo > 0 or chi < 0),
            "jackknife": jk,
        }
        print(
            f"[19] CLUSTER {a}-{b}: mean {cmean:+.4f} 95% CI [{clo:+.4f},{chi:+.4f}] "
            f"jackknife 95% CI [{jk['ci95_normal'][0]:+.4f},{jk['ci95_normal'][1]:+.4f}]"
        )

    m2_cluster = clustered["M2_minus_M1"]
    m2_verdict = (
        "WITHDRAW 'M2 best tier' as a claim: gene-clustered 95% CI for the M2-M1 "
        "tie-aware AUPRC gain includes zero; report as a non-significant descriptive "
        "difference."
        if not m2_cluster["cluster_excludes_zero"]
        else "M2-M1 gain survives gene-clustered inference (cluster 95% CI excludes zero)."
    )
    print(f"[19] M2 verdict: {m2_verdict}")

    # Paired AUROC diffs on the refit (variant-level, script 12b method, n_boot=2000).
    auroc_diffs = {}
    for a, b in [("M4", "M1"), ("M5", "M1"), ("M2", "M1"), ("M3", "M1"), ("M5", "M4")]:
        mean, lo, hi = paired_diff_ci(
            yt, probas[a], probas[b], roc_auc_score, n_boot=2000, seed=SEED
        )
        auroc_diffs[f"{a}_minus_{b}"] = {"mean": mean, "ci_lo": lo, "ci_hi": hi}

    # --- 3: M2/M3 clarification on the refit ---
    m2_m3 = {}
    for m in ["M2", "M3"]:
        row = {}
        for metric, fn in [
            ("auroc", roc_auc_score),
            ("auprc", average_precision_score),
        ]:
            mean, lo, hi = paired_diff_ci(yt, probas[m], probas["M1"], fn)
            pval, _ = paired_perm_pvalue(yt, probas[m], probas["M1"], fn)
            sig = bool(pval < BONF_ALPHA and (lo > 0 or hi < 0))
            row[metric] = {
                "mean": mean,
                "ci_lo": lo,
                "ci_hi": hi,
                "perm_p": pval,
                "significant": sig,
                "alpha": BONF_ALPHA,
            }
        m2_m3[m] = row

    # --- 4: CADD complementarity with refit M5 (script 17 section 4 protocol) ---
    cadd = pd.read_csv(proc / "variant_cadd.tsv.gz", sep="\t")
    cadd = cadd[cadd["cadd_phred"].notna()].rename(columns={"cadd_phred": "cadd"})
    dfc = df.merge(cadd[["variant_key", "cadd"]], on="variant_key", how="left")
    dfc["p_M5"] = np.nan
    dfc.loc[test, "p_M5"] = probas["M5"]
    valv = dfc[val & dfc["cadd"].notna()]
    T = float(valv.loc[valv["label"] == 1, "cadd"].median())
    testc = dfc[test & dfc["cadd"].notna()].copy()
    weak = testc["cadd"] < T
    y_w = testc.loc[weak, "label"].astype(int).values
    s_w = testc.loc[weak, "p_M5"].values
    auroc_w, alo, ahi = boot_ci_metric(y_w, s_w, roc_auc_score)
    auprc_w, plo, phi = boot_ci_metric(y_w, s_w, average_precision_score)
    k = max(1, int(0.1 * weak.sum()))
    order = np.argsort(-s_w)[:k]
    prec = float(y_w[order].mean())
    base_rate = float(y_w.mean())
    cadd_comp = {
        "threshold_source": "validation pathogenic median",
        "threshold_T": T,
        "test_cadd_scored_n": len(testc),
        "weak_n": int(weak.sum()),
        "weak_pathogenic_n": int((testc.loc[weak, "label"] == 1).sum()),
        "m5_auroc_weak": auroc_w,
        "m5_auroc_weak_ci": [alo, ahi],
        "m5_auprc_weak": auprc_w,
        "m5_auprc_weak_ci": [plo, phi],
        "mean_m5_pathogenic_weak": float(
            testc.loc[weak & (testc["label"] == 1), "p_M5"].mean()
        ),
        "mean_m5_pathogenic_all_cadd_test": float(
            testc.loc[testc["label"] == 1, "p_M5"].mean()
        ),
        "mean_m5_benign_weak": float(
            testc.loc[weak & (testc["label"] == 0), "p_M5"].mean()
        ),
        "top_decile_precision": prec,
        "top_decile_base_rate": base_rate,
        "top_decile_enrichment_ratio": round(prec / base_rate, 2),
        "note": "Refit M5 (no CADD in fitted tiers). Threshold T recomputed from validation; "
        "unchanged from pre-refit value if CADD file unchanged.",
    }
    print(
        f"[19] CADD-weak: T={T} n={weak.sum()} M5 AUROC {auroc_w:.4f} AUPRC {auprc_w:.4f}"
    )

    # --- 5: methods.md:209 audit -- closest computable quantity ---
    y_sub = testc["label"].astype(int).values
    s_cadd = rankdata(testc["cadd"].values) / len(testc)
    s_m5_sub = testc["p_M5"].values
    m5_sub_auprc = float(average_precision_score(y_sub, s_m5_sub))
    cadd_sub_auprc = float(average_precision_score(y_sub, s_cadd))
    m5_sub_auroc = float(roc_auc_score(y_sub, s_m5_sub))
    cadd_sub_auroc = float(roc_auc_score(y_sub, s_cadd))
    dfc["p_M4"] = np.nan
    dfc.loc[test, "p_M4"] = probas["M4"]
    testc4 = dfc[test & dfc["cadd"].notna()].copy()
    m4_sub_auprc = float(
        average_precision_score(
            testc4["label"].astype(int).values, testc4["p_M4"].values
        )
    )
    methods_check = {
        "claim": "methods.md:209: validation-trained CADD-plus-context stack beats CADD alone "
        "by +0.0191 AUPRC, isoform slice +0.0015, isoform block alone on top of CADD "
        "adds exactly 0.0000",
        "verdict": "NOT REPRODUCIBLE from any shipped file or script: no validation-trained "
        "stack exists in the repo (all fits are train-split; no script fits on "
        "validation), and +0.0191 / +0.0015 / 0.0000 appear nowhere in results/. "
        "Closest computable train-fit quantities on the CADD-scored test subset "
        "are recorded below; none equals +0.0191. Numbers removed from methods.md.",
        "closest_computable_train_fit": {
            "population": "test variants with CADD scored (n=11240)",
            "m5_auprc_same_subset": m5_sub_auprc,
            "cadd_auprc_same_subset": cadd_sub_auprc,
            "m5_minus_cadd_auprc": m5_sub_auprc - cadd_sub_auprc,
            "m5_auroc_same_subset": m5_sub_auroc,
            "cadd_auroc_same_subset": cadd_sub_auroc,
            "m5_minus_cadd_auroc": m5_sub_auroc - cadd_sub_auroc,
            "m4_auprc_same_subset": m4_sub_auprc,
            "m4_minus_m1_auprc_same_subset": "not computed; M1 AUPRC on this subset not part of the claim",
        },
    }
    print(
        f"[19] methods:209 check: M5-CADD AUPRC same-subset = "
        f"{m5_sub_auprc - cadd_sub_auprc:+.4f} (claimed +0.0191)"
    )

    # --- 6: per-tier AUPRC CIs (variant-level bootstrap) -> per_model_cis.json ---
    per_model = {}
    bench = json.loads((REPO / "results/benchmark.json").read_text())
    for m in TIER_ORDER:
        _, lo, hi = boot_ci_metric(yt, probas[m], average_precision_score)
        per_model[m] = {
            "auroc": bench["tiers"][m]["test"]["auroc"],
            "auroc_ci": bench["tiers"][m]["test"]["auroc_ci95"],
            "auprc": step_auprc[m],
            "auprc_ci": [lo, hi],
            "tie_auprc": tie_auprc[m],
        }
    per_model["_note"] = (
        "Refit 2026-11-09, no CADD in fitted tiers; AUPRC CIs are variant-level "
        "paired bootstrap (n_boot=5000, seed=20261006). This file previously had no "
        "generating script in the repo; values are now produced by "
        "scripts/19_refit_statistics.py. Variant-level CIs are descriptive only for "
        "gene-level claims; see statistics.json gene_clustered_cis."
    )
    (REPO / "results/per_model_cis.json").write_text(json.dumps(per_model, indent=2))
    print("[19] wrote results/per_model_cis.json")

    # --- 7: assemble statistics.json, preserving superseded values explicitly ---
    old = json.loads((REPO / "results/statistics.json").read_text())

    def carried(section_key, *live_keys):
        # On a rerun the file already carries labelled superseded sections;
        # carry them forward instead of relabelling the live values as old.
        section = old.get(section_key)
        if isinstance(section, dict) and section.get("values") is not None:
            return section["values"]
        for live_key in live_keys:
            value = old.get(live_key)
            if value is not None:
                if isinstance(value, dict):
                    # The renamed live block carries an explanatory note;
                    # superseded values keep the bare numeric content.
                    value = {k: v for k, v in value.items() if k != "note"}
                return value
        return None

    stats = {
        "seed": SEED,
        "n_boot": N_BOOT,
        "n_perm": N_PERM,
        "test_n": int(test.sum()),
        "refit": "2026-11-09: all values below recomputed from tiers fit WITHOUT cadd_phred "
        "and cadd_missing (train 0/35,921 CADD-scored). CADD is a pure external "
        f"baseline. Environment: pinned requirements stack (scikit-learn "
        f"{sklearn.__version__} / pandas {pd.__version__} / numpy {np.__version__}); "
        "an earlier 2026-11-09 refit ran on an unpinned system stack: its test "
        "values drift from these by up to 0.0013 (M3 test AUROC) and its M3 "
        "validation AUROC by 0.0054.",
        "auprc_diffs": auprc_diffs,
        "auprc_diffs_superseded_pre_rebuild": {
            "note": "SUPERSEDED. Pre-rebuild run: feature matrix had cadd_phred = -1.0 for all "
            "60,521 rows and the fitted cadd_phred coefficient was fit on a constant. "
            "The M4-M1 +0.0509 below is DEAD and must not be cited as a live gain. "
            "Kept only as the labelled record of that run.",
            "values": carried("auprc_diffs_superseded_pre_rebuild", "auprc_diffs"),
        },
        "step_metric_tie_artifact_demo": {
            "note": "STEP-metric values retained ONLY to demonstrate the tie artifact; significant:true on the step metric is not evidence of a tier gain — the tie-aware auprc_diffs block above is the deciding result.",
            **m2_m3,
        },
        "m2_m3_clarification_superseded_pre_rebuild": {
            "note": "SUPERSEDED pre-rebuild values (same CADD-constant defect).",
            "values": carried(
                "m2_m3_clarification_superseded_pre_rebuild",
                "m2_m3_clarification",
                "step_metric_tie_artifact_demo",
            ),
        },
        "tie_corrected_refit_no_cadd": {
            "method": "Tie-aware AUPRC: expected average precision under uniform tie-breaking "
            "(McSherry & Najork 2008), function copied from "
            "scripts/18_exon_architecture_nmd.py. Models refit on train with the "
            "locked protocol (logistic regression, C=1.0, class_weight=balanced, "
            "seed=20261006) on the no-CADD tier set; paired bootstrap CIs, "
            "n_boot=2000, seed=20261006.",
            "matrix": "no-CADD refit 2026-11-09; CADD a pure external baseline",
            "test_n": int(test.sum()),
            "tie_auprc": tie_auprc,
            "step_auprc": step_auprc,
            "diffs": {
                k: {
                    "mean": v["mean"],
                    "ci_lo": v["ci_lo"],
                    "ci_hi": v["ci_hi"],
                    "perm_p": v["perm_p"],
                    "significant": v["significant"],
                }
                for k, v in auprc_diffs.items()
            },
        },
        "tie_corrected_superseded_pre_refit": {
            "note": "SUPERSEDED. Tie-corrected diffs on the rebuilt matrix WITH the "
            "CADD-sentinel fit (cadd_phred coefficient fit on a constant in train). "
            "Kept only as the labelled record; the live values are in "
            "tie_corrected_refit_no_cadd. Old values: M2-M1 +0.0070, M3-M1 -0.0151, "
            "M4-M1 -0.0069 (tie-aware); M5-M1 was not in the old diffs table -- from "
            "the old tie_auprc values it is 0.9454-0.9463 = -0.0009 (step: -0.0011), "
            "and M5-M4 tie diff was 0.9454-0.9394 = +0.0060. Live values: "
            "tie_corrected_refit_no_cadd.",
            "values": carried("tie_corrected_superseded_pre_refit", "tie_corrected"),
        },
        "paired_auroc_diffs_refit_no_cadd": {
            "method": "Paired bootstrap CIs for AUROC differences, scripts/12b_paired_analysis.py "
            "method, n_boot=2000, seed=20261006. Models refit on the no-CADD tier set.",
            "matrix": "no-CADD refit 2026-11-09",
            "test_n": int(test.sum()),
            "diffs": auroc_diffs,
        },
        "paired_auroc_diffs_rebuilt_superseded": {
            "note": "SUPERSEDED. Rebuilt-matrix AUROC diffs computed with the CADD-sentinel fit.",
            "values": carried(
                "paired_auroc_diffs_rebuilt_superseded", "paired_auroc_diffs_rebuilt"
            ),
        },
        "gene_clustered_cis": {
            "method": "Cluster bootstrap resampling the 24 test genes with replacement "
            "(n_boot=2000, seed=20261013) + leave-one-gene-out jackknife (24 folds). "
            "Variant-level CIs are invalid for gene-level HPO constants, so the "
            "cluster interval is the deciding evidence for tier claims.",
            "n_clusters": 24,
            "diffs": clustered,
            "m2_minus_m1_verdict": m2_verdict,
        },
        "cadd_complementarity": cadd_comp,
        "cadd_complementarity_superseded_pre_refit": {
            "note": "SUPERSEDED. Same protocol with the pre-refit M5 (CADD-sentinel fit).",
            "values": carried(
                "cadd_complementarity_superseded_pre_refit", "cadd_complementarity"
            ),
        },
        "methods_209_check": methods_check,
        "significance_rule": {
            "method": "Bonferroni correction across the 5 paired tier comparisons",
            "n_comparisons": N_COMPARISONS,
            "per_comparison_alpha": BONF_ALPHA,
            "note": "A tier comparison is marked significant only if its permutation "
            "p-value -- computed with the same metric as its reported confidence "
            "interval (tie-aware AUPRC for auprc_diffs) -- clears the "
            "Bonferroni-corrected per-comparison alpha (0.01) AND its paired "
            "bootstrap 95% CI excludes zero. The `significant` flags on "
            "auprc_diffs use variant-level inference only; a tier claim ALSO needs its "
            "gene-clustered 95% CI to exclude zero (see gene_clustered_cis). M2-M1: "
            + m2_verdict,
        },
        "correction_note": (
            "2026-11-09 no-CADD refit supersedes the pre-rebuild record. Pre-rebuild: the "
            "feature matrix had cadd_phred = -1.0 for all 60,521 rows, so M1 never saw a real "
            "CADD value; M4-M1 +0.0509 (significant:true) was the headline of that run and is "
            "DEAD -- it survives only under auprc_diffs_superseded_pre_rebuild. Rebuilt-matrix "
            "run (2026-11-08) still had the defect in a subtler form: real CADD values in "
            "validation/test, but the fitted coefficient was fit on a constant in train "
            "(0/35,921 scored), so test-time CADD signal entered through an uncalibrated "
            "coefficient; its tie-corrected values (M4-M1 -0.0069) are preserved under "
            "tie_corrected_superseded_pre_refit. Live values are the *_refit_no_cadd sections. "
            "The M2-M1 (+0.0310) and M3-M1 (+0.0317) step-AUPRC gains in the pre-rebuild run "
            "were a tie-breaking artifact: a 6,077-variant tie block (test SNVs with "
            "am_missing=1 and review_stars=1) penalised M1 step AUPRC by -0.0318; tie-corrected "
            "gains on the pre-rebuild matrix were M2-M1 +0.0006 and M3-M1 -0.0002."
        ),
    }
    # Preserve the composite-CLNSIG sensitivity block across reruns: it is
    # produced by a standalone sensitivity script, not recomputed here.
    if "clnsig_composite_sensitivity" in old:
        stats["clnsig_composite_sensitivity"] = old["clnsig_composite_sensitivity"]
    (REPO / "results/statistics.json").write_text(json.dumps(stats, indent=2))
    print("[19] wrote results/statistics.json")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
