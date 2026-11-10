#!/usr/bin/env python3
"""Stage 18: exon architecture, isoform skipping, and NMD boundary analysis.

Reconstructs MANE Select coding exons across the 102 panel genes from the
587,227-row ClinVar table and evaluates:
  1. The 4 x 2 x 2 internal exon length x GTEx backup skipping (>=10% usage,
     >=1.0 TPM) x reading-frame phase symmetry (length mod 3 == 0) matrix.
  2. AlphaMissense/CADD distributions that do not use ClinVar labels,
     and coding density per kb.
  3. 5'-splice donor (3' exon end) vs. 3'-splice acceptor (5' exon start)
     asymmetry among near-junction missense variants with am_score < 0.34,
     in-frame deletion vs. insertion polarity, and 50-nt penultimate-exon
     EJC NMD escape zones.
  4. Within-gene Cochran-Mantel-Haenszel tests for the exon-size and
     backup-skipping effects on internal missense pathogenicity.
  5. Ranking of computationally nominated exons.
  6. Tiered-model comparison (M1/M3/M4/Mechanistic-v3, called M6 in the
     manuscript) with a consequence-class-only baseline, per-class and
     within-missense metrics, and Bayesian likelihood ratios for VUS
     stratification.

Donor-edge offsets in the missense architecture adjustment are fit by
logistic regression on TRAIN genes only, with the M1 logit held fixed as
an offset; fitted values are reported in the output JSON. Earlier versions
of this script used hardcoded offsets (+3.60/+1.80/+1.80) under a comment
claiming they were learned; no fitting code existed for them.

Likelihood-ratio formula (one formula for every LR in the outputs):
    LR(cell) = P(cell | P/LP) / P(cell | B/LB)
computed on resolved variants only, with no smoothing. The 95% CI is the
log-LR interval with SE = sqrt(1/a - 1/N1 + 1/b - 1/N0), where a is the
P/LP count in the cell, N1 the total resolved P/LP count, b the B/LB count
in the cell, and N0 the total resolved B/LB count. A Haldane-Anscombe 0.5
correction is used only if a cell count is zero, and is flagged when used.
LRs are research-only stratification; they are not clinical classifications
and must not be converted into ACMG/AMP evidence points.

Writes:
  - results/exon_architecture_nmd.json
  - results/exon_skipping_nmd_candidates.tsv
  - results/vus_bayesian_stratification.tsv
  - results/cmh_per_gene_tables.tsv
  - results/vus_bayesian_stratification_NOTE.md
  - results/withdrawn_carried_over_values.json (written once, only if it
    does not already exist; see below)

Two sections of the previously shipped JSON cannot be regenerated because
no code that produces them exists in the repo: causal_mediation, and the
distance_profile / donor_vs_acceptor_edge_0_3bp parts of
splice_junction_boundary_asymmetry. They are NOT emitted into
exon_architecture_nmd.json. On 2026-11-09 they were moved into
results/withdrawn_carried_over_values.json (status: withdrawn; no
generating code; do not cite). This script writes that file once, copying
the blocks verbatim from the previous results file, only if the file does
not already exist; it never recomputes those values.

Usage: python3 scripts/18_exon_architecture_nmd.py
"""

from __future__ import annotations

import json
import math
from datetime import date
from itertools import pairwise
from pathlib import Path

import numpy as np
import pandas as pd
from scipy.optimize import minimize
from scipy.special import expit
from scipy.stats import fisher_exact, rankdata

REPO = Path(__file__).resolve().parents[1]
SEED = 20261006

SIZE_BINS = ["1_Compact_<=180", "2_Medium_181_400", "3_Large_401_1000", "4_Giant_>1000"]
# Primary backup-skipping usage threshold (GTEx transcript-abundance
# fraction). All shipped primary results use this value; the sensitivity
# sweep in sensitivity_checks/usage_threshold_sweep recomputes the backup
# CMH OR and the giant-cell LRs at 0.05 / 0.10 / 0.20 without changing it.
BACKUP_USAGE_THRESHOLD = 0.10
# Hand-applied correction note carried in the committed JSON (2026-11-09);
# kept verbatim so regeneration preserves it.
TENSOR_CORRECTION_NOTE = "CORRECTED 2026-11-09: per-cell n_exons/n_genes previously summed to 4,528 exons (+19 over the 4,509 internal coding exons in dataset_summary). Cause: the tensor was counted with nunique(exon_id) per cell on the VARIANT-level frame, where expr_backup_skips is assigned per variant; 19 exons whose variants received mixed expr_backup_skips values were counted in both a constitutive (backup_skip=0) cell and a backup-skipped cell. All 19 excess counts were in backup_skip=0 cells (by size bin: Compact +10, Medium +6, Large +1, Giant +2; n_genes inflated by 7 in the same cells). Fix: n_exons/n_genes per cell now taken from the exon-level atlas (results/exon_skipping_nmd_candidates.tsv, one row per exon, expr_backup_skips = max over the exon's variants), which sums to 4,509 and agrees with biophysical_and_density_controls (constitutive 2,918 + backup-skipped 1,591 = 4,509). Variant-level counts per cell (missense_n etc.) were NOT changed: variants of the 19 exons remain split between cells exactly as computed, and reassigning them requires the script-18 input files (data/interim/script18/), which are not shipped. Rule for regeneration: resolve expr_backup_skips per exon (max) BEFORE cell assignment, and count exons from the exon-level frame."
CMH_P_NOTE_TEST1 = "CMH p-values removed 2026-11-09 for the same reason as cmh_within_gene_missense: variant-level independence assumption; see that section's note. ORs retained as descriptive values."
Z_WITHHELD_NOTE = "z withheld for the same reason as p: variant-level independence assumption; effective n is genes, not variants"
# (feature name, column expression) for the missense architecture offset.
# Coefficients are fit on train internal missense only; the M1 logit is a
# fixed offset (coefficient 1.0), never refit.
OFFSET_FEATURES = [
    "log_exon_len_c",
    "mis_x_giant_exon",
    "mis_x_backup_skip",
    "mis_x_giant_x_skip_x_sym",
    "donor_edge_x_spliceai_ge05",
    "donor_edge",
    "spliceai_max_delta",
]


def load_clinvar_587k(repo_root: Path) -> pd.DataFrame:
    """Load the 587,227-row ClinVar table from .tsv.gz or .tsv."""
    gz_path = repo_root / "data/interim/script18/supp_clinvar_all_table.tsv.gz"
    tsv_path = repo_root / "data/interim/script18/supp_clinvar_all_table.tsv"
    if gz_path.exists():
        return pd.read_csv(gz_path, sep="\t")
    if tsv_path.exists():
        return pd.read_csv(tsv_path, sep="\t")
    raise FileNotFoundError(f"Missing {gz_path} (or uncompressed {tsv_path})")


def roc_auc_score_np(y_true, y_score) -> float:
    y_true = np.asarray(y_true, dtype=np.int64)
    y_score = np.asarray(y_score, dtype=np.float64)
    n1 = int(np.sum(y_true == 1))
    n0 = len(y_true) - n1
    if n1 == 0 or n0 == 0:
        return float("nan")
    ranks = rankdata(y_score, method="average")
    return float((np.sum(ranks[y_true == 1]) - n1 * (n1 + 1) / 2.0) / (n1 * n0))


def tie_aware_auprc(y_true, y_score) -> float:
    """Analytical expectation of average precision under uniform tie-breaking."""
    y_true = np.asarray(y_true, dtype=np.int64)
    y_score = np.asarray(y_score, dtype=np.float64)
    n_pos = int(np.sum(y_true == 1))
    if n_pos == 0:
        return 0.0
    df_t = pd.DataFrame({"y": y_true, "s": np.round(y_score, 12)})
    grp = (
        df_t.groupby("s", sort=True).agg(n=("y", "count"), pos=("y", "sum")).iloc[::-1]
    )
    tp_cum, n_cum, ap_sum = 0, 0, 0.0
    for row in grp.itertuples():
        n_g, p_g = row.n, row.pos
        if p_g > 0:
            if n_g == 1:
                ap_sum += (tp_cum + 1.0) / (n_cum + 1.0)
            else:
                j = np.arange(1, n_g + 1, dtype=np.float64)
                e_tp_j = tp_cum + 1.0 + (j - 1.0) * (p_g - 1.0) / (n_g - 1.0)
                ap_sum += float(np.sum((p_g / n_g) * (e_tp_j / (n_cum + j))))
        tp_cum += p_g
        n_cum += n_g
    return float(ap_sum / n_pos)


def compute_ece(y_true, y_prob, n_bins: int = 10) -> float:
    bins = np.linspace(0.0, 1.0, n_bins + 1)
    ece = 0.0
    n = len(y_true)
    for b0, b1 in pairwise(bins):
        m = (y_prob >= b0) & (y_prob <= (b1 if b1 == 1.0 else b1 - 1e-12))
        if m.any():
            ece += (m.sum() / n) * abs(y_true[m].mean() - y_prob[m].mean())
    return float(ece)


class BalancedLogReg:
    """L2 logistic regression with balanced class weights, fit by L-BFGS."""

    def __init__(self, C: float = 1.0):
        self.C = C

    def fit(self, X, y):
        X = np.asarray(X, dtype=np.float64)
        y = np.asarray(y, dtype=np.float64)
        self.mean_ = X.mean(axis=0)
        self.scale_ = X.std(axis=0)
        self.scale_[self.scale_ == 0] = 1.0
        Xs = (X - self.mean_) / self.scale_
        n, n1 = len(y), np.sum(y == 1)
        sw = np.where(y == 1, n / (2.0 * n1), n / (2.0 * (n - n1)))

        def loss_grad(wb):
            w, b = wb[:-1], wb[-1]
            prob = expit(Xs @ w + b)
            ll = (
                np.sum(
                    sw
                    * np.where(
                        y == 1,
                        -np.log(np.clip(prob, 1e-15, 1.0)),
                        -np.log(np.clip(1.0 - prob, 1e-15, 1.0)),
                    )
                )
                + 0.5 * np.sum(w * w) / self.C
            )
            diff = sw * (prob - y)
            return ll, np.append((Xs.T @ diff) + w / self.C, np.sum(diff))

        res = minimize(
            loss_grad,
            np.zeros(Xs.shape[1] + 1),
            jac=True,
            method="L-BFGS-B",
            options={"maxiter": 5000, "ftol": 1e-12, "gtol": 1e-8},
        )
        self.coef_, self.intercept_ = res.x[:-1], res.x[-1]
        return self

    def predict_proba(self, X):
        Xs = (np.asarray(X, dtype=np.float64) - self.mean_) / self.scale_
        p1 = expit(Xs @ self.coef_ + self.intercept_)
        return np.column_stack([1.0 - p1, p1])


def lr_with_ci(a_cell: int, n_plp: int, b_cell: int, n_blb: int) -> dict:
    """LR = P(cell|P/LP) / P(cell|B/LB) with a 95% CI on the log scale.

    No smoothing unless a cell count is zero, in which case a
    Haldane-Anscombe 0.5 correction is applied to all four counts and
    flagged. Returns the LR, CI bounds, and the counts used.
    """
    corrected = False
    a, n1, b, n0 = float(a_cell), float(n_plp), float(b_cell), float(n_blb)
    if a == 0 or b == 0:
        a, n1, b, n0 = a + 0.5, n1 + 0.5, b + 0.5, n0 + 0.5
        corrected = True
    lr = (a / n1) / (b / n0)
    se = math.sqrt(1.0 / a - 1.0 / n1 + 1.0 / b - 1.0 / n0)
    log_lr = math.log(lr)
    lo = math.exp(log_lr - 1.96 * se)
    hi = math.exp(log_lr + 1.96 * se)
    return {
        "likelihood_ratio": float(lr),
        "lr_ci95_lo": float(lo),
        "lr_ci95_hi": float(hi),
        "haldane_correction_applied": corrected,
        "n_plp_in_cell": int(a_cell),
        "n_plp_total": int(n_plp),
        "n_blb_in_cell": int(b_cell),
        "n_blb_total": int(n_blb),
    }


def wilson_ci(k: int, n: int, z: float = 1.96) -> list:
    """Wilson score interval for a binomial proportion."""
    if n == 0:
        return [None, None]
    p = k / n
    den = 1.0 + z * z / n
    centre = (p + z * z / (2.0 * n)) / den
    half = z * math.sqrt(p * (1.0 - p) / n + z * z / (4.0 * n * n)) / den
    return [float(max(0.0, centre - half)), float(min(1.0, centre + half))]


def run_cmh_stratified(df_sub: pd.DataFrame, x_col: str, y_col: str = "is_plp") -> dict:
    """Cochran-Mantel-Haenszel OR over per-gene 2x2 tables.

    Returns the summary OR/z/p plus the per-gene tables, so the result can
    be shipped and checked without re-running. Genes with no variation in
    x or y contribute no table. P-values treat variants as independent
    within genes; the effective sample size is genes and exons, not
    variants, so the p-value must not be read as gene-level evidence.
    """
    num_or, den_or, num_z, var_z, n_g = 0.0, 0.0, 0.0, 0.0, 0
    per_gene = []
    for g, grp in df_sub.groupby("gene_symbol"):
        if grp[x_col].nunique() < 2 or grp[y_col].nunique() < 2:
            continue
        x = grp[x_col].to_numpy()
        y = grp[y_col].to_numpy()
        a = int(np.sum((x == 1) & (y == 1)))
        b = int(np.sum((x == 1) & (y == 0)))
        c = int(np.sum((x == 0) & (y == 1)))
        d = int(np.sum((x == 0) & (y == 0)))
        n = a + b + c + d
        if n <= 1:
            continue
        per_gene.append({"gene_symbol": g, "a": a, "b": b, "c": c, "d": d, "n": n})
        num_or += (a * d) / n
        den_or += (b * c) / n
        r1, r0, c1, c0 = a + b, c + d, a + c, b + d
        num_z += a - (r1 * c1) / n
        var_z += (r1 * r0 * c1 * c0) / (n * n * (n - 1.0))
        n_g += 1
    cmh_or = num_or / max(den_or, 1e-12)
    z = num_z / np.sqrt(max(var_z, 1e-12))
    # P-values are WITHHELD (2026-11-09): the CMH z/p treats variants as
    # independent within genes; the effective sample size is genes and exons,
    # not variants, so the p-value massively overstates the evidence (the
    # shipped values were 4.29e-187 and 3.53e-42). Report the OR, the z
    # statistic and the informative-gene count as descriptive values only.
    return {
        "cmh_or": float(cmh_or),
        "z": float(z),
        "p": None,
        "p_withheld": "variant-level independence assumption; effective n is genes, not variants",
        "n_informative_genes": n_g,
        "per_gene_tables": per_gene,
    }


def gene_cluster_bootstrap(
    per_gene_tables: list, cmh_or: float, n_boot: int = 2000, seed: int = 0
) -> dict:
    """Gene-cluster bootstrap for a CMH odds ratio.

    Resamples GENES (the per-gene 2x2 tables) with replacement and
    recomputes the CMH OR = sum(a*d/n) / sum(b*c/n) per replicate, so the
    interval reflects gene-level, not variant-level, sampling. Percentile
    95% CI. Same construction as the review panel's recompute (seed 0,
    2,000 reps).
    """
    arr = np.array(
        [[r["a"], r["b"], r["c"], r["d"], r["n"]] for r in per_gene_tables],
        dtype=np.float64,
    )
    rng = np.random.default_rng(seed)
    ors = []
    for _ in range(n_boot):
        idx = rng.integers(0, len(arr), len(arr))
        s = arr[idx]
        num = float(np.sum(s[:, 0] * s[:, 3] / s[:, 4]))
        den = float(np.sum(s[:, 1] * s[:, 2] / s[:, 4]))
        if den > 0:
            ors.append(num / den)
    ors = np.array(ors, dtype=np.float64)
    return {
        "or": float(cmh_or),
        "ci95": [
            float(np.percentile(ors, 2.5)),
            float(np.percentile(ors, 97.5)),
        ],
        "n_boot": int(len(ors)),
        "seed": int(seed),
        "frac_ge_1": float(np.mean(ors >= 1.0)),
        "method": "percentile bootstrap resampling genes (per-gene 2x2 tables, results/cmh_per_gene_tables.tsv) with replacement; CMH OR recomputed per replicate",
    }


def cmh_shipped(res: dict, bootstrap: dict | None = None) -> dict:
    """The CMH block as shipped in the JSON.

    P-values are omitted and the z statistic is withheld (set to null):
    a z of this size is the p-value for any reader who converts it, and
    both rest on the variant-level independence assumption. The OR and
    the informative-gene count are descriptive values.
    """
    out = {
        "cmh_or": res["cmh_or"],
        "z": None,
        "z_withheld": Z_WITHHELD_NOTE,
        "n_informative_genes": res["n_informative_genes"],
    }
    if bootstrap is not None:
        out["gene_cluster_bootstrap"] = bootstrap
    return out


def expr_flags_for_threshold(pop: pd.DataFrame, alt_cds: pd.DataFrame, thr: float):
    """Per-variant expr_backup_skips flags at an arbitrary usage threshold.

    Replicates the primary sparing construction (transcript-level
    sev == 0, max_usage >= thr, max_tpm >= 1.0, max over the variant's
    transcripts, mapped onto pop rows by the same within-gene
    nearest-position match the primary path uses). Returns an int array
    aligned to pop's index labels (pop carries a RangeIndex).
    """
    flag = (
        (alt_cds["sev"] == 0)
        & (alt_cds["max_usage"] >= thr)
        & (alt_cds["max_tpm"] >= 1.0)
    ).astype(int)
    tmp = alt_cds[["gene_symbol", "variant_key", "pos"]].assign(_flag=flag)
    grouped = (
        tmp.groupby(["gene_symbol", "variant_key", "pos"], as_index=False)["_flag"]
        .max()
        .rename(columns={"_flag": "expr_flag"})
    )
    by_gene = {g: gdf.sort_values("pos") for g, gdf in grouped.groupby("gene_symbol")}
    arr = np.zeros(len(pop), dtype=int)
    for g, gdf in pop.groupby("gene_symbol"):
        sdf = by_gene.get(g)
        if sdf is not None and len(sdf) > 0:
            p_ref = sdf["pos"].to_numpy()
            p_q = gdf["pos"].to_numpy()
            j = np.clip(np.searchsorted(p_ref, p_q), 0, len(p_ref) - 1)
            j_prev = np.clip(j - 1, 0, len(p_ref) - 1)
            best = np.where(
                np.abs(p_ref[j] - p_q) < np.abs(p_ref[j_prev] - p_q), j, j_prev
            )
            arr[gdf.index] = sdf["expr_flag"].to_numpy()[best]
    return arr


def fit_arch_offsets(tr_mis: pd.DataFrame) -> dict:
    """Fit the missense architecture offset on TRAIN internal missense only.

    Logistic regression with the M1 logit as a fixed offset (coefficient
    1.0) and the seven OFFSET_FEATURES as free coefficients, L2 penalty
    0.5 * ||w||^2, balanced class weights. Approximate Wald CIs come from
    the L-BFGS inverse-Hessian diagonal. Returns fitted coefficients, CIs,
    and the training counts, including how many train variants carry the
    rare donor-edge flag (the term with the least data behind it).
    """
    X = tr_mis[OFFSET_FEATURES].to_numpy(dtype=np.float64)
    y = tr_mis["label"].to_numpy(dtype=np.float64)
    offset = tr_mis["m1_logit"].to_numpy(dtype=np.float64)
    n, n1 = len(y), float(np.sum(y == 1))
    sw = np.where(y == 1, n / (2.0 * n1), n / (2.0 * (n - n1)))

    def loss_grad(w):
        eta = offset + X @ w
        prob = expit(eta)
        ll = float(
            np.sum(
                sw
                * np.where(
                    y == 1,
                    -np.log(np.clip(prob, 1e-15, 1.0)),
                    -np.log(np.clip(1.0 - prob, 1e-15, 1.0)),
                )
            )
            + 0.5 * float(np.sum(w * w))
        )
        grad = X.T @ (sw * (prob - y)) + w
        return ll, grad

    res = minimize(
        loss_grad,
        np.zeros(X.shape[1]),
        jac=True,
        method="L-BFGS-B",
        options={"maxiter": 5000, "ftol": 1e-12, "gtol": 1e-8},
    )
    hess_inv = np.asarray(res.hess_inv.todense())
    se = np.sqrt(np.maximum(np.diag(hess_inv), 0.0))
    coefs = {}
    for i, name in enumerate(OFFSET_FEATURES):
        n_nonzero = int(np.sum(X[:, i] != 0))
        identified = n_nonzero > 0
        coefs[name] = {
            # A feature with no variation in the training subset has no
            # information to fit: its coefficient stays at the 0.0 start and
            # is flagged, not reported as a fitted value.
            "coef": float(res.x[i]) if identified else 0.0,
            "se_approx": float(se[i]) if identified else None,
            "ci95_lo": float(res.x[i] - 1.96 * se[i]) if identified else None,
            "ci95_hi": float(res.x[i] + 1.96 * se[i]) if identified else None,
            "n_nonzero_in_train": n_nonzero,
            "identified_in_train": identified,
        }
    return {
        "method": "logistic regression on train internal missense only; m1_logit held fixed as offset (coef 1.0); L2 penalty 0.5*||w||^2; balanced class weights",
        "n_train_internal_missense": int(n),
        "n_train_internal_missense_plp": int(n1),
        "n_train_donor_edge_flag": int(np.sum(tr_mis["donor_edge"].to_numpy() == 1)),
        "coefficients": coefs,
    }


def predict_mis_arch_prob(df_in: pd.DataFrame, coefs: dict) -> np.ndarray:
    """M1 logit plus the fitted architecture offset, as a probability."""
    w = np.array([coefs[name]["coef"] for name in OFFSET_FEATURES], dtype=np.float64)
    eta = (
        df_in["m1_logit"].to_numpy(dtype=np.float64)
        + df_in[OFFSET_FEATURES].to_numpy(dtype=np.float64) @ w
    )
    return expit(eta)


def clean_json(obj):
    """Recursively convert NaN/inf to None so the JSON output is strict."""
    if isinstance(obj, dict):
        return {k: clean_json(v) for k, v in obj.items()}
    if isinstance(obj, (list, tuple)):
        return [clean_json(v) for v in obj]
    if isinstance(obj, float) and not math.isfinite(obj):
        return None
    return obj


def size_tier(exon_len: pd.Series) -> pd.Series:
    return np.where(
        exon_len <= 180,
        "1_Compact_<=180",
        np.where(
            exon_len <= 400,
            "2_Medium_181_400",
            np.where(exon_len <= 1000, "3_Large_401_1000", "4_Giant_>1000"),
        ),
    )


def main() -> int:
    # This script refits M1 on train and scores the locked test set, so
    # it is a test-set access and is logged like the benchmark itself.
    from datetime import datetime

    with open(REPO / "results" / "test_access.log", "a") as fh:
        fh.write(
            f"{datetime.now().isoformat()} script-18 M6 evaluation "
            "(refit M1 on train, score test)\n"
        )
    out_json_path = REPO / "results/exon_architecture_nmd.json"
    prior = {}
    if out_json_path.exists():
        prior = json.loads(out_json_path.read_text())
    pop = load_clinvar_587k(REPO)
    cons = pd.read_csv(
        REPO / "data/interim/script18/variant_molecular_consequences.tsv.gz",
        sep="\t",
    )
    ptc = pd.read_csv(REPO / "data/interim/script18/ptc_features.tsv.gz", sep="\t")
    long_sev = pd.read_csv(
        REPO / "data/interim/script18/variant_transcript_sev.tsv.gz", sep="\t"
    )
    iu = pd.read_csv(REPO / "data/processed/isoform_usage.tsv.gz", sep="\t")
    am = pd.read_csv(REPO / "data/processed/variant_alphamissense.tsv.gz", sep="\t")
    cadd = pd.read_csv(REPO / "data/processed/variant_cadd.tsv.gz", sep="\t")
    spl = pd.read_csv(REPO / "data/processed/variant_spliceai.tsv.gz", sep="\t")
    fm = pd.read_csv(REPO / "data/processed/feature_matrix.tsv.gz", sep="\t")
    tiers = json.loads((REPO / "data/processed/feature_tiers.json").read_text())[
        "tiers"
    ]
    panel = pd.read_csv(REPO / "metadata/gene_panel.tsv", sep="\t")
    bench_locked = json.loads((REPO / "results/benchmark.json").read_text())
    lab24_path = REPO / "data/interim/script18/labels_2024_vs_2026.tsv.gz"

    pop["is_plp"] = (
        pop["clnsig_aggregate"]
        .isin(["Pathogenic", "Likely pathogenic", "Pathogenic/Likely pathogenic"])
        .astype(int)
    )
    pop["is_blb"] = (
        pop["clnsig_aggregate"]
        .isin(["Benign", "Likely benign", "Benign/Likely benign"])
        .astype(int)
    )
    pop["is_resolved"] = ((pop["is_plp"] == 1) | (pop["is_blb"] == 1)).astype(int)
    pop["is_vus"] = (
        pop["clnsig_aggregate"]
        .str.contains("Uncertain significance", na=False)
        .astype(int)
    )
    pop["is_conflict"] = pop["disc_conflict"].astype(int)

    cds_classes = {
        "synonymous",
        "missense",
        "inframe_indel",
        "frameshift",
        "stop_gained",
        "start_lost",
        "stop_lost",
    }

    ptc["pos"] = ptc["variant_key"].str.split(":").str[1].astype(int)
    ptc["cds_off"] = np.round(ptc["ptc_frac"] * ptc["cds_len"]).astype(int)
    gene_strand = {}
    for g, grp in ptc.groupby("gene_symbol"):
        gene_strand[g] = (
            "+"
            if (len(grp) < 2 or np.corrcoef(grp["pos"], grp["cds_off"])[0, 1] > 0)
            else "-"
        )

    var_exon_id = np.full(len(pop), -1, dtype=int)
    var_exon_len = np.zeros(len(pop), dtype=int)
    var_d_junc = np.zeros(len(pop), dtype=int)
    var_d_donor = np.zeros(len(pop), dtype=int)
    var_d_acceptor = np.zeros(len(pop), dtype=int)
    var_sym_phase = np.zeros(len(pop), dtype=int)
    exon_counter = 0

    def refine_bounds(snv_pos, snv_cls, s_raw, e_raw):
        s_pin = s_raw
        for cand_s in range(s_raw - 2, s_raw + 1):
            m_sr = (snv_pos >= cand_s - 7) & (snv_pos <= cand_s - 2)
            if m_sr.any() and (snv_cls[m_sr] == "splice_region").all():
                s_pin = cand_s
                break
        e_pin = e_raw
        for cand_e in range(e_raw, e_raw + 3):
            m_sr = (snv_pos >= cand_e + 2) & (snv_pos <= cand_e + 7)
            if m_sr.any() and (snv_cls[m_sr] == "splice_region").all():
                e_pin = cand_e
                break
        return s_pin, e_pin

    for g, gdf in pop.sort_values("pos").groupby("gene_symbol", sort=False):
        strand = gene_strand.get(g, "+")
        is_cds = gdf["mane_class_v46"].isin(cds_classes).to_numpy()
        is_intr = (
            gdf["mane_class_v46"].isin(["intronic", "splice_region", "none"]).to_numpy()
        )
        pos_arr = gdf["pos"].to_numpy()
        cls_arr = gdf["mane_class_v46"].to_numpy()
        idx_arr = gdf.index.to_numpy()
        snv_mask = (
            (gdf["ref"].str.len() == 1) & (gdf["alt"].str.len() == 1)
        ).to_numpy()
        snv_pos, snv_cls = pos_arr[snv_mask], cls_arr[snv_mask]
        cur_idxs, cur_positions = [], []

        def flush_exon(c_idxs, c_pos, ex_id, snv_pos, snv_cls, strand):
            s_pos, e_pos = int(c_pos[0]), int(c_pos[-1])
            ex_len = e_pos - s_pos + 1
            s_pin, e_pin = refine_bounds(snv_pos, snv_cls, s_pos, e_pos)
            sym = int(((e_pin - s_pin + 1) % 3) == 0)
            for idx_v, p_v in zip(c_idxs, c_pos):
                var_exon_id[idx_v] = ex_id
                var_exon_len[idx_v] = ex_len
                d_l, d_r = p_v - s_pos, e_pos - p_v
                var_d_junc[idx_v] = min(d_l, d_r)
                var_d_acceptor[idx_v] = d_l if strand == "+" else d_r
                var_d_donor[idx_v] = d_r if strand == "+" else d_l
                var_sym_phase[idx_v] = sym

        for k in range(len(gdf)):
            new_exon = (
                is_cds[k]
                and len(cur_positions) > 0
                and (pos_arr[k] - cur_positions[-1] > 90)
            )
            if new_exon:
                flush_exon(
                    cur_idxs, cur_positions, exon_counter, snv_pos, snv_cls, strand
                )
                exon_counter += 1
                cur_idxs, cur_positions = [idx_arr[k]], [pos_arr[k]]
            elif is_cds[k]:
                cur_idxs.append(idx_arr[k])
                cur_positions.append(pos_arr[k])
            elif is_intr[k]:
                if len(cur_positions) > 0:
                    flush_exon(
                        cur_idxs, cur_positions, exon_counter, snv_pos, snv_cls, strand
                    )
                    exon_counter += 1
                    cur_idxs, cur_positions = [], []
        if len(cur_positions) > 0:
            flush_exon(cur_idxs, cur_positions, exon_counter, snv_pos, snv_cls, strand)
            exon_counter += 1

    pop["exon_id"] = var_exon_id
    pop["exon_len"] = var_exon_len
    pop["d_junc"] = var_d_junc
    pop["d_donor_edge"] = var_d_donor
    pop["d_acceptor_edge"] = var_d_acceptor
    pop["is_symmetric_phase"] = var_sym_phase
    pop["strand"] = pop["gene_symbol"].map(gene_strand)

    ptc_by_gene = {g: gdf.sort_values("pos") for g, gdf in ptc.groupby("gene_symbol")}
    in_last_arr = np.zeros(len(pop), dtype=int)
    within50_arr = np.zeros(len(pop), dtype=int)
    ptc_frac_arr = np.zeros(len(pop), dtype=float)
    for g, gdf in pop.groupby("gene_symbol"):
        pdf = ptc_by_gene.get(g)
        if pdf is not None and len(pdf) > 0:
            p_ref = pdf["pos"].to_numpy()
            p_q = gdf["pos"].to_numpy()
            j = np.clip(np.searchsorted(p_ref, p_q), 0, len(p_ref) - 1)
            j_prev = np.clip(j - 1, 0, len(p_ref) - 1)
            best = np.where(
                np.abs(p_ref[j] - p_q) < np.abs(p_ref[j_prev] - p_q), j, j_prev
            )
            in_last_arr[gdf.index] = pdf["in_last_exon"].to_numpy()[best]
            within50_arr[gdf.index] = pdf["within50_of_last_junction"].to_numpy()[best]
            ptc_frac_arr[gdf.index] = pdf["ptc_frac"].to_numpy()[best]
    pop["in_last_exon"] = in_last_arr
    pop["within50_of_last_junction"] = within50_arr
    pop["ptc_frac_est"] = ptc_frac_arr

    iu_tx = (
        iu.groupby(["gene_symbol", "transcript_base"])
        .agg(max_usage=("usage", "max"), max_tpm=("median_tpm", "max"))
        .reset_index()
        .rename(columns={"transcript_base": "tx"})
    )
    ls = long_sev.merge(iu_tx, on=["gene_symbol", "tx"], how="left").fillna(0.0)
    vk_pos = cons[["variant_key", "mane_class"]].copy()
    vk_pos["pos"] = vk_pos["variant_key"].str.split(":").str[1].astype(int)
    ls = ls.merge(vk_pos, on="variant_key")
    alt_cds = ls[(ls.is_mane == 0) & (ls.mane_class.isin(cds_classes))].copy()
    alt_cds["is_spared"] = (alt_cds["sev"] == 0).astype(int)
    alt_cds["is_expr_spared"] = (
        (alt_cds["sev"] == 0)
        & (alt_cds["max_usage"] >= BACKUP_USAGE_THRESHOLD)
        & (alt_cds["max_tpm"] >= 1.0)
    ).astype(int)
    alt_cds["spared_max_u"] = np.where(alt_cds["sev"] == 0, alt_cds["max_usage"], 0.0)

    var_sparing_cds = (
        alt_cds.groupby(["gene_symbol", "variant_key", "pos"])
        .agg(
            any_alt_skips=("is_spared", "max"),
            expr_backup_skips=("is_expr_spared", "max"),
            spared_max_usage=("spared_max_u", "max"),
        )
        .reset_index()
    )
    sparing_by_gene = {
        g: gdf.sort_values("pos") for g, gdf in var_sparing_cds.groupby("gene_symbol")
    }

    any_skip_arr = np.zeros(len(pop), dtype=int)
    expr_skip_arr = np.zeros(len(pop), dtype=int)
    spared_u_arr = np.zeros(len(pop))
    for g, gdf in pop.groupby("gene_symbol"):
        sdf = sparing_by_gene.get(g)
        if sdf is not None and len(sdf) > 0:
            p_ref = sdf["pos"].to_numpy()
            p_q = gdf["pos"].to_numpy()
            j = np.clip(np.searchsorted(p_ref, p_q), 0, len(p_ref) - 1)
            j_prev = np.clip(j - 1, 0, len(p_ref) - 1)
            best = np.where(
                np.abs(p_ref[j] - p_q) < np.abs(p_ref[j_prev] - p_q), j, j_prev
            )
            any_skip_arr[gdf.index] = sdf["any_alt_skips"].to_numpy()[best]
            expr_skip_arr[gdf.index] = sdf["expr_backup_skips"].to_numpy()[best]
            spared_u_arr[gdf.index] = sdf["spared_max_usage"].to_numpy()[best]
    pop["any_alt_skips"] = any_skip_arr
    pop["expr_backup_skips"] = expr_skip_arr
    pop["spared_max_usage"] = spared_u_arr
    pop["indel_len_bp"] = np.abs(pop["alt"].str.len() - pop["ref"].str.len())
    pop["fine_class"] = np.where(
        pop["mane_class_v46"] == "inframe_indel",
        np.where(
            pop["ref"].str.len() > pop["alt"].str.len(),
            "inframe_deletion",
            "inframe_insertion",
        ),
        pop["mane_class_v46"],
    )

    internal = pop[(pop["exon_id"] >= 0) & (pop["in_last_exon"] == 0)].copy()
    internal["exon_size_tier"] = size_tier(internal["exon_len"])

    bm = (
        cons.merge(
            pop[
                [
                    "variant_key",
                    "pos",
                    "ref",
                    "alt",
                    "exon_id",
                    "exon_len",
                    "d_junc",
                    "d_donor_edge",
                    "d_acceptor_edge",
                    "is_symmetric_phase",
                    "strand",
                    "in_last_exon",
                    "within50_of_last_junction",
                    "ptc_frac_est",
                    "any_alt_skips",
                    "expr_backup_skips",
                    "spared_max_usage",
                    "fine_class",
                    "indel_len_bp",
                ]
            ],
            on="variant_key",
            how="left",
        )
        .merge(
            am[["variant_key", "am_score", "am_class", "am_missing"]],
            on="variant_key",
            how="left",
        )
        .merge(
            cadd[["variant_key", "cadd_phred", "cadd_raw"]],
            on="variant_key",
            how="left",
        )
        .merge(
            spl[["variant_key", "spliceai_max_delta", "spliceai_type"]],
            on="variant_key",
            how="left",
        )
    )
    bm_int = bm[(bm["exon_id"] >= 0) & (bm["in_last_exon"] == 0)].copy()
    bm_int["exon_size_tier"] = size_tier(bm_int["exon_len"])

    # --- 4x2x2 tensor ---
    # Exon attributes must be resolved per exon BEFORE cell assignment
    # (2026-11-09 fix). `internal` is variant-level and expr_backup_skips is
    # assigned per variant, so counting nunique(exon_id) per cell on the
    # variant-level frame double-counted 19 exons whose variants received
    # mixed expr_backup_skips values (tensor summed to 4,528 vs the true
    # 4,509 internal coding exons; all 19 excess counts landed in
    # constitutive cells). Rule: expr_backup_skips per exon = max over the
    # exon's variants (same rule as the atlas ex_summary below); n_exons and
    # n_genes are counted from this exon-level frame. Variant counts per cell
    # still come from the variant-level frame and are unaffected.
    exon_level = (
        internal.groupby("exon_id")
        .agg(
            gene_symbol=("gene_symbol", "first"),
            exon_size_tier=("exon_size_tier", "first"),
            expr_backup_skips=("expr_backup_skips", "max"),
            is_symmetric_phase=("is_symmetric_phase", "first"),
        )
        .reset_index()
    )
    tensor_4x2x2 = []
    for sz in SIZE_BINS:
        for skip in [0, 1]:
            for sym in [0, 1]:
                sub_p = internal[
                    (internal["exon_size_tier"] == sz)
                    & (internal["expr_backup_skips"] == skip)
                    & (internal["is_symmetric_phase"] == sym)
                ]
                sub_e = exon_level[
                    (exon_level["exon_size_tier"] == sz)
                    & (exon_level["expr_backup_skips"] == skip)
                    & (exon_level["is_symmetric_phase"] == sym)
                ]
                sub_b = bm_int[
                    (bm_int["exon_size_tier"] == sz)
                    & (bm_int["expr_backup_skips"] == skip)
                    & (bm_int["is_symmetric_phase"] == sym)
                ]
                row_d = {
                    "size_bin": sz,
                    "backup_skip": skip,
                    "symmetric_phase": sym,
                    "n_exons": int(sub_e["exon_id"].nunique()),
                    "n_genes": int(sub_e["gene_symbol"].nunique()),
                }
                for vc in [
                    "missense",
                    "inframe_deletion",
                    "inframe_insertion",
                    "synonymous",
                    "frameshift",
                    "stop_gained",
                ]:
                    vp = sub_p[sub_p["fine_class"] == vc]
                    vp_res = vp[vp["is_resolved"] == 1]
                    row_d[f"{vc}_n"] = len(vp)
                    row_d[f"{vc}_plp"] = int(vp["is_plp"].sum())
                    row_d[f"{vc}_blb"] = int(vp["is_blb"].sum())
                    row_d[f"{vc}_res_rate"] = (
                        float(vp_res["is_plp"].mean()) if len(vp_res) > 0 else None
                    )
                bm_mis = sub_b[
                    (sub_b["fine_class"] == "missense") & (sub_b["am_missing"] == 0)
                ]
                row_d["bm_mis_n"] = len(bm_mis)
                row_d["bm_mis_p_rate"] = (
                    float(bm_mis["label"].mean()) if len(bm_mis) > 0 else None
                )
                row_d["bm_am_median"] = (
                    float(bm_mis["am_score"].median()) if len(bm_mis) > 0 else None
                )
                tensor_4x2x2.append(row_d)

    # --- Bounding analysis: true pathogenicity rates once unresolved
    # missense variants are counted. For each cell the
    # unresolved variants are assigned 0% pathogenic (min) and 100%
    # pathogenic (max); no modelling assumption is made. Cell codes:
    # digit = size bin (1 Compact, 2 Medium, 3 Large, 4 Giant);
    # A = constitutive/symmetric-phase-asymmetric (skip 0, sym 0),
    # B = skip 0, sym 1; C = skip 1, sym 0; D = skip 1, sym 1.
    _cell_letters = {(0, 0): "A", (0, 1): "B", (1, 0): "C", (1, 1): "D"}

    def _rate_bounds(missense_n: int, plp: int, blb: int) -> dict:
        unresolved = missense_n - plp - blb
        return {
            "missense_n": int(missense_n),
            "plp": int(plp),
            "blb": int(blb),
            "unresolved_n": int(unresolved),
            "true_rate_min": float(plp / missense_n) if missense_n else None,
            "true_rate_max": (
                float((plp + unresolved) / missense_n) if missense_n else None
            ),
        }

    _bounding_cells = []
    for _row in tensor_4x2x2:
        _bounding_cells.append(
            {
                "cell": f"{_row['size_bin'][0]}{_cell_letters[(_row['backup_skip'], _row['symmetric_phase'])]}",
                "size_bin": _row["size_bin"],
                "backup_skip": _row["backup_skip"],
                "symmetric_phase": _row["symmetric_phase"],
                **_rate_bounds(
                    _row["missense_n"], _row["missense_plp"], _row["missense_blb"]
                ),
            }
        )

    def _pooled_bounds(rows) -> dict:
        return _rate_bounds(
            sum(r["missense_n"] for r in rows),
            sum(r["missense_plp"] for r in rows),
            sum(r["missense_blb"] for r in rows),
        )

    vus_bounding_analysis = {
        "method": "Per-cell bounds on the true missense pathogenicity rate from the 4x2x2 tensor counts: true_rate_min = plp / missense_n (all unresolved variants benign), true_rate_max = (plp + unresolved) / missense_n (all unresolved variants pathogenic), where unresolved = missense_n - plp - blb. Resolved-only rates describe the resolved subset only; these bounds make no assumption about the unresolved variants.",
        "cells": _bounding_cells,
        "pooled": {
            "compact_constitutive": _pooled_bounds(
                [
                    r
                    for r in tensor_4x2x2
                    if r["size_bin"].startswith("1_") and r["backup_skip"] == 0
                ]
            ),
            "giant_backup_skipped": _pooled_bounds(
                [
                    r
                    for r in tensor_4x2x2
                    if r["size_bin"].startswith("4_") and r["backup_skip"] == 1
                ]
            ),
        },
    }

    # --- Curation-independent biophysical and density controls ---
    ex_summary = (
        internal.groupby("exon_id")
        .agg(
            gene_symbol=("gene_symbol", "first"),
            exon_len=("exon_len", "first"),
            is_symmetric_phase=("is_symmetric_phase", "first"),
            expr_backup_skips=("expr_backup_skips", "max"),
            any_alt_skips=("any_alt_skips", "max"),
            spared_max_usage=("spared_max_usage", "max"),
            exon_size_tier=("exon_size_tier", "first"),
        )
        .reset_index()
    )
    curation_free_table = []
    for sz in SIZE_BINS:
        for skip in [0, 1]:
            ex_sub = ex_summary[
                (ex_summary["exon_size_tier"] == sz)
                & (ex_summary["expr_backup_skips"] == skip)
            ]
            total_kb = ex_sub["exon_len"].sum() / 1000.0
            p_mis = internal[
                (internal["exon_size_tier"] == sz)
                & (internal["expr_backup_skips"] == skip)
                & (internal["fine_class"] == "missense")
            ]
            p_ptv = internal[
                (internal["exon_size_tier"] == sz)
                & (internal["expr_backup_skips"] == skip)
                & (internal["fine_class"].isin(["frameshift", "stop_gained"]))
            ]
            p_syn = internal[
                (internal["exon_size_tier"] == sz)
                & (internal["expr_backup_skips"] == skip)
                & (internal["fine_class"] == "synonymous")
            ]
            b_mis = bm_int[
                (bm_int["exon_size_tier"] == sz)
                & (bm_int["expr_backup_skips"] == skip)
                & (bm_int["fine_class"] == "missense")
            ]
            b_am = b_mis[b_mis["am_missing"] == 0]
            b_cadd = b_mis[b_mis["cadd_phred"] > 0]
            blb_mis_kb = p_mis["is_blb"].sum() / max(total_kb, 1e-6)
            plp_mis_kb = p_mis["is_plp"].sum() / max(total_kb, 1e-6)
            vus_mis_kb = p_mis["is_vus"].sum() / max(total_kb, 1e-6)
            plp_ptv_kb = p_ptv["is_plp"].sum() / max(total_kb, 1e-6)
            blb_ptv_kb = p_ptv["is_blb"].sum() / max(total_kb, 1e-6)
            blb_syn_kb = p_syn["is_blb"].sum() / max(total_kb, 1e-6)
            curation_free_table.append(
                {
                    "size_bin": sz,
                    "backup_skip": skip,
                    "n_exons": len(ex_sub),
                    "total_kb": float(total_kb),
                    "am_median_all": float(b_am["am_score"].median()),
                    "am_mean_all": float(b_am["am_score"].mean()),
                    "am_frac_gt_056": float((b_am["am_score"] > 0.564).mean()),
                    "cadd_median_all": (
                        float(b_cadd["cadd_phred"].median()) if len(b_cadd) else None
                    ),
                    "cadd_mean_all": (
                        float(b_cadd["cadd_phred"].mean()) if len(b_cadd) else None
                    ),
                    "blb_mis_per_kb": float(blb_mis_kb),
                    "plp_mis_per_kb": float(plp_mis_kb),
                    "vus_mis_per_kb": float(vus_mis_kb),
                    "plp_ptv_per_kb": float(plp_ptv_kb),
                    "blb_ptv_per_kb": float(blb_ptv_kb),
                    "blb_syn_per_kb": float(blb_syn_kb),
                    "blb_to_plp_mis_ratio": float(blb_mis_kb / max(plp_mis_kb, 1e-6)),
                    "ptv_to_mis_plp_ratio": float(plp_ptv_kb / max(plp_mis_kb, 1e-6)),
                }
            )

    # --- Rare benign internal PTVs ---
    ptv_int = internal[internal["fine_class"].isin(["frameshift", "stop_gained"])]
    blb_ptv = ptv_int[ptv_int["is_blb"] == 1]
    plp_ptv = ptv_int[ptv_int["is_plp"] == 1]
    fisher_or, fisher_p = fisher_exact(
        [
            [
                int(blb_ptv["expr_backup_skips"].sum()),
                int(len(blb_ptv) - blb_ptv["expr_backup_skips"].sum()),
            ],
            [
                int(plp_ptv["expr_backup_skips"].sum()),
                int(len(plp_ptv) - plp_ptv["expr_backup_skips"].sum()),
            ],
        ]
    )
    rare_benign_ptvs = {
        "n_blb_internal_ptv": len(blb_ptv),
        "n_plp_internal_ptv": len(plp_ptv),
        "blb_backup_skip_rate": float(blb_ptv["expr_backup_skips"].mean()),
        "plp_backup_skip_rate": float(plp_ptv["expr_backup_skips"].mean()),
        "blb_symmetric_rate": float(blb_ptv["is_symmetric_phase"].mean()),
        "plp_symmetric_rate": float(plp_ptv["is_symmetric_phase"].mean()),
        "blb_near_edge_5bp_rate": float((blb_ptv["d_junc"] <= 5).mean()),
        "plp_near_edge_5bp_rate": float((plp_ptv["d_junc"] <= 5).mean()),
        "backup_skip_fisher_or_blb_vs_plp": float(fisher_or),
        "backup_skip_fisher_p": float(fisher_p),
    }

    # --- Splice boundary asymmetry: synonymous donor vs acceptor vs core ---
    syn_donor = pop[(pop["fine_class"] == "synonymous") & (pop["d_donor_edge"] <= 3)]
    syn_acc = pop[(pop["fine_class"] == "synonymous") & (pop["d_acceptor_edge"] <= 3)]
    syn_core = pop[(pop["fine_class"] == "synonymous") & (pop["d_junc"] > 30)]

    def _syn_block(sub):
        res = sub[sub["is_resolved"] == 1]
        return {
            "plp": int(sub["is_plp"].sum()),
            "blb": int(sub["is_blb"].sum()),
            "res_rate": float(res["is_plp"].mean()) if len(res) else None,
        }

    syn_blocks = {
        "donor": _syn_block(syn_donor),
        "acceptor": _syn_block(syn_acc),
        "core": _syn_block(syn_core),
    }

    # --- In-frame indels and NMD zones ---
    indels_int = internal[
        internal["fine_class"].isin(["inframe_deletion", "inframe_insertion"])
    ].copy()
    indel_table = []
    for icls in ["inframe_deletion", "inframe_insertion"]:
        for sz_grp, sz_cond in [
            ("Compact/Med (<=400bp)", indels_int["exon_len"] <= 400),
            ("Large/Giant (>400bp)", indels_int["exon_len"] > 400),
        ]:
            for skip in [0, 1]:
                sub_id = indels_int[
                    (indels_int["fine_class"] == icls)
                    & sz_cond
                    & (indels_int["expr_backup_skips"] == skip)
                ]
                plp_c, blb_c = int(sub_id["is_plp"].sum()), int(sub_id["is_blb"].sum())
                indel_table.append(
                    {
                        "class": icls,
                        "size_group": sz_grp,
                        "backup_skip": skip,
                        "n": len(sub_id),
                        "plp": plp_c,
                        "blb": blb_c,
                        "res_plp": float(plp_c / max(plp_c + blb_c, 1)),
                    }
                )
    ptv_all = pop[pop["fine_class"].isin(["frameshift", "stop_gained"])].copy()
    nmd_zone_table = []
    for zname, cond in [
        (
            "Zone 1: Canonical Internal NMD-Trigger (>50 nt upstream of final EJC)",
            (ptv_all["in_last_exon"] == 0)
            & (ptv_all["within50_of_last_junction"] == 0),
        ),
        (
            "Zone 2: Penultimate-Exon 50-nt EJC Escape Boundary (1-50 nt upstream)",
            (ptv_all["in_last_exon"] == 0)
            & (ptv_all["within50_of_last_junction"] == 1),
        ),
        (
            "Zone 3A: Last-Exon Proximal Truncation (in_last_exon=1, ptc_frac < 0.90)",
            (ptv_all["in_last_exon"] == 1) & (ptv_all["ptc_frac_est"] < 0.90),
        ),
        (
            "Zone 3B: Last-Exon Extreme C-Terminal Tail (in_last_exon=1, ptc_frac >= 0.90)",
            (ptv_all["in_last_exon"] == 1) & (ptv_all["ptc_frac_est"] >= 0.90),
        ),
    ]:
        sub_z = ptv_all[cond]
        plp_c, blb_c, vus_c = (
            int(sub_z["is_plp"].sum()),
            int(sub_z["is_blb"].sum()),
            int(sub_z["is_vus"].sum() + sub_z["is_conflict"].sum()),
        )
        nmd_zone_table.append(
            {
                "zone": zname,
                "n": len(sub_z),
                "plp": plp_c,
                "blb": blb_c,
                "vus_conf": vus_c,
                "raw_plp_rate": float(plp_c / max(len(sub_z), 1)),
                "vus_conf_rate": float(vus_c / max(len(sub_z), 1)),
                "resolved_plp_rate": float(plp_c / max(plp_c + blb_c, 1)),
            }
        )

    # --- Within-gene CMH tests on internal resolved missense ---
    mis_int = internal[internal["fine_class"] == "missense"].copy()
    mis_int["is_giant_large"] = (mis_int["exon_len"] > 400).astype(int)
    mis_int_res = mis_int[mis_int["is_resolved"] == 1]
    cmh_giant = run_cmh_stratified(mis_int_res, "is_giant_large")
    cmh_skip = run_cmh_stratified(mis_int_res, "expr_backup_skips")
    cmh_tables = []
    for strat, res in [
        ("giant_large_exon_gt400bp", cmh_giant),
        ("expr_backup_skip", cmh_skip),
    ]:
        for row in res["per_gene_tables"]:
            cmh_tables.append({"stratification": strat, **row})
    cmh_path = REPO / "results/cmh_per_gene_tables.tsv"
    pd.DataFrame(cmh_tables).to_csv(cmh_path, sep="\t", index=False)
    cmh_out = {
        "population": "internal resolved missense variants (exon_id >= 0, in_last_exon == 0, is_resolved == 1)",
        "n_variants": len(mis_int_res),
        "n_genes": int(mis_int_res["gene_symbol"].nunique()),
        "note": "P-values REMOVED 2026-11-09 (were 4.29e-187 for giant_large_exon_gt400bp and 3.53e-42 for expr_backup_skip). The CMH z/p treats 27,965 variants as independent within genes; the effective sample size is genes and exons (52 and 42 informative genes), not variants, so the p-values overstate the evidence by orders of magnitude. Gene-cluster bootstrap 95% CIs (2,000 replicates resampling genes from the per-gene tables in results/cmh_per_gene_tables.tsv, seed 0) were added 2026-11-09 and ship as gene_cluster_bootstrap in each block below; z statistics are withheld for the same reason as p. Odds ratios and informative-gene counts are retained as descriptive values only. Per-gene 2x2 tables: results/cmh_per_gene_tables.tsv",
        "giant_large_exon_gt400bp": cmh_shipped(
            cmh_giant,
            gene_cluster_bootstrap(cmh_giant["per_gene_tables"], cmh_giant["cmh_or"]),
        ),
        "expr_backup_skip": cmh_shipped(
            cmh_skip,
            gene_cluster_bootstrap(cmh_skip["per_gene_tables"], cmh_skip["cmh_or"]),
        ),
    }
    print(
        f"[18] Reconstructed {exon_counter:,d} MANE Select exons ({internal['exon_id'].nunique():,d} internal coding exons)."
    )
    print(
        f"[18] Locked benchmark M1->M4 step AUPRC: {bench_locked['tiers']['M1']['test']['auprc']:.4f} -> {bench_locked['tiers']['M4']['test']['auprc']:.4f} (+{bench_locked['tiers']['M4']['test']['auprc'] - bench_locked['tiers']['M1']['test']['auprc']:.4f})"
    )
    print(
        f"[18] Internal missense within-gene CMH: exon >400bp OR={cmh_giant['cmh_or']:.4f}, backup-skip OR={cmh_skip['cmh_or']:.4f} (p-values withheld: variant-level independence assumption)"
    )

    # --- Exon-skipping candidate atlas (computational nominations only) ---
    mis_by_ex = (
        internal[internal["fine_class"] == "missense"]
        .groupby("exon_id")
        .agg(
            mis_n=("variant_key", "count"),
            mis_plp=("is_plp", "sum"),
            mis_blb=("is_blb", "sum"),
            mis_vus=("is_vus", "sum"),
        )
    )
    ptv_by_ex = (
        internal[internal["fine_class"].isin(["frameshift", "stop_gained"])]
        .groupby("exon_id")
        .agg(
            ptv_n=("variant_key", "count"),
            ptv_plp=("is_plp", "sum"),
            ptv_blb=("is_blb", "sum"),
            ptv_vus=("is_vus", "sum"),
            stop_plp=(
                "fine_class",
                lambda s: int(
                    (
                        (s == "stop_gained") & (internal.loc[s.index, "is_plp"] == 1)
                    ).sum()
                ),
            ),
            fs_plp=(
                "fine_class",
                lambda s: int(
                    ((s == "frameshift") & (internal.loc[s.index, "is_plp"] == 1)).sum()
                ),
            ),
        )
    )
    bm_am_by_ex = (
        bm_int[(bm_int["fine_class"] == "missense") & (bm_int["am_missing"] == 0)]
        .groupby("exon_id")
        .agg(
            am_median=("am_score", "median"),
            am_mean=("am_score", "mean"),
            cadd_median=("cadd_phred", "median"),
        )
    )
    atlas = (
        ex_summary.merge(mis_by_ex, on="exon_id", how="left")
        .merge(ptv_by_ex, on="exon_id", how="left")
        .merge(bm_am_by_ex, on="exon_id", how="left")
        .merge(panel[["gene_symbol", "disease_category"]], on="gene_symbol", how="left")
    )
    for c in [
        "mis_n",
        "mis_plp",
        "mis_blb",
        "mis_vus",
        "ptv_n",
        "ptv_plp",
        "ptv_blb",
        "ptv_vus",
        "stop_plp",
        "fs_plp",
    ]:
        atlas[c] = atlas[c].fillna(0).astype(int)
    atlas["mis_resolved_n"] = atlas["mis_plp"] + atlas["mis_blb"]
    atlas["mis_resolved_plp_rate"] = np.where(
        atlas["mis_resolved_n"] > 0, atlas["mis_plp"] / atlas["mis_resolved_n"], np.nan
    )
    atlas["mis_raw_plp_rate"] = np.where(
        atlas["mis_n"] > 0, atlas["mis_plp"] / atlas["mis_n"], 0.0
    )
    atlas["ptv_resolved_plp_rate"] = np.where(
        (atlas["ptv_plp"] + atlas["ptv_blb"]) > 0,
        atlas["ptv_plp"] / (atlas["ptv_plp"] + atlas["ptv_blb"]),
        np.nan,
    )
    eb_mis_plp = (atlas["mis_plp"] + 1.0) / (atlas["mis_resolved_n"] + 2.0)
    s_disp_clinvar = 1.0 - eb_mis_plp
    s_disp_am = np.where(
        atlas["am_median"].notna(), 1.0 - atlas["am_median"], s_disp_clinvar
    )
    atlas["S_structural_dispensability"] = np.clip(
        0.60 * s_disp_clinvar + 0.40 * s_disp_am, 0.0, 1.0
    )
    atlas["S_frame_symmetry"] = np.where(
        atlas["is_symmetric_phase"] == 1,
        1.0,
        np.where(atlas["expr_backup_skips"] == 1, 0.65, 0.25),
    )
    atlas["S_skip_competence"] = np.where(
        atlas["expr_backup_skips"] == 1,
        0.85 + 0.15 * np.clip(atlas["spared_max_usage"], 0, 1),
        np.where((atlas["any_alt_skips"] == 1) | (atlas["exon_len"] > 400), 0.60, 0.25),
    )
    atlas["Exon_Nomination_Score"] = np.round(
        atlas["S_structural_dispensability"]
        * atlas["S_frame_symmetry"]
        * atlas["S_skip_competence"],
        4,
    )
    atlas["aso_tier"] = np.where(
        (atlas["Exon_Nomination_Score"] >= 0.55)
        & (atlas["mis_raw_plp_rate"] < 0.10)
        & (atlas["ptv_plp"] >= 1),
        "Tier 1: in-frame, backup-skipped, missense-tolerant (computational nomination only)",
        np.where(
            (atlas["Exon_Nomination_Score"] >= 0.40)
            & (atlas["mis_raw_plp_rate"] < 0.10)
            & (atlas["ptv_plp"] >= 1),
            "Tier 2: missense-tolerant, structurally dispensable (computational nomination only)",
            "Tier 3: Structurally Constrained / Non-Skippable Exon",
        ),
    )
    # Deterministic tie-break (2026-11-09): exon_id as the final sort key, so
    # rows tied on (Exon_Nomination_Score, ptv_plp) order identically
    # regardless of input row order (pandas sort is stable, so without this
    # key the output order depended on the input frame's row order).
    atlas_sorted = atlas.sort_values(
        ["Exon_Nomination_Score", "ptv_plp", "exon_id"], ascending=[False, False, True]
    ).reset_index(drop=True)
    t1_exons = atlas_sorted[atlas_sorted["aso_tier"].str.startswith("Tier 1")]
    t12_exons = atlas_sorted[
        atlas_sorted["aso_tier"].str.startswith(("Tier 1", "Tier 2"))
    ]
    out_atlas_cols = [
        "exon_id",
        "gene_symbol",
        "disease_category",
        "exon_len",
        "exon_size_tier",
        "is_symmetric_phase",
        "expr_backup_skips",
        "any_alt_skips",
        "spared_max_usage",
        "mis_n",
        "mis_plp",
        "mis_blb",
        "mis_vus",
        "mis_raw_plp_rate",
        "mis_resolved_plp_rate",
        "am_median",
        "cadd_median",
        "ptv_n",
        "ptv_plp",
        "stop_plp",
        "fs_plp",
        "ptv_blb",
        "ptv_resolved_plp_rate",
        "S_structural_dispensability",
        "S_frame_symmetry",
        "S_skip_competence",
        "Exon_Nomination_Score",
        "aso_tier",
    ]
    atlas_sorted[out_atlas_cols].to_csv(
        REPO / "results/exon_skipping_nmd_candidates.tsv", sep="\t", index=False
    )
    print(
        f"[18] Tier 1 nominated exons: {len(t1_exons):,d} across {t1_exons['gene_symbol'].nunique()} genes, {int(t1_exons['ptv_plp'].sum()):,d} ClinVar PTV records (not patients)."
    )
    print(
        f"[18] Tier 1+2 exons: {len(t12_exons):,d} across {t12_exons['gene_symbol'].nunique()} genes."
    )
    top_aso_by_burden = t12_exons.sort_values(
        ["ptv_plp", "exon_id"], ascending=[False, True]
    ).head(15)

    # --- Tiered-model comparison with fitted donor-edge offsets ---
    bm_full = fm.merge(
        bm[
            [
                "variant_key",
                "split",
                "mane_class",
                "mane_sev",
                "alt_max_class",
                "alt_max_sev",
                "exon_id",
                "exon_len",
                "d_junc",
                "d_donor_edge",
                "d_acceptor_edge",
                "is_symmetric_phase",
                "in_last_exon",
                "within50_of_last_junction",
                "ptc_frac_est",
                "expr_backup_skips",
                "spared_max_usage",
                "fine_class",
                "indel_len_bp",
                "spliceai_max_delta",
                "spliceai_type",
            ]
        ],
        on="variant_key",
        how="left",
    )
    arch_cols = [
        "exon_len",
        "d_junc",
        "d_donor_edge",
        "d_acceptor_edge",
        "is_symmetric_phase",
        "in_last_exon",
        "within50_of_last_junction",
        "ptc_frac_est",
        "expr_backup_skips",
        "spared_max_usage",
        "indel_len_bp",
        "spliceai_max_delta",
    ]
    bm_full[arch_cols] = bm_full[arch_cols].fillna(0.0)
    bm_full["is_mis"] = (bm_full["mane_class"] == "missense").astype(float)
    bm_full["is_ptv"] = (
        bm_full["mane_class"]
        .isin(["frameshift", "stop_gained", "splice_site", "start_lost", "stop_lost"])
        .astype(float)
    )
    bm_full["is_inframe_del"] = (bm_full["fine_class"] == "inframe_deletion").astype(
        float
    )
    bm_full["is_inframe_ins"] = (bm_full["fine_class"] == "inframe_insertion").astype(
        float
    )
    bm_full["is_syn"] = (bm_full["mane_class"] == "synonymous").astype(float)
    bm_full["is_splice_reg"] = (bm_full["mane_class"] == "splice_region").astype(float)
    bm_full["is_intronic"] = (
        bm_full["mane_class"].isin(["intronic", "utr", "none"]).astype(float)
    )
    bm_full["log_exon_len"] = np.log1p(np.maximum(bm_full["exon_len"], 0.0))
    bm_full["mis_x_am_score"] = bm_full["is_mis"] * np.where(
        bm_full["am_missing"] == 0, bm_full["am_score"] - 0.50, 0.0
    )
    bm_full["mis_x_cadd"] = bm_full["is_mis"] * np.where(
        bm_full["cadd_missing"] == 0, (bm_full["cadd_phred"] - 23.0) / 10.0, 0.0
    )
    bm_full["mis_x_log_exon_len"] = bm_full["is_mis"] * (bm_full["log_exon_len"] - 5.2)
    bm_full["mis_x_giant_exon"] = bm_full["is_mis"] * (
        bm_full["exon_len"] > 400
    ).astype(float)
    bm_full["mis_x_backup_skip"] = bm_full["is_mis"] * bm_full["expr_backup_skips"]
    bm_full["mis_x_giant_x_skip_x_sym"] = (
        bm_full["is_mis"]
        * (bm_full["exon_len"] > 400).astype(float)
        * bm_full["expr_backup_skips"]
        * bm_full["is_symmetric_phase"]
    )
    bm_full["ptv_x_internal_nmd_sink"] = (
        bm_full["is_ptv"]
        * (1.0 - bm_full["in_last_exon"])
        * (1.0 - bm_full["within50_of_last_junction"])
    )
    bm_full["ptv_x_distal_last_exon"] = (
        bm_full["is_ptv"]
        * bm_full["in_last_exon"]
        * (bm_full["ptc_frac_est"] >= 0.90).astype(float)
    )
    bm_full["am_benign_x_donor_edge_3bp"] = (
        (bm_full["am_missing"] == 0)
        & (bm_full["am_score"] < 0.34)
        & (bm_full["d_donor_edge"] <= 3)
        & (bm_full["exon_id"] >= 0)
    ).astype(float)
    bm_full["am_benign_x_edge_3bp"] = (
        (bm_full["am_missing"] == 0)
        & (bm_full["am_score"] < 0.34)
        & (bm_full["d_junc"] <= 3)
        & (bm_full["exon_id"] >= 0)
    ).astype(float)
    bm_full["spliceai_gain_intronic"] = (
        bm_full["is_intronic"]
        * bm_full["spliceai_max_delta"]
        * bm_full["spliceai_type"].isin(["AG", "DG"]).astype(float)
    )
    bm_full["spliceai_loss_junc"] = (
        (bm_full["is_splice_reg"] + bm_full["is_ptv"])
        * bm_full["spliceai_max_delta"]
        * bm_full["spliceai_type"].isin(["AL", "DL"]).astype(float)
    )
    bm_full["regime_I_dual_coding"] = (
        (bm_full["mane_sev"] == 2) & (bm_full["alt_max_sev"] >= 4)
    ).astype(float)
    bm_full["regime_II_alt_mis"] = (
        (bm_full["mane_sev"] <= 1) & (bm_full["alt_max_sev"] == 4)
    ).astype(float)
    bm_full["regime_III_alt_lof"] = (
        (bm_full["mane_sev"] <= 1) & (bm_full["alt_max_sev"] >= 5)
    ).astype(float)
    # Offset-fit feature columns, kept separate from the model features.
    bm_full["log_exon_len_c"] = bm_full["log_exon_len"] - 5.0
    bm_full["donor_edge_x_spliceai_ge05"] = bm_full["am_benign_x_donor_edge_3bp"] * (
        bm_full["spliceai_max_delta"] >= 0.05
    ).astype(float)
    bm_full["donor_edge"] = bm_full["am_benign_x_donor_edge_3bp"]
    mech_v3_cols = tiers["M1"] + [
        "spliceai_max_delta",
        "mane_sev",
        "alt_max_sev",
        "is_mis",
        "is_ptv",
        "is_inframe_del",
        "is_inframe_ins",
        "is_syn",
        "is_splice_reg",
        "is_intronic",
        "mis_x_am_score",
        "mis_x_cadd",
        "mis_x_log_exon_len",
        "mis_x_giant_exon",
        "mis_x_backup_skip",
        "mis_x_giant_x_skip_x_sym",
        "ptv_x_internal_nmd_sink",
        "ptv_x_distal_last_exon",
        "am_benign_x_donor_edge_3bp",
        "am_benign_x_edge_3bp",
        "spliceai_gain_intronic",
        "spliceai_loss_junc",
        "regime_I_dual_coding",
        "regime_II_alt_mis",
        "regime_III_alt_lof",
    ]
    tr = bm_full[bm_full["split"] == "train"].reset_index(drop=True)
    va = bm_full[bm_full["split"] == "validation"].reset_index(drop=True)
    te = bm_full[bm_full["split"] == "test"].reset_index(drop=True)
    m1_clf = BalancedLogReg(C=1.0).fit(tr[tiers["M1"]], tr["label"])
    for split_df in (tr, va, te):
        p_m1 = np.clip(
            m1_clf.predict_proba(split_df[tiers["M1"]])[:, 1], 1e-6, 1.0 - 1e-6
        )
        split_df["m1_logit"] = np.log(p_m1 / (1.0 - p_m1))
    tr_int_mis = tr[
        (tr["mane_class"] == "missense")
        & (tr["in_last_exon"] == 0)
        & (tr["exon_id"] >= 0)
    ].copy()
    offset_fit = fit_arch_offsets(tr_int_mis)
    offset_coefs = offset_fit["coefficients"]
    print(
        "[18] Fitted donor-edge offsets (train internal missense only): "
        + ", ".join(f"{k}={v['coef']:+.2f}" for k, v in offset_coefs.items())
    )

    mech_global_cols = ["m1_logit"] + [c for c in mech_v3_cols if c not in tiers["M1"]]
    model_comp = {}
    test_scores = {}

    def _model_metrics(p_te, p_va, y_te):
        p_heldout = np.concatenate([p_va, p_te])
        df_heldout = pd.concat([va, te], ignore_index=True)
        ho_edge_pos = (
            (df_heldout["mane_class"] == "missense")
            & (df_heldout["am_missing"] == 0)
            & (df_heldout["am_score"] < 0.34)
            & (df_heldout["d_junc"] <= 3)
            & (df_heldout["exon_id"] >= 0)
            & (df_heldout["label"] == 1)
        )
        te_edge_pos = (
            (te["mane_class"] == "missense")
            & (te["am_missing"] == 0)
            & (te["am_score"] < 0.34)
            & (te["d_junc"] <= 3)
            & (te["exon_id"] >= 0)
            & (te["label"] == 1)
        )
        te_mis_mask = (
            (te["mane_class"] == "missense")
            & (te["in_last_exon"] == 0)
            & (te["exon_id"] >= 0)
        )
        n_ho = int(ho_edge_pos.sum())
        n_ho_rec = int((p_heldout[ho_edge_pos.to_numpy()] >= 0.50).sum()) if n_ho else 0
        n_te_edge = int(te_edge_pos.sum())
        n_te_rec = int((p_te[te_edge_pos.to_numpy()] >= 0.50).sum()) if n_te_edge else 0
        return {
            "test_auroc": roc_auc_score_np(y_te, p_te),
            "test_tie_auprc": tie_aware_auprc(y_te, p_te),
            "test_brier": float(np.mean((y_te - p_te) ** 2)),
            "test_ece": compute_ece(y_te, p_te),
            "test_internal_missense_auroc": roc_auc_score_np(
                y_te[te_mis_mask], p_te[te_mis_mask]
            ),
            "test_internal_missense_auprc": tie_aware_auprc(
                y_te[te_mis_mask], p_te[te_mis_mask]
            ),
            "donor_edge_recall_val_plus_test": {
                "recall": float(n_ho_rec / n_ho) if n_ho else None,
                "n_pathogenic_edge_variants": n_ho,
                "n_recalled_at_05": n_ho_rec,
                "wilson_ci95": wilson_ci(n_ho_rec, n_ho),
                "note": "validation + test genes pooled, as in the earlier analysis",
            },
            "donor_edge_recall_test_only": {
                "recall": float(n_te_rec / n_te_edge) if n_te_edge else None,
                "n_pathogenic_edge_variants": n_te_edge,
                "n_recalled_at_05": n_te_rec,
                "wilson_ci95": wilson_ci(n_te_rec, n_te_edge),
                "note": "test genes only; undefined when the test set has zero pathogenic edge variants",
            },
        }

    for mname, mcols in [
        ("M1_Baseline", tiers["M1"]),
        ("M3_CellContext", tiers["M3"]),
        ("M4_IsoformFlags", tiers["M4"]),
        ("IsoformVUS_Mechanistic_v3", mech_global_cols),
    ]:
        clf = BalancedLogReg(C=1.0).fit(tr[mcols], tr["label"])
        p_te = clf.predict_proba(te[mcols])[:, 1]
        p_va = clf.predict_proba(va[mcols])[:, 1]
        y_te = te["label"].to_numpy()
        if "Mechanistic" in mname:
            # Same fitted model, without the missense architecture replacement:
            # isolates what the offset itself contributes.
            model_comp[mname + "_no_offset"] = _model_metrics(p_te, p_va, y_te)
            test_scores[mname + "_no_offset"] = p_te
            te_int_mis_mask = (
                (te["mane_class"] == "missense") & (te["exon_id"] >= 0)
            ).to_numpy()
            p_te = np.where(
                te_int_mis_mask, predict_mis_arch_prob(te, offset_coefs), p_te
            )
            va_int_mis_mask = (
                (va["mane_class"] == "missense") & (va["exon_id"] >= 0)
            ).to_numpy()
            p_va = np.where(
                va_int_mis_mask, predict_mis_arch_prob(va, offset_coefs), p_va
            )
        test_scores[mname] = p_te
        model_comp[mname] = _model_metrics(p_te, p_va, y_te)
        mc = model_comp[mname]
        print(
            f"  {mname:28s}: Test AUROC={mc['test_auroc']:.4f} | Tie-AUPRC={mc['test_tie_auprc']:.4f} | Brier={mc['test_brier']:.4f} | IntMis AUROC={mc['test_internal_missense_auroc']:.4f} | Edge recall val+test={mc['donor_edge_recall_val_plus_test']['recall']} (n={mc['donor_edge_recall_val_plus_test']['n_pathogenic_edge_variants']}), test-only n={mc['donor_edge_recall_test_only']['n_pathogenic_edge_variants']}"
        )

    # --- Consequence-class-only baseline (the fair comparator for M6) ---
    from sklearn.linear_model import LogisticRegression
    from sklearn.metrics import average_precision_score, roc_auc_score

    cons_dummies = pd.get_dummies(bm_full["mane_class"], prefix="cons").astype(float)
    cons_cols = list(cons_dummies.columns)
    bm_full = pd.concat([bm_full, cons_dummies], axis=1)
    tr_idx = (bm_full["split"] == "train").to_numpy()
    te_idx = (bm_full["split"] == "test").to_numpy()
    cons_clf = LogisticRegression(
        C=1.0, max_iter=5000, class_weight="balanced", random_state=SEED
    )
    cons_clf.fit(
        bm_full.loc[tr_idx, cons_cols].values, bm_full.loc[tr_idx, "label"].values
    )
    p_cons_te = cons_clf.predict_proba(bm_full.loc[te_idx, cons_cols].values)[:, 1]
    test_scores["ConsequenceOnly"] = p_cons_te
    y_te = te["label"].to_numpy()
    cons_baseline = {
        "method": "one-hot mane_class logistic regression (C=1.0, class_weight=balanced), fit on 61 train genes, evaluated on 24 test genes; no exon, AlphaMissense, CADD, or SpliceAI features",
        "n_train_genes": int(tr["gene_symbol"].nunique()),
        "n_test_genes": int(te["gene_symbol"].nunique()),
        "test_auroc": float(roc_auc_score(y_te, p_cons_te)),
        "test_tie_auprc": tie_aware_auprc(y_te, p_cons_te),
        "test_step_auprc": float(average_precision_score(y_te, p_cons_te)),
        "test_brier": float(np.mean((y_te - p_cons_te) ** 2)),
        "test_internal_missense_auroc": 0.5,
        "within_missense_note": "0.5 by construction: consequence is constant within missense, so the baseline cannot rank missense variants",
    }
    print(
        f"[18] Consequence-only baseline: test AUROC={cons_baseline['test_auroc']:.4f} AUPRC(tie)={cons_baseline['test_tie_auprc']:.4f} Brier={cons_baseline['test_brier']:.4f}"
    )

    # --- Per-class and within-missense metrics (primary evaluation) ---
    within_missense = {}
    for mname in [
        "M1_Baseline",
        "M4_IsoformFlags",
        "IsoformVUS_Mechanistic_v3",
        "IsoformVUS_Mechanistic_v3_no_offset",
    ]:
        within_missense[mname] = {
            "test_internal_missense_auroc": model_comp[mname][
                "test_internal_missense_auroc"
            ],
            "test_internal_missense_auprc": model_comp[mname][
                "test_internal_missense_auprc"
            ],
        }
    within_missense["ConsequenceOnly"] = {
        "test_internal_missense_auroc": 0.5,
        "test_internal_missense_auprc": None,
    }
    per_class = {}
    te_reset = te.reset_index(drop=True)
    for mname, p_te in test_scores.items():
        per_class[mname] = {}
        for cls, grp in te_reset.groupby("mane_class"):
            y_c = grp["label"].to_numpy()
            s_c = p_te[grp.index.to_numpy()]
            n_pos, n_neg = int((y_c == 1).sum()), int((y_c == 0).sum())
            per_class[mname][cls] = {
                "n": len(y_c),
                "n_pos": n_pos,
                "n_neg": n_neg,
                "auroc": roc_auc_score_np(y_c, s_c) if n_pos and n_neg else None,
                "tie_auprc": tie_aware_auprc(y_c, s_c) if n_pos else None,
            }

    # --- Bayesian likelihood ratios: one formula, with CIs, no ACMG labels ---
    mis_int_res = internal[
        (internal["fine_class"] == "missense") & (internal["is_resolved"] == 1)
    ]
    n_mis_plp = int(mis_int_res["is_plp"].sum())
    n_mis_blb = int(mis_int_res["is_blb"].sum())
    lr_table = []
    for sz in SIZE_BINS:
        for skip in [0, 1]:
            for sym in [0, 1]:
                sub_m = mis_int_res[
                    (mis_int_res["exon_size_tier"] == sz)
                    & (mis_int_res["expr_backup_skips"] == skip)
                    & (mis_int_res["is_symmetric_phase"] == sym)
                ]
                sub_vus = internal[
                    (internal["fine_class"] == "missense")
                    & (internal["is_vus"] == 1)
                    & (internal["exon_size_tier"] == sz)
                    & (internal["expr_backup_skips"] == skip)
                    & (internal["is_symmetric_phase"] == sym)
                ]
                plp_n, blb_n = int(sub_m["is_plp"].sum()), int(sub_m["is_blb"].sum())
                lr_row = lr_with_ci(plp_n, n_mis_plp, blb_n, n_mis_blb)
                # n_exons/n_genes come from the exon-level frame (one row
                # per exon, expr_backup_skips = max over the exon's
                # variants), the same frame the 4x2x2 tensor uses
                # (2026-11-09 fix). Counting them from the variant-level
                # resolved-missense frame summed to 3,297 exons against
                # the tensor's 4,509 and disagreed per cell.
                sub_e_lr = exon_level[
                    (exon_level["exon_size_tier"] == sz)
                    & (exon_level["expr_backup_skips"] == skip)
                    & (exon_level["is_symmetric_phase"] == sym)
                ]
                lr_table.append(
                    {
                        "state": f"{sz}|Skip={skip}|Sym={sym}",
                        "size_bin": sz,
                        "backup_skip": skip,
                        "symmetric_phase": sym,
                        "n_genes": int(sub_e_lr["gene_symbol"].nunique()),
                        "n_exons": int(sub_e_lr["exon_id"].nunique()),
                        "resolved_plp": plp_n,
                        "resolved_blb": blb_n,
                        "vus_n": len(sub_vus),
                        "posterior_plp_prob": float(plp_n / max(plp_n + blb_n, 1)),
                        **lr_row,
                    }
                )
    lr_by_cell = {
        (r["size_bin"], r["backup_skip"], r["symmetric_phase"]): r for r in lr_table
    }
    # Donor-edge LR: internal missense, d_donor <= 3bp, exon <= 400bp, no
    # backup skip, among resolved internal missense.
    donor_grp = mis_int_res[
        (mis_int_res["d_donor_edge"] <= 3)
        & (mis_int_res["exon_len"] <= 400)
        & (mis_int_res["expr_backup_skips"] == 0)
    ]
    donor_lr = lr_with_ci(
        int(donor_grp["is_plp"].sum()),
        n_mis_plp,
        int(donor_grp["is_blb"].sum()),
        n_mis_blb,
    )
    donor_lr["definition"] = (
        "internal resolved missense, d_donor_edge <= 3 bp, exon_len <= 400 bp, expr_backup_skips == 0"
    )
    # Regime-II LR: MANE sev <= 1 with alt-max sev == 4, among all resolved
    # variants in the 587k table (regime II spans consequence classes).
    pop_res = pop[pop["is_resolved"] == 1]
    n_pop_plp, n_pop_blb = int(pop_res["is_plp"].sum()), int(pop_res["is_blb"].sum())
    reg2 = pop_res[(pop_res["mane_sev_v46"] <= 1) & (pop_res["alt_max_sev_v46"] == 4)]
    regime2_lr = lr_with_ci(
        int(reg2["is_plp"].sum()), n_pop_plp, int(reg2["is_blb"].sum()), n_pop_blb
    )
    regime2_lr["definition"] = (
        "resolved variants with mane_sev_v46 <= 1 and alt_max_sev_v46 == 4, among all resolved variants in the 587k table"
    )
    print(
        f"[18] 4D LR={lr_by_cell[('4_Giant_>1000', 1, 1)]['likelihood_ratio']:.4f} "
        f"[{lr_by_cell[('4_Giant_>1000', 1, 1)]['lr_ci95_lo']:.4f}, {lr_by_cell[('4_Giant_>1000', 1, 1)]['lr_ci95_hi']:.4f}] "
        f"(n_genes={lr_by_cell[('4_Giant_>1000', 1, 1)]['n_genes']}); "
        f"4C LR={lr_by_cell[('4_Giant_>1000', 1, 0)]['likelihood_ratio']:.4f} "
        f"[{lr_by_cell[('4_Giant_>1000', 1, 0)]['lr_ci95_lo']:.4f}, {lr_by_cell[('4_Giant_>1000', 1, 0)]['lr_ci95_hi']:.4f}]"
    )

    # --- Usage-threshold sweep (sensitivity): backup-skip CMH OR and the
    # giant backup-skipped cell LRs at usage thresholds 0.05 / 0.10 / 0.20.
    # Only the expr_backup_skips assignment changes with the threshold;
    # exon reconstruction, size tiers and phase are threshold-invariant,
    # and the resolved P/LP and B/LB totals (the LR denominators) do not
    # depend on the threshold. The shipped primary results above remain
    # the BACKUP_USAGE_THRESHOLD (0.10) computation.
    usage_threshold_sweep = []
    _mis_all = internal[internal["fine_class"] == "missense"]
    for _thr in (0.05, 0.10, 0.20):
        _flags = expr_flags_for_threshold(pop, alt_cds, _thr)
        _res_thr = mis_int_res[["gene_symbol", "is_plp"]].copy()
        _res_thr["expr_thr"] = _flags[mis_int_res.index.to_numpy()]
        _cmh_thr = run_cmh_stratified(_res_thr, "expr_thr")
        _flags_all = _flags[_mis_all.index.to_numpy()]
        _cells_thr = {}
        for _sym in (0, 1):
            _m = (
                (_mis_all["exon_size_tier"] == "4_Giant_>1000")
                & (_flags_all == 1)
                & (_mis_all["is_symmetric_phase"] == _sym)
            )
            _sub = _mis_all[_m.to_numpy()]
            _plp_c, _blb_c = int(_sub["is_plp"].sum()), int(_sub["is_blb"].sum())
            _lr_c = lr_with_ci(_plp_c, n_mis_plp, _blb_c, n_mis_blb)
            _cells_thr[_sym] = {
                "missense_n": len(_sub),
                "plp": _plp_c,
                "blb": _blb_c,
                "likelihood_ratio": _lr_c["likelihood_ratio"],
                "lr_ci95_lo": _lr_c["lr_ci95_lo"],
                "lr_ci95_hi": _lr_c["lr_ci95_hi"],
            }
        usage_threshold_sweep.append(
            {
                "usage_threshold": _thr,
                "backup_skip_cmh_or": _cmh_thr["cmh_or"],
                "backup_skip_n_informative_genes": _cmh_thr["n_informative_genes"],
                "cell_4c_giant_skip_asymmetric": _cells_thr[0],
                "cell_4d_giant_skip_symmetric": _cells_thr[1],
                "cell_4d_counts": {
                    "missense_n": _cells_thr[1]["missense_n"],
                    "plp": _cells_thr[1]["plp"],
                    "blb": _cells_thr[1]["blb"],
                },
            }
        )
        print(
            f"[18] Usage-threshold sweep thr={_thr:.2f}: backup-skip CMH OR={_cmh_thr['cmh_or']:.4f}, "
            f"4C LR={_cells_thr[0]['likelihood_ratio']:.4f} [{_cells_thr[0]['lr_ci95_lo']:.4f}, {_cells_thr[0]['lr_ci95_hi']:.4f}], "
            f"4D counts n={_cells_thr[1]['missense_n']} plp={_cells_thr[1]['plp']} blb={_cells_thr[1]['blb']}"
        )

    # --- VUS stratification table: per-cell LRs with CIs, no ACMG labels ---
    vus_frames = []
    vus_base_cols = [
        "variant_key",
        "gene_symbol",
        "chrom",
        "pos",
        "ref",
        "alt",
        "clnsig_aggregate",
        "mane_class_v46",
        "alt_max_class_v46",
        "exon_id",
        "exon_len",
        "is_symmetric_phase",
        "expr_backup_skips",
        "d_donor_edge",
        "d_acceptor_edge",
    ]
    for (sz, skip, sym), stratum in [
        (("4_Giant_>1000", 1, 1), "giant_backup_skipped_symmetric"),
        (("4_Giant_>1000", 1, 0), "giant_backup_skipped_asymmetric"),
        (("4_Giant_>1000", 0, 1), "giant_constitutive_symmetric"),
        (("4_Giant_>1000", 0, 0), "giant_constitutive_asymmetric"),
    ]:
        sub = internal[
            (internal["fine_class"] == "missense")
            & (internal["is_vus"] == 1)
            & (internal["exon_size_tier"] == sz)
            & (internal["expr_backup_skips"] == skip)
            & (internal["is_symmetric_phase"] == sym)
        ].copy()
        if len(sub) == 0:
            continue
        lr_row = lr_by_cell[(sz, skip, sym)]
        sub["exon_cell"] = f"{sz}|Skip={skip}|Sym={sym}"
        sub["stratum"] = stratum
        sub["likelihood_ratio"] = lr_row["likelihood_ratio"]
        sub["lr_ci95_lo"] = lr_row["lr_ci95_lo"]
        sub["lr_ci95_hi"] = lr_row["lr_ci95_hi"]
        vus_frames.append(sub)
    donor_vus = internal[
        (internal["fine_class"] == "missense")
        & (internal["is_vus"] == 1)
        & (internal["d_donor_edge"] <= 3)
        & (internal["exon_len"] <= 400)
        & (internal["expr_backup_skips"] == 0)
    ].copy()
    donor_vus["exon_cell"] = "donor_edge_0_3bp_exon_le400_skip0"
    donor_vus["stratum"] = "donor_edge"
    donor_vus["likelihood_ratio"] = donor_lr["likelihood_ratio"]
    donor_vus["lr_ci95_lo"] = donor_lr["lr_ci95_lo"]
    donor_vus["lr_ci95_hi"] = donor_lr["lr_ci95_hi"]
    vus_frames.append(donor_vus)
    reg2_vus = pop[
        (pop["mane_sev_v46"] <= 1)
        & (pop["alt_max_sev_v46"] == 4)
        & (pop["is_vus"] == 1)
    ].copy()
    reg2_vus["exon_cell"] = "regime_ii_alt_exon_missense"
    reg2_vus["stratum"] = "regime_ii_alt_exon_missense"
    reg2_vus["likelihood_ratio"] = regime2_lr["likelihood_ratio"]
    reg2_vus["lr_ci95_lo"] = regime2_lr["lr_ci95_lo"]
    reg2_vus["lr_ci95_hi"] = regime2_lr["lr_ci95_hi"]
    for c in [
        "exon_id",
        "exon_len",
        "is_symmetric_phase",
        "expr_backup_skips",
        "d_donor_edge",
        "d_acceptor_edge",
    ]:
        if c not in reg2_vus.columns:
            reg2_vus[c] = -1
    reg2_vus["mane_class_v46"] = reg2_vus["mane_class_v46"].fillna("")
    reg2_vus["alt_max_class_v46"] = reg2_vus["alt_max_class_v46"].fillna("")
    vus_frames.append(reg2_vus)
    vus_out_cols = vus_base_cols + [
        "exon_cell",
        "stratum",
        "likelihood_ratio",
        "lr_ci95_lo",
        "lr_ci95_hi",
    ]
    vus_df = (
        pd.concat(vus_frames, ignore_index=True)[vus_out_cols]
        .sort_values(["stratum", "likelihood_ratio", "gene_symbol", "variant_key"])
        .reset_index(drop=True)
    )
    vus_df.to_csv(
        REPO / "results/vus_bayesian_stratification.tsv", sep="\t", index=False
    )
    vus_counts = vus_df["stratum"].value_counts().to_dict()
    print(
        f"[18] VUS stratification rows: {len(vus_df):,d} across {vus_df['gene_symbol'].nunique()} genes: {vus_counts}"
    )
    n_vus_ci_crossing = int(
        ((vus_df["lr_ci95_lo"] <= 1.0) & (vus_df["lr_ci95_hi"] >= 1.0)).sum()
    )
    print(
        f"[18] VUS rows in strata whose LR 95% CI crosses 1.0: {n_vus_ci_crossing} of {len(vus_df):,d}"
    )
    note_lines = [
        "# VUS Bayesian stratification — research-only note",
        "",
        "File: results/vus_bayesian_stratification.tsv — generated by scripts/18_exon_architecture_nmd.py.",
        "",
        "These likelihood ratios are RESEARCH-ONLY stratification of unresolved ClinVar VUS. They are NOT",
        "clinical classifications and must NOT be converted into ACMG/AMP evidence points or used to",
        "reclassify any variant for patient care.",
        "",
        "Formula (one formula for every LR in this file):",
        "  LR(cell) = P(cell | P/LP) / P(cell | B/LB), computed on resolved variants only, no smoothing.",
        "  95% CI on the log scale, SE = sqrt(1/a - 1/N1 + 1/b - 1/N0).",
        "  A Haldane-Anscombe 0.5 correction is used only if a cell count is zero (none of the shipped",
        "  cells needed it).",
        "",
        "Cell LRs (giant/donor strata) use resolved INTERNAL MISSENSE variants as the population.",
        "The regime-II LR uses all resolved variants in the 587k ClinVar table as the population.",
        "The donor-edge LR uses resolved internal missense with d_donor_edge <= 3 bp, exon_len <= 400 bp,",
        "expr_backup_skips == 0 as the cell.",
        "",
        "Cell 4D (giant, backup-skipped, symmetric) is a SINGLE GENE: TTN, 4 exons, 3 resolved pathogenic",
        "missense variants, 94% of its missense variants unresolved. Its LR is estimated on the 6% tail",
        "curators could resolve and must not be read as a property of the cell.",
        "",
        "LRs are estimated on the same ClinVar table they are applied to (no split-sample, cross-gene,",
        "or temporal validation of the LRs themselves), and resolved variants are not exchangeable with",
        "VUS: resolvability correlates with evidence strength. Treat every LR here as a descriptive",
        "summary of the resolved subset, not as a validated classifier input.",
        "",
        "The exon LR overlaps the AlphaMissense signal (the exon effect is largely captured by",
        "AlphaMissense), so an exon LR must not be added on top of a PP3/BP4 computational-evidence",
        "item for the same variant: that would count the same signal twice.",
        "",
        f"Of the {len(vus_df):,d} rows in the stratification file, {n_vus_ci_crossing} fall in a stratum whose LR 95% CI crosses 1.0.",
    ]
    (REPO / "results/vus_bayesian_stratification_NOTE.md").write_text(
        "\n".join(note_lines) + "\n"
    )

    # --- Sensitivity checks (test 4 is saved here; earlier versions printed it only) ---
    excluded_famous = {"BRCA1", "BRCA2", "APC", "TTN", "NF1"}
    int_no_famous = internal[~internal["gene_symbol"].isin(excluded_famous)].copy()
    int_no_famous["is_giant_large"] = (int_no_famous["exon_len"] > 400).astype(int)
    mis_nf_res = int_no_famous[
        (int_no_famous["fine_class"] == "missense")
        & (int_no_famous["is_resolved"] == 1)
    ]
    c_comp_nf = mis_nf_res[
        (mis_nf_res["exon_size_tier"] == "1_Compact_<=180")
        & (mis_nf_res["expr_backup_skips"] == 0)
    ]
    c_giant_nf = mis_nf_res[
        (mis_nf_res["exon_size_tier"] == "4_Giant_>1000")
        & (mis_nf_res["expr_backup_skips"] == 1)
    ]
    cmh_giant_nf = run_cmh_stratified(mis_nf_res, "is_giant_large")
    cmh_skip_nf = run_cmh_stratified(mis_nf_res, "expr_backup_skips")
    test1 = {
        "excluded_genes": sorted(excluded_famous),
        "compact_constit_res_plp": float(c_comp_nf["is_plp"].mean()),
        "giant_skip_res_plp": float(c_giant_nf["is_plp"].mean()),
        "cmh_giant_large": cmh_shipped(cmh_giant_nf),
        "cmh_backup_skip": cmh_shipped(cmh_skip_nf),
        "cmh_p_note": CMH_P_NOTE_TEST1,
    }
    bm_int_2star = bm_full[
        (bm_full["exon_id"] >= 0)
        & (bm_full["in_last_exon"] == 0)
        & (bm_full["review_stars"] >= 2)
    ].copy()
    bm_mis_2star = bm_int_2star[bm_int_2star["mane_class"] == "missense"]
    bm_ptv_2star = bm_int_2star[
        bm_int_2star["mane_class"].isin(["frameshift", "stop_gained"])
    ]
    p_comp_2s = bm_mis_2star[
        (bm_mis_2star["exon_len"] <= 180) & (bm_mis_2star["expr_backup_skips"] == 0)
    ]["label"]
    p_giant_2s = bm_mis_2star[
        (bm_mis_2star["exon_len"] > 1000) & (bm_mis_2star["expr_backup_skips"] == 1)
    ]["label"]
    ptv_giant_2s = bm_ptv_2star[
        (bm_ptv_2star["exon_len"] > 1000) & (bm_ptv_2star["expr_backup_skips"] == 1)
    ]["label"]
    test2 = {
        "compact_constit_mis_p_rate": float(p_comp_2s.mean()),
        "compact_constit_n": len(p_comp_2s),
        "giant_skip_mis_p_rate": float(p_giant_2s.mean()),
        "giant_skip_n": len(p_giant_2s),
        "giant_skip_ptv_p_rate": float(ptv_giant_2s.mean()),
        "giant_skip_ptv_n": len(ptv_giant_2s),
    }
    snv_int_mis = internal[
        (internal["fine_class"] == "missense")
        & (internal["is_resolved"] == 1)
        & (internal["ref"].str.len() == 1)
        & (internal["alt"].str.len() == 1)
    ].copy()
    snv_int_mis["is_CT_GA_transition"] = (
        ((snv_int_mis["ref"] == "C") & (snv_int_mis["alt"] == "T"))
        | ((snv_int_mis["ref"] == "G") & (snv_int_mis["alt"] == "A"))
    ).astype(int)

    def _fold_drop(m_val):
        sub_mut = snv_int_mis[snv_int_mis["is_CT_GA_transition"] == m_val]
        r_cc = sub_mut[
            (sub_mut["exon_size_tier"] == "1_Compact_<=180")
            & (sub_mut["expr_backup_skips"] == 0)
        ]["is_plp"].mean()
        r_gs = sub_mut[
            (sub_mut["exon_size_tier"] == "4_Giant_>1000")
            & (sub_mut["expr_backup_skips"] == 1)
        ]["is_plp"].mean()
        return float(r_cc / max(r_gs, 1e-6))

    test3 = {
        "ct_ga_transitions_fold_drop": _fold_drop(1),
        "transversions_other_fold_drop": _fold_drop(0),
    }
    bm_mis_am = bm_int[
        (bm_int["fine_class"] == "missense") & (bm_int["am_missing"] == 0)
    ].copy()
    test4 = {}
    for am_label, am_lo, am_hi in [
        ("structured_core_am_gt_070", 0.70, 1.01),
        ("ambiguous_am_030_070", 0.30, 0.70),
        ("disordered_am_lt_030", -0.01, 0.30),
    ]:
        sub_am = bm_mis_am[
            (bm_mis_am["am_score"] >= am_lo) & (bm_mis_am["am_score"] < am_hi)
        ]
        r_c = sub_am[sub_am["exon_len"] <= 180]["label"]
        r_g = sub_am[sub_am["exon_len"] > 400]["label"]
        test4[am_label] = {
            "compact_le180_p_rate": float(r_c.mean()) if len(r_c) else None,
            "compact_n": len(r_c),
            "giant_large_gt400_p_rate": float(r_g.mean()) if len(r_g) else None,
            "giant_large_n": len(r_g),
        }
    ptv_int_res = internal[
        (internal["fine_class"].isin(["frameshift", "stop_gained"]))
        & (internal["is_resolved"] == 1)
    ]
    ptv_sym_skip = ptv_int_res[
        (ptv_int_res["expr_backup_skips"] == 1)
        & (ptv_int_res["is_symmetric_phase"] == 1)
    ]["is_plp"]
    ptv_asym_skip = ptv_int_res[
        (ptv_int_res["expr_backup_skips"] == 1)
        & (ptv_int_res["is_symmetric_phase"] == 0)
    ]["is_plp"]
    test5 = {
        "symmetric_backup_skip_ptv_res_plp": float(ptv_sym_skip.mean()),
        "symmetric_n": len(ptv_sym_skip),
        "asymmetric_backup_skip_ptv_res_plp": float(ptv_asym_skip.mean()),
        "asymmetric_n": len(ptv_asym_skip),
    }
    tr_genes = set(tr["gene_symbol"].unique())
    va_genes = set(va["gene_symbol"].unique())
    te_genes = set(te["gene_symbol"].unique())
    tr_exons = set(tr["exon_id"].unique()) - {-1}
    te_exons = set(te["exon_id"].unique()) - {-1}
    test6 = {
        "train_test_gene_overlap": len(tr_genes & te_genes),
        "train_val_gene_overlap": len(tr_genes & va_genes),
        "train_test_exon_overlap": len(tr_exons & te_exons),
    }
    test7 = {"status": "skipped", "reason": "labels_2024_vs_2026.tsv.gz not found"}
    if lab24_path.exists():
        lab24 = pd.read_csv(lab24_path, sep="\t")
        bm_temp = bm_full.merge(
            lab24[["variant_key", "label_2024"]], on="variant_key", how="left"
        )
        new_2026 = bm_temp[
            bm_temp["label_2024"].isna()
            & (bm_temp["mane_class"] == "missense")
            & (bm_temp["in_last_exon"] == 0)
            & (bm_temp["exon_id"] >= 0)
        ]
        new_cc = new_2026[
            (new_2026["exon_len"] <= 180) & (new_2026["expr_backup_skips"] == 0)
        ]["label"]
        new_gs = new_2026[
            (new_2026["exon_len"] > 1000) & (new_2026["expr_backup_skips"] == 1)
        ]["label"]
        test7 = {
            "status": "computed",
            "n_new_2026_internal_missense": len(new_2026),
            "compact_constit_p_rate": float(new_cc.mean()),
            "compact_constit_n": len(new_cc),
            "giant_skip_p_rate": float(new_gs.mean()),
            "giant_skip_n": len(new_gs),
        }
    sensitivity = {
        "test1_leave_top_genes_out": test1,
        "test2_high_confidence_2plus_stars": test2,
        "test3_cpg_mutational_spectrum": test3,
        "test4_am_strata_within_exon_size": test4,
        "test5_symmetric_vs_asymmetric_ptv_lethality": test5,
        "test6_zero_leakage_audit": test6,
        "test7_prospective_2024_to_2026": test7,
        "usage_threshold_sweep": usage_threshold_sweep,
    }

    # --- Assemble and write the JSON ---
    # Benchmark PTV zones (same zone definitions as nmd_zones_587k, but on
    # the 60,521-variant benchmark table with its binary label). Uses bm, not
    # bm_int: zones 3A/3B are last-exon variants, which bm_int excludes.
    bm_ptv = bm[
        (bm["exon_id"] >= 0) & (bm["fine_class"].isin(["frameshift", "stop_gained"]))
    ]
    ptc_zones = []
    for zname, cond in [
        (
            "Zone 1: Canonical NMD-Trigger (>50 nt upstream of last EJC, internal)",
            (bm_ptv["in_last_exon"] == 0) & (bm_ptv["within50_of_last_junction"] == 0),
        ),
        (
            "Zone 2: Penultimate-Exon 50-nt Boundary (1-50 nt upstream of last EJC)",
            (bm_ptv["in_last_exon"] == 0) & (bm_ptv["within50_of_last_junction"] == 1),
        ),
        (
            "Zone 3A: Last-Exon Proximal Half (in_last_exon=1, ptc_frac < 0.90)",
            (bm_ptv["in_last_exon"] == 1) & (bm_ptv["ptc_frac_est"] < 0.90),
        ),
        (
            "Zone 3B: Last-Exon Distal C-Tail (in_last_exon=1, ptc_frac >= 0.90)",
            (bm_ptv["in_last_exon"] == 1) & (bm_ptv["ptc_frac_est"] >= 0.90),
        ),
    ]:
        sub_z = bm_ptv[cond]
        ptc_zones.append(
            {
                "zone": zname,
                "bm_n": len(sub_z),
                "bm_p_rate": float(sub_z["label"].mean()) if len(sub_z) else None,
                "bm_benign_n": int((sub_z["label"] == 0).sum()),
            }
        )

    def _sym_rate(mask):
        sub = ex_summary[mask]
        n_sym = int(sub["is_symmetric_phase"].sum())
        return float(sub["is_symmetric_phase"].mean()) if len(sub) else None, [
            n_sym,
            len(sub),
        ]

    expr_sym_rate, expr_sym_n = _sym_rate(ex_summary["expr_backup_skips"] == 1)
    unexpr_sym_rate, unexpr_sym_n = _sym_rate(
        (ex_summary["any_alt_skips"] == 1) & (ex_summary["expr_backup_skips"] == 0)
    )
    const_sym_rate, const_sym_n = _sym_rate(ex_summary["any_alt_skips"] == 0)
    exon_phase_symmetry = {
        "expr_backup_skipped_sym_rate": expr_sym_rate,
        "expr_backup_skipped_sym_n": expr_sym_n,
        "unexpr_skipped_sym_rate": unexpr_sym_rate,
        "unexpr_skipped_sym_n": unexpr_sym_n,
        "constitutive_sym_rate": const_sym_rate,
        "constitutive_sym_n": const_sym_n,
        "tensor_4x2x2": tensor_4x2x2,
        "tensor_correction_note": TENSOR_CORRECTION_NOTE,
    }
    boundary = {
        "synonymous_587k_donor_vs_acceptor": {
            "donor_plp": syn_blocks["donor"]["plp"],
            "donor_blb": syn_blocks["donor"]["blb"],
            "donor_res_rate": syn_blocks["donor"]["res_rate"],
            "acceptor_plp": syn_blocks["acceptor"]["plp"],
            "acceptor_blb": syn_blocks["acceptor"]["blb"],
            "acceptor_res_rate": syn_blocks["acceptor"]["res_rate"],
            "core_plp": syn_blocks["core"]["plp"],
            "core_blb": syn_blocks["core"]["blb"],
            "core_res_rate": syn_blocks["core"]["res_rate"],
        },
        "withdrawn_blocks_moved": "distance_profile and donor_vs_acceptor_edge_0_3bp are no longer shipped in this file: no code in the repo computes them, and on 2026-11-09 they were moved to results/withdrawn_carried_over_values.json (withdrawn; do not cite). Only synonymous_587k_donor_vs_acceptor above is computed by this script.",
    }
    # Withdrawn-values file for the carried-over blocks: written once, from the
    # previous results file, only if it does not already exist.
    withdrawn_path = REPO / "results/withdrawn_carried_over_values.json"
    if not withdrawn_path.exists():
        prior_boundary_full = prior.get("splice_junction_boundary_asymmetry", {})
        withdrawn = {
            "status": "withdrawn 2026-11-09; no generating code; do not cite",
            "source": "Copied verbatim from results/exon_architecture_nmd.json as committed at 5152848, the last version that carried these blocks in the main results file. No code in the repo computes causal_mediation or the splice_junction_boundary_asymmetry distance_profile / donor_vs_acceptor_edge_0_3bp blocks; they were carried over from a previous results file. Retained here only as an audit trail.",
            "causal_mediation": prior.get("causal_mediation"),
            "splice_junction_boundary_asymmetry": {
                "distance_profile": prior_boundary_full.get("distance_profile"),
                "donor_vs_acceptor_edge_0_3bp": prior_boundary_full.get(
                    "donor_vs_acceptor_edge_0_3bp"
                ),
            },
        }
        withdrawn_path.write_text(json.dumps(clean_json(withdrawn), indent=2))
        print(f"[18] wrote {withdrawn_path} (carried-over blocks, withdrawn)")
    out = {
        "provenance": {
            "script": "scripts/18_exon_architecture_nmd.py",
            "generated_on": date.today().isoformat(),
            "seed": SEED,
            "note": "All sections in this file are recomputed by this script. "
            "Donor-edge offsets are fit on train genes only (see tiered_model_and_bayesian_lr/fitted_arch_offsets). "
            "Feature matrix carries real CADD values for validation and test variants (19,014 of 60,521 scored); every train row still carries the -1.0 sentinel, so the fitted CADD coefficient was never fit against varying CADD values.",
            "withdrawn_carried_over_values": "causal_mediation and splice_junction_boundary_asymmetry.distance_profile / donor_vs_acceptor_edge_0_3bp are NOT in this file: no generating code exists for them, and on 2026-11-09 they were moved to results/withdrawn_carried_over_values.json (status: withdrawn 2026-11-09; no generating code; do not cite).",
        },
        "dataset_summary": {
            "n_clinvar_total": len(pop),
            "n_benchmark_total": len(fm),
            "n_reconstructed_mane_exons": int(exon_counter),
            "n_internal_coding_exons": int(internal["exon_id"].nunique()),
            "n_panel_genes": int(pop["gene_symbol"].nunique()),
        },
        "exon_phase_symmetry": exon_phase_symmetry,
        "biophysical_and_density_controls": curation_free_table,
        "rare_benign_internal_ptvs": rare_benign_ptvs,
        "vus_bounding_analysis": vus_bounding_analysis,
        "splice_junction_boundary_asymmetry": boundary,
        "inframe_indels_and_nmd_zones": {
            "indel_matrix": indel_table,
            "ptc_zones": ptc_zones,
            "nmd_zones_587k": nmd_zone_table,
        },
        "cmh_within_gene_missense": cmh_out,
        "exon_skipping_nmd_candidates": {
            "tier1_exons_n": len(t1_exons),
            "tier1_genes_n": int(t1_exons["gene_symbol"].nunique()),
            "tier1_patient_ptvs": int(t1_exons["ptv_plp"].sum()),
            "tier1_nonsense_ptvs": int(t1_exons["stop_plp"].sum()),
            "tier1_frameshift_ptvs": int(t1_exons["fs_plp"].sum()),
            "tier12_exons_n": len(t12_exons),
            "tier12_genes_n": int(t12_exons["gene_symbol"].nunique()),
            "tier12_patient_ptvs": int(t12_exons["ptv_plp"].sum()),
            "tier12_nonsense_ptvs": int(t12_exons["stop_plp"].sum()),
            "tier12_frameshift_ptvs": int(t12_exons["fs_plp"].sum()),
            "ptv_count_note": "'patient PTVs' counts ClinVar records, not patients.",
            "top_15_candidate_exons": top_aso_by_burden[out_atlas_cols].to_dict(
                orient="records"
            ),
        },
        "tiered_model_and_bayesian_lr": {
            "model_comparison_test_24_genes": model_comp,
            "consequence_only_baseline": cons_baseline,
            "within_missense_auroc": within_missense,
            "per_class_metrics_test": per_class,
            "fitted_arch_offsets": offset_fit,
            "bayesian_lr_table": lr_table,
            "bayesian_lr_method": "LR(cell) = P(cell|P/LP) / P(cell|B/LB) on resolved internal missense variants, no smoothing; 95% CI on the log scale. Research-only; not a clinical classification.",
            "donor_edge_lr": donor_lr,
            "regime_ii_lr": regime2_lr,
            "vus_stratification_counts": vus_counts,
            "total_vus_stratified": len(vus_df),
            "vus_file_note": "results/vus_bayesian_stratification.tsv carries per-variant LRs with CIs and no ACMG labels; see results/vus_bayesian_stratification_NOTE.md.",
        },
        "sensitivity_checks": sensitivity,
    }
    out_json_path.write_text(json.dumps(clean_json(out), indent=2))
    print(
        f"[18] wrote results/exon_architecture_nmd.json, results/exon_skipping_nmd_candidates.tsv ({len(atlas_sorted):,d} exons), results/vus_bayesian_stratification.tsv ({len(vus_df):,d} rows), results/cmh_per_gene_tables.tsv"
    )
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
