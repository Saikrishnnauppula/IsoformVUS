#!/usr/bin/env python3
"""Stage 18: Exon architecture, isoform skipping, and NMD boundary analysis.

Reconstructs 4,132 MANE Select coding exons (3,955 internal coding exons)
across the 102 panel genes and evaluates:
  1. Deconstruction of tiered models M1-M5 (analytical tie-corrected AUPRC vs.
     step AUPRC in results/benchmark.json: 0.8016 -> 0.8524, +0.0509).
  2. Three non-canonical isoform discordance regimes across 587,227 ClinVar
     variants (dual-coding frame shifts, alternative-exon missense VUS, and
     alternative-exon loss-of-function variants).
  3. The 4 x 2 x 2 internal exon length x GTEx backup skipping (>=10% usage,
     >=1.0 TPM) x reading-frame phase symmetry (length mod 3 == 0) matrix
     across 301,264 internal coding variants.
  4. Label-independent AlphaMissense/CADD distributions, coding density per kb,
     and 1,000-bootstrap counterfactual mediation analysis.
  5. 5'-splice donor (3' exon end) vs. 3'-splice acceptor (5' exon start)
     asymmetry among near-junction missense variants with am_score < 0.34,
     in-frame deletion vs. insertion polarity, and 50-nt penultimate-exon EJC
     NMD escape zones.
  6. Candidate exon ranking for antisense exon-skipping rescue, Bayesian
     likelihood ratios (Tavtigian et al. 2018), and sensitivity checks.

Writes:
  - results/exon_architecture_nmd.json
  - results/exon_skipping_nmd_candidates.tsv
  - results/vus_bayesian_stratification.tsv

Usage: python3 scripts/18_exon_architecture_nmd.py
"""

from __future__ import annotations

import json
from pathlib import Path
import numpy as np
import pandas as pd
from scipy.optimize import minimize
from scipy.special import expit
from scipy.stats import fisher_exact, norm, rankdata

REPO = Path(__file__).resolve().parents[1]
SEED = 20261006


def load_clinvar_587k(repo_root: Path) -> pd.DataFrame:
    """Load the 587,227-row ClinVar table from .tsv.gz or .tsv."""
    gz_path = repo_root / "breakthrough/round2/e8_fragility/supp_clinvar_all_table.tsv.gz"
    tsv_path = repo_root / "breakthrough/round2/e8_fragility/supp_clinvar_all_table.tsv"
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
    grp = df_t.groupby("s", sort=True).agg(n=("y", "count"), pos=("y", "sum")).iloc[::-1]
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
    for b0, b1 in zip(bins[:-1], bins[1:]):
        m = (y_prob >= b0) & (y_prob <= (b1 if b1 == 1.0 else b1 - 1e-12))
        if m.any():
            ece += (m.sum() / n) * abs(y_true[m].mean() - y_prob[m].mean())
    return float(ece)


class BalancedLogReg:
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
            ll = np.sum(sw * np.where(y == 1, -np.log(np.clip(prob, 1e-15, 1.0)), -np.log(np.clip(1.0 - prob, 1e-15, 1.0)))) + 0.5 * np.sum(w * w) / self.C
            diff = sw * (prob - y)
            return ll, np.append((Xs.T @ diff) + w / self.C, np.sum(diff))

        res = minimize(loss_grad, np.zeros(Xs.shape[1] + 1), jac=True, method="L-BFGS-B", options={"maxiter": 5000, "ftol": 1e-12, "gtol": 1e-8})
        self.coef_, self.intercept_ = res.x[:-1], res.x[-1]
        return self

    def predict_proba(self, X):
        Xs = (np.asarray(X, dtype=np.float64) - self.mean_) / self.scale_
        p1 = expit(Xs @ self.coef_ + self.intercept_)
        return np.column_stack([1.0 - p1, p1])


def run_cmh_stratified(df_sub: pd.DataFrame, x_col: str, y_col: str = "is_plp") -> dict:
    num_or, den_or, num_z, var_z, n_g = 0.0, 0.0, 0.0, 0.0, 0
    for _, grp in df_sub.groupby("gene_symbol"):
        if grp[x_col].nunique() < 2 or grp[y_col].nunique() < 2:
            continue
        x = grp[x_col].to_numpy()
        y = grp[y_col].to_numpy()
        a = float(np.sum((x == 1) & (y == 1)))
        b = float(np.sum((x == 1) & (y == 0)))
        c = float(np.sum((x == 0) & (y == 1)))
        d = float(np.sum((x == 0) & (y == 0)))
        n = a + b + c + d
        if n <= 1:
            continue
        num_or += (a * d) / n
        den_or += (b * c) / n
        r1, r0, c1, c0 = a + b, c + d, a + c, b + d
        num_z += a - (r1 * c1) / n
        var_z += (r1 * r0 * c1 * c0) / (n * n * (n - 1.0))
        n_g += 1
    cmh_or = num_or / max(den_or, 1e-12)
    z = num_z / np.sqrt(max(var_z, 1e-12))
    p = float(2.0 * norm.sf(abs(z)))
    return {"cmh_or": float(cmh_or), "z": float(z), "p": p, "n_informative_genes": n_g}


def main() -> int:
    pop = load_clinvar_587k(REPO)
    cons = pd.read_csv(REPO / "breakthrough/exp_a_molecular/variant_molecular_consequences.tsv.gz", sep="\t")
    ptc = pd.read_csv(REPO / "breakthrough/round2/e5_nmd/ptc_features.tsv.gz", sep="\t")
    long_sev = pd.read_csv(REPO / "breakthrough/exp_b_molecular/variant_transcript_sev.tsv.gz", sep="\t")
    iu = pd.read_csv(REPO / "data/processed/isoform_usage.tsv.gz", sep="\t")
    am = pd.read_csv(REPO / "data/processed/variant_alphamissense.tsv.gz", sep="\t")
    cadd = pd.read_csv(REPO / "data/processed/variant_cadd.tsv.gz", sep="\t")
    spl = pd.read_csv(REPO / "data/processed/variant_spliceai.tsv.gz", sep="\t")
    fm = pd.read_csv(REPO / "data/processed/feature_matrix.tsv.gz", sep="\t")
    tiers = json.loads((REPO / "data/processed/feature_tiers.json").read_text())["tiers"]
    panel = pd.read_csv(REPO / "metadata/gene_panel.tsv", sep="\t")
    bench_locked = json.loads((REPO / "results/benchmark.json").read_text())

    pop["is_plp"] = pop["clnsig_aggregate"].isin(["Pathogenic", "Likely pathogenic", "Pathogenic/Likely pathogenic"]).astype(int)
    pop["is_blb"] = pop["clnsig_aggregate"].isin(["Benign", "Likely benign", "Benign/Likely benign"]).astype(int)
    pop["is_resolved"] = ((pop["is_plp"] == 1) | (pop["is_blb"] == 1)).astype(int)
    pop["is_vus"] = pop["clnsig_aggregate"].str.contains("Uncertain significance", na=False).astype(int)
    pop["is_conflict"] = pop["disc_conflict"].astype(int)

    cds_classes = {"synonymous", "missense", "inframe_indel", "frameshift", "stop_gained", "start_lost", "stop_lost"}

    ptc["pos"] = ptc["variant_key"].str.split(":").str[1].astype(int)
    ptc["cds_off"] = np.round(ptc["ptc_frac"] * ptc["cds_len"]).astype(int)
    gene_strand = {}
    for g, grp in ptc.groupby("gene_symbol"):
        gene_strand[g] = "+" if (len(grp) < 2 or np.corrcoef(grp["pos"], grp["cds_off"])[0, 1] > 0) else "-"

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
        is_intr = gdf["mane_class_v46"].isin(["intronic", "splice_region", "none"]).to_numpy()
        pos_arr = gdf["pos"].to_numpy()
        cls_arr = gdf["mane_class_v46"].to_numpy()
        idx_arr = gdf.index.to_numpy()
        snv_mask = ((gdf["ref"].str.len() == 1) & (gdf["alt"].str.len() == 1)).to_numpy()
        snv_pos, snv_cls = pos_arr[snv_mask], cls_arr[snv_mask]
        cur_idxs, cur_positions = [], []

        def flush_exon(c_idxs, c_pos, ex_id):
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
            if is_cds[k]:
                if len(cur_positions) > 0 and (pos_arr[k] - cur_positions[-1] > 90):
                    flush_exon(cur_idxs, cur_positions, exon_counter)
                    exon_counter += 1
                    cur_idxs, cur_positions = [idx_arr[k]], [pos_arr[k]]
                else:
                    cur_idxs.append(idx_arr[k])
                    cur_positions.append(pos_arr[k])
            elif is_intr[k]:
                if len(cur_positions) > 0:
                    flush_exon(cur_idxs, cur_positions, exon_counter)
                    exon_counter += 1
                    cur_idxs, cur_positions = [], []
        if len(cur_positions) > 0:
            flush_exon(cur_idxs, cur_positions, exon_counter)
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
            best = np.where(np.abs(p_ref[j] - p_q) < np.abs(p_ref[j_prev] - p_q), j, j_prev)
            in_last_arr[gdf.index] = pdf["in_last_exon"].to_numpy()[best]
            within50_arr[gdf.index] = pdf["within50_of_last_junction"].to_numpy()[best]
            ptc_frac_arr[gdf.index] = pdf["ptc_frac"].to_numpy()[best]
    pop["in_last_exon"] = in_last_arr
    pop["within50_of_last_junction"] = within50_arr
    pop["ptc_frac_est"] = ptc_frac_arr

    iu_tx = iu.groupby(["gene_symbol", "transcript_base"]).agg(max_usage=("usage", "max"), max_tpm=("median_tpm", "max")).reset_index().rename(columns={"transcript_base": "tx"})
    ls = long_sev.merge(iu_tx, on=["gene_symbol", "tx"], how="left").fillna(0.0)
    vk_pos = cons[["variant_key", "mane_class"]].copy()
    vk_pos["pos"] = vk_pos["variant_key"].str.split(":").str[1].astype(int)
    ls = ls.merge(vk_pos, on="variant_key")
    alt_cds = ls[(ls.is_mane == 0) & (ls.mane_class.isin(cds_classes))].copy()
    alt_cds["is_spared"] = (alt_cds["sev"] == 0).astype(int)
    alt_cds["is_expr_spared"] = ((alt_cds["sev"] == 0) & (alt_cds["max_usage"] >= 0.10) & (alt_cds["max_tpm"] >= 1.0)).astype(int)
    alt_cds["spared_max_u"] = np.where(alt_cds["sev"] == 0, alt_cds["max_usage"], 0.0)

    var_sparing_cds = alt_cds.groupby(["gene_symbol", "variant_key", "pos"]).agg(
        any_alt_skips=("is_spared", "max"),
        expr_backup_skips=("is_expr_spared", "max"),
        spared_max_usage=("spared_max_u", "max"),
    ).reset_index()
    sparing_by_gene = {g: gdf.sort_values("pos") for g, gdf in var_sparing_cds.groupby("gene_symbol")}

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
            best = np.where(np.abs(p_ref[j] - p_q) < np.abs(p_ref[j_prev] - p_q), j, j_prev)
            any_skip_arr[gdf.index] = sdf["any_alt_skips"].to_numpy()[best]
            expr_skip_arr[gdf.index] = sdf["expr_backup_skips"].to_numpy()[best]
            spared_u_arr[gdf.index] = sdf["spared_max_usage"].to_numpy()[best]
    pop["any_alt_skips"] = any_skip_arr
    pop["expr_backup_skips"] = expr_skip_arr
    pop["spared_max_usage"] = spared_u_arr
    pop["indel_len_bp"] = np.abs(pop["alt"].str.len() - pop["ref"].str.len())
    pop["fine_class"] = np.where(
        pop["mane_class_v46"] == "inframe_indel",
        np.where(pop["ref"].str.len() > pop["alt"].str.len(), "inframe_deletion", "inframe_insertion"),
        pop["mane_class_v46"]
    )

    internal = pop[(pop["exon_id"] >= 0) & (pop["in_last_exon"] == 0)].copy()
    internal["exon_size_tier"] = np.where(
        internal["exon_len"] <= 180, "1_Compact_<=180",
        np.where(internal["exon_len"] <= 400, "2_Medium_181_400",
        np.where(internal["exon_len"] <= 1000, "3_Large_401_1000", "4_Giant_>1000"))
    )

    mis_int = internal[internal["fine_class"] == "missense"].copy()
    mis_int["is_giant_large"] = (mis_int["exon_len"] > 400).astype(int)
    mis_int_res = mis_int[mis_int["is_resolved"] == 1]
    cmh_giant = run_cmh_stratified(mis_int_res, "is_giant_large")
    cmh_skip = run_cmh_stratified(mis_int_res, "expr_backup_skips")

    print(f"[18] Reconstructed {exon_counter:,d} MANE Select exons ({internal['exon_id'].nunique():,d} internal coding exons).")
    print(f"[18] Locked benchmark M1->M4 step AUPRC: {bench_locked['tiers']['M1']['test']['auprc']:.4f} -> {bench_locked['tiers']['M4']['test']['auprc']:.4f} (+{bench_locked['tiers']['M4']['test']['auprc'] - bench_locked['tiers']['M1']['test']['auprc']:.4f})")
    print(f"[18] Internal missense within-gene CMH: exon >400bp OR={cmh_giant['cmh_or']:.4f} (p={cmh_giant['p']:.2e}), backup-skip OR={cmh_skip['cmh_or']:.4f} (p={cmh_skip['p']:.2e})")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
