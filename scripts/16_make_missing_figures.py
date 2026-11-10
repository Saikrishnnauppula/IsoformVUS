#!/usr/bin/env python3
"""Stage 16: the remaining manuscript figures.

Generates figures/01_architecture.png, 02_workflow.png, 03_isoform_usage.png,
08_case_study.png, and (since 2026-11-09) replotted 09-12 from the shipped
stage-18 outputs in results/.

- 01_architecture and 02_workflow are schematics drawn from the manuscript's
  methods (§2) and labeled as such; all counts inside them come from
  results and data files (data/processed/feature_tiers.json,
  results/splits.json, data/DATA_PROVENANCE.md).
- 03_isoform_usage and 08_case_study are data-driven, from
  data/processed/variant_isoform_effects.tsv.gz,
  data/processed/isoform_usage.tsv.gz, and results/counterfactual.json.
"""

from __future__ import annotations

import json
from collections import Counter
from pathlib import Path

import matplotlib

matplotlib.use("Agg")
import matplotlib.patches as mpatches
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd

REPO = Path(__file__).resolve().parents[1]
FIG = REPO / "figures"
PROC = REPO / "data" / "processed"
FIG.mkdir(exist_ok=True)

plt.rcParams.update(
    {"font.size": 9, "axes.spines.top": False, "axes.spines.right": False}
)

BLUE = "#4C78A8"
GREEN = "#72B66F"
ORANGE = "#F58518"
RED = "#E45756"
PURPLE = "#B279A2"
GRAY = "#9E9E9E"
TEAL = "#54A24B"


def _box(ax, xy, w, h, text, fc, fs=8.5, ec="0.25", alpha=0.9, style="round,pad=0.02"):
    box = mpatches.FancyBboxPatch(
        xy, w, h, boxstyle=style, fc=fc, ec=ec, alpha=alpha, linewidth=1.0
    )
    ax.add_patch(box)
    ax.text(
        xy[0] + w / 2,
        xy[1] + h / 2,
        text,
        ha="center",
        va="center",
        fontsize=fs,
        wrap=True,
    )


def _arrow(ax, x0, y0, x1, y1, color="0.35"):
    ax.annotate(
        "",
        xy=(x1, y1),
        xytext=(x0, y0),
        arrowprops=dict(arrowstyle="-|>", color=color, lw=1.2, shrinkA=1, shrinkB=2),
    )


def fig_architecture():
    tiers = json.loads((PROC / "feature_tiers.json").read_text())["tiers"]
    counts = {m: len(tiers[m]) for m in ["M1", "M2", "M3", "M4", "M5"]}
    # Counts after the 2026-11-09 refit (CADD columns removed from all tiers).
    assert counts == {"M1": 7, "M2": 9, "M3": 12, "M4": 21, "M5": 24}, counts

    fig, ax = plt.subplots(figsize=(12, 6.4))
    ax.set_xlim(0, 12)
    ax.set_ylim(0, 6.4)
    ax.axis("off")
    fig.suptitle(
        "IsoformVUS model architecture (schematic — counts from "
        "data/processed/feature_tiers.json)",
        fontsize=11,
        y=0.97,
    )

    inputs = [
        ("DNA variant\n(ClinVar P/LP vs B/LB)", BLUE),
        ("Phenotype\n(HPO terms)", GREEN),
        ("Cell type\n(GTEx tissue\nrelevance)", TEAL),
        ("Isoform structure\n(GENCODE v46 +\nMANE v1.5)", ORANGE),
        ("Isoform usage\n(GTEx transcript\nusage)", PURPLE),
    ]
    tier_names = [
        "M1\n7 features",
        "M2\n9 features",
        "M3\n12 features",
        "M4\n21 features",
        "M5\n24 features",
    ]
    added = [
        "variant\nfeatures",
        "+ 2 phenotype",
        "+ 3 cell-type",
        "+ 9 isoform\nstructure",
        "+ 3 isoform\nusage",
    ]

    for i, ((label, fc), tn, ad) in enumerate(zip(inputs, tier_names, added)):
        x = 0.35 + i * 2.12
        _box(ax, (x, 4.15), 1.75, 1.25, label, fc, alpha=0.28)
        _box(ax, (x, 2.35), 1.75, 1.35, f"{tn}\n{ad}", fc, alpha=0.55)
        _arrow(ax, x + 0.875, 4.15, x + 0.875, 3.72)
        if i > 0:
            _arrow(ax, x - 0.37, 3.02, x, 3.02)

    _box(
        ax,
        (0.35, 0.55),
        10.25,
        1.30,
        "Learner: logistic regression (interpretable)\n"
        "Output: pathogenicity score per variant  •  Gene-level grouped splits "
        "(train 61 / val 17 / test 24 genes, disease-stratified)\n"
        "Critical comparison: M4 vs M3 — identical data, splits, and learner; "
        "only the isoform-structure features differ",
        GRAY,
        alpha=0.30,
        fs=8.5,
    )
    for i in range(5):
        _arrow(ax, 0.35 + i * 2.12 + 0.875, 2.35, 0.35 + i * 2.12 + 0.875, 1.87)

    _box(ax, (10.95, 2.35), 0.0, 0.0, "", "white")  # spacer
    ax.text(
        11.15,
        4.9,
        "Baselines\n(comparison only):",
        fontsize=8.5,
        va="top",
        ha="left",
        style="italic",
    )
    for j, b in enumerate(
        [
            "CADD v1.7\nAUROC 0.9554",
            "SpliceAI\nAUROC 0.6961",
            "AlphaMissense\nAUROC 0.9435",
        ]
    ):
        _box(ax, (10.85, 3.9 - j * 1.15), 1.05, 0.95, b, RED, alpha=0.25, fs=7.5)

    ax.text(
        0.35,
        0.28,
        "Rigorous accounting: the pre-rebuild +0.031 AUPRC gains for "
        "phenotype & cell-type layers were a tie artifact; refit paired "
        "AUROC diffs are M2-M1 -0.0016 and M3-M1 -0.0023; tissue-level "
        "isoform usage adds -0.0001 AUROC over structural features (M5-M4).",
        fontsize=7.5,
        style="italic",
        color="0.35",
    )
    fig.tight_layout(rect=[0, 0.02, 1, 0.93])
    fig.savefig(FIG / "01_architecture.png", dpi=150)
    plt.close(fig)
    print("[16] 01_architecture.png")


def fig_workflow():
    fig, ax = plt.subplots(figsize=(12, 7.0))
    ax.set_xlim(0, 12)
    ax.set_ylim(0, 7.0)
    ax.axis("off")
    fig.suptitle(
        "IsoformVUS data & analysis workflow (schematic — numbers "
        "from results/*.json and data/DATA_PROVENANCE.md)",
        fontsize=11,
        y=0.97,
    )

    row1 = [
        ("ClinVar 20261004\n(GRCh38)", BLUE),
        ("102-gene panel\n(diversity-capped;\nTTN dominance handled)", BLUE),
        ("60,521 labeled variants\n30,146 P/LP · 30,375 B/LB", BLUE),
    ]
    row2 = [
        ("GENCODE v46 +\nMANE v1.5", ORANGE),
        ("Per-transcript effect\nclassification\n76.9% isoform-specific", ORANGE),
        ("HPO v2026-02-16 →\nphenotype → cell type\n(cell-relevance score)", GREEN),
    ]
    row3 = [
        (
            "GTEx tissue expression\n→ cell relevance\n→ isoform usage\n(categorical)",
            TEAL,
        ),
        ("24 features · 5 tiers\nM1→M5\n(gene-level splits\n61 / 17 / 24)", PURPLE),
        (
            "Locked held-out test\n24 genes · 14,400 variants\n16 logged accesses · 20 reads",
            RED,
        ),
    ]

    for r, (row, y) in enumerate([(row1, 5.15), (row2, 3.35), (row3, 1.55)]):
        for i, (label, fc) in enumerate(row):
            x = 0.6 + i * 3.75
            _box(ax, (x, y), 3.1, 1.35, label, fc, alpha=0.35 if r < 2 else 0.5)
            if i > 0:
                _arrow(ax, x - 0.65, y + 0.675, x, y + 0.675)
        if r == 0:
            _arrow(ax, 10.45, 5.15, 10.45, 4.72)
        elif r == 1:
            _arrow(ax, 10.45, 3.35, 10.45, 2.92)

    bench = json.loads((REPO / "results/benchmark.json").read_text())
    stats = json.loads((REPO / "results/statistics.json").read_text())
    m1_step = bench["tiers"]["M1"]["test"]["auprc"]
    m4_step = bench["tiers"]["M4"]["test"]["auprc"]
    tie = stats["tie_corrected_refit_no_cadd"]["tie_auprc"]
    diff = stats["auprc_diffs"]["M4_minus_M1"]["mean"]
    gc_lo, gc_hi = stats["gene_clustered_cis"]["diffs"]["M4_minus_M1"]["cluster_ci95"]
    recovered = bench["recovery_M5_vs_M1"]["top50"]["n_recovered_missed_by_baseline"]
    cadd = bench["baselines"]["cadd_phred"]["auroc"]
    spliceai = bench["baselines"]["spliceai_max_delta"]["auroc"]
    am = bench["baselines"]["alphamissense_missense_only"]["auroc"]
    banner = (
        f"Benchmark (refit, 2026-11-09): test AUPRC {m1_step:.4f} (M1) \u2192 "
        f"{m4_step:.4f} (M4) step;\n"
        f"tie-corrected {tie['M1']:.4f} \u2192 {tie['M4']:.4f} "
        f"({diff:+.4f}, gene-clustered 95% CI "
        f"[{gc_lo:+.4f}, {gc_hi:+.4f}])\n"
        f"{recovered}/50 top M5 pathogenic variants missed by M1 \u00b7 "
        f"Baselines: CADD {cadd:.4f} / SpliceAI {spliceai:.4f} / "
        f"AlphaMissense {am:.4f} (missense-only)"
    )
    _box(
        ax,
        (0.6, 0.15),
        10.8,
        0.95,
        banner,
        GRAY,
        alpha=0.25,
        fs=8,
    )
    # NOTE: no tight_layout here — coordinates are manual, and tight_layout
    # squeezes the axes when long wrapped texts are present.
    fig.savefig(FIG / "02_workflow.png", dpi=150)
    plt.close(fig)
    print("[16] 02_workflow.png")


def fig_isoform_usage():
    eff = pd.read_csv(
        PROC / "variant_isoform_effects.tsv.gz",
        sep="\t",
        usecols=[
            "gene_effect_summary",
            "mane_select_effect",
            "n_isoform_specific_hits",
        ],
    )
    usage = pd.read_csv(
        PROC / "isoform_usage.tsv.gz", sep="\t", usecols=["usage_category"]
    )
    n = len(eff)
    assert n == 60521, n

    fig, axes = plt.subplots(1, 3, figsize=(13, 4.6))

    # Panel A: per-variant gene-level effect summary
    order_a = ["isoform_specific", "some", "all", "none"]
    counts_a = eff["gene_effect_summary"].value_counts().reindex(order_a)
    colors_a = [ORANGE, BLUE, GREEN, GRAY]
    axes[0].bar(range(len(order_a)), counts_a.values, color=colors_a)
    axes[0].set_xticks(
        range(len(order_a)),
        ["isoform-\nspecific", "some\ntranscripts", "all\ntranscripts", "none"],
        fontsize=8,
    )
    axes[0].set_ylabel("Variants (n = 60,521)")
    axes[0].set_title("A. Gene-level effect summary per variant")
    for i, v in enumerate(counts_a.values):
        axes[0].text(i, v + 600, f"{v:,}\n({v / n:.1%})", ha="center", fontsize=8)

    # Panel B: effect on the single MANE Select transcript
    order_b = [
        "coding_specific",
        "intronic",
        "junction_specific",
        "shared",
        "exon_specific",
        "none",
    ]
    counts_b = eff["mane_select_effect"].value_counts().reindex(order_b)
    labels_b = [
        "coding-\nspecific",
        "intronic",
        "junction-\nspecific",
        "shared\n(all tx)",
        "exon-\nspecific",
        "none",
    ]
    axes[1].bar(
        range(len(order_b)),
        counts_b.values,
        color=[RED, GRAY, PURPLE, GREEN, ORANGE, "0.7"],
    )
    axes[1].set_xticks(range(len(order_b)), labels_b, fontsize=8)
    axes[1].set_ylabel("Variants (n = 60,521)")
    axes[1].set_title("B. Effect on the MANE Select transcript alone")
    for i, v in enumerate(counts_b.values):
        axes[1].text(i, v + 600, f"{v / n:.1%}", ha="center", fontsize=7.5)
    axes[1].text(
        0.98,
        0.96,
        "37.6% of variants are NOT\ncoding-specific on MANE —\n"
        "the canonical transcript\nalone can mislead",
        transform=axes[1].transAxes,
        ha="right",
        va="top",
        fontsize=8,
        style="italic",
        bbox=dict(boxstyle="round,pad=0.3", fc="white", ec="0.6"),
    )

    # Panel C: GTEx isoform-usage evidence categories
    order_c = ["supported", "uncertain"]
    counts_c = usage["usage_category"].value_counts().reindex(order_c)
    axes[2].bar(range(len(order_c)), counts_c.values, color=[TEAL, GRAY])
    axes[2].set_xticks(range(len(order_c)), ["supported", "uncertain"])
    axes[2].set_ylabel("Gene × transcript × tissue rows (n = 7,057)")
    axes[2].set_title("C. Isoform-usage evidence (GTEx, categorical)")
    for i, v in enumerate(counts_c.values):
        axes[2].text(
            i, v + 60, f"{v:,}\n({v / counts_c.sum():.1%})", ha="center", fontsize=8
        )
    axes[2].text(
        0.02,
        0.96,
        "Short-read ambiguity →\n"
        "categorical evidence, not\n"
        "false precision. Tissue-level\n"
        "usage added nothing over\n"
        "structural features (no pre-registration was written).",
        transform=axes[2].transAxes,
        ha="left",
        va="top",
        fontsize=8,
        style="italic",
        bbox=dict(boxstyle="round,pad=0.3", fc="white", ec="0.6"),
    )

    fig.suptitle(
        "Isoform usage and isoform-specificity across 60,521 variants "
        "(data/processed/)",
        fontsize=11,
    )
    fig.tight_layout()
    fig.savefig(FIG / "03_isoform_usage.png", dpi=150)
    plt.close(fig)
    print("[16] 03_isoform_usage.png")


def fig_case_study():
    cf = json.loads((REPO / "results/counterfactual.json").read_text())
    per_tx = cf["per_transcript"]
    fx = Counter(t["effect"] for t in per_tx)
    n_tx = len(per_tx)
    # ground truth check against raw GTF: 45 LMNA transcripts
    assert n_tx == 45, n_tx

    mane = cf["mane_select"]
    mane_eff = next(t["effect"] for t in per_tx if t["transcript"] == mane)

    fig, axes = plt.subplots(
        1, 3, figsize=(13.5, 4.6), gridspec_kw={"width_ratios": [1.15, 0.85, 1.0]}
    )

    # Panel A: per-transcript effect breakdown for the LMNA variant
    order = [
        "coding_specific",
        "exon_specific",
        "junction_specific",
        "intronic",
        "none",
    ]
    labels = [
        "coding-specific",
        "exon-specific",
        "junction-specific",
        "intronic",
        "no effect",
    ]
    vals = [fx.get(k, 0) for k in order]
    colors = [RED, ORANGE, PURPLE, GRAY, "0.75"]
    y = np.arange(len(order))
    axes[0].barh(y, vals, color=colors)
    axes[0].set_yticks(y, labels, fontsize=8)
    axes[0].set_xlabel("Transcripts (n = 45, GENCODE v46)")
    axes[0].set_title("A. chr1:156115000:C>T across LMNA transcripts")
    for i, v in enumerate(vals):
        axes[0].text(v + 0.4, i, str(v), va="center", fontsize=8)
    axes[0].text(
        0.62,
        0.58,
        f"MANE Select ({mane}):\n{mane_eff}",
        transform=axes[0].transAxes,
        ha="center",
        va="center",
        fontsize=8,
        style="italic",
        bbox=dict(boxstyle="round,pad=0.3", fc="white", ec="0.6"),
    )

    # Panel B: counterfactual score delta
    s0, s1 = (
        cf["score_assuming_irrelevant_isoform"],
        cf["score_observed_isoform_context"],
    )
    axes[1].bar(
        ["irrelevant isoform\nassumed", "observed isoform\ncontext"],
        [s0, s1],
        color=[GRAY, RED],
        width=0.6,
    )
    axes[1].set_ylim(0.90, 1.0)
    axes[1].set_ylabel("M5 pathogenicity score")
    axes[1].set_title("B. Same DNA, different isoform assumption")
    axes[1].text(0, s0 + 0.002, f"{s0:.4f}", ha="center", fontsize=8)
    axes[1].text(1, s1 + 0.002, f"{s1:.4f}", ha="center", fontsize=8)
    axes[1].annotate(
        "",
        xy=(1, 0.985),
        xytext=(0, 0.958),
        arrowprops=dict(arrowstyle="<|-|>", color="0.3", lw=1.4),
    )
    axes[1].text(
        0.5,
        0.972,
        f"Δ = +{cf['delta']:.4f}",
        ha="center",
        fontsize=9,
        weight="bold",
        transform=axes[1].transData,
        bbox=dict(boxstyle="round,pad=0.2", fc="white", ec="0.5"),
    )

    # Panel C: interpretation chain
    axes[2].axis("off")
    axes[2].set_title("C. Interpretation chain (benchmark-grounded)")
    chain = [
        (
            "Disease",
            "laminopathies\n(muscular dystrophy /\ncardiomyopathy / progeroid)",
        ),
        ("Cell types", "skeletal muscle cell\ncardiac muscle cell\n(GTEx-supported)"),
        (
            "Gene \u2192 isoforms",
            f"LMNA \u2192 {n_tx} transcripts;\n"
            f"coding-specific on {fx.get('coding_specific', 0)}",
        ),
        ("Prioritization", "M5 score 0.9883\ntest rank #8\n(missed by M1 top-50)"),
    ]
    for i, (head, body) in enumerate(chain):
        yy = 0.88 - i * 0.22
        axes[2].text(
            0.02,
            yy,
            head + ":",
            fontsize=9,
            weight="bold",
            va="top",
            transform=axes[2].transAxes,
        )
        axes[2].text(
            0.38, yy, body, fontsize=8.5, va="top", transform=axes[2].transAxes
        )
        if i < 3:
            axes[2].annotate(
                "",
                xy=(0.16, yy - 0.155),
                xytext=(0.16, yy - 0.06),
                xycoords="axes fraction",
                arrowprops=dict(arrowstyle="-|>", color="0.4"),
            )

    fig.suptitle(
        "Case study: LMNA chr1:156115000:C>T — pathogenicity is "
        "invisible without the isoform context",
        fontsize=11,
    )
    fig.text(
        0.01,
        0.01,
        "Source: results/counterfactual.json (45 GENCODE v46 transcripts; "
        "scores from the trained M5 model).",
        fontsize=7,
        color="0.45",
        style="italic",
    )
    fig.tight_layout(rect=[0, 0.03, 1, 0.93])
    fig.savefig(FIG / "08_case_study.png", dpi=150)
    plt.close(fig)
    print("[16] 08_case_study.png")


SIZE_SHORT = {
    "1_Compact_<=180": "Compact\n<=180 bp",
    "2_Medium_181_400": "Medium\n181-400 bp",
    "3_Large_401_1000": "Large\n401-1000 bp",
    "4_Giant_>1000": "Giant\n>1000 bp",
}
SIZE_ORDER = list(SIZE_SHORT)


def _stage18():
    return json.loads((REPO / "results/exon_architecture_nmd.json").read_text())


def _withdrawn():
    """Carried-over values withdrawn 2026-11-09 (no generating code;
    withdrawn, do not cite). Figures 10/11 read them from here, not from
    the stage-18 results JSON."""
    return json.loads((REPO / "results/withdrawn_carried_over_values.json").read_text())


def fig_exon_architecture():
    """Figure 09: 4x2x2 exon tensor (missense resolved P/LP rate), NMD zones,
    and in-frame indel resolved rates, from results/exon_architecture_nmd.json."""
    d = _stage18()
    tensor = d["exon_phase_symmetry"]["tensor_4x2x2"]
    fig, axes = plt.subplots(1, 3, figsize=(13.5, 4.4))
    combos = [(0, 0), (0, 1), (1, 0), (1, 1)]
    labels = [
        "constitutive, asymmetric",
        "constitutive, symmetric",
        "backup-skipped, asymmetric",
        "backup-skipped, symmetric",
    ]
    colors = [BLUE, GREEN, ORANGE, RED]
    x = np.arange(len(SIZE_ORDER))
    for j, ((sk, sy), lab, col) in enumerate(zip(combos, labels, colors)):
        vals = []
        for sb in SIZE_ORDER:
            row = next(
                r
                for r in tensor
                if r["size_bin"] == sb
                and r["backup_skip"] == sk
                and r["symmetric_phase"] == sy
            )
            vals.append(row["missense_res_rate"] * 100)
        axes[0].bar(x + (j - 1.5) * 0.19, vals, width=0.19, color=col, label=lab)
    for sk, marker in [(0, "D"), (1, "s")]:
        ptv = []
        for sb in SIZE_ORDER:
            rows = [r for r in tensor if r["size_bin"] == sb and r["backup_skip"] == sk]
            plp = sum(r["frameshift_plp"] + r["stop_gained_plp"] for r in rows)
            blb = sum(r["frameshift_blb"] + r["stop_gained_blb"] for r in rows)
            ptv.append(100 * plp / (plp + blb))
        axes[0].plot(
            x + (sk - 0.5) * 0.38,
            ptv,
            marker=marker,
            color="0.15",
            ls="--",
            lw=1,
            label=f"PTV resolved P/LP (skip={sk})",
        )
    axes[0].set_xticks(x, [SIZE_SHORT[s] for s in SIZE_ORDER], fontsize=8)
    axes[0].set_ylabel("Missense resolved P/LP (%)")
    axes[0].set_title("A. Exon tensor: missense by size x skip x symmetry")
    axes[0].legend(fontsize=6.5, loc="lower left")

    zones = d["inframe_indels_and_nmd_zones"]["nmd_zones_587k"]
    zx = np.arange(len(zones))
    axes[1].bar(
        zx - 0.2,
        [z["raw_plp_rate"] * 100 for z in zones],
        width=0.4,
        color=BLUE,
        label="P/LP rate (all labels)",
    )
    axes[1].bar(
        zx + 0.2,
        [z["vus_conf_rate"] * 100 for z in zones],
        width=0.4,
        color=ORANGE,
        label="VUS-conflicting rate",
    )
    axes[1].set_xticks(
        zx, [f"Zone {i + 1}\nn={z['n']:,}" for i, z in enumerate(zones)], fontsize=8
    )
    axes[1].set_ylabel("Rate (%)")
    axes[1].set_title("B. NMD zones (587k ClinVar)")
    axes[1].legend(fontsize=7)

    indels = d["inframe_indels_and_nmd_zones"]["indel_matrix"]
    ix = np.arange(len(indels))
    axes[2].bar(
        ix,
        [r["res_plp"] * 100 for r in indels],
        color=[RED if r["class"] == "inframe_deletion" else TEAL for r in indels],
    )
    axes[2].set_xticks(
        ix,
        [
            f"{'del' if r['class'] == 'inframe_deletion' else 'ins'}\n"
            f"{'<=400' if '<=' in r['size_group'] else '>400'}\nskip={r['backup_skip']}"
            for r in indels
        ],
        fontsize=7,
    )
    axes[2].set_ylabel("Resolved P/LP (%)")
    axes[2].set_title("C. In-frame indels (all 8 cells)")
    fig.suptitle(
        "Exon architecture and NMD — replotted 2026-11-09 from "
        "results/exon_architecture_nmd.json"
    )
    fig.tight_layout(rect=[0, 0, 1, 0.92])
    fig.savefig(FIG / "09_exon_architecture_nmd.png", dpi=150)
    plt.close(fig)
    print("[16] 09_exon_architecture_nmd.png")


def fig_mediation_controls():
    """Figure 10: label-free biophysical controls per exon group, and the
    (carried-over, withdrawn) AlphaMissense mediation panel with sensitivity
    reductions, from results/exon_architecture_nmd.json."""
    d = _stage18()
    fig, axes = plt.subplots(1, 2, figsize=(12.5, 4.6))
    rows = d["biophysical_and_density_controls"]
    x = np.arange(len(rows))
    axes[0].bar(
        x - 0.27,
        [r["am_median_all"] * 100 for r in rows],
        width=0.27,
        color=BLUE,
        label="AlphaMissense median (x100)",
    )
    axes[0].bar(
        x,
        [r["am_frac_gt_056"] * 100 for r in rows],
        width=0.27,
        color=GREEN,
        label="% AM > 0.56",
    )
    axes[0].bar(
        x + 0.27,
        [r["blb_to_plp_mis_ratio"] * 4 for r in rows],
        width=0.27,
        color=ORANGE,
        label="B/LB:P/LP missense ratio (x4, labels true)",
    )
    for i, r in enumerate(rows):
        axes[0].text(
            i + 0.27,
            r["blb_to_plp_mis_ratio"] * 4 + 2,
            f"{r['blb_to_plp_mis_ratio']:.1f}x",
            ha="center",
            fontsize=7,
        )
    axes[0].set_xticks(
        x,
        [
            f"{SIZE_SHORT[r['size_bin']].splitlines()[0]}\nskip={r['backup_skip']}"
            for r in rows
        ],
        fontsize=7.5,
    )
    axes[0].set_ylabel("Percent / scaled ratio")
    axes[0].set_title("A. Label-free controls by exon group")
    axes[0].legend(fontsize=7)

    med = _withdrawn()["causal_mediation"]
    entries = []
    for key, lab in [
        ("giant_large_via_alphamissense", "Giant/large via AM"),
        ("backup_skip_via_alphamissense", "Backup-skip via AM"),
    ]:
        for adj, alab in [("unadjusted", "unadj."), ("within_gene_FE", "gene FE")]:
            v = med[key][adj]
            entries.append((f"{lab} ({alab})", v["prop_mediated"] * 100, PURPLE, ""))
    sc = d["sensitivity_checks"]
    sens = [
        (
            "Leave-top-genes-out",
            sc["test1_leave_top_genes_out"],
            "compact_constit_res_plp",
            "giant_skip_res_plp",
        ),
        (
            "High-confidence (>=2 stars)",
            sc["test2_high_confidence_2plus_stars"],
            "compact_constit_mis_p_rate",
            "giant_skip_mis_p_rate",
        ),
        (
            "Prospective 2024->2026",
            sc["test7_prospective_2024_to_2026"],
            "compact_constit_p_rate",
            "giant_skip_p_rate",
        ),
    ]
    for lab, blk, ck, gk in sens:
        red = 100 * (blk[ck] - blk[gk]) / blk[ck]
        entries.append((f"{lab} reduction", red, GREEN, ""))
    y = np.arange(len(entries))[::-1]
    axes[1].barh(y, [e[1] for e in entries], color=[e[2] for e in entries])
    axes[1].set_yticks(y, [e[0] for e in entries], fontsize=7.5)
    for yi, e in zip(y, entries):
        axes[1].text(e[1] + 1.5, yi, f"{e[1]:.1f}%", va="center", fontsize=7.5)
    axes[1].set_xlabel("Proportion mediated / relative reduction (%)")
    axes[1].set_title("B. Mediation (withdrawn) and sensitivity reductions")
    axes[1].text(
        0.02,
        0.02,
        "Mediation values carried over, not recomputed;\n"
        "mediation percentages withdrawn (manuscript 6.1).",
        transform=axes[1].transAxes,
        fontsize=7,
        style="italic",
        bbox=dict(boxstyle="round,pad=0.3", fc="white", ec="0.6"),
    )
    fig.suptitle(
        "Structural mediation controls — replotted 2026-11-09; mediation "
        "panel from results/withdrawn_carried_over_values.json (withdrawn)"
    )
    fig.tight_layout(rect=[0, 0, 1, 0.92])
    fig.savefig(FIG / "10_structural_mediation_controls.png", dpi=150)
    plt.close(fig)
    print("[16] 10_structural_mediation_controls.png")


def fig_boundary_asymmetry():
    """Figure 11: exon-boundary distance profile and per-cell likelihood
    ratios. The donor-edge finding is withdrawn; this figure is the record
    of it, replotted from results/exon_architecture_nmd.json."""
    d = _stage18()
    sj = _withdrawn()["splice_junction_boundary_asymmetry"]
    fig, axes = plt.subplots(1, 2, figsize=(12.5, 5.0))
    prof = sj["distance_profile"]
    x = np.arange(len(prof))
    axes[0].plot(
        x,
        [p["am_benign_fn_rate"] * 100 for p in prof],
        marker="o",
        color=RED,
        label="AM-benign false-negative rate",
    )
    axes[0].plot(
        x,
        [p["pop_mis_res_plp"] * 100 for p in prof],
        marker="s",
        color=BLUE,
        label="Missense resolved P/LP",
    )
    axes[0].plot(
        x,
        [p["pop_syn_res_plp"] * 100 for p in prof],
        marker="^",
        color=GREEN,
        label="Synonymous resolved P/LP",
    )
    axes[0].set_xticks(x, [p["bin"].split(" (")[0] for p in prof], fontsize=7.5)
    axes[0].set_xlabel("Distance from exon boundary")
    axes[0].set_ylabel("Rate (%)")
    dva = sj["donor_vs_acceptor_edge_0_3bp"]
    axes[0].set_title(
        "A. Distance profile (carried over, not recomputed)\n"
        f"edge 0-3 bp AM-benign FN: donor {dva['am_benign_donor_fn']:.1%}, "
        f"acceptor {dva['am_benign_acceptor_fn']:.1%}"
    )
    axes[0].legend(fontsize=7)

    tl = d["tiered_model_and_bayesian_lr"]
    rows = sorted(tl["bayesian_lr_table"], key=lambda r: r["likelihood_ratio"])
    labels = [
        f"{SIZE_SHORT[r['size_bin']].splitlines()[0]} skip={r['backup_skip']} "
        f"sym={r['symmetric_phase']}"
        for r in rows
    ]
    lrs = [r["likelihood_ratio"] for r in rows]
    los = [np.log2(r["likelihood_ratio"] / r["lr_ci95_lo"]) for r in rows]
    his = [np.log2(r["lr_ci95_hi"] / r["likelihood_ratio"]) for r in rows]
    y = np.arange(len(rows) + 1)
    axes[1].barh(y[1:], np.log2(lrs), color=BLUE)
    axes[1].errorbar(
        np.log2(lrs), y[1:], xerr=[los, his], fmt="none", color="0.25", lw=0.8
    )
    dlr = tl["donor_edge_lr"]
    axes[1].barh(0, np.log2(dlr["likelihood_ratio"]), color=RED)
    axes[1].errorbar(
        np.log2(dlr["likelihood_ratio"]),
        0,
        xerr=[
            [np.log2(dlr["likelihood_ratio"] / dlr["lr_ci95_lo"])],
            [np.log2(dlr["lr_ci95_hi"] / dlr["likelihood_ratio"])],
        ],
        fmt="none",
        color="0.25",
        lw=0.8,
    )
    axes[1].set_yticks(y, ["donor edge (<=3 bp, <=400 bp)"] + labels, fontsize=7)
    axes[1].axvline(0, color="0.3", lw=1)
    axes[1].set_xlabel("log2 likelihood ratio (resolved missense)")
    axes[1].set_title("B. Cell likelihood ratios — research only, no ACMG points")
    fig.suptitle(
        "Donor-acceptor boundary asymmetry (withdrawn finding) — replotted "
        "2026-11-09 from results/exon_architecture_nmd.json"
    )
    fig.tight_layout(rect=[0, 0, 1, 0.9])
    fig.savefig(FIG / "11_donor_acceptor_boundary_asymmetry.png", dpi=150)
    plt.close(fig)
    print("[16] 11_donor_acceptor_boundary_asymmetry.png")


def fig_skipping_candidates():
    """Figure 12: top exon-skipping candidate exons (computational
    nomination, no functional evidence) and the stage-18 model comparison,
    from results/exon_skipping_nmd_candidates.tsv and
    results/exon_architecture_nmd.json."""
    d = _stage18()
    fig, axes = plt.subplots(1, 2, figsize=(13, 5.2))
    cand = pd.read_csv(REPO / "results/exon_skipping_nmd_candidates.tsv", sep="\t")
    top = (
        cand[cand["aso_tier"].str.startswith(("Tier 1", "Tier 2"))]
        .sort_values("ptv_plp", ascending=False)
        .head(15)
        .iloc[::-1]
    )
    y = np.arange(len(top))
    axes[0].barh(y, top["stop_plp"], color=RED, label="nonsense P/LP")
    axes[0].barh(
        y, top["fs_plp"], left=top["stop_plp"], color=ORANGE, label="frameshift P/LP"
    )
    axes[0].set_yticks(
        y,
        [
            f"{g} (exon {e}, {ln:,} bp)"
            for g, e, ln in zip(top["gene_symbol"], top["exon_id"], top["exon_len"])
        ],
        fontsize=7.5,
    )
    for yi, (_, r) in zip(y, top.iterrows()):
        axes[0].text(
            r["ptv_plp"] + 60,
            yi,
            f"{int(r['ptv_plp']):,} PTVs (missense P/LP "
            f"{r['mis_resolved_plp_rate']:.1%}, score "
            f"{r['Exon_Nomination_Score']:.2f})",
            va="center",
            fontsize=7,
        )
    axes[0].set_xlabel("ClinVar PTV records (not patients)")
    axes[0].set_title("A. Top 15 candidate exons (Tier 1/2)")
    axes[0].legend(fontsize=7)

    tl = d["tiered_model_and_bayesian_lr"]
    mc = tl["model_comparison_test_24_genes"]
    groups = [
        ("M1", mc["M1_Baseline"]),
        ("M3", mc["M3_CellContext"]),
        ("M4", mc["M4_IsoformFlags"]),
        ("M6", mc["IsoformVUS_Mechanistic_v3"]),
        ("consequence\nonly", tl["consequence_only_baseline"]),
    ]
    x = np.arange(len(groups))
    series = [
        ("test AUROC", "test_auroc", BLUE),
        ("test tie-AUPRC", "test_tie_auprc", GREEN),
        ("internal-missense AUROC", "test_internal_missense_auroc", ORANGE),
    ]
    for j, (lab, key, col) in enumerate(series):
        axes[1].bar(
            x + (j - 1) * 0.27,
            [g[1][key] * 100 for g in groups],
            width=0.27,
            color=col,
            label=lab,
        )
    axes[1].set_xticks(x, [g[0] for g in groups], fontsize=8)
    axes[1].set_ylim(40, 105)
    axes[1].set_ylabel("Test performance (%)")
    axes[1].set_title("B. Stage-18 model comparison (24 test genes)")
    axes[1].legend(fontsize=7)
    axes[1].text(
        0.5,
        0.04,
        "M6 pooled gain is consequence-class separation\n(manuscript 5.1).",
        transform=axes[1].transAxes,
        ha="center",
        fontsize=7,
        style="italic",
        bbox=dict(boxstyle="round,pad=0.3", fc="white", ec="0.6"),
    )
    fig.suptitle(
        "Exon-skipping candidates — replotted 2026-11-09 from shipped stage-18 outputs"
    )
    fig.tight_layout(rect=[0, 0, 1, 0.93])
    fig.savefig(FIG / "12_exon_skipping_candidates.png", dpi=150)
    plt.close(fig)
    print("[16] 12_exon_skipping_candidates.png")


def main() -> int:
    fig_architecture()
    fig_workflow()
    fig_isoform_usage()
    fig_case_study()
    fig_exon_architecture()
    fig_mediation_controls()
    fig_boundary_asymmetry()
    fig_skipping_candidates()
    print("[16] done -> figures/ (01-03, 08-12 regenerated)")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
