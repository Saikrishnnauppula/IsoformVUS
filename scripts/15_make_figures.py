#!/usr/bin/env python3
"""Stage 15: publication figures (all from real artifacts; never hand-drawn).

Generates figures/04_benchmark.png, 05_ablation.png, 06_topk_recovery.png,
07_counterfactual.png from results/*.json. Architecture/workflow/usage/case
figures are schematic and labeled as such.
"""

from __future__ import annotations

import json
from pathlib import Path

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np

REPO = Path(__file__).resolve().parents[1]
FIG = REPO / "figures"
FIG.mkdir(exist_ok=True)

plt.rcParams.update({"font.size": 10, "axes.spines.top": False,
                     "axes.spines.right": False})


def fig_benchmark():
    b = json.loads((REPO / "results/benchmark.json").read_text())
    models = ["M1", "M2", "M3", "M4", "M5"]
    auroc = [b["tiers"][m]["test"]["auroc"] for m in models]
    lo = [b["tiers"][m]["test"]["auroc_ci95"][0] for m in models]
    hi = [b["tiers"][m]["test"]["auroc_ci95"][1] for m in models]
    auprc = [b["tiers"][m]["test"]["auprc"] for m in models]

    fig, axes = plt.subplots(1, 2, figsize=(11, 4.2))
    x = np.arange(len(models))
    axes[0].bar(x, auroc, yerr=[np.array(auroc) - np.array(lo),
                                np.array(hi) - np.array(auroc)],
                capsize=4, color="#4C78A8")
    axes[0].set_xticks(x, models)
    axes[0].set_ylabel("Test AUROC [95% CI]")
    axes[0].set_ylim(0.75, 0.82)
    axes[0].set_title("Held-out test AUROC by tier")
    axes[1].bar(x, auprc, color="#72B66F")
    axes[1].set_xticks(x, models)
    axes[1].set_ylabel("Test AUPRC")
    axes[1].set_ylim(0.78, 0.87)
    axes[1].set_title("Held-out test AUPRC by tier")
    for ax in axes:
        for i, v in enumerate(auroc if ax == axes[0] else auprc):
            ax.text(i, v + 0.001, f"{v:.4f}", ha="center", fontsize=8)
    fig.suptitle("IsoformVUS held-out benchmark (24 genes, 14,400 variants, "
                 "gene-level splits)", fontsize=11)
    fig.tight_layout()
    fig.savefig(FIG / "04_benchmark.png", dpi=150)
    plt.close(fig)
    print("[15] 04_benchmark.png")


def fig_ablation():
    # Paired differences computed in 12b; hard-code from verified output with
    # source note (values from scripts/12b_paired_analysis.py run 2026-10-06).
    comps = [("M4-M1\n(+isoform struct.)", 0.0043, -0.0002, 0.0086),
             ("M5-M4\n(+usage)", -0.0000, -0.0002, 0.0001),
             ("M3-M1\n(+cell type)", -0.0021, -0.0048, 0.0007),
             ("M2-M1\n(+phenotype)", -0.0018, -0.0047, 0.0010)]
    labels = [c[0] for c in comps]
    means = [c[1] for c in comps]
    lo = [c[2] for c in comps]
    hi = [c[3] for c in comps]
    colors = ["#4C78A8" if l > 0 else "#E45756" for l in lo]

    fig, ax = plt.subplots(figsize=(9, 4.2))
    x = np.arange(len(labels))
    ax.bar(x, means, yerr=[np.array(means) - np.array(lo),
                           np.array(hi) - np.array(means)],
           capsize=5, color=colors)
    ax.axhline(0, color="black", linewidth=0.8)
    ax.set_xticks(x, labels)
    ax.set_ylabel("Paired AUROC difference [95% bootstrap CI]")
    ax.set_title("Ablation: does each biological layer add signal? (held-out test)")
    for i, (m, l, h) in enumerate(zip(means, lo, hi)):
        ax.text(i, h + 0.0004, f"{m:+.4f}", ha="center", fontsize=8)
    fig.tight_layout()
    fig.savefig(FIG / "05_ablation.png", dpi=150)
    plt.close(fig)
    print("[15] 05_ablation.png")


def fig_recovery():
    b = json.loads((REPO / "results/benchmark.json").read_text())
    rec = b["recovery_M5_vs_M1"]
    ks = [10, 20, 50]
    m5 = [rec[f"top{k}"]["n_new_topk_pathogenic"] for k in ks]
    m1 = [rec[f"top{k}"]["n_base_topk_pathogenic"] for k in ks]
    missed = [rec[f"top{k}"]["n_recovered_missed_by_baseline"] for k in ks]

    x = np.arange(len(ks))
    w = 0.35
    fig, ax = plt.subplots(figsize=(8, 4.2))
    ax.bar(x - w / 2, m1, w, label="M1 top-k pathogenic", color="#BAB0AC")
    ax.bar(x + w / 2, m5, w, label="M5 top-k pathogenic", color="#4C78A8")
    for i, (k, m) in enumerate(zip(ks, missed)):
        ax.text(i + w / 2, m5[i] + 0.6, f"+{m} missed\nby M1", ha="center",
                fontsize=8, color="#4C78A8")
    ax.set_xticks(x, [f"top-{k}" for k in ks])
    ax.set_ylabel("Pathogenic variants in top-k")
    ax.set_title("M5 recovers pathogenic variants that M1 misses (held-out test)")
    ax.legend()
    fig.tight_layout()
    fig.savefig(FIG / "06_topk_recovery.png", dpi=150)
    plt.close(fig)
    print("[15] 06_topk_recovery.png")


def fig_counterfactual():
    cf = json.loads((REPO / "results/counterfactual.json").read_text())
    effects = {}
    for p in cf["per_transcript"]:
        effects[p["effect"]] = effects.get(p["effect"], 0) + 1

    fig, axes = plt.subplots(1, 2, figsize=(11, 4.2))
    labels = list(effects)
    axes[0].bar(labels, [effects[l] for l in labels], color="#F58518")
    axes[0].set_ylabel("Transcripts")
    axes[0].set_title(f"Per-transcript effect of {cf['variant']}\n({cf['gene']}, "
                      f"{sum(effects.values())} GENCODE v46 transcripts)")
    axes[0].tick_params(axis="x", rotation=20)
    scores = [cf["score_assuming_irrelevant_isoform"],
              cf["score_observed_isoform_context"]]
    axes[1].bar(["Irrelevant\nisoform assumed", "Observed\nisoform context"],
                scores, color=["#BAB0AC", "#4C78A8"])
    axes[1].set_ylabel("M5 pathogenicity score")
    axes[1].set_ylim(0.9, 1.0)
    axes[1].set_title(f"Same DNA, different isoform assumption\n"
                      f"(Δ = {cf['delta']:+.4f})")
    for i, s in enumerate(scores):
        axes[1].text(i, s + 0.002, f"{s:.4f}", ha="center", fontsize=9)
    fig.tight_layout()
    fig.savefig(FIG / "07_counterfactual.png", dpi=150)
    plt.close(fig)
    print("[15] 07_counterfactual.png")


def main() -> int:
    fig_benchmark()
    fig_ablation()
    fig_recovery()
    fig_counterfactual()
    print("[15] done -> figures/")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
