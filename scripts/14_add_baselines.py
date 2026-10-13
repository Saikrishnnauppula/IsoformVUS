#!/usr/bin/env python3
"""Add SpliceAI + CADD baselines to the held-out benchmark.

Reads:
    data/processed/variant_spliceai.tsv.gz  (variant_key, spliceai_max_delta)
    data/processed/variant_cadd.tsv.gz      (variant_key, phred)
Joins to the locked test set, computes the full metric panel for each
baseline on its scored subset, and appends to results/benchmark.json +
results/benchmark.md. Also compares each baseline head-to-head with M5
on the same subset (paired AUROC difference).

Refuses to run until BOTH enrichment files exist.
"""

from __future__ import annotations

import gzip
import json
import sys
from pathlib import Path

import numpy as np
import pandas as pd
from sklearn.metrics import roc_auc_score

REPO = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(REPO / "src"))
from isoformvus.evaluation.metrics import full_report  # noqa: E402

ACCESS_LOG = REPO / "results" / "test_access.log"


def main() -> int:
    proc = REPO / "data" / "processed"
    sp_path = proc / "variant_spliceai.tsv.gz"
    cadd_path = proc / "variant_cadd.tsv.gz"
    missing = [p for p in (sp_path, cadd_path) if not p.exists()]
    if missing:
        print(f"[baselines] waiting on: {[p.name for p in missing]}",
              file=sys.stderr)
        return 2

    from datetime import datetime
    with open(ACCESS_LOG, "a") as fh:
        fh.write(f"{datetime.now().isoformat()} baseline enrichment "
                 "evaluation (spliceai+cadd)\n")

    df = pd.read_csv(proc / "feature_matrix.tsv.gz", sep="\t")
    tiers = json.loads((proc / "feature_tiers.json").read_text())["tiers"]
    splits = {a["gene"]: a["split"]
              for a in json.loads((REPO / "results/splits.json").read_text())["assignments"]}
    df["split"] = df["gene_symbol"].map(splits)
    test = df["split"] == "test"
    y = df.loc[test, "label"].astype(int).values
    dft = df[test].reset_index(drop=True)

    # M5 scores for head-to-head (same protocol as benchmark).
    from sklearn.linear_model import LogisticRegression
    train = df["split"] == "train"
    cols = tiers["M5"]
    clf = LogisticRegression(C=1.0, max_iter=5000, class_weight="balanced",
                             random_state=20261006)
    clf.fit(df.loc[train, cols].values, df.loc[train, "label"].values)
    m5 = clf.predict_proba(dft[cols].values)[:, 1]
    # Carry labels + M5 scores as columns so post-merge alignment is exact
    # (merge resets the index; positional indexing would misalign).
    dft["_y"] = y
    dft["_m5"] = m5

    bench = json.loads((REPO / "results/benchmark.json").read_text())
    lines_extra = []

    # --- SpliceAI ---
    sp = pd.read_csv(sp_path, sep="\t")
    m = dft.merge(sp, on="variant_key", how="inner")
    if len(m):
        yl = m["_y"].values.astype(int)
        m5_sub = m["_m5"].values
        s = m["spliceai_max_delta"].values
        rep = full_report(yl, s)
        rep["n_subset"] = int(len(m))
        # Head-to-head vs M5 on the same subset.
        d = (roc_auc_score(yl, m5_sub)
             - roc_auc_score(yl, s))
        rep["m5_minus_baseline_auroc"] = float(d)
        bench["baselines"]["spliceai_max_delta"] = rep
        lines_extra.append(
            f"- SpliceAI max delta (n={len(m)}): AUROC={rep['auroc']:.4f} "
            f"AUPRC={rep['auprc']:.4f} | M5 minus SpliceAI AUROC={d:+.4f}")
        print(f"[baselines] SpliceAI n={len(m)} AUROC={rep['auroc']:.4f} "
              f"(M5-{d:+.4f})")

    # --- CADD ---
    cd = pd.read_csv(cadd_path, sep="\t")
    # Feature matrix already carries stale partial-CADD columns; rename the
    # fresh enrichment columns to avoid _x/_y collision on merge.
    cd = cd.rename(columns={"cadd_phred": "cadd_phred_enrich",
                            "cadd_raw": "cadd_raw_enrich"})
    cd = cd[cd["cadd_phred_enrich"].notna()]
    m2 = dft.merge(cd, on="variant_key", how="inner")
    if len(m2) > 100:
        from scipy.stats import rankdata
        yl2 = m2["_y"].values.astype(int)
        m5_sub2 = m2["_m5"].values
        s_raw = m2["cadd_phred_enrich"].values
        # PHRED is not a probability: rank-uniformize to [0,1] so the full
        # metric panel (incl. Brier) is valid. Rank-based AUROC/AUPRC unchanged.
        s = rankdata(s_raw) / len(s_raw)
        rep = full_report(yl2, s)
        rep["n_subset"] = int(len(m2))
        d = (roc_auc_score(yl2, m5_sub2)
             - roc_auc_score(yl2, s))
        rep["m5_minus_baseline_auroc"] = float(d)
        bench["baselines"]["cadd_phred"] = rep
        lines_extra.append(
            f"- CADD PHRED (n={len(m2)}): AUROC={rep['auroc']:.4f} "
            f"AUPRC={rep['auprc']:.4f} | M5 minus CADD AUROC={d:+.4f}")
        print(f"[baselines] CADD n={len(m2)} AUROC={rep['auroc']:.4f} "
              f"(M5-{d:+.4f})")

    (REPO / "results/benchmark.json").write_text(json.dumps(bench, indent=2))
    with open(REPO / "results/benchmark.md", "a") as fh:
        fh.write("\n## Baselines added post-hoc (enrichment completed)\n\n")
        fh.write("\n".join(lines_extra) + "\n")
    print("[baselines] benchmark.json + benchmark.md updated")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
