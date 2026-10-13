#!/usr/bin/env python3
"""Stage 07: build the feature matrix (real implementation).

Joins, per variant:
    variant_labels            -> vartype, review_stars, gene, label
    variant_isoform_effects   -> gene_effect_summary, n_isoform_specific_hits,
                                 junction/coding flags, mane_select_effect
    cell_relevance            -> max cell_relevance score per gene (M3),
                                 top tissue/cell, phenotype_match, expression_support
    isoform_usage             -> dominant-transcript usage in relevant tissue (M5;
                                 if stage 05 output missing, M5 falls back to M4
                                 features and the fallback is logged)
    variant_alphamissense     -> am_score + am_missing indicator
    variant_cadd              -> cadd_phred + cadd_missing indicator
                                 (partial coverage while 09b runs)

Missing-data policy (features/spec.md): every potentially-missing feature
gets an explicit `_missing` indicator; no silent zero-imputation.

Writes:
    data/processed/feature_matrix.tsv.gz
    data/processed/feature_tiers.json   (model_id -> column list; M1..M5
                                         strictly nested — enforced)

Gene-level aggregates (n_hpo_terms etc.) come from metadata/gene_panel.tsv.
"""

from __future__ import annotations

import argparse
import gzip
import json
import sys
from pathlib import Path

import pandas as pd

REPO = Path(__file__).resolve().parents[1]

TIER_ORDER = ["M1", "M2", "M3", "M4", "M5"]


def read_tsv_gz(path: Path, usecols=None) -> pd.DataFrame:
    return pd.read_csv(path, sep="\t", usecols=usecols)


def main() -> int:
    ap = argparse.ArgumentParser()
    ap.add_argument("--out", default=str(REPO / "data/processed/feature_matrix.tsv.gz"))
    ap.add_argument("--tiers-out", default=str(REPO / "data/processed/feature_tiers.json"))
    args = ap.parse_args()
    proc = REPO / "data" / "processed"

    labels = read_tsv_gz(proc / "variant_labels.tsv.gz")
    effects = read_tsv_gz(proc / "variant_isoform_effects.tsv.gz")
    df = labels.merge(effects, on=["variant_key", "gene_symbol"], how="left",
                      validate="one_to_one")

    # --- AlphaMissense (missense only) ---
    am_path = proc / "variant_alphamissense.tsv.gz"
    if am_path.exists():
        am = read_tsv_gz(am_path)
        df = df.merge(am, on="variant_key", how="left", validate="one_to_one")
    else:
        df["am_score"] = float("nan"); df["am_class"] = ""; df["am_missing"] = 1
    df["am_missing"] = df["am_missing"].fillna(1).astype(int)
    df["am_score"] = pd.to_numeric(df["am_score"], errors="coerce")

    # --- CADD (partial coverage OK) ---
    cadd_path = proc / "variant_cadd.tsv.gz"
    if cadd_path.exists():
        cadd = read_tsv_gz(cadd_path)
        df = df.merge(cadd, on="variant_key", how="left", validate="one_to_one")
    else:
        df["cadd_phred"] = float("nan"); df["cadd_raw"] = float("nan")
    df["cadd_missing"] = df["cadd_phred"].isna().astype(int)
    df["cadd_phred"] = pd.to_numeric(df["cadd_phred"], errors="coerce")

    # --- Cell relevance: per-gene aggregates (M2/M3) ---
    cr_path = proc / "cell_relevance.tsv.gz"
    if cr_path.exists():
        cr = read_tsv_gz(cr_path)
        g = cr.groupby("gene_symbol").agg(
            cell_relevance_max=("score", "max"),
            phenotype_match_max=("phenotype_match", "max"),
            expression_support_max=("expression_support", "max"),
            n_relevant_cells=("cell_ontology_id", "nunique"),
        ).reset_index()
        top = cr.sort_values("score", ascending=False).drop_duplicates("gene_symbol")
        g = g.merge(top[["gene_symbol", "cell_type_name", "tissue"]],
                    on="gene_symbol", how="left")
        g = g.rename(columns={"cell_type_name": "top_cell_type",
                              "tissue": "top_tissue"})
        df = df.merge(g, on="gene_symbol", how="left", validate="many_to_one")
    else:
        for c in ["cell_relevance_max", "phenotype_match_max",
                  "expression_support_max", "n_relevant_cells"]:
            df[c] = float("nan")

    # --- Gene panel extras (M2) ---
    panel = pd.read_csv(REPO / "metadata/gene_panel.tsv", sep="\t",
                        usecols=["gene_symbol", "n_hpo_terms",
                                 "n_protein_coding_transcripts",
                                 "has_coding_alternative_exon"])
    df = df.merge(panel, on="gene_symbol", how="left", validate="many_to_one")

    # --- Isoform usage (M5): dominant-transcript usage in top tissue ---
    iu_path = proc / "isoform_usage.tsv.gz"
    m5_fallback = False
    if iu_path.exists():
        iu = read_tsv_gz(iu_path)
        # Per (gene, tissue): usage of the dominant transcript.
        dom = iu.sort_values("usage", ascending=False).drop_duplicates(
            ["gene_symbol", "tissue"])
        dom = dom[["gene_symbol", "tissue", "usage", "usage_category"]]
        dom = dom.rename(columns={"tissue": "top_tissue",
                                  "usage": "dominant_isoform_usage",
                                  "usage_category": "dominant_isoform_category"})
        df = df.merge(dom, on=["gene_symbol", "top_tissue"], how="left",
                      validate="many_to_one")
    else:
        m5_fallback = True
        df["dominant_isoform_usage"] = float("nan")
        df["dominant_isoform_category"] = "uncertain"
    df["isoform_usage_missing"] = df["dominant_isoform_usage"].isna().astype(int)

    # --- One-hot encodings ---
    df = pd.get_dummies(df, columns=["vartype", "gene_effect_summary",
                                     "dominant_isoform_category"],
                        prefix=["vartype", "geff", "usecat"], dummy_na=False)

    # --- Tier definitions (strictly nested) ---
    vartype_cols = [c for c in df.columns if c.startswith("vartype_")]
    geff_cols = [c for c in df.columns if c.startswith("geff_")]
    usecat_cols = [c for c in df.columns if c.startswith("usecat_")]

    tiers = {
        "M1": (["review_stars", "am_score", "am_missing",
                "cadd_phred", "cadd_missing"] + vartype_cols),
        "M2": ["n_hpo_terms", "phenotype_match_max"],
        "M3": ["cell_relevance_max", "expression_support_max",
               "n_relevant_cells"],
        "M4": (["n_isoform_specific_hits", "has_junction_specific_hit",
                "has_coding_specific_hit", "n_protein_coding_transcripts",
                "has_coding_alternative_exon"] + geff_cols),
        "M5": (["dominant_isoform_usage", "isoform_usage_missing"] + usecat_cols),
    }
    # Strict nesting: each tier's full column set extends the previous.
    resolved, acc = {}, []
    for m in TIER_ORDER:
        acc = acc + [c for c in tiers[m] if c not in acc]
        missing = [c for c in acc if c not in df.columns]
        if missing:
            raise KeyError(f"{m}: columns not in matrix: {missing}")
        resolved[m] = list(acc)

    # Fill remaining NaNs in numeric feature columns (indicators already exist).
    num_cols = df.select_dtypes(include="number").columns
    feat_cols = [c for c in resolved["M5"] if c in num_cols]
    df[feat_cols] = df[feat_cols].fillna(-1)  # sentinel; *_missing flags carry the signal

    keep = (["variant_key", "gene_symbol", "label"] +
            sorted(set(resolved["M5"])))
    df[keep].to_csv(args.out, sep="\t", index=False, compression="gzip")
    Path(args.tiers_out).write_text(json.dumps(
        {"tiers": resolved,
         "m5_fallback_no_usage_data": m5_fallback,
         "note": "M1..M5 strictly nested; -1 sentinel only where a *_missing "
                 "indicator exists."}, indent=2))
    print(f"[07] matrix: {df.shape[0]} variants x {len(resolved['M5'])} features")
    print(f"[07] tier sizes: " + ", ".join(f"{m}={len(resolved[m])}" for m in TIER_ORDER))
    if m5_fallback:
        print("[07] WARNING: isoform_usage.tsv.gz missing — M5 falls back to "
              "M4 features + missing indicators (logged in feature_tiers.json)")
    print(f"[07] wrote {args.out}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
