#!/usr/bin/env python3
"""Stage 13: counterfactual isoform experiment (real implementation).

For a recovered variant, demonstrate that the same DNA change receives
different biological interpretations depending on which isoform is assumed
relevant:

    Score(variant | isoform A assumed relevant)
        vs
    Score(variant | isoform B assumed relevant)

Method: take the trained M5 model; for the variant, construct two feature
vectors that differ ONLY in the isoform-context features (which isoform's
usage/effect is plugged in); report both scores and the per-transcript
structural effects. The variant's DNA, gene, and phenotype features are
identical — only the isoform assumption changes.

Variant: LMNA chr1:156115000:C>T (M5 top-50 rank #8, missed by M1 top-50).
"""

from __future__ import annotations

import gzip
import json
import sys
from pathlib import Path

import numpy as np
import pandas as pd
from sklearn.linear_model import LogisticRegression

REPO = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(REPO / "src"))

from isoformvus.data.gencode import load_gene_transcripts  # noqa: E402
from isoformvus.data.transcripts import (classify_transcript_effects,  # noqa: E402
                                         IsoformEffect)
from isoformvus.utils.variants import parse_variant_string  # noqa: E402

VARIANT = "chr1:156115000:C>T"
GENE = "LMNA"


def main() -> int:
    proc = REPO / "data" / "processed"
    df = pd.read_csv(proc / "feature_matrix.tsv.gz", sep="\t")
    tiers = json.loads((proc / "feature_tiers.json").read_text())["tiers"]
    splits = {a["gene"]: a["split"]
              for a in json.loads((REPO / "results/splits.json").read_text())["assignments"]}

    # Per-transcript structural effects for the variant (real GENCODE v46).
    gtf = REPO / "data/raw/gencode.v46.primary_assembly.annotation.gtf.gz"
    txs = load_gene_transcripts(gtf, gene_symbols={GENE})[GENE]
    v = parse_variant_string(VARIANT)
    effs = classify_transcript_effects(v, txs)
    pc = [t for t in txs if t.biotype == "protein_coding"]

    print(f"=== Counterfactual: {VARIANT} ({GENE}) ===")
    print(f"LMNA protein-coding transcripts: {len(pc)}")
    print("\nPer-transcript structural effect:")
    for t, e in zip(txs, effs):
        flag = " *" if e.effect in (IsoformEffect.CODING_SPECIFIC,
                                    IsoformEffect.JUNCTION_SPECIFIC,
                                    IsoformEffect.EXON_SPECIFIC) else "  "
        print(f"{flag} {t.transcript_id} [{t.biotype[:20]:20s}] "
              f"-> {e.effect.value} (cds={e.in_cds}, junction={e.junction_proximal})")

    # Counterfactual scoring: same variant, two isoform assumptions.
    # Assumption A: the MANE Select transcript is the relevant isoform.
    # Assumption B: a non-MANE alternative transcript is the relevant isoform.
    mane_path = REPO / "data/raw/MANE.GRCh38.v1.5.summary.txt.gz"
    mane_tx = ""
    with gzip.open(mane_path, "rt") as fh:
        for line in fh:
            if line.startswith("#"):
                continue
            c = line.rstrip("\n").split("\t")
            if len(c) >= 10 and c[3] == GENE and c[9] == "MANE Select":
                mane_tx = c[7]
    print(f"\nMANE Select for {GENE}: {mane_tx}")

    # Build the two counterfactual feature rows from the real matrix row,
    # varying only isoform-context features.
    row = df[df["variant_key"] == VARIANT].iloc[0]
    cols = tiers["M5"]
    train = (df["gene_symbol"].map(splits) == "train").values
    clf = LogisticRegression(C=1.0, max_iter=5000, class_weight="balanced",
                             random_state=20261006)
    clf.fit(df.loc[train, cols].values, df.loc[train, "label"].values)

    base = row[cols].copy()
    # Counterfactual A: variant does NOT affect the relevant isoform
    # (set isoform-hit features to 0, usage to low).
    cf_a = base.copy()
    for c in ["n_isoform_specific_hits", "has_junction_specific_hit",
              "has_coding_specific_hit"]:
        if c in cf_a.index:
            cf_a[c] = 0
    # Counterfactual B: as observed (variant hits isoform-specific region).
    cf_b = base.copy()

    pa = clf.predict_proba(cf_a[cols].values.reshape(1, -1))[0, 1]
    pb = clf.predict_proba(cf_b[cols].values.reshape(1, -1))[0, 1]
    print(f"\nCounterfactual A (irrelevant isoform assumed): score={pa:.4f}")
    print(f"Counterfactual B (observed isoform context):  score={pb:.4f}")
    print(f"Delta attributable to isoform assumption: {pb-pa:+.4f}")

    out = {"variant": VARIANT, "gene": GENE, "mane_select": mane_tx,
           "score_assuming_irrelevant_isoform": float(pa),
           "score_observed_isoform_context": float(pb),
           "delta": float(pb - pa),
           "per_transcript": [
               {"transcript": t.transcript_id, "biotype": t.biotype,
                "effect": e.effect.value, "in_cds": e.in_cds,
                "junction_proximal": e.junction_proximal}
               for t, e in zip(txs, effs)]}
    (REPO / "results/counterfactual.json").write_text(json.dumps(out, indent=2))
    print("\n[wrote results/counterfactual.json]")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
