# Prototype tiered experiment — validation results (2026-10-06)

**Status:** prototype on VALIDATION genes; the held-out TEST set (14,400
variants, 24 genes) has not been touched. Full tiered models with
CADD/SpliceAI/AlphaMissense features are pending (stages 07/10).

## Setup

- Train: 35,921 variants (61 genes) · Validation: 10,200 variants (17 genes)
- Learner: L2 logistic regression (C=1.0, class_weight=balanced, seed 20261006)
- M1 (variant-only): vartype one-hot + ClinVar review stars
- M4 (+isoform): M1 + gene_effect_summary one-hot + n_isoform_specific_hits
  + has_junction_specific_hit + has_coding_specific_hit (+ mane_affected)

## Results (validation)

| Model | AUROC | AUPRC |
|---|---|---|
| M1 variant-only | 0.7484 | 0.7403 |
| M4 +isoform | **0.7814** | **0.8152** |
| M4-noMANE +isoform | 0.7813 | 0.8174 |

Delta (M4 − M1): **+0.033 AUROC, +0.075 AUPRC**.

## Ascertainment check

MANE Select transcripts are chosen partly for clinical relevance, so
`mane_affected` could leak ascertainment. Excluding all MANE-derived
features changes AUROC by −0.0001 — the gain is driven by structural
features (junction-specific hits, coding-specific hits), not MANE.

## Top M4 coefficients (log-odds, validation fit)

| Feature | Coef | Reading |
|---|---|---|
| vartype_SNV | −2.18 | SNVs less likely P than indels here (frameshift enrichment among P indels) |
| has_junction_specific_hit | +1.41 | splice-junction hits in a subset of isoforms → pathogenic |
| mane_affected | +1.30 | (ascertainment-checked: contributes ~nothing — see above) |
| has_coding_specific_hit | +1.05 | coding hits unique to isoform subsets → pathogenic |
| vartype_deletion | +0.88 | |

## Honest caveats

1. **Weak M1 baseline.** vartype + review stars is deliberately minimal; the
   isoform delta may shrink once CADD/SpliceAI/AlphaMissense join M1.
   The critical test is the full stage-10 comparison, not this prototype.
2. **review_stars as a feature** is a label-quality proxy that correlates
   with ascertainment; it is disclosed, and stage 10 will test dropping it.
3. **Structural, not yet cell-type-aware.** These isoform features use no
   cell-type or phenotype information — the M2/M3/M5 tiers are pending.
4. **Validation only.** Test-set numbers come from stage 11; these
   validation numbers must not be presented as final.
