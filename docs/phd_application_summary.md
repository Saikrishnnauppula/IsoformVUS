# IsoformVUS — PhD application summary

*For Bioinformatics / Computational Genomics PhD applications (Fall 2027).*

## The problem

Exome sequencing for rare disease routinely returns **variants of uncertain
significance (VUSs)** — rare variants that cannot currently be called
pathogenic or benign. Computational prioritization (CADD, SpliceAI,
AlphaMissense) triages these candidates, but scores them from sequence alone,
agnostic to the patient's phenotype and to the cell type where disease
manifests. Meanwhile, ~95% of human multi-exon genes are alternatively
spliced, and cell types use different isoforms — so a variant's consequence
can depend on which isoform the disease-relevant cell expresses. That
information is systematically unused by current prioritizers.

## Why it matters

A variant in an exon skipped by the brain's dominant isoform is silent in
brain no matter how "deleterious" its sequence score; conversely, a
moderate-scoring variant in a brain-specific exon can be devastating. Published
evidence bounds the scale: ~2.5% of pathogenic/likely-pathogenic ClinVar
variants map exclusively to non-MANE sequence (MANE Plus Clinical preprint),
and case reports (DNM1 brain-specific exon 10a, rescued via AlphaGenome
Atlas; TRNT1 isoform-selective alleles, Ly et al. 2025) show single-transcript
pipelines miss real diagnoses. Even a small recovery rate matters when each
recovery is a diagnosed patient.

## Research hypothesis

Incorporating disease-relevant cell-type information and isoform usage
improves prioritization of pathogenic rare variants over conventional
variant-effect predictors — specifically by recovering pathogenic variants
that gene-level or tissue-agnostic approaches rank poorly. Tested
experimentally, with a pre-registered honest null.

## Methodology

1. **Literature-first**: 24-method review (2024–2026) identifying IMPPROVE
   (2025) as closest competitor — phenotype + single-cell expression but
   gene-level, no isoform layer — and a conservative novelty statement.
2. **~100 Mendelian gene panel** with documented inclusion criteria (ClinVar
   P/LP counts, benign controls, multi-isoform genes, expression + cell-type
   data availability).
3. **Phenotype → cell-type relevance pipeline**: HPO terms → tissues → cell
   types via measurable evidence (not hard-coded), producing a
   cell-relevance score.
4. **Cell-type → isoform usage**: categorical usage (dominant / supported /
   low / uncertain) with uncertainty widths; short-read ambiguity handled
   explicitly, long-read evidence tracked.
5. **Isoform-specific variant-effect classification** per transcript, with
   SpliceAI re-interpreted in isoform context.
6. **Tiered interpretable models** (logistic regression → LightGBM): M1
   variant-only … M5 all + isoform-effect; the critical M4-vs-M3 comparison
   isolates the isoform layer's contribution on identical data and splits.
7. **Leakage prevention**: gene-level grouped splits (never variant-level),
   locked held-out test, paralog sensitivity analysis.
8. **Evaluation**: AUROC/AUPRC/top-k recall with bootstrap CIs; headlined by
   *recovered pathogenic variants that baselines miss*; separate
   isoform-selective subgroup analysis; counterfactual isoform experiments;
   SHAP explainability.

## Datasets (all public)

ClinVar (labels), gnomAD (frequencies/controls), GENCODE v46 + MANE
(transcripts), HPO (phenotypes), GTEx (tissue expression), CELLxGENE/HCA
(single-cell), SpliceAI / CADD / AlphaMissense (features + baselines).
Full provenance (versions, builds, licenses, dates) in
`data/DATA_PROVENANCE.md`.

## Results

Held-out test set (24 genes, 14,400 variants, gene-level splits — test set
never touched during development):

| Model | Test AUROC [95% CI] | Test AUPRC |
|---|---|---|
| M1 variant-only | 0.7985 [0.7917, 0.8061] | 0.8016 |
| M2 +phenotype | 0.7967 | 0.8327 |
| M3 +cell type | 0.7964 | 0.8330 |
| M4 +isoform structure | 0.8029 [0.7955, 0.8105] | 0.8524 |
| M5 +isoform usage | 0.8029 | 0.8523 |

Paired bootstrap: M4−M1 = +0.0043 [−0.0002, +0.0086] (n.s. on AUROC);
phenotype, cell-type, and usage layers null as implemented. Headline
finding: M5's top-50 test variants (all pathogenic) include **23 missed by
M1's top-50** — isoform features change *which* pathogenic variants get
prioritized. Counterfactual: same LMNA variant scores +0.042 higher under
its observed isoform context vs an irrelevant-isoform assumption.
AlphaMissense baseline (missense-only, n=1,441): AUROC 0.9435 — stronger
but circular (ClinVar-trained) and variant-class-restricted.

Honest verdict: isoform structure carries modest, real signal (strongest
in ranking behavior); phenotype/cell-type layers need single-cell isoform
resolution to be fairly tested — the precise next experiment.

## Comparison with existing approaches

- **vs CADD/SpliceAI/AlphaMissense**: they are features and baselines, not
  competitors to beat on their own terms; IsoformVUS adds the
  phenotype/cell-type/isoform conditioning they omit.
- **vs Exomiser/LIRICAL/V2P**: phenotype-aware but gene-level; no
  cell-type or isoform modeling.
- **vs IMPPROVE**: the direct extension — same phenotype + single-cell +
  interpretable design, plus the isoform-usage and isoform-specific
  consequence layers, with the ablation IMPPROVE never ran.

## Limitations (selection)

Transcript annotation incompleteness (~80% of disease-variant-bearing
alternative transcripts absent from references); short-read isoform
ambiguity (handled via uncertainty categories); ClinVar ascertainment bias;
~100-gene scope; single-cell coverage gaps. Full list with bias directions:
`docs/limitations.md`.

## Future research

Single-cell long-read isoform quantification (Isosceles/FLAMES) as ground
truth; developmental-stage-specific isoform usage; extending beyond
Mendelian germline to de novo neurodevelopmental cohorts; prospective
validation with clinical collaborators — ranking held-out VUSs for
functional follow-up.
