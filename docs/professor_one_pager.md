# IsoformVUS — one-pager for professors

**Candidate:** Faizan Hussain Syed · CS + ML background · Fall 2027 Bioinformatics/Computational Genomics PhD applicant

## One sentence

IsoformVUS tests whether knowing *which isoform the disease-relevant cell
type uses* improves computational prioritization of rare-disease VUSs —
extending phenotype-aware models (IMPPROVE-class) with the isoform layer
they omit.

## Why a professor should care

- **Real clinical gap.** VUS triage is the bottleneck in rare-disease
  sequencing; current scorers (CADD/SpliceAI/AlphaMissense) are
  phenotype- and cell-type-blind.
- **Biology is on our side.** Alternative splicing is pervasive, cell-type
  isoform preferences are documented, and published cases (DNM1 exon 10a;
  TRNT1 isoform-selective alleles) prove single-transcript pipelines miss
  diagnoses. The question is whether this converts into *measurable*
  prioritization gain — asked here as a controlled experiment, not an
  assumption.
- **Clean experimental design.** Tiered ablation (variant → +phenotype →
  +cell type → +isoform → +isoform effect) isolates exactly what each
  biological layer contributes; gene-level splits prevent the standard
  leakage; the headline metric is *pathogenic variants recovered that
  baselines miss*, with CIs.

## What's built (repo: IsoformVUS/)

Literature review (24 methods, 2024–2026) + conservative novelty statement;
~100-gene panel design; phenotype→cell-type relevance module;
cell-type→isoform usage module with uncertainty handling; per-transcript
variant-effect classifier; 5-tier interpretable modeling scaffold
(logreg → LightGBM); full metric panel incl. top-k recovery; 31 passing
tests; provenance-tracked public data only (ClinVar, gnomAD, GENCODE v46,
HPO, GTEx, CELLxGENE, SpliceAI/CADD/AlphaMissense).

## Intellectual contribution

Not a bigger black box: the smallest testable version of the question,
with the null result pre-registered as publishable. The design directly
answers IMPPROVE (2025) — the closest published method — by running the
ablation it never ran.

## Ask

Feedback on the experimental design (especially the isoform-usage
evidence model and the leakage strategy), and whether the lab's
disease area offers a validation cohort or long-read data for the
isoform-selective subgroup.

## Status

Benchmark complete on the locked held-out test set: isoform-structure
features give +0.05 AUPRC and recover 23/50 top-ranked pathogenic variants
missed by the variant-only model (AUROC +0.0043, n.s.); phenotype/cell-type
layers null as implemented — honest nulls with a diagnosed cause
(tissue-level resolution). Design docs, code, 37 tests, and figures are
complete and reproducible; CADD baseline enrichment still running.
