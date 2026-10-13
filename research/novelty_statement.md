# IsoformVUS — Novelty Statement

*(One page. Conservative, committee-grade. Literature as of 2026-10-06. Full citations and comparison table: `literature_review.md`.)*

## The claim

IsoformVUS is, to the best of our knowledge, the first **interpretable, benchmarked variant-prioritization model that jointly combines HPO phenotype conditioning, single-cell-derived cell-type relevance, cell-type-specific isoform usage, and isoform-specific functional consequence** — trained on ClinVar labels with gene-level splits, evaluated with an ablation that isolates the isoform layer's contribution, and benchmarked against CADD, SpliceAI, AlphaMissense, Exomiser, and IMPPROVE.

## What we do NOT claim

We do not claim that nobody has considered isoforms in variant interpretation. The biology is established: Ly et al. (*Molecular Cell*, 2025) demonstrated isoform-selective alleles via alternative start codons (TRNT1/SIFD); Weykopf et al. (*Nature Communications*, 2026) showed ~80% of alternative transcripts carrying disease variants are absent from reference annotations (DPP9 lung-isoform case); AlphaMissense (2023) and ESM1b (Brandes et al., 2023) already score variants in alternative-isoform contexts; and the AlphaGenome Atlas preprint (2026) rescued a real DNM1 noncoding VUS through a brain-specific isoform. We also do not claim novelty for phenotype-aware prioritization (Exomiser, 2015; LIRICAL, 2022; V2P, 2025) or for adding single-cell expression to prioritization models (VARPP, 2019; **IMPPROVE**, Chapman & Lassmann, 2025).

## The specific gap

No published method occupies the intersection of all four layers in a single trained, benchmarked model:

- **IMPPROVE (2025)**, our closest competitor, shares the phenotype + cell-type + interpretable-model design — but its expression features are **gene-level**. It cannot distinguish a variant hitting a brain-dominant isoform from one hitting a blood-dominant isoform of the same gene. Our ablation (cell-type model vs cell-type + isoform model) tests exactly the value IMPPROVE leaves unmeasured.
- **AlphaGenome Atlas/AVI (2026 preprint)** has biosample-aware splicing and demonstrated an isoform-selective VUS rescue — but it is phenotype-free, trained on allele-frequency proxy labels rather than disease labels, and performs no with/without-isoform ablation. It is a molecular-effects engine; IsoformVUS is the phenotype-conditioned clinical-prioritization layer.
- **Weykopf et al. (2026)** and **Ly et al. (2025)** establish the biology but build no prioritization model; **tissue-specific splicing predictors** (AbSplice, MTSplice, SpliceTransformer) and **CI-SpliceAI** predict splicing effects without phenotype, cell-type, or disease-label training.

## Why the experiment still matters if the hypothesis is weak

The literature bounds the expected effect size: an estimated ~2.5% of pathogenic/likely-pathogenic variants map exclusively to MANE Plus Clinical-specific sequence (preprint, Exeter group) — so the isoform-selective subgroup is small but real. Our design therefore evaluates this subgroup separately rather than relying on genome-wide AUROC, and the ablation is pre-registered to report a null result honestly if the isoform layer adds nothing. Either outcome is publishable: a positive result extends IMPPROVE-class models one biological layer deeper; a null result is the first direct experimental test of whether isoform conditioning improves clinical variant prioritization.

## Verification status

Two citations need final confirmation before manuscript submission (Weykopf et al. exact DOI; V2P author list); the MANE Plus Clinical estimate is from a preprint; CI-SpliceAI is cited as an open-source tool pending paper verification. Three directly relevant works appeared within nine months of this review — the novelty check must be re-run at submission time.
