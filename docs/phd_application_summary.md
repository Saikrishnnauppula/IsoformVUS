# IsoformVUS — PhD Application Research Summary

**Applicant:** Sai Krishna Uppula · Computational Genomics & Bioinformatics  
**Repository:** [github.com/Saikrishnnauppula/IsoformVUS](https://github.com/Saikrishnnauppula/IsoformVUS)

## Research Problem

Clinical exome and genome sequencing routinely return variants of uncertain significance (VUS). Standard variant-effect predictors (CADD, SpliceAI, AlphaMissense) score coding variants on a single canonical transcript without modeling how alternative isoform expression, internal exon length, and nonsense-mediated decay (NMD) interact across variant classes.

## Experimental Design & Self-Audit Discovery Progression

1. **Baseline Tiered Benchmark (Stages 01–17; 60,521 Variants Across 102 Mendelian Genes)**:
   - Built a leakage-safe, gene-stratified benchmark (61 train / 17 validation / 24 held-out test genes; n = 14,400 test variants) evaluating sequential logistic-regression models (`M1` variant scores → `M2` +HPO phenotype → `M3` +cell-type relevance → `M4` +isoform structural flags → `M5` +isoform usage).
   - On standard step AUPRC (`results/benchmark.json`), `M4` improved over `M1` from **0.8016 to 0.8524 (+0.0509, p < 0.001; AUROC 0.7985 → 0.8028)** and recovered 23 of the top-50 pathogenic variants missed by `M1`.
2. **Follow-Up Falsification of the M2/M3 Tie-Discretization Artifact (`breakthrough/round2/e1_e2`)**:
   - Rather than accepting the intermediate +0.031 step-AUPRC gain of `M2`/`M3` at face value, I designed gene-level shuffle controls and found that shuffling HPO and GTEx tables still produced the same +0.031 step-AUPRC jump.
   - Mathematical diagnosis revealed that 6,077 test SNVs with missing AlphaMissense scores were tied in `M1`; standard right-endpoint step AUPRC penalized that tie block by −0.0321. Implementing the exact analytical tie-corrected AUPRC estimator (McSherry & Najork 2008) showed `M1 = 0.8334`, `M2 = 0.8340`, `M3 = 0.8332`, and `M4 = 0.8517` (+0.0183), proving `M2`/`M3` were null and 94.5% of `M4`'s gain reflected canonical coding-sequence membership.
3. **Exon Architecture, Backup Skipping, and Co-Transcriptional NMD (`scripts/18_exon_architecture_nmd.py`)**:
   - Expanding to all **587,227 unselected ClinVar variants** across the 102 genes (4,132 reconstructed MANE Select exons), I identified three non-canonical isoform regimes when MANE Select is unhit: benign dual-coding frame shifts (n = 10,800, 99.93% benign), VUS-trapped alternative-exon missense variants (n = 2,173, 64.15% VUS, 89.2% benign when resolved), and expression-dependent alternative-exon loss-of-function variants (n = 877, 28.6%–35.0% resolved pathogenic).
   - Across **301,264 strictly internal canonical coding variants** (`in_last_exon = 0`), stratifying exons by length (≤180 bp vs. >1,000 bp), GTEx-expressed backup skipping (≥10% usage), and reading-frame symmetry (`L mod 3 == 0`) revealed a **62.5-fold reduction** in resolved missense pathogenicity (**73.13% → 1.17%**, within-gene CMH p = 4.29 × 10⁻¹⁸⁷ and 3.53 × 10⁻⁴², mediated 77.9%–94.7% by AlphaMissense 3D structural constraint), while internal frameshift and nonsense variants in those exact same exons remained **99.85%–100.00% pathogenic** due to nuclear co-transcriptional splicing into EJC-dependent NMD.
   - Resolving strand orientation showed that **88.2% (15/17)** of AlphaMissense's 0–3 bp exon-boundary false negatives (`am_score < 0.34`) occur at the **3' end of the exon (5'-splice donor consensus `MAG|GTRAGT`)**, where 65.22% (15/23) are pathogenic vs. 13.33% (2/15) at the 5' exon start (OR = 12.19, p = 0.0024).
   - Combining these mechanistic features (`M6`) increased held-out test AUROC to **0.9930** (tie-corrected AUPRC **0.9934**), calibrated Tavtigian et al. (2018) Bayesian likelihood ratios across 23,477 ClinVar VUS, and prioritized 304 internal exons harboring 10,910 pathogenic patient PTVs for antisense oligonucleotide (ASO) exon-skipping rescue.
