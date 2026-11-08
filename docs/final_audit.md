# Final Verification & Reproducibility Audit

**Author & Reviewer:** Sai Krishna Uppula  
**Last Verified:** 2026-11-08  

## 1. Benchmark & Statistical Integrity
- [x] **Locked Gene-Level Splits (`results/splits.json`)**: 61 train genes (36,600 variants), 17 validation genes (9,521 variants), and 24 test genes (14,400 variants). Zero shared genes or exons across splits.
- [x] **Exact Reproduction of `results/benchmark.json`**:
  - `M1` (variant-only): Test AUROC = **0.7985**, step AUPRC = **0.8016**, tie-corrected AUPRC = **0.8334**.
  - `M2` (+HPO phenotype): Test AUROC = **0.7967**, step AUPRC = **0.8327**, tie-corrected AUPRC = **0.8340**.
  - `M3` (+cell-type expression): Test AUROC = **0.7953**, step AUPRC = **0.8318**, tie-corrected AUPRC = **0.8332**.
  - `M4` (+isoform structural flags): Test AUROC = **0.8028**, step AUPRC = **0.8524** (ΔAUPRC vs. M1 = **+0.0509**, 95% CI [+0.0473, +0.0546]), tie-corrected AUPRC = **0.8517** (+0.0183).
  - `M5` (+dominant isoform usage): Test AUROC = **0.8027**, step AUPRC = **0.8523**, tie-corrected AUPRC = **0.8516**.
  - `M6` (+exon architecture, NMD boundary, and 5'-donor edge): Test AUROC = **0.9930**, tie-corrected AUPRC = **0.9934**, Brier = **0.0307**.
- [x] **Resolution of the M2/M3 Step-AUPRC Artifact (`breakthrough/round2/e1_e2`)**: Documented across `README.md`, `paper/manuscript.md`, and `paper/FINAL_REPORT.md` that the +0.031 step-AUPRC jump in `M2`/`M3` was caused by breaking a 6,077-SNV tie block in `M1` under right-endpoint step evaluation.

## 2. Stage 18 Exon Architecture & NMD Pipeline Verification
- [x] **Standalone Execution (`scripts/18_exon_architecture_nmd.py`)**: Uses only repository-relative paths (`REPO = Path(__file__).resolve().parents[1]`), supports both compressed (`.tsv.gz`) and uncompressed (`.tsv`) ClinVar 587k tables, and is invoked in `scripts/run_all.sh`.
- [x] **Gene-Stratified & Label-Independent Controls (`results/exon_architecture_nmd.json`)**:
  - Within-gene Cochran–Mantel–Haenszel tests confirm protective effects of exon length >400 bp (CMH OR = 0.2005, p = 4.29 × 10⁻¹⁸⁷) and GTEx backup exon skipping ≥10% (CMH OR = 0.3461, p = 3.53 × 10⁻⁴²) for internal missense variants, contrasted with 99.85%–100.00% lethality for internal PTVs.
  - Label-independent AlphaMissense median drops 9.1-fold (0.943 → 0.104) and benign-to-pathogenic missense density per kilobase shifts 63.4-fold (0.38× → 24.11×).
  - Leave-top-genes-out jackknife (excluding *BRCA1*, *BRCA2*, *APC*, *TTN*, *NF1*), ≥2-star ClinVar filtering, CpG transition controls, and post-2024 prospective variant holdouts all replicate the primary findings.
