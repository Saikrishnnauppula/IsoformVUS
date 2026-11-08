# IsoformVUS — Research Overview for Faculty

**Candidate:** Sai Krishna Uppula · Computational Genomics & Bioinformatics  
**Repository:** [github.com/Saikrishnnauppula/IsoformVUS](https://github.com/Saikrishnnauppula/IsoformVUS)

## Core Question

How do alternative transcript structure, internal exon architecture, and nonsense-mediated decay (NMD) interact to determine the clinical pathogenicity of rare coding and splice-region variants across Mendelian disease genes?

## Summary of Contributions

1. **Rigorous Negative-Control Falsification of Tie-Breaking Artifacts**:
   - On a gene-stratified benchmark of 60,521 ClinVar variants across 102 Mendelian genes (24 held-out test genes, n = 14,400), adding isoform structural features (`M4`) improved step AUPRC from **0.8016 to 0.8524 (+0.0509)**.
   - Using permutation and analytical tie-corrected AUPRC (McSherry & Najork 2008), I demonstrated that an apparent +0.031 step-AUPRC gain in intermediate HPO/GTEx gene-level tiers (`M2`/`M3`) was an artifact of breaking a 6,077-SNV tie block in `M1` (`0.8334` tie-corrected), and that 94.5% of `M4`'s gain over `M3` was driven by a 1-bit canonical coding indicator.
2. **Internal Exon-Architecture × Backup-Skipping × NMD Decoupling Across 587,227 ClinVar Variants**:
   - Reconstructing all 4,132 MANE Select exons (3,955 internal coding exons) revealed that while resolved missense pathogenicity drops **62.5-fold** from **73.13%** in compact (≤180 bp), constitutive, phase-symmetric exons to **1.17%** in giant (>1,000 bp), GTEx backup-skipped (≥10% usage), symmetric exons (within-gene CMH p = 4.29 × 10⁻¹⁸⁷), internal nonsense and frameshift variants in those same dispensable exons remain **99.85%–100.00% pathogenic** because pre-mRNA splicing commits the transcript to NMD in the nucleus before translation.
3. **U1 snRNP 5'-Splice Donor Asymmetry in AlphaMissense Exon-Edge False Negatives**:
   - Because AlphaMissense takes only the spliced amino-acid sequence as input, it misses splice-disrupting missense variants within 3 bp of exon junctions (44.74% false-negative rate vs. 15.23% in the exon core, p = 4.9 × 10⁻⁵). Separating strand orientation proved that **88.2% (15/17)** of these false negatives sit at the **3' end of the exon (5'-splice donor `MAG|GTRAGT`)**, where **65.22% (15/23)** are pathogenic versus **13.33% (2/15)** at the 5' exon start (OR = 12.19, p = 0.0024).
4. **Translational Outputs (`results/`)**:
   - Integrating these features (`M6`) reaches **0.9930 AUROC / 0.9934 tie-corrected AUPRC** on the 24 unseen test genes, provides Bayesian ACMG/AMP likelihood ratios across 23,477 ClinVar VUS (`results/vus_bayesian_stratification.tsv`), and catalogs 304 phase-symmetric/backup-skipped internal exons containing 10,910 pathogenic patient PTVs amenable to antisense exon skipping (`results/exon_skipping_nmd_candidates.tsv`).
