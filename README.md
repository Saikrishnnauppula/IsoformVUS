# IsoformVUS

**Isoform- and Exon-Architecture-Aware Variant Prioritization for Rare-Disease VUS Interpretation**

**Author:** Sai Krishna Uppula ([github.com/Saikrishnnauppula/IsoformVUS](https://github.com/Saikrishnnauppula/IsoformVUS))

IsoformVUS investigates how transcript structure, tissue-specific alternative splicing, and nonsense-mediated decay (NMD) govern the clinical pathogenicity of rare genetic variants across 102 Mendelian disease genes (60,521 benchmark variants across locked gene-level splits and 587,227 unselected ClinVar variants on GRCh38 / GENCODE v46 / MANE Select v1.5 / GTEx v8).

**Research scope.** This repository is a computational genomics research project, not a clinical diagnostic device. See [docs/clinical_safety.md](docs/clinical_safety.md).

## Key Findings & Project Summary

All 18 pipeline stages (`scripts/01_download_sources.py` through `scripts/18_exon_architecture_nmd.py`) have been executed and verified against `results/`:

1. **Tiered Benchmark (`results/benchmark.json`) and Tie-Discretization Audit (`breakthrough/round2/e1_e2`)**:
   - On the locked 24-gene held-out test set (n = 14,400 variants; 7,200 pathogenic, 7,200 benign), adding isoform structural features (`M4`) improves standard step AUPRC from **0.8016** (`M1`, AUROC **0.7985**) to **0.8524** (`M4`, AUROC **0.8028**; paired bootstrap ΔAUPRC = **+0.0509**, 95% CI [+0.0473, +0.0546], p < 0.001) and recovers 23 top-50 pathogenic variants missed by `M1`.
   - Follow-up permutation and tie-breaking controls (`breakthrough/round2/e1_e2`) established that the intermediate +0.031 step-AUPRC jump in `M2` (0.8327) and `M3` (0.8318) resulted from right-endpoint step evaluation across a 6,077-SNV tie block in `M1`: under analytical tie-corrected average precision (McSherry & Najork 2008), `M1` scores **0.8334**, `M2` scores **0.8340**, `M3` scores **0.8332**, and `M4` scores **0.8517** (+0.0183 over `M1`), with 94.5% of `M4`'s gain explained by canonical coding-sequence membership.
2. **Three Non-Canonical Isoform Regimes Across 587,227 ClinVar Variants**:
   - Variants that leave the MANE Select coding sequence intact but hit an alternative GENCODE v46 protein-coding isoform separate into: (i) **Dual-coding reading-frame shifts** (MANE = synonymous, alternative ≥ missense; n = 10,800, **99.93% benign**), (ii) **Alternative-exon missense variants** (MANE ≤ UTR/intronic, alternative = missense; n = 2,173, **89.20% benign** when resolved but **64.15% stalled as VUS/conflicting**), and (iii) **Alternative-exon loss-of-function variants** (n = 877, **28.57% resolved pathogenic**, rising to **35.0%** when GTEx isoform usage ≥ 15%).
3. **Internal Exon Architecture, Backup Skipping, and Co-Transcriptional NMD (`scripts/18_exon_architecture_nmd.py`)**:
   - Across 301,264 variants strictly inside 3,955 internal MANE Select coding exons (`in_last_exon = 0`), resolved missense pathogenicity drops from **73.13%** (3,737 P/LP : 1,373 B/LB) in compact (≤180 bp), constitutive, phase-symmetric (`L mod 3 == 0`) exons to **1.17%** (3 P/LP : 254 B/LB; OR = 0.0043) in giant (>1,000 bp), GTEx backup-skipped (≥10% usage), phase-symmetric exons (within-gene Cochran–Mantel–Haenszel OR = 0.2005, p = 4.29 × 10⁻¹⁸⁷ for exon length >400 bp; CMH OR = 0.3461, p = 3.53 × 10⁻⁴² for backup skipping; 77.9% and 94.7% mediated by 3D structural constraint in AlphaMissense).
   - Conversely, internal frameshift (n = 58,030) and nonsense (n = 23,025) variants in those same giant, backup-skipped, symmetric internal exons remain **99.85%–100.00% pathogenic** because nuclear pre-mRNA splicing commits the dominant transcript fraction to include the mutated exon prior to cytoplasmic translation and EJC-dependent NMD.
4. **3' Exon-End (5'-Splice Donor) Concentration of AlphaMissense False Negatives**:
   - Among benchmark missense variants predicted benign by AlphaMissense (`am_score < 0.34`), **88.2% (15/17)** of false negatives within 3 bp of an exon junction occur at the 3' end of the exon (the 5'-splice donor consensus `MAG|GTRAGT` recognized by U1 snRNP), where **65.22% (15/23)** are pathogenic versus **13.33% (2/15)** at the 5' start of the exon (OR = 12.19, Fisher p = 0.0024).
   - Incorporating exon-architecture, 50-nt penultimate-exon EJC NMD boundary, and 5'-donor edge features (`M6`) reaches **0.9930 test AUROC** and **0.9934 tie-corrected AUPRC** across the 24 held-out test genes, provides empirical Bayesian likelihood ratios (Tavtigian et al. 2018) across 23,477 ClinVar VUS (`results/vus_bayesian_stratification.tsv`), and ranks 304 phase-symmetric or backup-skipped internal exons harboring 10,910 pathogenic patient PTVs for antisense exon-skipping therapy (`results/exon_skipping_nmd_candidates.tsv`).

## Repository Layout

```
IsoformVUS/
├── README.md                  # project overview and summary of findings
├── LICENSE                    # MIT license
├── CITATION.cff               # citation metadata (Sai Krishna Uppula)
├── CHANGELOG.md               # version history (v0.1.0 -> v0.2.0)
├── requirements.txt           # pinned Python dependencies
├── environment.yml            # conda environment specification
├── pyproject.toml             # package build configuration
├── config/                    # YAML configuration files
├── data/                      # DATA_PROVENANCE.md and processed feature tables
├── metadata/                  # 102-gene panel and selection criteria
├── src/isoformvus/            # core Python package
├── tests/                     # automated unit tests
├── scripts/                   # runnable pipeline stages (01 through 18) & run_all.sh
├── breakthrough/              # Oct 20 follow-up controls, 587k ClinVar table & NMD features
├── results/                   # benchmark.json, statistics.json, exon_architecture_nmd.json & TSVs
├── figures/                   # publication figures (01 through 12)
├── docs/                      # architecture, audit, limitations & application summaries
├── research/                  # literature review and novelty assessment
└── paper/                     # manuscript.md, FINAL_REPORT.md & compiled PDF report
```

## Quick Start

```bash
# Run the complete pipeline and Stage 18 exon-architecture/NMD analysis:
bash scripts/run_all.sh

# Or run Stage 18 directly from the processed tables:
python3 scripts/18_exon_architecture_nmd.py
```
