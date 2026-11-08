# IsoformVUS: Exon Architecture, Alternative Isoform Expression, and Co-Transcriptional NMD Govern Rare-Disease Variant Interpretation

**Author:** Sai Krishna Uppula  
**Repository:** [github.com/Saikrishnnauppula/IsoformVUS](https://github.com/Saikrishnnauppula/IsoformVUS)  
**Genome Build & Annotations:** GRCh38, GENCODE v46 / v44, [MANE Select v1.5](https://www.nature.com/articles/s41586-022-04558-8), GTEx v8  

## Abstract

Clinical genome sequencing routinely returns variants of uncertain significance (VUS). Standard variant-effect predictors—such as CADD, [SpliceAI](https://www.cell.com/cell/fulltext/S0092-8674(18)31629-5), and [AlphaMissense](https://www.science.org/doi/10.1126/science.adg7492)—evaluate coding variants primarily on the canonical transcript without modeling how exon length, alternative exon skipping in human tissues, and nonsense-mediated decay (NMD) interact across variant classes. Here we analyze 587,227 ClinVar variants and a gene-stratified benchmark of 60,521 labeled variants (30,146 pathogenic/likely pathogenic; 30,375 benign/likely benign) across 102 Mendelian disease genes (4,132 reconstructed MANE Select exons, including 3,955 internal coding exons).

First, we audit our baseline tiered models (M1–M5) and show that right-endpoint step evaluation of average precision penalized the discrete M1 baseline by −0.0321 across a 6,077-variant tie block ([McSherry & Najork 2008](https://www.microsoft.com/en-us/research/publication/computing-information-retrieval-performance-measures-efficiently-in-the-presence-of-tied-scores/)), accounting for the apparent M2/M3 ranking gain, while 94.5% of M4's remaining gain over M3 reflects canonical coding-sequence membership. Second, when a variant leaves the MANE Select coding sequence intact, alternative-transcript hits separate into benign dual-coding reading-frame shifts (n = 10,800, 99.93% benign), alternative-exon missense variants disproportionately stalled as VUS (n = 2,173, 64.15% VUS/conflicting, resolved benign rate 89.2%), and alternative-exon loss-of-function variants whose clinical pathogenicity scales with GTEx isoform usage (28.6%–35.0% resolved pathogenic when isoform usage ≥ 15%).

Third, examining 301,264 variants located strictly inside internal MANE Select coding exons, we find a pronounced divergence between amino-acid substitutions and premature termination codons (PTCs). Across internal exons stratified by nucleotide length (compact ≤180 bp vs. giant >1,000 bp), GTEx-expressed exon skipping by an alternative protein-coding isoform (≥10% usage, ≥1.0 TPM), and exon reading-frame symmetry (length divisible by 3), resolved missense pathogenicity drops from 73.13% (3,737 pathogenic : 1,373 benign) in compact, constitutive, symmetric exons to 1.17% (3 pathogenic : 254 benign; odds ratio = 0.0043) in giant, backup-skipped, symmetric exons. This pattern holds within individual genes (Cochran–Mantel–Haenszel odds ratio = 0.2005, p = 4.29 × 10⁻¹⁸⁷ for exon length >400 bp; odds ratio = 0.3461, p = 3.53 × 10⁻⁴² for expressed backup skipping) and is mediated primarily by 3D structural constraint rather than nucleotide conservation (77.9% and 94.7% mediated by AlphaMissense in counterfactual mediation analysis). In contrast, nonsense and frameshift variants within those same giant, backup-skipped, phase-symmetric internal exons remain 99.85%–100.00% pathogenic (n = 3,498). Because pre-mRNA splicing commits the dominant transcript pool to include the mutated exon in the nucleus before cytoplasmic translation, an internal PTC triggers nonsense-mediated decay of the canonical mRNA regardless of whether the encoded peptide segment is structurally dispensable.

Fourth, by resolving transcriptional strand orientation across all 102 genes, we show that among missense variants predicted to be structurally benign by AlphaMissense (score < 0.34), false negatives within 3 bp of an exon boundary are concentrated at the 3' end of the exon (the 5'-splice donor site recognized by U1 snRNP base-pairing across positions −3 to −1), where 65.22% (15/23) are pathogenic compared with 13.33% (2/15) at the 5' end of the exon (odds ratio = 12.19, Fisher p = 0.0024). Incorporating these exon-architecture, NMD-zone ([Nagy & Maquat 1998](https://pubmed.ncbi.nlm.nih.gov/9644970/)), and splice-boundary features into a tiered model (M6) improves held-out test AUROC across 24 unseen genes from 0.8028 (M4) to 0.9930 (tie-corrected AUPRC 0.9934), yields empirical Bayesian likelihood ratios ([Tavtigian et al. 2018](https://pmc.ncbi.nlm.nih.gov/articles/PMC6336098/)) that stratify 23,477 unresolved ClinVar VUS, and identifies 304 phase-symmetric or naturally skipped internal exons harboring 10,910 pathogenic patient PTCs that are candidates for antisense oligonucleotide (ASO) exon-skipping therapy ([Lauffer et al. 2024](https://pmc.ncbi.nlm.nih.gov/articles/PMC10770028/)).

## 1. Introduction

Exome and genome sequencing for suspected Mendelian disorders frequently identify rare variants whose clinical significance remains uncertain. Computational tools such as CADD, [SpliceAI](https://www.cell.com/cell/fulltext/S0092-8674(18)31629-5), and [AlphaMissense](https://www.science.org/doi/10.1126/science.adg7492) provide essential sequence- and structure-based priors, yet most pipelines evaluate coding consequences on a single canonical transcript—typically [MANE Select](https://www.nature.com/articles/s41586-022-04558-8)—without accounting for the exon architecture of the gene or the expression of alternative splice isoforms across human tissues ([Cummings et al. 2020](https://pmc.ncbi.nlm.nih.gov/articles/PMC7334041/)).

In the initial stages of the IsoformVUS project (Stages 01–17), we assembled 60,521 ClinVar-labeled variants across 102 Mendelian disease genes and trained five sequential logistic-regression models (M1: variant scores → M2: +HPO phenotype → M3: +cell-type expression → M4: +isoform structural flags → M5: +tissue isoform usage). While M4 improved standard step AUPRC from 0.8016 to 0.8524 (+0.0509, 95% CI [+0.0473, +0.0546]; tie-corrected AUPRC 0.8334 to 0.8517), subsequent diagnostic audits raised two fundamental questions: why did shuffling HPO and GTEx tables in M2 and M3 still produce an apparent +0.031 AUPRC gain over M1, and why did continuous isoform-usage weighting in M5 add no discrimination over binary structural flags in M4?

Beyond resolving those modeling artifacts, a broader biological question remained unaddressed. Previous isoform-aware studies have focused on variants that fall outside the canonical transcript but hit an alternative isoform. Far less is known about the inverse setting: when a variant directly alters an internal coding exon of the canonical MANE Select transcript, how do exon length, alternative exon skipping in normal human tissues, and reading-frame phase jointly shape the tolerance of missense and in-frame indel variants versus premature termination codons?

## 2. Deconstruction of Tiered Models M1–M5 and Non-Canonical Isoform Regimes

### 2.1 Resolution of the M2/M3 Tie-Breaking Artifact and M4 Shared-Exon Proxy

In the 14,400-variant test split (24 held-out genes), 6,077 single-nucleotide variants (SNVs)—primarily intronic and synonymous variants—lack AlphaMissense scores and share an identical CADD missingness pattern, receiving an identical predicted probability under the 4-feature M1 baseline. Standard right-endpoint step evaluation of average precision (`sklearn.metrics.average_precision_score`) evaluates tied blocks at their worst-case rank, penalizing M1 by −0.0321 relative to the exact analytical expectation under uniform random tie-breaking ([McSherry & Najork 2008](https://www.microsoft.com/en-us/research/publication/computing-information-retrieval-performance-measures-efficiently-in-the-presence-of-tied-scores/)):

| Model Tier | Test AUROC | Standard Step AUPRC | Tie-Corrected AUPRC (McSherry & Najork 2008) | Within-Gene Mean Tie-Corrected AUPRC |
| :--- | :---: | :---: | :---: | :---: |
| **M1 (CADD + AlphaMissense)** | 0.7985 | 0.8016 | **0.8334** | 0.8299 |
| **M2 (+HPO Phenotype Counts)** | 0.7967 | 0.8327 | **0.8340** (+0.0006) | 0.8299 |
| **M3 (+Cell-Type Relevance)** | 0.7953 | 0.8318 | **0.8332** (−0.0002) | 0.8297 |
| **M3 + Canonical Coding Indicator (1-bit)** | 0.8005 | 0.8498 | **0.8506** (+0.0174) | 0.8454 |
| **M4 (+Isoform Structural Flags)** | 0.8028 | 0.8524 | **0.8517** (+0.0183) | 0.8448 |
| **M5 (+Dominant Isoform Usage)** | 0.8027 | 0.8523 | **0.8516** (+0.0182) | 0.8447 |

Because M2 and M3 add gene-level constants, they break the 6,077-variant inter-gene tie block without changing within-gene variant ordering (within-gene AUPRC = 0.8299 across M1, M2, and M3). Furthermore, in Stage 06 (`06_variant_isoform_effects.py`), coding and splice-junction hits were recorded whenever a variant overlapped at least one protein-coding transcript but fewer than the total number of annotated transcripts (`n_overlapped < n_transcripts`). Because 99 of the 102 panel genes have non-coding or partial GENCODE transcripts that omit canonical exons, M4's coding/junction flags were 96.57% concordant with a single 1-bit indicator of whether the variant lies in the MANE Select coding sequence (`M3 + Canonical Coding Indicator` tie-corrected AUPRC = 0.8506 vs. M4 = 0.8517).

### 2.2 Three Regimes of Canonical-Idle Isoform Discordance Across 587,227 ClinVar Variants

Across the unselected 587,227 ClinVar variants in our 102 genes (`supp_clinvar_all_table.tsv.gz`), 16,660 variants have a more severe predicted consequence on an alternative GENCODE v46 protein-coding transcript than on MANE Select. Rather than forming a uniform class, these discordant variants separate into three distinct biological regimes:

1. **Regime I — Dual-Coding Reading-Frame Shifts (MANE = Synonymous, Alternative ≥ Missense; n = 10,800)**: These variants lie inside canonical MANE exons where an alternative isoform uses a shifted reading frame. Once resolved in ClinVar, **99.93% (7,128 / 7,133)** are Benign or Likely Benign.
2. **Regime II — Alternative-Exon Missense Variants (MANE = Intronic/UTR, Alternative = Missense; n = 2,173)**: Among resolved records, 89.20% (570 / 639) are Benign/Likely Benign, consistent with weak purifying selection on minor alternative exons ([Pozo et al. 2022](https://www.nature.com/articles/s42003-022-03382-w)). However, **64.15% (1,394 / 2,173)** remain trapped as VUS or conflicting in ClinVar (compared with 12.68% of canonical intronic variants; odds ratio = 11.16, p < 10⁻³⁰⁰) because standard annotation pipelines flag the alternative-isoform missense change without weighting by transcript expression.
3. **Regime III — Alternative-Exon Loss-of-Function Variants (MANE = Intronic/UTR, Alternative = Nonsense/Frameshift/Canonical Splice; n = 877)**: Resolved pathogenicity is **28.57% (60 / 210)**—a 16.5-fold enrichment over canonical intronic variants (p = 2.95 × 10⁻⁴²)—and rises to **35.0%** when the affected alternative isoform has ≥15% GTEx usage and ≥1.0 median TPM in disease-relevant tissues (such as alternative 3'-terminal exons of *PCDH15*, *SCN1A*, and *USH2A*).

## 3. Exon Architecture, Alternative Exon Skipping, and Co-Transcriptional NMD in Internal Coding Exons

### 3.1 The 4 × 2 × 2 Internal Exon Matrix Across 301,264 Coding Variants

We next examined variants that directly alter internal coding exons of the canonical MANE Select transcript (excluding 3'-terminal exons to remove classical last-exon NMD escape; [Abou Tayoun et al. 2018](https://pmc.ncbi.nlm.nih.gov/articles/PMC6185798/)). Across all 102 genes, we reconstructed 3,955 internal coding exons containing 301,264 ClinVar coding variants (220,209 missense, 9,985 in-frame indels, 58,030 frameshift, and 23,025 nonsense variants) and classified each exon along three coordinates:

* **Exon nucleotide length**: compact splice-defined exons (≤180 bp), medium exons (181–400 bp), large exons (401–1,000 bp), and giant exons (>1,000 bp).
* **GTEx-expressed backup exon skipping**: whether an alternative protein-coding transcript with ≥10% isoform usage and ≥1.0 median TPM in GTEx v8 skips the exon entirely (`Skip = 1`) versus constitutive inclusion (`Skip = 0`).
* **Exon reading-frame phase symmetry**: whether the exon nucleotide length is an exact multiple of 3 (`L mod 3 == 0`, symmetric phase, so skipping preserves the downstream reading frame) versus asymmetric (`L mod 3 != 0`; [Light & Elofsson 2013](https://academic.oup.com/nar/article-abstract/41/8/4409/2408986)).

Across the 3,955 internal exons, GTEx-expressed backup exon skipping is significantly enriched in phase-symmetric exons (**46.18%, 580/1,256** symmetric among backup-skipped exons vs. **39.12%, 451/1,153** among constitutive exons; Fisher p = 5.0 × 10⁻⁴). Stratifying all internal coding variants across the resulting 4 × 2 × 2 matrix reveals an 18.1-fold reduction in resolved missense pathogenicity across the 4 × 2 length-by-skipping classes—and a **62.5-fold reduction** when including phase symmetry—while internal frameshift and nonsense variants remain 99.85%–100.00% pathogenic in every state:

| Internal Exon Architecture State (Excluding 3'-Terminal Exons) | Phase Symmetry | 587k Missense Total (P/LP : B/LB) | 587k Missense Resolved P/LP % | 587k Internal Frameshift Resolved P/LP % (n) | 587k Internal Nonsense Resolved P/LP % (n) | Bayesian Likelihood Ratio ([Tavtigian et al. 2018](https://pmc.ncbi.nlm.nih.gov/articles/PMC6336098/)) |
| :--- | :---: | :---: | :---: | :---: | :---: | :--- |
| **1A. Compact (≤180 bp) · Constitutive** | Asymmetric | 53,875 (5,586 : 2,205) | **71.70%** | **99.90%** (13,341) | **99.92%** (5,522) | 1.880 (Indeterminate) |
| **1B. Compact (≤180 bp) · Constitutive** | Symmetric | 33,428 (3,737 : 1,373) | **73.13%** | **99.86%** (8,572) | **99.97%** (3,372) | 2.020 (Near Supporting Pathogenic) |
| **1C. Compact (≤180 bp) · Backup-skipped** | Asymmetric | 16,877 (1,000 : 1,055) | **48.66%** | **99.97%** (5,414) | **99.95%** (2,071) | 0.704 (Indeterminate) |
| **1D. Compact (≤180 bp) · Backup-skipped** | Symmetric | 9,689 (845 : 616) | **57.84%** | **100.00%** (3,121) | **100.00%** (1,180) | 1.018 (Indeterminate) |
| **2A. Medium (181–400 bp) · Constitutive** | Asymmetric | 25,480 (2,175 : 1,239) | **63.71%** | **99.86%** (6,845) | **99.96%** (2,837) | 1.303 (Indeterminate) |
| **2B. Medium (181–400 bp) · Constitutive** | Symmetric | 11,787 (1,145 : 564) | **67.00%** | **99.90%** (3,122) | **100.00%** (1,264) | 1.506 (Indeterminate) |
| **2C. Medium (181–400 bp) · Backup-skipped** | Asymmetric | 11,351 (503 : 769) | **39.54%** | **100.00%** (2,740) | **99.92%** (1,186) | 0.486 (Near Supporting Benign) |
| **2D. Medium (181–400 bp) · Backup-skipped** | Symmetric | 6,482 (248 : 391) | **38.81%** | **100.00%** (1,481) | **100.00%** (631) | **0.471 (BP Supporting Benign, −1 pt)** |
| **3A. Large (401–1,000 bp) · Constitutive** | Asymmetric | 9,687 (280 : 338) | **45.31%** | **99.95%** (2,418) | **100.00%** (947) | 0.615 (Indeterminate) |
| **3B. Large (401–1,000 bp) · Constitutive** | Symmetric | 2,557 (67 : 110) | **37.85%** | **100.00%** (613) | **100.00%** (241) | **0.453 (BP Supporting Benign, −1 pt)** |
| **3C. Large (401–1,000 bp) · Backup-skipped** | Asymmetric | 4,612 (122 : 319) | **27.66%** | **99.77%** (965) | **100.00%** (415) | **0.285 (BP Supporting Benign, −1 pt)** |
| **3D. Large (401–1,000 bp) · Backup-skipped** | Symmetric | 1,239 (43 : 124) | **25.75%** | **100.00%** (275) | **100.00%** (90) | **0.259 (BP Supporting Benign, −1 pt)** |
| **4A. Giant (>1,000 bp) · Constitutive** | Asymmetric | 15,609 (192 : 1,224) | **13.56%** | **100.00%** (3,975) | **99.93%** (1,497) | **0.117 (BM Moderate Benign, −2 pts)** |
| **4B. Giant (>1,000 bp) · Constitutive** | Symmetric | 5,807 (62 : 503) | **10.97%** | **100.00%** (2,501) | **99.85%** (796) | **0.092 (BM Moderate Benign, −2 pts)** |
| **4C. Giant (>1,000 bp) · Backup-skipped** | Asymmetric | 7,431 (42 : 831) | **4.81%** | **100.00%** (1,440) | **100.00%** (738) | **0.038 (BS Strong Benign, −4 pts)** |
| **4D. Giant (>1,000 bp) · Backup-skipped** | Symmetric | 4,298 (3 : 254) | **1.17%** | **100.00%** (2,707) | **99.85%** (791) | **0.010 (BS Strong Benign, −4 pts)** |

To verify that this divergence operates within individual genes rather than reflecting differences between genes, we performed gene-stratified Cochran–Mantel–Haenszel (CMH) tests across all informative panel genes. Controlling for gene identity, large/giant internal exon size (>400 bp vs. ≤400 bp) is strongly protective for missense variants (**CMH odds ratio = 0.2005, Z = −29.17, p = 4.29 × 10⁻¹⁸⁷** across 52 genes; paired Wilcoxon lower in 40/54 genes, p = 5.74 × 10⁻⁷), as is GTEx-expressed backup exon skipping (**CMH odds ratio = 0.3461, Z = −13.61, p = 3.53 × 10⁻⁴²** across 42 genes). For internal PTVs in those same genes, neither large exon size (**CMH odds ratio = 1.6093, Z = +9.01**) nor backup skipping (**CMH odds ratio = 1.7367, Z = +8.25**) reduces pathogenicity.

### 3.2 Label-Independent Biophysical Controls, Coding Density per Kilobase, and Causal Mediation

Because clinical curations use ACMG criteria, we tested whether the 4 × 2 pattern holds independently of ClinVar pathogenic/benign labels by examining continuous [AlphaMissense](https://www.science.org/doi/10.1126/science.adg7492) scores, CADD v1.7 PHRED scores, and variant density per kilobase of coding sequence across the 4 × 2 matrix:

| Internal Exon Architecture State | Coding Sequence (kb) | Median AlphaMissense (% > 0.564) | Median CADD PHRED | Benign Missense / kb | Pathogenic Missense / kb | Benign : Pathogenic Missense Ratio | Pathogenic PTV / kb |
| :--- | :---: | :---: | :---: | :---: | :---: | :---: | :---: |
| **1A. Compact (≤180 bp) · Constitutive** | 208.9 kb | **0.943** (74.0%) | **26.60** | 17.1 / kb | **44.6 / kb** | **0.38×** | **146.5 / kb** |
| **1B. Compact (≤180 bp) · Backup-skipped** | 109.1 kb | **0.807** (63.1%) | 24.70 | 15.3 / kb | 16.9 / kb | **0.91×** | 84.9 / kb |
| **2A. Medium (181–400 bp) · Constitutive** | 95.6 kb | **0.946** (70.4%) | 25.60 | 18.9 / kb | 34.7 / kb | **0.54×** | 138.5 / kb |
| **2B. Medium (181–400 bp) · Backup-skipped** | 92.6 kb | **0.666** (54.8%) | 25.90 | 12.5 / kb | 8.1 / kb | **1.54×** | 60.5 / kb |
| **3A. Large (401–1,000 bp) · Constitutive** | 28.0 kb | **0.644** (52.5%) | 23.30 | 16.0 / kb | 12.4 / kb | **1.29×** | 139.0 / kb |
| **3B. Large (401–1,000 bp) · Backup-skipped** | 23.1 kb | **0.220** (37.4%) | 22.90 | 19.1 / kb | 7.1 / kb | **2.68×** | 71.5 / kb |
| **4A. Giant (>1,000 bp) · Constitutive** | 57.0 kb | **0.138** (18.1%) | 14.45 | **30.3 / kb** | 4.5 / kb | **6.80×** | **148.6 / kb** |
| **4B. Giant (>1,000 bp) · Backup-skipped** | 54.1 kb | **0.104** (**4.6%**) | 20.90 | **20.0 / kb** | **0.8 / kb** | **24.11×** | **102.7 / kb** |

Without using ClinVar labels, median AlphaMissense drops 9.1-fold from 0.943 in compact constitutive exons to 0.104 in giant backup-skipped exons. In parallel, the ratio of benign to pathogenic missense variants per kilobase shifts 63.4-fold from 0.38× (17.1 benign vs. 44.6 pathogenic per kb) to 24.11× (20.0 benign vs. 0.8 pathogenic per kb), whereas pathogenic PTV density remains high across all eight states (60.5–148.6 per kb). Across the 8 internal exon states, the spatial distribution of pathogenic PTVs (n = 78,284) diverges from pathogenic missense variants (n = 16,050) with a Kullback–Leibler divergence of 0.4892 bits (and 0.8084 bits between benign and pathogenic missense variants).

In counterfactual mediation analysis (1,000 bootstrap resamples), 3D structural constraint measured by AlphaMissense mediates **77.9%** (95% CI [69.2%, 89.9%]; within-gene fixed-effects **70.9%**, 95% CI [57.1%, 90.8%]) of the large/giant exon effect (total effect = −0.3820, indirect effect = −0.2976, direct effect = −0.0844) and **94.7%** (within-gene fixed-effects **72.6%**) of the expressed backup-skipping effect. By contrast, CADD PHRED sequence conservation fails to mediate the within-gene backup-skipping effect (−22.1%). Biophysical selection explains both findings: compact exons (≤180 bp) bounded by spliceosomal exon definition typically encode tightly packed globular cores, whereas exons expanded beyond >400–1,000 bp or naturally spliced out in ≥10% of normal transcripts predominantly encode peripheral or intrinsically disordered segments ([Buljan et al. 2012](https://www.cell.com/molecular-cell/fulltext/S1097-2765(12)00302-4); [Light & Elofsson 2015](https://pmc.ncbi.nlm.nih.gov/articles/PMC4477647/)). However, when a nonsense or frameshift mutation occurs inside such an exon, nuclear pre-mRNA splicing still incorporates the exon into 60%–90% of nascent canonical transcripts before cytoplasmic translation, triggering EJC-dependent NMD and depleting functional mRNA. Consistent with this mechanism, among the 53 rare Benign/Likely Benign internal PTVs out of 78,337 resolved internal PTVs (0.068%), **39.62% (21/53)** lie within 5 bp of an exon-intron junction compared with 9.56% of pathogenic internal PTVs (odds ratio = 6.21, Fisher p = 2.9 × 10⁻⁸), where splice-junction proximity can induce exon skipping or alternative splice-site usage ([medRxiv 2022](https://www.medrxiv.org/content/10.1101/2022.03.03.22271834v1.full-text)).

## 4. 5'-Splice Donor Exon-End Asymmetry at Exon Boundaries and the 50-nt EJC NMD Boundary

### 4.1 3' Exon-End (5'-Donor) Concentration of AlphaMissense Exon-Boundary False Negatives

Because protein structure and language models ([AlphaMissense](https://www.science.org/doi/10.1126/science.adg7492)) take only the spliced amino-acid sequence as input, they lack coordinates for exon-intron boundaries ([medRxiv Oct 2026](https://www.medrxiv.org/content/10.64898/2026.10.09.26365070v1)). Using our strand-resolved exon coordinates, we separated each variant's distance to the **5' start of the exon (the 3'-splice acceptor site `YAG|G`, recognized by U2AF)** from its distance to the **3' end of the exon (the 5'-splice donor site `MAG|GTRAGT`, where U1 snRNP base-pairs directly across positions −3, −2, and −1 of the exon)**:

| Exonic Boundary Distance Bin | Benchmark Missense with AM < 0.34 (n) | Pathogenic Fraction (False Negative Rate %) | 587k ClinVar Synonymous Resolved (P/LP : B/LB) | 587k Synonymous Resolved P/LP % | Primary Splicing / Structural Mechanism |
| :--- | :---: | :---: | :---: | :---: | :--- |
| **3' Exon End (d ≤ 3 bp, 5'-splice donor edge)** | 23 | **65.22%** (15 / 23) | 2,456 (10 : 2,446) | **0.407%** | **U1 snRNP exonic `MAG` base-pairing disruption (OR = 12.19, p = 0.0024 vs. 5' exon start)** |
| **5' Exon Start (d ≤ 3 bp, 3'-splice acceptor edge)** | 15 | **13.33%** (2 / 15) | 2,934 (3 : 2,931) | **0.102%** | Intronic polypyrimidine/`YAG` dominates U2AF recognition (only +1 G is exonic) |
| **Combined Exon Edge (d ≤ 3 bp)** | 38 | **44.74%** (17 / 38) | 5,390 (13 : 5,377) | **0.241%** | **3.0-fold increase in false negatives vs. exon interior (p = 4.9 × 10⁻⁵)** |
| **Proximal Flank (d = 4–10 bp)** | 168 | **26.79%** (45 / 168) | 6,412 (9 : 6,403) | **0.140%** | Proximal exonic splicing enhancer (ESE) flank |
| **Near Interior (d = 11–50 bp)** | 694 | **22.48%** (156 / 694) | 31,420 (32 : 31,388) | **0.102%** | Compact and medium exon interior |
| **Deep Exon Core (d = 51–300 bp)** | 696 | **20.11%** (140 / 696) | 23,262 (18 : 23,244) | **0.077%** | Medium and large exon interior |
| **Giant Exon Core (d > 300 bp)** | 434 | **3.46%** (15 / 434) | 18,894 (0 : 18,894) | **0.000%** | **Disordered linker interior (0 / 18,894 pathogenic synonymous)** |

Among benchmark missense variants predicted benign by AlphaMissense (`am_score < 0.34`), **88.2% (15/17)** of all false negatives within 3 bp of an exon boundary occur at the 3' end of the exon (5'-splice donor site), where the pathogenic rate reaches **65.22% (15/23)** versus **13.33% (2/15)** at the 5' end of the exon (odds ratio = 12.19, Fisher p = 0.0024). Across 118,877 ClinVar synonymous variants, resolved pathogenicity at the 3' exon end (0.407%, 10:2,446) is likewise 4.0-fold higher than at the 5' exon start (0.102%, 3:2,931) and drops to 0.000% (0/18,894) at >300 bp inside giant exons.

### 4.2 In-Frame Indel Polarity and the 50–55 nt Penultimate-Exon EJC Boundary

Across 9,985 internal in-frame indels in ClinVar, in-frame deletions are substantially more pathogenic than in-frame insertions in compact/medium constitutive exons (**97.7% resolved P/LP**, 1,156:27 for deletions vs. **82.4%**, 566:121 for insertions; odds ratio = 9.15, p < 10⁻³⁰), reflecting register disruption in packed secondary structures. In large/giant backup-skipped exons (>400 bp), resolved pathogenicity drops 4.0-fold for in-frame deletions to **24.2%** (22:69; odds ratio = 0.0074, p < 10⁻⁵⁵) and to **35.2%** (43:79) for in-frame insertions.

For premature termination codons (n = 87,779 across the 102 genes), mapping the exact nucleotide distance to the 3'-most exon-exon junction ([Nagy & Maquat 1998](https://pubmed.ncbi.nlm.nih.gov/9644970/)) separates four clinical regimes:

* **Canonical internal NMD-trigger zone (>50 nt upstream of the final exon-exon junction; n = 80,127)**: **96.66% raw P/LP**, **1.87% VUS/conflicting**, and **99.93% resolved P/LP** (77,452:51).
* **Penultimate-exon 50-nt EJC escape zone (1–50 nt upstream of the final exon-exon junction; n = 928)**: raw P/LP drops to **89.66%** while VUS/conflicting classifications rise **4.7-fold to 8.73%** (81/928, p < 10⁻²⁵), showing that NMD escape in the distal 50 nt of the penultimate exon increases clinical classification uncertainty even within an internal exon.
* **Last-exon proximal truncation (`in_last_exon = 1`, <90% of CDS; n = 4,298)**: **94.00% raw P/LP** and **4.77% VUS/conflicting** (99.90% resolved P/LP, 4,040:4).
* **Last-exon distal C-terminal tail (`in_last_exon = 1`, ≥90% of CDS; n = 2,426)**: raw P/LP drops to **70.20%** while VUS/conflicting classifications surge **14.6-fold to 27.37%** (664/2,426; resolved P/LP 97.71%, 1,703:40; note that in C-terminal degron genes such as *ASXL3*, NMD-escaping truncations can act via gain-of-function; [ASXL3 Collaborative Group 2026](https://www.medrxiv.org/content/10.64898/2026.07.20.26358515v1.full)).

## 5. Candidate Exons for Antisense Exon-Skipping Therapy and Tiered Model Evaluation

### 5.1 Prioritizing Internal Exons for Antisense Oligonucleotide (ASO) Exon Skipping

Because nonsense and frameshift variants in structurally dispensable, phase-symmetric or naturally skipped internal exons cause disease by routing the nascent transcript into NMD rather than by removing an essential protein domain, steric-blocking ASOs that induce exon skipping can restore an in-frame, functional transcript ([Lauffer et al. 2024](https://pmc.ncbi.nlm.nih.gov/articles/PMC10770028/)). We scored all 3,955 internal exons across three criteria—empirical-Bayes missense/AlphaMissense tolerance, reading-frame symmetry (`L mod 3 == 0`), and GTEx backup exon skipping—identifying **304 Tier 1 exons** across 36 genes (harboring **10,910 pathogenic patient PTVs**: 3,057 nonsense and 7,853 frameshift) and **929 Tier 1+2 exons** across 81 genes (harboring **20,823 pathogenic patient PTVs**: 6,277 nonsense and 14,546 frameshift; full table in `results/exon_skipping_nmd_candidates.tsv`):

| Exon ID & Gene | Disease Category | Exon Length | Phase Symmetric (`L mod 3 == 0`) | GTEx Backup Skip Usage | Internal Missense Resolved P/LP % (Median AM) | Pathogenic Patient PTVs (Nonsense : Frameshift) | Candidate Score |
| :--- | :--- | :---: | :---: | :---: | :---: | :---: | :---: |
| **Exon 955 (*BRCA2*)** | Hereditary cancer / Fanconi | 4,932 bp | **Yes** | 0.00 (giant linker) | **0.7%** (0.113) | **2,666** (479 : 2,187) | **0.569 (Tier 1)** |
| **Exon 1331 (*BRCA1* Exon 11)** | Hereditary breast/ovarian | 3,423 bp | No (Δ11q donor) | **0.69** | **0.3%** (0.115) | **2,266** (366 : 1,900) | **0.589 (Tier 1)** |
| **Exon 3975 (*TTN*)** | Cardiomyopathy | 17,097 bp | **Yes** | **0.20** | **1.3%** (0.583) | **1,197** (417 : 780) | **0.664 (Tier 1)** |
| **Exon 2588 (*ANKRD11*)** | KBG syndrome | 6,541 bp | No | **0.14** | **4.4%** (0.098) | **612** (179 : 433) | **0.527 (Tier 2)** |
| **Exon 3941 (*TTN*)** | Cardiomyopathy | 5,602 bp | **Yes** | **0.20** | **2.6%** (0.241) | **299** (115 : 184) | **0.769 (Tier 1)** |
| **Exon 4400 (*BARD1*)** | Hereditary cancer | 949 bp | **Yes** | 0.25 | **0.0%** (0.076) | **240** (65 : 175) | **0.571 (Tier 1)** |
| **Exon 226 (*DSP*)** | Arrhythmogenic cardiomyopathy | 2,289 bp | No (DSPIa) | **0.56** | **5.9%** (0.135) | **239** (77 : 162) | **0.544 (Tier 2)** |
| **Exon 1820 (*KMT2D*)** | Kabuki syndrome | 2,789 bp | **Yes** | 0.00 (giant linker) | **1.5%** (0.164) | **217** (109 : 108) | **0.553 (Tier 1)** |
| **Exon 3912 (*NSD1*)** | Sotos syndrome | 2,552 bp | No | **0.19** | **0.0%** (0.096) | **211** (55 : 156) | **0.546 (Tier 2)** |
| **Exon 4268 (*TTN*)** | Cardiomyopathy | 2,773 bp | **Yes** | **0.20** | **0.0%** (0.161) | **156** (71 : 85) | **0.814 (Tier 1)** |
| **Exon 523 (*CDKL5*)** | DEE2 epilepsy | 959 bp | **Yes** | 0.09 | **0.0%** (0.112) | **130** (30 : 100) | **0.565 (Tier 1)** |
| **Exon 3537 (*TSC1*)** | Tuberous sclerosis | 555 bp | No | **0.13** | **0.0%** (0.104) | **126** (26 : 100) | **0.530 (Tier 2)** |

### 5.2 Out-of-Fold Model Evaluation (M6) and Sensitivity Checks

Adding the exon-architecture, NMD-zone, and 5'-donor boundary features to our tiered benchmark (M6, evaluated strictly on the 24 held-out test genes, n = 14,400, with zero gene or exon overlap) improves global test AUROC from 0.8028 (M4) to **0.9930** and tie-corrected AUPRC from 0.8517 to **0.9934** (Brier score 0.0307 vs. 0.1654), while preserving internal missense discrimination (AUROC 0.9506, AUPRC 0.9784) and increasing recall on held-out 0–3 bp exon-edge variants with low AlphaMissense scores (`am_score < 0.34`) from 0.0% (M1/M3/M4) to **42.9%**. Under the [Tavtigian et al. (2018)](https://pmc.ncbi.nlm.nih.gov/articles/PMC6336098/) Bayesian framework, giant backup-skipped symmetric internal exons achieve Strong Benign support (LR = 0.010, −4 points; n = 3,279 VUS) and asymmetric exons achieve LR = 0.038 (−4 points; n = 3,630 VUS), stratifying 23,477 ClinVar VUS (`results/vus_bayesian_stratification.tsv`).

We verified robustness across seven sensitivity checks:
1. **Leave-top-genes-out jackknife**: Excluding *BRCA1*, *BRCA2*, *APC*, *TTN*, and *NF1* preserves the internal missense gradient (73.39% in compact constitutive vs. 7.05% in giant backup-skipped exons; CMH odds ratio = 0.2348, p = 7.91 × 10⁻¹³¹ for exon size >400 bp; CMH odds ratio = 0.5043, p = 5.01 × 10⁻¹⁵ for backup skipping).
2. **Review-status filtering (≥2 stars)**: Restricting to ≥2-star ClinVar variants (expert panel or multiple submitters with no conflicts) yields 85.48% pathogenic missense in compact constitutive exons (n = 1,315) versus 5.38% in giant backup-skipped exons (n = 93), while PTVs in giant backup-skipped exons remain 100.00% pathogenic (n = 272).
3. **Mutational spectrum control**: C>T and G>A transitions (which include hypermutable CpG sites) show an 18.7-fold drop (67.94% → 3.64%), matching the 17.7-fold drop across all other substitutions (74.47% → 4.20%).
4. **Within-exon structural islands**: Stratifying by AlphaMissense (<0.30, 0.30–0.70, >0.70) confirms that rare structured domain islands (`am_score > 0.70`) inside large/giant exons retain high pathogenicity (90.05%, n = 201), while disordered segments (`am_score < 0.30`) in large/giant exons have a 3.6-fold lower pathogenic rate than in compact exons (6.16% vs. 22.18%).
5. **Phase symmetry and internal PTV lethality**: Resolved PTV pathogenicity in phase-symmetric backup-skipped internal exons is 99.986% (7,350 / 7,351) versus 99.966% (14,730 / 14,735) in asymmetric exons, confirming that endogenous splicing does not automatically skip PTC-bearing exons without antisense intervention.
6. **Split integrity**: Train (61 genes), validation (17 genes), and test (24 genes) splits share 0 genes and 0 exons.
7. **Temporal holdout (post-2024 classifications)**: Among 2,408 internal missense variants newly resolved between 2024 and 2026 (`labels_2024_vs_2026.tsv.gz`), pathogenicity drops 12.2-fold from 78.61% (n = 1,183) in compact constitutive exons to 6.45% (n = 93) in giant backup-skipped exons.

## 6. Limitations and Reproducibility

This study has several limitations. First, our panel covers 102 Mendelian disease genes enriched for haploinsufficient and loss-of-function mechanisms; extending the 4 × 2 × 2 exon matrix to all ~20,000 human protein-coding genes using gnomAD v4 population allele frequencies and multiplexed assays of variant effect (MAVE) is a natural next step. Second, GTEx v8 measures bulk tissue transcript usage from short-read RNA-seq; single-cell long-read transcriptomics will refine isoform usage in developmental cell types. Third, while our exon-skipping candidate table identifies exons where amino-acid sequence is tolerant to variation and reading frame is preserved, experimental validation in patient-derived cells (measuring allele-specific mRNA rescue under NMD inhibition and ASO treatment) is required before clinical translation. All code, data provenance, and reproducible scripts (`scripts/01_download_sources.py` through `scripts/18_exon_architecture_nmd.py`) are available in the repository.

## 7. References

* **Abou Tayoun, A. N. et al. (2018)**. Recommendations for interpreting the loss of function PVS1 ACMG/AMP variant criterion. [Human Mutation 39, 1517–1524](https://pmc.ncbi.nlm.nih.gov/articles/PMC6185798/).
* **ASXL3 Collaborative Group (2026)**. ASXL3 truncating patient variants mediate transcriptional gain-of-function via NMD escape and C-terminal degron loss. [medRxiv 2026.07.20.26358515v1](https://www.medrxiv.org/content/10.64898/2026.07.20.26358515v1.full).
* **Buljan, M. et al. (2012)**. Tissue-specific splicing of disordered segments that embed binding motifs directs interaction networks. [Molecular Cell 46, 871–883](https://www.cell.com/molecular-cell/fulltext/S1097-2765(12)00302-4).
* **Cheng, J. et al. (2023)**. Accurate proteome-wide missense variant effect prediction with AlphaMissense. [Science 381, eadg7492](https://www.science.org/doi/10.1126/science.adg7492).
* **Cummings, B. B. et al. (2020)**. Transcript expression-aware annotation improves rare variant interpretation (pext). [Nature 581, 452–458](https://pmc.ncbi.nlm.nih.gov/articles/PMC7334041/).
* **Jaganathan, K. et al. (2019)**. Predicting splicing from primary sequence with deep learning (SpliceAI). [Cell 176, 535–548](https://www.cell.com/cell/fulltext/S0092-8674(18)31629-5).
* **Lauffer, M. C. et al. (2024)**. Possibilities and limitations of antisense oligonucleotide therapies for genetic diseases. [Nature Reviews Genetics / PMC10770028](https://pmc.ncbi.nlm.nih.gov/articles/PMC10770028/).
* **Light, S. & Elofsson, A. (2013)**. The impact of splice isoforms on protein domain architecture and exon-phase symmetry. [Nucleic Acids Research 41, 4409–4419](https://academic.oup.com/nar/article-abstract/41/8/4409/2408986).
* **Light, S. & Elofsson, A. (2015)**. Splice junctions are constrained by protein disorder. [Nucleic Acids Research 43, 4814–4822](https://pmc.ncbi.nlm.nih.gov/articles/PMC4477647/).
* **McSherry, F. & Najork, M. (2008)**. Computing information retrieval performance measures efficiently in the presence of tied scores. [ECIR 2008, LNCS 4956, 414–421](https://www.microsoft.com/en-us/research/publication/computing-information-retrieval-performance-measures-efficiently-in-the-presence-of-tied-scores/).
* **Morales, J. et al. (2022)**. A joint NCBI and EMBL-EBI transcript set for clinical genomics and research (MANE). [Nature 604, 310–315](https://www.nature.com/articles/s41586-022-04558-8).
* **Nagy, E. & Maquat, L. E. (1998)**. A rule for termination-codon position within intron-containing genes: when nonsense affects RNA abundance. [Trends in Biochemical Sciences 23, 198–199](https://pubmed.ncbi.nlm.nih.gov/9644970/).
* **Pozo, F. et al. (2022)**. Assessing the functional relevance of splice isoforms. [Communications Biology 5, 459](https://www.nature.com/articles/s42003-022-03382-w).
* **Spliceogenic Missense Consortium (2026)**. Investigating the pathogenicity of missense spliceogenic variants with AlphaMissense and SpliceAI. [medRxiv 2026.10.09.26365070v1](https://www.medrxiv.org/content/10.64898/2026.10.09.26365070v1).
* **Tavtigian, S. V. et al. (2018)**. Modeling the ACMG/AMP variant classification guidelines as a Bayesian classification framework. [Genetics in Medicine 20, 1054–1060](https://pmc.ncbi.nlm.nih.gov/articles/PMC6336098/).
