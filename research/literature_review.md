# Literature Review: Isoform- and Cell-Type-Aware Variant Prioritization for Rare-Disease VUS Interpretation

**Project:** IsoformVUS
**Compiled:** 2026-10-06 (review window: 2024–2026, plus foundational earlier work)
**Compiler note:** Every entry below was found via live web search on 2026-10-06. Entries marked "citation pending verification" need a final DOI/title check before the manuscript. No paper is cited from training memory alone.

---

## 1. Scope and method of this review

This review covers: (a) isoform-aware and transcript-aware variant interpretation; (b) cell-type- and single-cell-aware variant prioritization; (c) phenotype-aware variant prioritization; (d) tissue/cell-specific splicing-effect predictors; (e) the baseline predictors CADD, SpliceAI, and AlphaMissense; (f) long-read and single-cell isoform biology relevant to variant interpretation. It does not claim completeness: the literature is moving fast (at least three directly relevant papers/preprints appeared in 2025–2026), and any "gap" claim must be re-checked at manuscript time.

---

## 2. Comparison table

| # | Method / paper | Year | Variant prioritization? | Phenotype usage | Cell-type / single-cell usage | Isoform / transcript usage | Key limitation |
|---|---|---|---|---|---|---|---|
| 1 | **IMPPROVE** — Chapman & Lassmann | 2025 | Yes (1,866 HPO-specific Random Forest models) | **Yes** — models built per HPO term via Phen2Gene gene lists | **Yes** — GTEx bulk + Tabula Sapiens single-cell expression features | **None** — gene-level expression only | No isoform layer; gene-level features cannot distinguish isoform-selective effects |
| 2 | **AlphaGenome Atlas / AVI** — Cheng J et al. (DeepMind/Broad/Exeter) | 2026 (preprint) | Yes — AVI score, PHRED-scaled, SHAP-attributed | None | Biosample-level (RNA-seq/DNase tracks incl. glutamatergic neurons) | **Yes** — splicing feature flags tissue-specific exon usage (e.g. brain-specific DNM1 exon 10a) | Not phenotype-aware; trained on gnomAD AF proxy labels, not ClinVar disease labels; noncoding-focused |
| 3 | **Weykopf et al.** (Edinburgh) | 2026 | No — analysis + case study, not a model | Disease-variant mapping (common + rare) | Tissue-level (tissue-specific isoforms) | **Yes** — long-read-resolved alternative isoforms; DPP9 lung-specific isoform case | No prioritization model; no phenotype→cell mapping; healthy-tissue data |
| 4 | **V2P** | 2025 | Yes — multi-task predictor (23 top-level HPO) | **Yes** — phenotype-specific prediction heads | None stated | None stated | No cell-type or isoform context |
| 5 | **Exomiser** — Smedley et al.; Jacobsen et al. | 2015–2025 | Yes — clinical gold standard, phenotype-driven | **Yes** — HPO similarity (PhenIX/PHIVE/hiPHIVE) + ACMG auto-classifier | None | MANE/Plus Clinical transcripts for reporting, not modeled | Gene-level; no expression/cell-type/isoform-usage features |
| 6 | **LIRICAL** — Kelly et al. | 2022 | Yes — likelihood-ratio genotype–phenotype model | **Yes** — HPO likelihood ratios | None | None | Gene-level; no molecular cell/isoform context |
| 7 | **VARPP** — Chapman & Lassmann | 2019 | Yes | **Yes** — phenotype-specific (Phenolyzer gene panels) | Tissue-level (GTEx expression specificity) | None | Predecessor of IMPPROVE; no single-cell, no isoform |
| 8 | **popEVE** — Orenbuch et al. | 2025 | Yes — missense severity ranking, cross-gene comparable | None | None | None | Protein-level only; no phenotype, cell type, or isoform |
| 9 | **AlphaMissense** — Cheng et al. | 2023 | Yes (missense pathogenicity) | None | None | Partial — canonical transcript scored; separate 60k-isoform resource released "for future research" | Missense-only; no phenotype/cell context; threshold calibrated on bulk ClinVar |
| 10 | **CADD** — Rentzsch et al.; CADD-Splice 2021 | 2019/2021 | Yes — deleteriousness score (all variant classes) | None | None | None | Tissue/cell-agnostic; VUS scores fall between benign and pathogenic (region-aware CADD thresholds, 2025) |
| 11 | **SpliceAI** — Jaganathan et al. | 2019 | Splicing-effect prediction, not disease prioritization | None | None | Implicit — scores all annotated junctions; not isoform-aware | Sequence-only; no expression, phenotype, or cell context |
| 12 | **AbSplice** — Wagner, Çelik et al. | 2023 | Aberrant-splicing probability per tissue | None | **Tissue-level** — 49 GTEx tissues via SpliceMaps | Implicit via tissue-specific splice-site annotations | Tissue, not cell-type; splicing-only |
| 13 | **MTSplice** — Cheng et al. | 2021 | Tissue-specific Δψ splicing prediction | None | Tissue-level | Transcript-context via splicing quantification | Splicing-only; tissue, not cell-type |
| 14 | **SpliceTransformer** — You et al. | 2024 | Tissue-specific splicing-effect score | None | Tissue-level (15 tissues) | None | Splicing-only; tissue, not cell-type |
| 15 | **CI-SpliceAI** — Strauch | 2022 | Splicing prediction retrained on all GENCODE isoforms | None | None | **Yes** — transcript-aware SpliceAI | Splicing-only; no phenotype/cell context; preprint status needs verification |
| 16 | **ESM1b isoform analysis** — Brandes et al. | 2023 | Missense effect prediction (protein language model) | None | None | **Yes** — per-isoform scores; ~1.8M isoform-sensitive variants | No phenotype/cell context; no prioritization benchmark on isoform-sensitive set |
| 17 | **Ly et al.** (Whitehead/MIT) + SwissIsoform (in development) | 2025 | Isoform-selective allele identification | Disease (rare disease patients) | Localization-level (mitochondrial vs cytosolic) | **Yes** — alternative start-codon N-terminal isoforms; TRNT1 SIFD | N-terminal isoforms only; tool not yet released |
| 18 | **Liu et al.** | 2025 | No — genetic architecture of isoform variation (WGS+RNA-seq) | None | None | **Yes** — isoform QTLs | Not a variant prioritizer |
| 19 | **cellSTAAR** — Lin/Natarajan labs | 2025 | Association testing (complex traits), not Mendelian VUS | None | **Yes** — scATAC-seq cell-type annotations | None | Complex-trait association, not rare-disease prioritization |
| 20 | **Phen2Gene** — Zhao et al. (cited via IMPPROVE) | 2020 | Gene-phenotype mapping (used by IMPPROVE) | **Yes** | None | None | Gene-level mapping tool, not a variant scorer |

---

## 3. Per-method entries

### 3.1 Baseline variant-effect predictors (no phenotype, cell type, or isoform modeling)

#### CADD (Rentzsch et al., 2019; CADD-Splice, 2021)
- **Publication/year:** Rentzsch P, Witten D, Cooper GM, Shendure J, Kircher M. "CADD: predicting the deleteriousness of variants throughout the human genome." *Nucleic Acids Research* 47(D1):D886–D894, 2019. CADD-Splice: Rentzsch P, Schubach M, Shendure J, Kircher M. "CADD-Splice—improving genome-wide variant effect prediction using deep learning-derived splice scores." *Genome Medicine* 13:31, 2021.
- **Datasets:** Observed vs simulated variants (human–chimp ancestral contrast); proxy-benign from high-frequency human variants; CADD-Splice integrates SpliceAI and MMSplice features.
- **Phenotype/cell-type/isoform:** none. Integrative, tissue-agnostic deleteriousness score; v1.6 added splice features via CADD-Splice.
- **Variant prioritization capability:** genome-wide, all variant classes; widely used baseline; CADD 1.6 chosen by IMPPROVE specifically because it is independent of ClinVar (avoids training circularity).
- **Limitations:** region-specific performance varies; a 2025 study on region-aware CADD thresholds (NAR Genomics & Bioinformatics, lqaf157) found CADD VUS distributions fall between benign and pathogenic, with AUC < 0.7 outside splicing/3′UTR/intronic regions — i.e., CADD alone is a weak VUS resolver. Tissue-agnostic.
- **How IsoformVUS differs:** CADD is an input feature, not the model; IsoformVUS adds phenotype, cell-type, and isoform layers that CADD deliberately omits.

#### SpliceAI (Jaganathan et al., 2019)
- **Publication/year:** Jaganathan K et al. "Predicting splicing from primary sequence with deep learning." *Cell* 176:535–548.e24, 2019.
- **Datasets:** trained on pre-mRNA sequence → splice junctions genome-wide; validated on RNA-seq and autism/ID de novo mutations (estimated 9–11% of pathogenic rare-disease mutations are cryptic-splice variants).
- **Phenotype/cell-type/isoform:** none. Scores all annotated junctions identically regardless of transcript, tissue, or cell type.
- **Variant prioritization capability:** the standard splicing-effect predictor; not itself a disease-prioritization model.
- **Limitations:** tissue- and isoform-agnostic; tissue-specific splicing extensions exist separately (AbSplice, MTSplice, SpliceTransformer, CI-SpliceAI) but none combine splicing predictions with phenotype or cell-type priors in a VUS classifier.
- **How IsoformVUS differs:** uses SpliceAI as a molecular feature but re-interprets it against the isoform used by the disease-relevant cell type (the counterfactual-isoform experiment).

#### AlphaMissense (Cheng et al., 2023)
- **Publication/year:** Cheng J et al. "Accurate proteome-wide missense variant effect prediction with AlphaMissense." *Science* 381(6664):eadg7492, 2023.
- **Datasets:** 71M missense predictions across 19,233 canonical proteins; plus a released expanded resource of predictions across 60,000 alternative transcript isoforms, explicitly described as "for future research and evaluation of isoform-specific effects."
- **Phenotype/cell-type/isoform:** scores are per canonical transcript in the main resource; isoform predictions exist as an unexploited side resource. No phenotype or cell-type usage.
- **Variant prioritization capability:** strong missense pathogenicity classifier (90% precision on ClinVar at chosen cutoffs).
- **Limitations:** missense-only; no indels, splice, noncoding, or structural variants; single 0.56 threshold calibrated on bulk ClinVar; no mechanistic output; isoform resource published but not used in the paper's clinical evaluations.
- **How IsoformVUS differs:** would be a *consumer* of the AlphaMissense isoform resource — selecting the isoform-specific score matching the disease-relevant cell type rather than the canonical one — inside a phenotype-aware model.

#### popEVE (Orenbuch et al., 2025)
- **Publication/year:** Orenbuch R, Shearer CA, Kollasch AW, Spinner AD, Hopf T, van Niekerk L, Franceschi D, Dias M, Frazer J, Marks DS. "Proteome-wide model for human disease genetics." *Nature Genetics*, 2025. DOI: 10.1038/s41588-025-02400-1.
- **Datasets:** cross-species evolutionary variation + UK Biobank + gnomAD population calibration; validated on ~31,058 severe developmental disorder (SDD) trios; 123 novel candidate genes.
- **Phenotype/cell-type/isoform:** none — deliberately phenotype-free; severity scale comparable across genes.
- **Variant prioritization capability:** best-in-class de novo missense prioritization in SDD (98% of known causal variants ranked top); severity scale distinguishes childhood-lethal from adult-onset.
- **Limitations:** missense only; no phenotype, tissue, cell-type, or isoform dimension.
- **How IsoformVUS differs:** orthogonal axis — popEVE asks "how severe is this variant anywhere"; IsoformVUS asks "is this variant pathogenic *in this phenotype and cell type via this isoform*."

### 3.2 Phenotype-aware variant prioritization (no isoform layer)

#### Exomiser (Smedley et al., 2015; active through 2025)
- **Publication/year:** Smedley D et al. "Next-generation diagnostics and disease-gene discovery with the Exomiser." *Nature Protocols* 10:2004–2015, 2015; reanalysis/benchmark updates through 2025 (Vestito et al., *npj Genomic Medicine* 2024/2025; Genome Medicine 2025 optimization study).
- **Datasets:** ClinVar, HPO, model-organism phenotypes, PPI networks, 100,000 Genomes Project cohorts.
- **Phenotype:** yes — hiPHIVE/PhenIX/PHIVE phenotype-similarity scoring; ~66% of diagnostic variants in top 10 with default hiPHIVE; automated ACMG/AMP classifier.
- **Cell-type/isoform:** none. MANE Select/Plus Clinical transcripts used for consistent reporting, not as model features.
- **Limitations:** gene-level; cannot resolve variants whose effect depends on the isoform used in the affected cell type; independent benchmarks (e.g., GIM Open 2024) show sensitivity drops on P/LP-only sets and note HPO sparsity as a failure mode.
- **How IsoformVUS differs:** Exomiser is the clinical reference to beat on phenotype-driven ranking; IsoformVUS adds cell-type and isoform layers Exomiser lacks, and targets VUS where phenotype alone is ambiguous.

#### LIRICAL (Kelly et al., 2022)
- **Publication/year:** Kelly C et al. *Genome Medicine*, 2022.
- **Datasets:** HPO, ClinVar, gnomAD frequencies.
- **Phenotype:** yes — genotype–phenotype likelihood-ratio model.
- **Cell-type/isoform:** none. Complements Exomiser in benchmarks (DDD + in-house cohorts).
- **How IsoformVUS differs:** same gap as Exomiser — no molecular cell/isoform context.

#### VARPP (Chapman & Lassmann, 2019) and IMPPROVE (Chapman & Lassmann, 2025)
- **VARPP publication/year:** "Personalised analytics for rare disease diagnostics." *Nature Communications*, 2019. DOI: 10.1038/s41467-019-13345-5 (cited as ref. 14 in IMPPROVE; authors confirmed via IMPPROVE's peer-review record).
- **IMPPROVE publication/year:** Chapman T, Lassmann T. "Single-cell data combined with phenotypes improves variant interpretation." *BMC Genomics*, 2025. DOI: 10.1186/s12864-025-11711-w. PMC12117811.
- **Datasets:** IMPPROVE: Phen2Gene HPO→gene mapping, ClinVar pathogenic variants, gnomAD benign variants, CADD v1.6 scores, GENCODE v26 exons, GTEx bulk RNA-seq, Tabula Sapiens single-cell RNA-seq. 1,866 HPO-term Random Forest models.
- **Phenotype:** yes — per-HPO-term models; expression-informed models beat CADD alone in 90% of phenotypes.
- **Cell-type:** yes — Tabula Sapiens expression aggregated by cell type and tissue; single-cell features beat bulk in 86% of models vs CADD, and beat GTEx in 67% of models. The authors explicitly position this as the first phenotype-specific tool enriched with single-cell data.
- **Isoform:** **none** — features are gene-level expression values (TPM by tissue/cell type). A variant is scored by the gene's expression in relevant cells, with no notion of which isoform the cell uses or whether the variant touches that isoform.
- **Limitations:** gene-level expression cannot distinguish a variant that hits a brain-dominant isoform from one that hits a blood-dominant isoform of the same gene; no isoform-usage or isoform-specific consequence features; no ablation of cell type vs isoform (isoform never modeled).
- **How IsoformVUS differs:** IMPPROVE is the closest published competitor — it shares the phenotype + cell-type + CADD + interpretable-model design. IsoformVUS adds exactly the missing layer: isoform usage in the relevant cell type and isoform-specific variant consequence. The ablation (cell-type model vs cell-type + isoform model) is designed to test precisely the value IMPPROVE leaves on the table.

#### V2P (2025)
- **Publication/year:** "Expanding the utility of variant effect predictions with phenotype-specific models." *Nature Communications*, 2025. DOI: 10.1038/s41467-025-66607-w. (Author list pending verification.)
- **Datasets:** HGMD, ClinVar, CADD, Capice comparisons; multi-task across 23 top-level HPO phenotypes.
- **Phenotype:** yes — phenotype-specific prediction heads improve causal-variant ranks over CADD/Capice across 411/690 dataset-phenotype pairs.
- **Cell-type/isoform:** none stated.
- **How IsoformVUS differs:** phenotype-aware baseline to compare against; lacks the cell-type and isoform axes.

### 3.3 Isoform-aware biology and methods (no phenotype-aware prioritization model)

#### Weykopf et al. (2026) — "Disease-associated genetic variants can cause missense effects in tissue-specific protein isoforms"
- **Publication/year:** Weykopf G, Badonyi M, Friman ET, Livesey BJ, Hird EF, Stanton CM, Vitart V, Su J, Murphy L, Baillie JK, Gorrell M, Marsh JA, Bickmore WA, Biddie S, Ioannou A, Coutts A, Wham M. *Nature Communications*, published 16 June 2026. (Exact DOI pending verification; coverage: https://conexiant.com/pathology/articles/alternative-isoforms-may-refine-variant-interpretation/ and Nature Communications subject listings.)
- **Datasets:** population variant catalogs + disease-associated variants + long-read RNA-seq + targeted full-length capture sequencing; functional validation of a lung-epithelium-specific DPP9 isoform (Leu8Pro in an alternative first exon, altered enzymatic activity, altered NLRP1 interaction).
- **Findings:** ~80% of alternative transcripts carrying disease-associated variants are absent from current reference annotations; variant-bearing alternative isoforms are more tissue-specific than their reference counterparts; structural modeling predicted the largest effects in ClinVar pathogenic variants, intermediate in VUS.
- **Phenotype:** disease variants mapped, but no phenotype→cell mapping model. **Cell-type:** tissue-level (lung epithelium). **Prioritization model:** none — an analysis and case-study framework.
- **Limitations (authors' own):** mostly computational predictions; one functional example (DPP9); healthy-tissue transcriptomics; missense-focused; many alternative transcripts may be NMD targets.
- **How IsoformVUS differs:** Weykopf et al. provide the empirical motivation and the kind of evidence (tissue-specific alternative isoforms with hidden coding consequences) IsoformVUS would operationalize. IsoformVUS turns this insight into a trained, benchmarked, phenotype- and cell-type-aware prioritization model with ablation — which Weykopf et al. do not attempt.

#### ESM1b isoform analysis (Brandes et al., 2023)
- **Publication/year:** Brandes N et al. "Genome-wide prediction of disease variant effects with a deep protein language model." *Nature Genetics*, 2023. DOI: 10.1038/s41588-023-01465-0.
- **Isoform content:** scored variants across primary vs alternative isoforms; ~1.8M isoform-sensitive variants across ~9,000 genes (85% of genes with curated alternative isoforms); enriched near splice junctions and in splicing-disrupted domains; e.g., 170 variants with LLR difference > 4 between p53 primary and Δ133p53β isoforms, including 3 ClinVar VUS.
- **Limitations:** no phenotype or cell-type context; isoform-sensitivity defined on sequence scores alone, not benchmarked as a clinical prioritization task.
- **How IsoformVUS differs:** ESM1b demonstrates the *signal* (isoform-dependent scores) at scale; IsoformVUS tests whether conditioning on the disease-relevant cell type's isoform usage converts that signal into better clinical prioritization.

#### Ly et al. (2025) — isoform-selective alleles via alternative start codons; SwissIsoform
- **Publication/year:** Ly J, Di Bernardo M, Tao YF, Khalizeva E, Giuliano CJ, Lourido S, Fleming MD, Cheeseman IM. "Alternative start codon selection shapes mitochondrial function and rare human diseases." *Molecular Cell* 85(22):4198–4214.e11, 2025. DOI: 10.1016/j.molcel.2025.10.013. PubMed: https://pubmed.ncbi.nlm.nih.gov/41205602/
- **Content:** nonsense/frameshift/missense variants between two translation-initiation sites selectively eliminate the longer protein isoform; isoform-specific TRNT1 mutations explain heterogeneous SIFD phenotypes; even a structurally benign conservative substitution (TRNT1 M30L) is pathogenic because it disrupts an alternative start codon. The Cheeseman lab is building "SwissIsoform," a clinician-facing tool to flag isoform-affecting mutations (in development, not released).
- **Relevance:** the canonical "isoform-selective allele" demonstration and the direct justification for a TRNT1 case study in IsoformVUS.
- **How IsoformVUS differs:** Ly et al. cover N-terminal (alternative-start) isoforms in mitochondrial biology; IsoformVUS is variant-class-general and phenotype/cell-type-general, with a benchmarked model rather than a per-gene analysis.

#### Liu et al. (2025) — genetic architecture of isoform variation
- **Publication/year:** Liu C, Joehanes R, Ma J et al. "Integrating whole genome and transcriptome sequencing to characterize the genetic architecture of isoform variation." *Nature Communications*, 2025. DOI: 10.1038/s41467-025-64336-8.
- **Content:** WGS + RNA-seq integration mapping genetic determinants of isoform diversity. Not a variant prioritizer; useful background on how common variation shapes isoform usage.
- **How IsoformVUS differs:** no clinical prioritization component.

### 3.4 Tissue/cell-specific splicing predictors

#### AbSplice (Wagner, Çelik et al., 2023)
- **Publication/year:** *Nature Genetics*, 2023. DOI: 10.1038/s41588-023-01373-3. Precomputed genome-wide scores (hg38/hg19) on Zenodo; UCSC track released April 2024.
- **Content:** probability that a variant causes aberrant splicing in each of 49 GTEx tissues, built on tissue-specific splice-site annotations (SpliceMaps) + MMSplice/SpliceAI features. AbSplice-RNA can additionally use patient RNA from accessible tissues.
- **Limitations:** tissue (not cell-type) resolution; splicing-only; no phenotype; not a disease-prioritization model.
- **How IsoformVUS differs:** AbSplice is a feature source (tissue-specific splicing) that IsoformVUS would consume and re-weight by cell-type relevance and isoform usage.

#### MTSplice (Cheng et al., 2021); SpliceTransformer (You et al., 2024); CI-SpliceAI (Strauch, 2022)
- **MTSplice:** Cheng J et al., *Genome Biology* 22:94, 2021 — tissue-specific Δψ splicing prediction.
- **SpliceTransformer:** You et al., *Nature Communications*, 2024 (DOI 10.1038/s41467-024-53088-6) — per-variant splicing scores plus Y/N tissue-usage-change flags for 15 tissues.
- **CI-SpliceAI:** open-source tool (`cispliceai`; Strauch, 2022) retraining SpliceAI on all GENCODE isoforms — the only explicitly *transcript-aware* splicing predictor found. Paper/preprint status needs verification.
- **Shared gap:** all are splicing-effect predictors; none is a phenotype-aware, cell-type-aware, isoform-aware disease-variant classifier.

### 3.5 The genome-scale molecular models (closest 2025–2026 entrants)

#### AlphaGenome (Linder J et al., 2025)
- **Publication/year:** "Advancing regulatory variant effect prediction with AlphaGenome." *Nature*, 2025. DOI: 10.1038/s41586-025-10014-0. (First author per Nature page; confirm at citation time.)
- **Content:** 1-Mb-context sequence-to-function model predicting TF binding, chromatin accessibility, histone marks, RNA abundance, and splicing across biosamples; improved eQTL magnitude/direction prediction over Borzoi/Enformer (Spearman ρ 0.39→0.49).
- **Limitations:** sequence-to-function foundation model, not a clinical prioritization model; no phenotype input.

#### AlphaGenome Atlas / AVI (Cheng J et al., 2026 preprint)
- **Publication/year:** *medRxiv*, 2026-09-16. DOI: 10.64898/2026.09.16.26363192. https://www.medrxiv.org/content/10.64898/2026.09.16.26363192v1 (Lead authors: Jun Cheng, Kyle R. Taylor, Lauren Nicolaisen, Joshua Pan, Clare Bycroft et al.; supervisors incl. Heidi Rehm, Anne O'Donnell-Luria, Pushmeet Kohli, Žiga Avsec.)
- **Content:** in-silico saturation mutagenesis of the human genome; unified interpretable AVI score (PHRED-scaled) combining AlphaGenome ISM predictions, AlphaMissense, VEP protein LOF annotations, PhastCons 470-way, Zoonomia Cactus 241-way; SHAP feature attributions per variant; trained on gnomAD v4.1 allele-frequency proxy labels (AF > 0.001 benign proxy; AF < 0.001 impactful proxy).
- **The isoform-relevant case:** an epileptic-encephalopathy proband (GREGoR cohort) whose top-ranked variant was a heterozygous noncoding VUS in intron 10 of DNM1 affecting a *brain-specific transcript isoform* — exon 10a absent from the MANE Select transcript; 69% of the AVI score attributed to the splicing feature; predictions pointed to altered splice-site usage in glutamatergic-neuron biosamples; blood RNA-seq had been inconclusive because exon 10a is not expressed in blood. This is the closest published demonstration of the IsoformVUS hypothesis in a real case.
- **Phenotype:** none — AVI is phenotype-free. **Cell-type:** biosample-level predictions (305 DNase / 371 RNA-seq tracks), but no explicit phenotype→cell-type relevance model.
- **Limitations:** preprint (not peer-reviewed at review time); proxy-label training (AF-based, not disease labels); noncoding-focused; no HPO phenotype input; no explicit isoform-usage weighting — the isoform insight emerges from track-level splicing predictions, not from modeling which isoform the relevant cell type uses.
- **How IsoformVUS differs:** AVI is the molecular-effects engine; IsoformVUS is the clinical-prioritization layer — phenotype-conditioned, ClinVar-labeled, cell-type-relevance-weighted, with an explicit isoform-usage feature and the critical with/without-isoform ablation AVI never runs.

### 3.6 Long-read and single-cell isoform infrastructure

- **LORALS (Glinos et al., 2022):** "Transcriptome variation in human tissues revealed by long-read sequencing." *Nature*, 2022 (DOI 10.1038/s41586-022-05035-y). Showed that supplementing GENCODE with novel long-read transcripts reassigns VEP consequence categories for rare variants (Extended Data Fig. 15) — direct evidence that annotation completeness changes variant interpretation.
- **Isosceles (2024):** "Accurate long-read transcript discovery and quantification at single-cell, pseudo-bulk and bulk resolution." *Nature Communications*, 2024 (DOI 10.1038/s41467-024-51584-3). Benchmarked best for single-cell/pseudobulk isoform quantification vs IsoQuant/FLAMES.
- **FLAMES (Holmqvist et al., 2021):** single-cell full-length isoform discovery/quantification (Bioconductor).
- **ISAHESVAL-NL (2025):** allele-level isoform analysis of heterozygous splicing variants from nanopore data (*Scientific Reports*, DOI 10.1038/s41598-025-14566-z) — PYGM/McArdle disease example.
- **RAPID (2026):** targeted long-read RNA workflow for functional resolution of splicing VUS in rare disease (*Genome Medicine*, DOI 10.1186/s13073-026-01754-3) — clinical workflow, not a predictor.
- **Relevance to IsoformVUS:** these are the data-generation and quantification methods that make cell-type isoform usage measurable; none is a variant prioritization model.

### 3.7 Standards and transcript infrastructure

- **MANE (Morales et al., 2022):** "A joint NCBI and EMBL-EBI transcript set for clinical genomics and research." *Nature*, 2022 (DOI 10.1038/s41586-022-04558-8). The MANE authors themselves warn that "selection of one transcript does not imply that the rich biology of the human genome can be reduced to one transcript at each locus" and that tissue-specific/developmental mechanisms are "not within the scope of the MANE Select set" — the standardization that creates the isoform blind spot.
- **MANE Plus Clinical evaluation (preprint, Exeter group):** systematic evaluation of 65 genes with both MS and MPC transcripts; ~2.5% of (likely) pathogenic variants mapped *exclusively* to MPC-specific coding sequence, rising to ≥25% for 7/65 genes; >4M ClinVar variants contain clinically relevant exons outside current MANE annotations (e.g., PCDH15 alternative cytoplasmic-domain exon with 14 P/LP variants). Preprint: researchsquare.com/article/rs-9680068. **This is the strongest quantitative estimate of how many pathogenic variants a single-transcript pipeline misses.**

### 3.8 Cell-type-aware rare-variant methods (complex traits — scope note)

- **cellSTAAR (2025):** *Nature Methods* (DOI 10.1038/s41592-025-02919-5). Integrates WGS with scATAC-seq cell-type annotations to boost rare-variant association testing for noncoding regions (TOPMed ~60k, UK Biobank ~190k; lipid traits). Cell-type-aware, but for complex-trait association testing, not Mendelian VUS prioritization.

### 3.9 Phenotype-benchmark literature (for evaluation design)

- **Jacobsen et al. (2022):** "Phenotype-driven approaches to enhance variant prioritization and diagnosis of rare disease." (PMC9288531) — Exomiser recalls 82.6–93.6% of 4,877 diagnoses in top 1–10; LIRICAL close behind.
- **Genome Medicine (2025):** "An optimized variant prioritization process for rare disease diagnostics: recommendations for Exomiser and Genomiser" (DOI 10.1186/s13073-025-01546-1) — phenotype algorithm choice (hiPHIVE best: 66.6% top-10) and HPO quality materially affect rankings.
- **Evaluation of phenotype-driven gene prioritization (2022):** LIRICAL and AMELIE complementary on 305 DDD + 209 in-house cases (PMC8921623).
- **"A phenotype centric benchmark of variant prioritisation tools" (2018, npj Genomic Medicine):** HPO-term filtering methodology (≥25 ClinVar pathogenic variants per term) that IMPPROVE and any HPO-stratified benchmark should follow.

---

## 4. Closest competitors — ranked

1. **IMPPROVE (2025).** Shares the most design DNA with IsoformVUS: HPO-specific, interpretable (Random Forest), CADD-based, bulk + single-cell expression features, ClinVar/gnomAD labels. The single cleanest contrast: IMPPROVE's expression features are gene-level; IsoformVUS adds the isoform-usage and isoform-specific-consequence layers IMPPROVE never models. The ablation "cell-type model vs cell-type + isoform model" directly extends IMPPROVE's "bulk vs single-cell" experiment one level deeper.
2. **AlphaGenome Atlas / AVI (2026 preprint).** The strongest molecular-effects competitor: genome-wide, biosample-aware splicing, SHAP-interpretable, and already demonstrated an isoform-selective rare-disease VUS rescue (DNM1 brain-specific exon 10a). It is phenotype-free, proxy-label-trained, and noncoding-focused — IsoformVUS is its phenotype-conditioned, ClinVar-labeled, coding-inclusive complement, not a replacement.
3. **Weykopf et al. (2026).** The empirical manifesto for the hypothesis: tissue-specific alternative isoforms hide coding consequences of disease variants; DPP9 lung-isoform case. Not a model — the evidence base IsoformVUS operationalizes.
4. **V2P (2025).** The other phenotype-specific deep predictor; no cell-type or isoform axes. Useful as a phenotype-aware baseline alongside Exomiser.
5. **Exomiser (2015–2025).** The clinical reference standard for phenotype-driven ranking. Gene-level; the benchmark IsoformVUS must report top-k recovery against.
6. **Ly et al. (2025) / SwissIsoform.** The isoform-selective-allele proof (TRNT1 SIFD); scope-limited to alternative start codons; the planned tool is not released, leaving the general isoform-selective prioritization space open.

---

## 5. What already exists vs what is new

**Already exists (do not re-claim):**
- Phenotype-aware variant prioritization (Exomiser 2015, LIRICAL 2022, VARPP 2019, V2P 2025).
- Adding single-cell/tissue expression to prioritization models (VARPP 2019; IMPPROVE 2025).
- Tissue-specific splicing-effect prediction (MTSplice 2021, AbSplice 2023, SpliceTransformer 2024).
- Transcript-aware splicing scoring (CI-SpliceAI 2022).
- Isoform-dependent missense scoring (ESM1b/Brandes 2023; AlphaMissense 60k-isoform resource 2023).
- Isoform-selective disease alleles as a concept (Ly et al. 2025; Weykopf et al. 2026; DNM1 case in AVI 2026).
- Quantitative evidence that single-transcript pipelines miss pathogenic variants (~2.5% of P/LP variants exclusively in MANE Plus Clinical-specific sequence; ≥25% for 7/65 genes — preprint).

**Not found (as of 2026-10-06) — the defensible gap:**
No published method combines, in one interpretable, benchmarked model: (i) HPO phenotype conditioning, (ii) cell-type relevance derived from single-cell data, (iii) isoform usage in the disease-relevant cell type, and (iv) isoform-specific functional consequence — with gene-level splits, an ablation isolating the isoform layer's contribution, and benchmarking against CADD/SpliceAI/AlphaMissense plus phenotype-aware tools (Exomiser, IMPPROVE). IMPPROVE has (i)+(ii) but no isoform layer; AVI has a form of (ii)+(iii) at the molecular level but no phenotype conditioning or disease-label training; Weykopf et al. and Ly et al. demonstrate the biology but build no prioritization model.

---

## 6. Implications for IsoformVUS design (from the literature)

1. **Use gene-level splits, not variant-level splits.** IMPPROVE trains per-HPO-term models on variants pooled across genes; a follow-up PhD-grade project must go further (gene- or disease-level splits) to avoid memorizing gene-specific expression patterns.
2. **Do not train on ClinVar if the baseline was.** IMPPROVE's choice of CADD (ClinVar-independent) as the base score is the right circularity hygiene; replicate it.
3. **Expect the isoform layer to matter most in a subgroup.** The MANE Plus Clinical estimate (~2.5% of P/LP variants) and the Weykopf 80%-unannotated figure suggest the isoform-selective subgroup is small but real — design the benchmark to evaluate this subgroup separately rather than expecting genome-wide AUROC jumps.
4. **Short-read ambiguity is a documented weakness.** Weykopf et al. and the isoform-quantification literature agree that full-length long-read data (Isosceles/FLAMES) is the ground truth for isoform usage; short-read estimates need uncertainty handling.
5. **The space is moving fast.** Three directly relevant works (V2P, Weykopf, AVI) appeared within ~9 months of this review. Re-run the novelty check before manuscript submission, and cite the preprints as preprints.

---

## 7. Source URLs

- IMPPROVE: https://pmc.ncbi.nlm.nih.gov/articles/PMC12117811/ (DOI 10.1186/s12864-025-11711-w)
- AlphaGenome Atlas / AVI: https://www.medrxiv.org/content/10.64898/2026.09.16.26363192v1
- AlphaGenome (Nature 2025): https://www.nature.com/articles/s41586-025-10014-0
- Weykopf et al. coverage: https://conexiant.com/pathology/articles/alternative-isoforms-may-refine-variant-interpretation/ ; Nature Communications listing: https://www.nature.com/subjects/genetic-variation/ncomms
- V2P: https://www.nature.com/articles/s41467-025-66607-w
- Ly et al.: https://pubmed.ncbi.nlm.nih.gov/41205602/ ; https://www.cell.com/molecular-cell/fulltext/S1097-2765(25)00854-8
- popEVE: https://www.nature.com/articles/s41588-025-02400-1
- AlphaMissense: https://www.science.org/doi/10.1126/science.adg7492
- ESM1b/Brandes: https://www.nature.com/articles/s41588-023-01465-0
- SpliceAI: https://doi.org/10.1016/j.cell.2018.12.015 ; https://github.com/Illumina/SpliceAI
- AbSplice: https://www.nature.com/articles/s41588-023-01373-3 ; https://github.com/gagneurlab/absplice
- SpliceTransformer: https://www.nature.com/articles/s41467-024-53088-6
- CADD: https://shendurelab.github.io/documents/rentzsch_nucleic_2019.pdf
- CADD-Splice: https://link.springer.com/article/10.1186/s13073-021-00835-9
- Region-aware CADD thresholds (2025): https://academic.oup.com/nargab/article/7/4/lqaf157/8328387
- Exomiser reanalysis (npj Genomic Medicine): https://www.researchgate.net/publication/387176821_Efficient_reinterpretation_of_rare_disease_cases_using_Exomiser
- Exomiser/Genomiser optimization (Genome Medicine 2025): https://link.springer.com/article/10.1186/s13073-025-01546-1
- Jacobsen et al. 2022 (phenotype-driven approaches): https://pmc.ncbi.nlm.nih.gov/articles/PMC9288531/
- MANE (Nature 2022): https://pmc.ncbi.nlm.nih.gov/articles/PMC9007741/
- MANE Plus Clinical evaluation (preprint): https://www.researchsquare.com/article/rs-9680068/v1.pdf
- LORALS long-read transcriptome (Nature 2022): https://www.nature.com/articles/s41586-022-05035-y
- Isosceles (2024): https://www.nature.com/articles/s41467-024-51584-3
- FLAMES (Bioconductor): https://master.bioconductor.org/packages/3.19/bioc/manuals/FLAMES/man/FLAMES.pdf
- ISAHESVAL-NL (2025): https://www.nature.com/articles/s41598-025-14566-z
- RAPID (2026): https://link.springer.com/article/10.1186/s13073-026-01754-3
- Liu et al. isoform variation (2025): https://doi.org/10.1038/s41467-025-64336-8
- cellSTAAR (Nature Methods 2025): https://link.springer.com/article/10.1038/s41592-025-02919-5
- ClinVar: https://www.ncbi.nlm.nih.gov/clinvar/ ; news/updates: https://www.ncbi.nlm.nih.gov/clinvar/docs/news/
- IMPPROVE peer-review record (VARPP lineage): https://static-content.springer.com/openpeerreview/art%3A10.1186%2Fs12864-025-11711-w/12864_2025_11711_ReviewerReport_V2_R1.pdf

---

*End of literature review. Companion file: novelty_statement.md*
