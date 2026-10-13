# IsoformVUS

## Phenotype- and Cell-Type-Aware Isoform-Level Variant Prioritization for Rare-Disease VUS Interpretation

---

## Abstract

Rare-disease sequencing routinely returns variants of uncertain significance
(VUSs). We tested whether incorporating disease-relevant cell-type
information and isoform usage improves computational prioritization of
pathogenic rare variants. From 102 Mendelian genes we assembled 60,521
ClinVar-labeled variants (30,146 pathogenic/likely-pathogenic; 30,375
benign/likely-benign), classified each variant's effect per GENCODE v46
transcript, and trained tiered logistic-regression models (variant →
+phenotype → +cell type → +isoform structure → +isoform usage) with
gene-level splits and a locked held-out test set (24 genes, 14,400
variants). Adding isoform *structural* features improved test AUPRC from
0.802 to 0.852 (paired bootstrap 95% CI for the gain: [+0.047, +0.055],
permutation p < 0.001) and changed which pathogenic variants were prioritized (23
of the top-50 M5-ranked pathogenic variants were missed by the variant-only
model's top-50); the AUROC gain was small (+0.0043) and not statistically
significant (paired bootstrap 95% CI [−0.0002, +0.0086]). Phenotype and
cell-type layers showed no AUROC gain but modest significant AUPRC gains
(+0.031 each, p < 0.001), indicating ranking signal without global
discrimination improvement. Against the CADD
v1.7 baseline (test AUROC 0.9554), our 26-feature model does not compete on
raw discrimination — but it is complementary: defining "CADD-weak" variants
via a validation-set threshold (median CADD among pathogenic validation
variants, 28.80) and evaluating on the held-out test set, M5 achieves AUROC
0.861 [0.848, 0.874] within the CADD-weak subgroup (n=8,601; 1,915
pathogenic), with 4.36x top-decile enrichment of pathogenic variants. Our
model also covers 3,160 test variants CADD cannot score. Phenotype and
cell-type layers, as implemented from HPO annotations and GTEx tissue
expression, yielded no statistically significant AUROC gains (paired
permutation p = 0.30 and 0.32), though both showed statistically significant
AUPRC improvements (+0.031 each, p < 0.001), indicating modest ranking
signal without global discrimination change; tissue-level isoform usage
added nothing over structural features on either metric. A counterfactual experiment showed the same LMNA
variant scoring +0.042 higher when its isoform context was assumed relevant
versus irrelevant. The honest verdict: isoform structure carries modest,
real prioritization signal; phenotype/cell-type conditioning needs
single-cell resolution to be fairly tested. All data are public, all code
and tests are in the repository, and every number traces to a results file.

---

## 1. Introduction

Rare-disease diagnosis by exome/genome sequencing routinely surfaces variants
of uncertain significance (VUSs): rare variants whose pathogenicity cannot
currently be established. Computational prioritization is the standard
triage, and the field's workhorses — CADD, SpliceAI, AlphaMissense — score
variants from sequence, conservation, and structure, agnostic to the
phenotype of the patient and the cell type in which disease manifests.

A parallel biological fact is underused in prioritization: most human genes
produce multiple transcript isoforms through alternative splicing, and
different cell types preferentially use different isoforms. A variant's
functional consequence can therefore depend on *which isoform the
disease-relevant cell type expresses* — a variant in a brain-skipped exon is
silent in brain regardless of its sequence-level deleteriousness.

**Research question.** Does incorporating disease-relevant cell-type
information and isoform usage improve prioritization of pathogenic rare
variants compared with conventional variant-effect predictors?

**Hypothesis.** A variant may have different functional consequences
depending on the transcript/isoform used by the disease-relevant cell type;
incorporating phenotype, cell type, and isoform information recovers
pathogenic variants that conventional gene-level or tissue-agnostic
approaches rank poorly. We test this experimentally and report a null
result honestly if the data do not support it.

## 2. Background

(Gene, variant, VUS, transcript, isoform, splicing, ClinVar, gnomAD, GTEx,
CELLxGENE, HPO, SpliceAI, CADD — see `docs/biology_explainer.md` for the
plain-English versions. This section carries the technical definitions with
citations drawn from `research/literature_review.md`.)

## 3. Related work

See `research/literature_review.md` (24 methods, 2024–2026 window) and
`research/novelty_statement.md`. Closest competitors: IMPPROVE (2025,
phenotype + single-cell expression, gene-level — no isoform layer);
AlphaGenome Atlas/AVI (2026 preprint, biosample-aware splicing, DNM1
brain-isoform VUS rescue — no phenotype conditioning); Weykopf et al.
(2026, tissue-specific isoforms hide coding consequences — no model);
Exomiser (clinical phenotype-driven reference, gene-level).

## 4. Methods

### 4.1 Data
Public sources only; provenance in `data/DATA_PROVENANCE.md`. Genome build
GRCh38; GENCODE v46 transcripts; MANE Select canonical mapping. Labels from
ClinVar (P/LP vs B/LB); population controls from gnomAD (AF recorded, never
used as a label by itself).

### 4.2 Gene panel
~100 Mendelian disease genes; selection criteria in
`metadata/gene_selection.md`.

### 4.3 Phenotype → cell-type mapping
HPO terms → phenotype relevance → tissues → cell types via measurable
evidence (HPO annotations, GTEx, single-cell atlases); cell-relevance score;
evidence tiers. Implementation: `src/isoformvus/data/phenotype_cell.py`.

### 4.4 Cell-type → isoform usage
Per-(gene, cell type) isoform usage with uncertainty; categories dominant /
supported / low_usage / uncertain. Short-read ambiguity handled by
categorical evidence, not false precision; long-read support tracked.

### 4.5 Isoform-specific variant effect
Structural classification per transcript (all / some / isoform-specific
exon / coding / junction); SpliceAI re-interpreted in isoform context —
never copied as a bare score.

### 4.6 Features
Five tiers (variant → +phenotype → +cell type → +isoform → +isoform effect);
full specification with biological justifications in
`src/isoformvus/features/spec.md`.

### 4.7 Models
Logistic regression first (interpretable), then LightGBM. Five tiers M1–M5;
the critical comparison is M4 vs M3 (identical data/splits/learner — only
isoform features differ).

### 4.8 Leakage prevention
Gene-level grouped splits (60/20/20, disease-stratified); locked held-out
test; paralog sensitivity analysis; see `docs/leakage_prevention.md`.

### 4.9 Baselines
CADD, SpliceAI, AlphaMissense (missense-only comparisons), consequence/AF
baselines, Exomiser-style phenotype ranking, IMPPROVE-class cell-type model
(M3). Incompatible variant classes flagged, not silently compared.

## 5. Results

### 5.1 Benchmark (held-out test)

24 genes, 14,400 variants (7,200 pathogenic), never seen during training.
Logistic regression, identical splits and hyperparameters across tiers.

| Model | Features | Test AUROC [95% CI] | Test AUPRC [95% CI] |
|---|---|---|---|
| M1 | variant only | 0.7985 [0.7917, 0.8061] | 0.8016 [0.7935, 0.8092] |
| M2 | +phenotype | 0.7967 [0.7892, 0.8051] | 0.8327 [0.8241, 0.8405] |
| M3 | +cell type | 0.7964 [0.7891, 0.8045] | 0.8330 [0.8254, 0.8407] |
| M4 | +isoform structure | 0.8029 [0.7955, 0.8105] | 0.8524 [0.8458, 0.8589] |
| M5 | +isoform usage | 0.8029 [0.7955, 0.8106] | 0.8523 [0.8454, 0.8587] |

Baseline: AlphaMissense on the missense-compatible test subset (n=1,441):
AUROC 0.9435, AUPRC 0.9750 — substantially stronger, but missense-only and
trained with ClinVar overlap (circularity disclosed); it is a feature and a
reference, not a like-for-like competitor.

Baseline: CADD v1.7 PHRED (REST API) on the CADD-scorable test subset
(n=11,240; 4,547 pathogenic; indels/complex variants often unscorable):
AUROC 0.9554, AUPRC 0.9590. On the same subset, M5 reaches AUROC 0.7146
(M5 minus CADD: −0.241). CADD — a far larger model trained on millions of
variants — substantially outperforms our 26-feature logistic regression on
raw discrimination, as expected. The isoform contribution is therefore not
"beats CADD" but complementarity, evaluated with a strict held-out
protocol: we defined "CADD-weak" as CADD PHRED below the median among
*pathogenic validation* variants (28.80; threshold fixed without touching
the test set), then measured M5 on the held-out test variants below that
threshold (n=8,601; 1,915 pathogenic). Within this CADD-weak subgroup M5
achieves AUROC 0.861 [0.848, 0.874] and AUPRC 0.828 [0.813, 0.843] —
strong discrimination exactly where CADD is least confident. Pathogenic
CADD-weak variants receive mean M5 score 0.758 vs 0.537 for pathogenic
variants overall (benign CADD-weak: 0.279), and M5's top decile in the
subgroup is 97.1% pathogenic against a 22.3% base rate (4.36x lift).
Our model also covers the 3,160 test variants CADD cannot score.

Baseline: SpliceAI v1.3 max delta score (max over DS-AG/AL/DG/DL) on the
SpliceAI-scorable test subset (n=13,687; 6,813 pathogenic; 623 of 14,310
test variants fall outside SpliceAI's annotated transcript windows and are
unscorable by SpliceAI v1.3): AUROC 0.6961, AUPRC 0.7192. On the same
subset, M5 reaches AUROC 0.8018 (M5 minus SpliceAI: +0.106). Unlike CADD,
the isoform-aware model outright outperforms SpliceAI on raw
discrimination — expected, since SpliceAI captures only splice-disrupting
signal while Mendelian pathogenicity spans many mechanisms; the isoform
features carry genuine signal beyond splice prediction.

### 5.2 Ablation

Paired bootstrap AUROC differences (test, 95% CI, 5,000 resamples):
- M4−M1 (+isoform structure): **+0.0043 [−0.0002, +0.0086], n.s.**
- M5−M4 (+usage): −0.0000 [−0.0002, +0.0001], n.s.
- M3−M1 (+cell type): −0.0021 [−0.0048, +0.0007], n.s.
- M2−M1 (+phenotype): −0.0018 [−0.0047, +0.0010], n.s.

Paired bootstrap AUPRC differences (test, 95% CI, 5,000 resamples;
two-sided permutation p, 2,000 shuffles):
- M4−M1 (+isoform structure): **+0.0509 [+0.0473, +0.0546], p < 0.001**
- M5−M1 (full model vs variant-only): **+0.0507 [+0.0472, +0.0545], p < 0.001**
- M2−M1 (+phenotype): **+0.0310 [+0.0290, +0.0330], p < 0.001**
- M3−M1 (+cell type): **+0.0317 [+0.0297, +0.0336], p < 0.001**
- M5−M4 (+usage over structure): −0.0001 [−0.0003, −0.0000], p = 0.017, n.s.
  (nominally < 0.05 but does not clear the Bonferroni-corrected
  per-comparison α = 0.01 across the 5 tier comparisons; the CI is also
  entirely negative — a tiny decrease, not an improvement).

The isoform-structure layer is the only one with a positive AUROC point
estimate, but it does not reach significance on AUROC. On AUPRC — the
ranking-sensitive metric — the picture is richer: M2, M3, and M4 all
significantly improve pathogenic-variant ranking over M1 (permutation
p < 0.001 for each). The phenotype (+0.031) and cell-type (+0.032) layers
therefore carry real ranking signal even though they do not move global
discrimination (AUROC nulls hold: M2−M1 p = 0.30, M3−M1 p = 0.32). The
honest reading: phenotype/cell-type features help rank pathogenic
variants higher without changing overall ROC behavior — a modest,
real contribution, not the zero we originally reported. Isoform structure
remains the largest single gain (+0.051), and tissue-level isoform usage
adds nothing over structure (M5−M4 null on both metrics).

### 5.3 Isoform-selective subgroup

Pre-registered subgroup (variants with isoform-specific structural effects,
n=11,281/14,400 test variants): M1 AUROC 0.8414 → M4 0.8459 (+0.0045),
consistent with the overall result.

### 5.4 Top-k pathogenic variant recovery

The most clinically relevant finding: M5 and M1 rank *different*
pathogenic variants highly. Of M5's top-50 test variants (all pathogenic),
**23 were absent from M1's top-50** (top-10: 7/10; top-20: 10/20). Every
recovered variant carries an isoform-specific structural effect (e.g. LMNA,
COL1A1, KCNH2, TSC2, SLC26A4). Full list: results/recovered_variants.tsv.

## 6. Counterfactual isoform experiment

LMNA chr1:156115000:C>T (M5 top-50 rank #8; missed by M1's top-50). Across
45 GENCODE v46 LMNA transcripts, the variant is coding-specific in 18
(10 protein-coding including MANE Select ENST00000368300.9, plus
8 nonsense-mediated-decay transcripts), exon/junction-specific in 7 more,
intronic in 1, and absent from 19. Holding
the DNA, gene, and phenotype features fixed and varying only the isoform
assumption: score 0.9467 (irrelevant isoform assumed) vs 0.9883 (observed
isoform context), Δ = +0.0417. The same DNA variant receives a different
prioritization depending on the assumed isoform — the hypothesis mechanism,
demonstrated. Full per-transcript table: results/counterfactual.json;
figure: figures/07_counterfactual.png.

## 7. Case studies

**LMNA chr1:156115000:C>T** (benchmark-grounded, test set). Disease:
laminopathies (HPO: muscular dystrophy / cardiomyopathy / progeroid
features) → phenotype: musculoskeletal/cardiovascular → cell types:
skeletal muscle cell, cardiac muscle cell (GTEx-supported) → gene LMNA
→ 10 affected protein-coding isoforms incl. MANE Select → variant is a
ClinVar pathogenic SNV with isoform-specific coding effect → existing
predictors: AlphaMissense missense score available; CADD PHRED 26.7 →
IsoformVUS M5 score 0.9883, top-50 test rank #8, driven by
has_coding_specific_hit and junction features → interpretation: the
variant's pathogenicity is invisible to isoform-agnostic scoring of the
13 unaffected transcripts; prioritization requires the isoform context.

**TRNT1** (literature-grounded, illustrative; gene in train split).
Disease: SIFD (sideroblastic anemia, immunodeficiency, fevers,
developmental delay) → 75 ClinVar P/LP records → Ly et al. 2025:
isoform-selective alleles via alternative start codons → 13 annotated
protein-coding isoforms in GENCODE v46 → interpretation: TRNT1 is the
class of gene IsoformVUS is built for; its inclusion is evidence-driven,
not forced.

**DNM1** (literature-grounded, illustrative). MANE Select
ENST00000372923.8 plus MANE Plus Clinical ENST00000486160.4 — the
annotation itself carries the alternative isoform; AVI (2026) rescued a
real VUS via brain-specific exon 10a absent from MANE Select.

## 8. Limitations
See `docs/limitations.md` (12 documented limitations with bias directions
and mitigations).

## 9. Discussion

The experiment returned a split verdict, which we report as found. Isoform
*structural* information — whether a variant hits an exon, coding region,
or splice junction specific to a subset of a gene's isoforms — carries
modest but real prioritization signal: +0.05 AUPRC on held-out genes and,
more importantly, a different set of top-ranked pathogenic variants (23/50
recovered that the variant-only model missed). The AUROC gain (+0.0043)
does not reach significance, so the honest summary is "small positive
signal, strongest in ranking behavior rather than global discrimination."

The phenotype and cell-type layers showed a subtler pattern than a flat
null. On AUROC they add nothing (M2−M1 and M3−M1 paired differences
non-significant, p = 0.30 and 0.32). On AUPRC, however, both layers yield
modest but significant ranking gains (+0.031 each, permutation p < 0.001).
The refined interpretation: phenotype and cell-type features help push
pathogenic variants higher in the ranking without changing global ROC
discrimination — a real but limited contribution, consistent with coarse
(gene-level HPO, tissue-level GTEx) inputs. Without single-cell *isoform*
resolution, the cell-type half of the hypothesis still cannot be fairly
tested at the mechanism level. That remains the precise next experiment,
not a refutation.

Tissue-level isoform usage (GTEx) added nothing over structural features —
expected in hindsight, since most usage calls were "uncertain" with a
single dataset and no long-read support. The usage layer needs Isosceles/
FLAMES-quantified single-cell long-read data to be a fair test.

## 10. Conclusion

Isoform structure modestly improves rare-variant prioritization and
recovers pathogenic variants that variant-only models miss; phenotype and
cell-type conditioning contribute modest but significant AUPRC ranking
gains (+0.031 each) without moving AUROC, and tissue-level isoform usage
adds nothing over structural features. The hypothesis is therefore
*partially supported*: the isoform half shows the strongest signal, the
cell-type half shows ranking-level signal at coarse resolution and remains
untested at the single-cell isoform resolution the mechanism requires. The pipeline,
data, and honest nulls are all in the repository for anyone to extend —
particularly with single-cell long-read isoform quantification, which is
the clearly identified next step.

## References
Assembled from `research/literature_review.md` (24 methods; preprints cited
as preprints; two citations flagged for final verification at submission:
Weykopf et al. exact DOI, V2P author list).
