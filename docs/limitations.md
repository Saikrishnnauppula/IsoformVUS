# Scientific limitations

No weaknesses hidden. Each limitation below names the failure mode, its
expected direction of bias, and the mitigation or diagnostic in place.

## 1. Transcript annotation errors and incompleteness
GENCODE v46 is the best available annotation and still incomplete: Weykopf et
al. (2026) estimate ~80% of alternative transcripts carrying disease variants
are absent from reference annotations. Direction: we *underestimate* the
isoform-selective subgroup. Mitigation: long-read datasets as supporting
evidence; MANE Plus Clinical regions flagged; limitation quantified, not
hand-waved.

## 2. Uncertainty in isoform expression
Short-read scRNA-seq cannot reliably resolve full-length isoforms (3' bias,
multi-mapping ambiguity). Direction: noisy usage estimates attenuate the
isoform signal toward zero (bias toward the null). Mitigation: categorical
usage with explicit uncertainty widths (see `data/isoform.py`); low expression
is never equated with irrelevance; long-read support tracked per call.

## 3. Short-read limitations
Covered in (2); additionally, junction-spanning reads are sparse at scRNA-seq
depth, so splice-junction usage is the least certain quantity we use.

## 4. Incomplete single-cell coverage
Not every phenotype-relevant tissue has a high-quality annotated single-cell
dataset; rare cell types are under-sampled. Direction: cell-relevance scores
for under-sampled types are conservative (low evidence, not low relevance).
Dataset-of-origin is kept in provenance; genes lacking cell-type data are
excluded from the panel rather than imputed.

## 5. Cell-type annotation errors
"Excitatory neuron" labels differ across atlases. Mitigation: Cell Ontology
mapping, multi-dataset aggregation, evidence tiers (curated > expression >
inferred).

## 6. Phenotype ambiguity
HPO terms are broad ("Seizure" HP:0001250); patients are described with
variable granularity. Direction: noisy phenotype input dilutes the phenotype
layer. Mitigation: information-content weighting of HPO terms; ablation
(Model 2 vs Model 1) measures whether the phenotype layer survives the noise.

## 7. ClinVar label bias
Well-studied genes, European ancestries, and submitters' practices dominate.
VUSs are not a random sample. Direction: optimistic performance on
well-characterized biology; poor generalization to understudied genes.
Mitigation: gene-level splits, disease-stratified sampling, explicit
discussion; VUSs never used as training labels.

## 8. Population-control limitations
gnomAD AF>1% variants are *controls*, not certified benign; rare ≠ pathogenic
and common ≠ benign are respected throughout. Control type is a recorded
column, and ClinVar B/LB vs population-common controls are analyzed
separately.

## 9. Model bias
Logistic regression assumes linear log-odds; LightGBM can overfit small
gene panels. Mitigation: gene-level splits, validation-based model choice,
calibration reporting (Brier score), and a pre-registered decision rule for
whether a neural network is even attempted (only if both learners underfit
validation).

## 10. Dataset leakage
Covered in `docs/leakage_prevention.md`: gene-level grouped splits, locked
held-out test, paralog sensitivity analysis, predictor-score hygiene
(CADD chosen for ClinVar-independence per IMPPROVE's precedent).

## 11. Generalization limitations
~100 Mendelian genes cannot represent all rare-disease biology; somatic,
mitochondrial, and repeat-expansion mechanisms are out of scope by design.
Claims are scoped to the evaluated panel and variant classes.

## 12. Moving literature
Three directly relevant works appeared within nine months of the 2026-10-06
review (V2P, Weykopf et al., AlphaGenome Atlas/AVI). The novelty check must be
re-run before manuscript submission; preprints cited as preprints.
