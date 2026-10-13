# Leakage prevention

## Threat model

The central leakage risk: variants from the *same gene* share gene-level
properties (constraint, expression, disease mechanism). A random variant-level
split lets the model memorize "variants in gene X are usually pathogenic" and
rewards it at test time — inflating metrics without teaching variant-level
biology. Related risks: same-disease genes sharing pathways; ClinVar submitters
batch-labeling; paralogs.

## Strategy (mandatory)

1. **Gene-level grouped splits.** Genes — never individual variants — are assigned
   to TRAIN / VALIDATION / HELD-OUT TEST (60/20/20 by gene count, stratified by
   disease category). All variants of a gene stay in one split. Implementation:
   `src/isoformvus/data/splits.py`, seeds fixed in `config/base.yaml`.
2. **Held-out test is locked.** Split assignment is written once
   (`results/splits.json`) before any model training. Test variants are never
   used for feature selection, hyperparameter tuning, early stopping, or
   threshold choice. Any accidental access is logged as a protocol deviation.
3. **Disease-level sensitivity split.** As a secondary analysis, re-split by
   disease category (all genes of one disease family held out) to test
   cross-disease generalization. Reported alongside, not instead of, the primary.
4. **No VUS in training labels.** VUSs are never a training class; they are the
   *application* target. (A separate descriptive analysis may rank held-out VUSs;
   those rankings are not evaluation metrics.)
5. **Predictor-score hygiene.** CADD/SpliceAI/AlphaMissense were trained on data
   that overlaps ClinVar. We do not retrain them; we use their scores as
   features/baselines and *disclose* that baseline comparisons inherit their
   training-data overlap. The critical comparison (same model ± isoform features)
   is leakage-free by construction since both arms see identical splits.
6. **Paralog check.** Held-out genes' close paralogs (Ensembl paralogy) in train
   are flagged; a sensitivity analysis drops them to test memorization via
   homology.

## Verification

- `tests/test_splits.py` asserts: no gene appears in two splits; no variant key
  appears in two splits; stratification balance within tolerance; determinism
  across runs with the same seed.
- The final audit (`docs/final_audit.md`) re-verifies split integrity from
  `results/splits.json` independently of the training code.
