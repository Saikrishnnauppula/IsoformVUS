# Clinical safety

IsoformVUS is a **research project**, not a clinical diagnostic tool.

## Hard rules

1. **No diagnosis.** Nothing in this repository diagnoses a patient or determines
   clinical care. Scores are research prioritization signals.
2. **No automatic ACMG classification.** Model output must never be converted into
   ACMG/AMP pathogenicity classifications (Pathogenic / Likely Pathogenic / VUS /
   Likely Benign / Benign) without independent expert review and orthogonal evidence.
3. **Three distinct things, never conflated:**
   - *computational evidence* — what the model computed and why (features, scores);
   - *research prioritization* — a ranked shortlist for follow-up experiments;
   - *clinical interpretation* — a judgment made by qualified professionals with
     patient context, family data, and functional validation.
4. **No clinical language in outputs.** Explanations say "prioritized for follow-up"
   or "consistent with further investigation", never "diagnostic of" or "causal".
5. **Uncertainty is surfaced, not hidden.** Every ranked variant carries its
   isoform-usage uncertainty category and data-support level (see
   `docs/isoform_analysis.md`).

## Intended use

- Hypothesis generation for rare-disease variant interpretation research.
- Benchmarking whether isoform/cell-type context adds signal over conventional
  predictors.
- Teaching material for computational genomics (with the biology explainer).

## Misuse

Do not deploy IsoformVUS scores in a clinical pipeline, do not present them to
patients or families as results, and do not cite them as standalone evidence in a
diagnostic report.
