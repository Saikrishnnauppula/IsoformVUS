# Final scientific audit

**Audit completed:** 2026-10-12. Each check performed independently where
feasible; evidence recorded as file + command, not vibes.
Items marked [x] verified; details below.

## Biology
- [x] Every biological assumption in `docs/biology_explainer.md` re-checked
      against the literature review; no folk-biology claims.
      (Verified 2026-10-12: LMNA counterfactual transcript counts corrected
      22 → 45 against GENCODE v46 in results/counterfactual.json.)
- [x] Isoform-usage categories reflect measurement reality (short-read
      ambiguity not overstated; long-read support honestly labeled).
      (Figure 03 panel C explicitly labels 77.4% of usage calls "uncertain";
      manuscript §9 states the usage layer needs long-read data for a fair test.)
- [x] Case studies' claimed mechanisms verified against primary sources.
      (TRNT1: Ly et al. 2025 isoform-selective alleles; DNM1: AVI 2026 preprint
      brain-specific exon 10a rescue; both cited as literature-grounded, TRNT1
      explicitly "not forced".)

## Data
- [x] `data/DATA_PROVENANCE.md` complete: every input has source, URL,
      version, download date, genome build, license, preprocessing.
      (Updated 2026-10-12: gnomAD marked DEFERRED/unreachable; CELLxGENE
      marked deferred; SpliceAI v1.3.1 + CADD v1.7 completion recorded with
      dates. No [PENDING] rows remain.)
- [x] Genome build consistent (GRCh38) across all inputs; SpliceAI's
      bundled GENCODE V24-era annotation documented as the sole exception.
- [x] Transcript annotation consistent (GENCODE v46; MANE v1.5); no silent
      mixing of annotation versions between stages.

## Machine learning
- [x] `results/splits.json` regenerated from scratch and matches the used
      splits (determinism check). (Original audit 2026-10-06: exact match.)
- [x] No gene appears in more than one split; no test variant key in train.
      (Gene-level splits: 61/17/24 genes; verified in leakage_prevention.md.)
- [x] Test set accessed once, logged in `results/test_access.log`.
- [x] All tier comparisons use identical variants, splits, learner, and
      hyperparameters — only feature tiers differ. (LogisticRegression
      C=1.0, balanced, seed 20261006; re-fit verified 2026-10-12 in
      scripts/17_statistical_strengthening.py reproduces benchmark.json
      point estimates to 4 decimals.)
- [x] Baselines scored on compatible variant classes only; incompatible
      classes flagged (AlphaMissense missense-only n=1,441 with ClinVar
      circularity caveat; CADD n=11,240 with indel-unscorable note;
      SpliceAI n=13,687 with 623 unscorable disclosed).

## Statistics
- [x] Metrics appropriate for imbalanced labels (AUPRC + top-k alongside
      AUROC, not AUROC alone).
- [x] Confidence intervals computed (bootstrap, seeded) and reported;
      model comparisons use paired tests, not eyeballing.
      (2026-10-12 strengthening: per-model AUPRC 95% CIs added to Table 5.1;
      paired AUPRC-difference CIs + permutation p-values in §5.2 —
      M4−M1 +0.0509 [+0.0473,+0.0546] p<0.001; M2−M1 +0.0310 p<0.001;
      M3−M1 +0.0317 p<0.001; M5−M4 n.s. All in results/statistics.json.)
- [x] Subgroup analyses labeled exploratory with CIs; no subgroup
      cherry-picking (isoform-selective subgroup pre-registered;
      CADD-weak subgroup threshold defined on VALIDATION, applied to TEST —
      results/statistics.json §cadd_complementarity).

## Reproducibility
- [x] `bash scripts/run_all.sh` exists; pipeline stages are numbered
      scripts/01–17 with fixed seeds.
- [x] Random seeds fixed and recorded (20261006 throughout); package
      versions pinned in `requirements.txt`.
- [x] `pytest tests/` passes: 37/37 (verified 2026-10-12).

## Novelty
- [x] Literature review covers 2024–2026 (24 methods); novelty statement
      conservative and defensible (closest competitor IMPPROVE is gene-level;
      no equivalent phenotype+cell+isoform tiered design found).
      Two citations flagged for submission-time verification (Weykopf DOI,
      V2P author list) — noted in manuscript References.
- [x] No claim of "nobody has ever done this" anywhere in paper or docs.
      (grep verified 2026-10-12.)

## Scientific integrity
- [x] Every number in `paper/manuscript.md` traces to a file in `results/`.
      (2026-10-12: 12-number end-to-end spot-check — all consistent, including
      the new AUPRC CIs, CADD-weak AUROC 0.861, and threshold 28.80.)
- [x] All `[RESULT PENDING]` markers resolved (grep: 0 remaining).
- [x] Null results reported as null results (M4−M1 AUROC n.s.; M5−M4 n.s.;
      M2/M3 AUROC n.s. with p-values; gnomAD deferred).

## PhD value
- [ ] The finished project demonstrates independent research ability:
      literature-grounded question, controlled experiment, honest statistics,
      documented limitations — not merely software engineering.

**Audit outcome:** PASS with noted limitations (see below)
**Auditor:** Muse (automated checks) + human review required before submission
**Date:** 2026-10-06

## Findings

- Biology: assumptions checked against literature review; isoform-usage
  uncertainty honestly labeled (5,465/7,057 calls "uncertain"); TRNT1/DNM1
  case studies evidence-driven.
- Data: GRCh38 consistent; GENCODE v46 + MANE v1.5 pinned; provenance table
  complete except pending CADD/SpliceAI rows.
- ML: splits re-derived deterministically (MATCH); test accessed once
  (logged); M4-vs-M1 uses identical data/splits/learner.
- Statistics: bootstrap CIs + paired tests reported; subgroup
  pre-registered; nulls reported as nulls.
- Reproducibility: 37/37 tests pass; run_all.sh documents stage order;
  seeds pinned; CADD enrichment still running (deferred baseline).
- Novelty: literature re-run required at submission (field moving fast);
  no "nobody has done this" claims; preprints cited as preprints.
- Integrity: 5/5 spot-checked numbers trace to results files; no
  [RESULT PENDING] remains in results sections; nulls reported.
- PhD value: independent research design demonstrated — literature-grounded
  question, controlled ablation, honest statistics, documented limitations.

## Open items (not blockers for the repo, required before any submission)

1. CADD baseline comparison pending 09b completion (~7h API fetch).
2. SpliceAI local run not yet executed (needs reference FASTA + TF env).
3. Single-cell isoform resolution for the cell-type layer (the identified
   next experiment).
4. Human re-read of manuscript prose.
