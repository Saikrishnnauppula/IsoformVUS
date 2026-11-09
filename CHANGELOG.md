# Changelog

All notable changes to IsoformVUS. Versions follow semantic versioning.
Current version: 0.1.0.

## [Unreleased] - 2026-11-09 (review fix round)

- Gene-cluster bootstrap for the CMH odds ratios (script 18, 2,000
  resamples of genes, seed 0), shipped in
  results/exon_architecture_nmd.json: length OR 0.2005, 95% CI
  [0.1119, 0.3384]; backup-skip OR 0.3461, 95% CI [0.1280, 1.0299]
  (includes 1 — backup skipping is now reported as descriptive). The
  shipped CMH z statistics are withheld at source, as the p-values
  already were.
- Bounding analysis over unresolved variants shipped in the same JSON
  (vus_bounding_analysis): cell 1B true rate in [11.07%, 95.91%], cell
  4D in [0.07%, 94.09%]; the intervals overlap completely.
- Annotation exon atlas (scripts/20_build_exon_atlas.py,
  data/processed/mane_exon_atlas.tsv.gz) from GENCODE v46 + MANE v1.5,
  and validation + annotated re-derivation of the 4x2x2 tensor
  (scripts/21_validate_exon_reconstruction.py,
  results/exon_reconstruction_validation.json, docs/exon_atlas_validation.md):
  strand agreement 102/102; 14.65% of matched reconstructed exons differ
  in length from annotation by >10%; 625 reconstructed exons have no
  annotated counterpart; only 3,761 strict-internal CDS exons exist
  against 4,509 reconstructed. Annotated tensor is now primary in the
  manuscript (1B 78.11%, 4D 1.38%; annotated CMH length OR 0.1860
  [0.0946, 0.3247], backup OR 0.2720 [0.1090, 0.9218]).
- Usage-threshold sweep shipped in the stage-18 JSON: backup CMH OR
  0.6063 / 0.3461 / 0.3650 at thresholds 0.05 / 0.10 / 0.20.
- Composite-CLNSIG sensitivity in results/statistics.json
  (clnsig_composite_sensitivity): excluding the 17 test rows with
  drug_response / low_penetrance composite significances changes M1–M5
  test metrics by at most 0.000862.
- Carried-over mediation and distance-profile values moved out of
  results/exon_architecture_nmd.json into
  results/withdrawn_carried_over_values.json (status: withdrawn
  2026-11-09; no generating code; do not cite). Script 18 no longer
  emits them; figures 10/11 read them from the quarantine file.
- Candidate TSV: ASO_Rescuability_Score renamed to Exon_Nomination_Score;
  Tier 1/2 labels now say "computational nomination only". The
  bayesian_lr_table is recomputed on the corrected exon-level frame and
  its exon counts sum to the tensor's 4,509.
- Figure 02's banner is generated from results/benchmark.json +
  results/statistics.json; Figure 08 panel C reads its counts from
  results/counterfactual.json; Figure 01's footnote carries the refit
  paired AUROC diffs. statistics.json's step-metric block is renamed to
  step_metric_tie_artifact_demo with a not-evidence note.
- run_all.sh refuses to run stages 10-11 under an unpinned interpreter
  (exit 3 with instructions). Test-access logging in scripts 10_11 and
  14 now appends on completion, not at start. Degenerate top-10 recall
  CIs in benchmark.json are null with an explanatory note key.
- Manuscript retitled ("…— A Benchmark Audit"; the "and NMD" phrase is
  gone); Section 4.1 reports the annotated tensor, the clustered CIs,
  and the bounding analysis; Section 4.2 rows recoded C1–C8; McSherry &
  Najork DOI corrected to 10.1007/978-3-540-78646-7_38; Buljan 2012 and
  Weykopf 2026 now cited in the text. New data-grounded tests pin the
  tensor sums, CMH values, clustered CIs, bounding block, LR table, and
  the annotated-frame headline cells.

## [Unreleased] - 2026-11-09 (professor-audit fixes)

- Refit regenerated on the pinned `requirements.txt` stack
  (scikit-learn 1.5.1 / pandas 2.2.2 / numpy 1.26.4). The earlier
  2026-11-09 refit ran on an unpinned system stack; its M3 test values
  drifted (AUROC 0.7975 -> 0.7962, step AUPRC 0.8335 -> 0.8328) and its
  M3 validation AUROC by 0.0054. Live values now come from
  scripts 10_11 + 14 + 19 run on the pinned stack, in that order.
- `scripts/19_refit_statistics.py`: the `auprc_diffs` permutation p-value
  now uses the same tie-aware AUPRC as the reported CI (it used step
  AUPRC before), and `significant` requires perm p < 0.01 (Bonferroni)
  AND the bootstrap CI to exclude zero. M2-M1 and M3-M1 variant-level
  flags are now `false` (they were `true` while their CIs included
  zero); M4-M1 and M5-M1 stay `true`. Script 19 also carries its
  labelled superseded sections forward on reruns instead of
  relabelling live values as superseded.
- Stage-18 inputs ship: the five tables in `data/interim/script18/`
  (the 587k ClinVar table gzipped to 7.9 MB) are committed under a
  `.gitignore` negation. Stage 18 runs from a fresh clone and reproduced
  its committed TSVs byte-identically (JSON values within 2e-6).
- Figures 09-12 are replotted by `scripts/16_make_missing_figures.py`
  from the shipped stage-18 outputs; they are new renderings, not the
  original images, and stale values in the originals (donor-edge LR
  4.52; M6 internal-missense AUROC 95.1%) are corrected to the shipped
  values (1.817; 0.9321).
- `scripts/run_all.sh` has two modes: default verify mode (stages
  runnable from shipped data, runs to completion from a fresh clone)
  and `--full` (original 19-stage rebuild behind a fail-fast preflight).
  Bit-exact full reproduction remains impossible: the locked ClinVar
  weekly returns HTTP 404 and `data/raw/` is not shipped.
- Test-set access count updated to the log as committed: 16 logged
  accesses, 20 reads in all. The 2026-11-09 entries are a script-18
  verification rerun, the pinned-stack benchmark refit rerun, and a
  baseline-enrichment run.
- `docs/limitations.md` and `docs/methods.md`: the dead rebuilt run's
  step values now read M1 0.9465 / M4 0.9395 everywhere (two places said
  0.9463 / 0.9394, which are that run's tie-corrected values).
- The access-log test now pins the committed 16-entry prefix and
  permits later appended entries: the log is append-only, and a
  verify-mode run appends 4 real entries before its own pytest stage,
  which the exact-count version of the test rejected.

## [0.1.0] - 2026-10-13

First public release.

### Results (locked held-out benchmark, 24 unseen genes / 14,400 test variants)

- Isoform-structure features (M4 vs M1): test AUPRC 0.802 → 0.852 (+0.0509, 95% CI [+0.0473, +0.0546], permutation p < 0.001).
- Phenotype (M2) and cell-type (M3) layers: no AUROC gains (permutation p = 0.30, 0.32). The +0.031 step-AUPRC gains listed here were later shown to be a tie artifact (see 2026-11-08 entry); tie-corrected gains are +0.0006 and -0.0002.
- Tissue-level isoform usage (M5 vs M4): no gain on either metric (AUPRC Δ ≈ -0.0001, permutation p = 0.017, n.s. under Bonferroni α = 0.01).
- CADD v1.7 baseline (test AUROC 0.9554) outperforms the 26-feature model outright. Inside the held-out CADD-weak subgroup (n = 8,601), M5 reaches AUROC 0.861 [0.848, 0.874] with 4.36x top-decile pathogenic enrichment.
- SpliceAI v1.3.1 baseline (test AUROC 0.6961) is beaten by M5 (+0.106).
- AlphaMissense (missense-only, AUROC 0.9435; ClinVar-training caveat disclosed).

### Reproducibility notes

- Full raw-data pipeline is **not** independently reproducible from this release: `data/raw/` (~5.1 GB) and `data/interim/` (~3.1 GB) are excluded; see `data/DATA_PROVENANCE.md` for exact source versions and re-download instructions.
- Reviewers can verify headline results without raw data via `results/benchmark.json`, `results/statistics.json`, `results/per_model_cis.json`, and `data/processed/feature_matrix.tsv.gz`.
- Pinned dependencies tested in a clean Python 3.11 environment (`environment.yml`); 37/37 tests passing.

### Data checksums (v0.1.0, sha256)

- `data/processed/feature_matrix.tsv.gz`: 42926422601cd7c808e1e4d8551c469355b82bfda6779eae50f788cb65210bad
- `data/processed/variant_cadd.tsv.gz`: 12c99ed606711a02241787298138187f90a9796a461bee22115280ecba9e750b
- `data/processed/variant_spliceai.tsv.gz`: f9dbaacf1af62c870d8a2f2f2e4395b9c5cbf912bed26869db3b997eda78e0b4
- `results/benchmark.json`: 14eabd7b7bae16bc54a8ca5a52e5a6194b029bfbfae67bff3ca8866fe49f956e
- `results/statistics.json`: 8d3c6a085b106631ac602cd9fce8955215f70bc21ff1700a8890713b96cc591b

The feature-matrix checksum above is the pre-CADD build (all `cadd_phred` values -1.0); the rebuild landed on 2026-11-08 (see the entry below).

## 2026-11-09 - No-CADD refit and statistics rebuild

- `cadd_phred` and `cadd_missing` removed from every fitted tier in `data/processed/feature_tiers.json`: CADD was scored for 0 of 35,921 train rows, so the fitted coefficient was fit on a constant and test-time CADD signal inflated the 2026-11-08 rebuilt tiers (M1 step AUPRC 0.9465). Benchmark refit (`scripts/10_11_benchmark.py` + `scripts/14_add_baselines.py`): step AUPRC M1 0.8015, M2 0.8328, M3 0.8335, M4 0.8524, M5 0.8523; tie-corrected M1 0.8336, M4 0.8526. CADD is a pure external baseline now (AUROC 0.9554 / AUPRC 0.9590 on its 11,240-variant subset; it beats every tier on AUPRC). The rebuilt run's values are superseded, not deleted: they sit under labelled `*_superseded_*` keys in `results/statistics.json`.
- `results/statistics.json` rebuilt by the new `scripts/19_refit_statistics.py`: top-level `auprc_diffs` are the refit tie-aware diffs (M4-M1 +0.0189 [+0.0154, +0.0228]); the dead pre-rebuild +0.0509 survives only under `auprc_diffs_superseded_pre_rebuild`, labelled SUPERSEDED/DEAD. Gene-clustered CIs added (cluster bootstrap over the 24 test genes): M2-M1 +0.0007 [-0.0029, +0.0045] — "M2 best tier" withdrawn; M4-M1 +0.0189 [+0.0107, +0.0278]. `cadd_complementarity` recomputed on the refit (T = 28.8, weak n = 8,601, M5 AUROC 0.8610). `results/per_model_cis.json` is now generated by script 19.
- `docs/methods.md` §209 claim removed: the "validation-trained CADD-plus-context stack beats CADD alone by +0.0191 AUPRC (isoform slice +0.0015)" numbers appear in no results file and no shipped script fits a validation-trained stack; the closest computable train-fit quantity is M5 minus CADD AUPRC -0.2367 on the CADD-scored test subset (audit trail in `results/statistics.json`, key `methods_209_check`).
- Script 18: CMH p-values (4.29e-187, 3.53e-42) withheld at source — they treat variants as independent within genes; ORs and informative-gene counts remain as descriptive values. The 19-exon tensor mismatch is resolved: per-cell exon counts are now taken from the exon-level frame (backup skipping resolved per exon as the maximum over its variants) and sum to 4,509, agreeing with the exon atlas; the variant-level counts of those 19 exons remain split between cells (reassigning them needs the script-18 inputs). `results/exon_architecture_nmd.json` carries a `tensor_correction_note`.

## 2026-11-09 - Manuscript corrections

- Abstract and Sections 1, 2.1, 3.1, 5.1, 5.3, 6.1 (new item 15), and 6.2 rewritten around the refit: the rebuilt run's "M2 best tier" and "M1-M4 gain gone" claims are withdrawn as a CADD artifact; the M4 gain is reported tie-corrected (+0.0189, gene-clustered CI excluding zero). Figure callouts for Figures 01–08 added; the LMNA counterfactual and the TRNT1/DNM1 literature case studies are reported in Section 5.3 with their caveats (the counterfactual is built on M5, whose usage features add nothing over M4).
- Section 3.2 (three isoform regimes) cut from the abstract and body, and the AlphaMissense mediation percentages (77.9% / 94.7%) retracted in Section 4.2: neither has a shipped source (inputs in `data/interim/script18/`, builder code deleted). Both withdrawals are recorded in Section 6.1 (items 8 and 10). CMH p-values deleted from Section 4.1. Indel counts corrected to 1,155 of 1,182; tie penalty corrected to -0.0318; the M6 comparison now uses tie-corrected AUPRC for both models (0.9928 against 0.9856).
- Figures 09–12 marked as unregenerable archival figures in `figures/README.md` and Section 6.2; Figure 11 is retained only as a record of the withdrawn donor-edge finding. The AI-use statement is completed (AI tools used for language editing and code review; the author takes responsibility for all results). References for Chapman & Lassmann 2025, Buljan et al. 2012, Ly et al. 2025, Weykopf et al. 2026, and Cheng et al. 2026 (AlphaGenome Atlas) verified and added to Section 7.

## 2026-11-09 - Professor-check fixes (code, README, infrastructure)

- `results/test_access.log` is now tracked: the original timestamped log (13 entries, 2026-10-06 to 2026-11-08) sat untracked in the working tree because `.gitignore` carried a bare `*.log` rule. `.gitignore` now exempts it (`!results/test_access.log`), and the stale "9 prior accesses" comment in `scripts/run_all.sh` reads 13 logged / 17 in all. README, `docs/methods.md`, and `docs/limitations.md` use the same 13-logged / 17-reads wording (the 4 unlogged script-18 reads are disclosed in the manuscript and docs); figure 02's hard-coded "accessed 13 times, logged" now reads "13 logged accesses · 17 reads".
- Stage 01 ClinVar pin moved from the dead `clinvar_20261004.vcf.gz` (HTTP 404) to the live weekly `clinvar_20261108.vcf.gz` (HEAD-verified 2026-11-09), with a fallback to the unversioned latest file; `config/data_sources.yaml` and the `scripts/03_build_variant_labels.py` default VCF path match. A rebuild from the new pin is a different ClinVar version; `docs/data.md` states what still blocks a full fresh-clone run (live GTEx/CADD/SpliceAI services, `data/raw/` and `data/interim/` absent, script-18 inputs absent).
- README fixes: `pyproject.toml` pins with `==` (the "unpinned" claim was wrong, also fixed in `docs/methods.md`); the runner is 19 numbered stages in 23 script invocations (stage 04 runs two scripts, stage 09 runs five, stages 10–11 share one); a "Verify without re-downloading raw data" section was added (the section `scripts/01_download_sources.py`'s error text already pointed at). `scripts/03b_gnomad_af.py` remains unwired and no gnomAD AF feature exists; the README now says so.
- Dead code removed: `src/isoformvus/models/tiered.py` (zero references; the whole `models/` package went with it) and the empty `src/isoformvus/explainability/` package. `config/base.yaml` no longer points at the nonexistent `docs/leakage_prevention.md` or `src/isoformvus/features/spec.py` (the file is `spec.md`), and its dead LightGBM block is gone.
- Dependencies pruned to what the code imports (numpy, pandas, scipy, scikit-learn, matplotlib, pytest): pysam, pyBigWig, biopython, shap, click, tqdm, joblib, seaborn, jupyter, nbformat, lightgbm, and pyyaml removed from `requirements.txt` and `pyproject.toml` after grep showed zero imports of each.
- New data-grounded tests in `tests/test_shipped_evidence.py`: feature-matrix shape/NaNs, CADD scored counts by split recomputed from the shipped files (train 0/35,921, validation 7,774/10,200, test 11,240/14,400), zero gene overlap, benchmark/statistics cross-checks, a subsample refit check, and a guard that `statistics.json` carries no live significant +0.0509 (that guard fails until the statistics fix lands).

## 2026-11-08 - Sync from Drive folder and review fixes

- Synced the repository to the Drive-folder project state. This added `scripts/18_exon_architecture_nmd.py` (4,593 reconstructed MANE Select exons; 4x2x2 exon length x backup skipping x phase matrix; donor/acceptor boundary analysis; mediation; Tavtigian likelihood ratios), `results/exon_architecture_nmd.json`, `results/exon_skipping_nmd_candidates.tsv`, `results/vus_bayesian_stratification.tsv`, and figures 09–12. The version stays 0.1.0; the version heading from the synced changelog was removed.
- An internal pre-submission review recommended against the manuscript in its original form. In response: the M6 model, VUS likelihood ratios, ASO candidate list, donor-edge finding, and "co-transcriptional NMD" mechanism were withdrawn; the test set was disclosed as accessed 13 times as logged (17 reads in all); the reference list was rebuilt (6 wrong, 3 partly wrong, 2 missing in the reviewed version; the review's proposed SpliceAI and Pozo fixes were both wrong on re-check); the stub pipeline stages were deleted or implemented and the runner rewired; and the docs were rewritten (`docs/methods.md`, `docs/data.md`, `docs/limitations.md`). The review report itself was not copied into the repository.
- Landed 2026-11-08: feature matrix rebuilt with real CADD values (19,014 of 60,521 scored) and benchmark re-run from the shipped matrix. The old M1-M4 gain (+0.0509) is gone (M1 AUPRC 0.9465, M4 0.9395, M2 best at 0.9533); `statistics.json` tie flags corrected and tie-corrected values added; Fig 05 reads values from `statistics.json`; `results/splits.json` pointer fixed; M6 donor-edge offsets fit on train genes only; CMH odds ratios shipped with per-gene tables; sensitivity check 4 saved; script 18 writes its outputs and logs test accesses; `breakthrough/` deleted and its five script-18 inputs moved to `data/interim/script18/`; Makefile targets (`test`, `benchmark`, `figures`, `pipeline`, `lint`, `format`) and CI (pytest + ruff/black) added; the manuscript was rewritten around the rebuilt numbers and the reference list replaced with 11 verified entries.
- New checksums (2026-11-08): `data/processed/feature_matrix.tsv.gz` a67f15d71eba05e31a84ef60d36d8f14177cf04959aa8033360d8a09b847b6a6; `results/benchmark.json` a5b8c7d3d8d19063400ae754f4f0f5949a4ef833075f0beaa8b61739d7e0cad3; `results/statistics.json` 0272ee0f2203ca5894af89fd71adfee44f202f47f3bab378ed13337e20ea1393; `results/exon_architecture_nmd.json` 532f0e877ee2f2fc083cd0521e8430dd68b6189674d7c2052b35dce5aecfca5d.
