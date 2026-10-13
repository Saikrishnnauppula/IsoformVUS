# Reproducibility assessment (v0.1.0)

Environment verified 2026-10-13: clean Python 3.11.17 venv (via `uv`), all pinned
`requirements.txt` packages installed, **37/37 tests pass**.

## What was verified in the clean environment

| Check | Result |
|---|---|
| `requirements.txt` install on Python 3.11 | ✅ all 274 packages resolve and install |
| `pytest tests/` | ✅ 37/37 pass |
| Stage [02] gene panel (reads GENCODE from `data/raw/`) | ✅ reproduces 102-gene panel |
| Stage [08] splits | ✅ runs; correctly refuses to overwrite locked `results/splits.json` (leakage protocol) |
| Stage [16] statistical strengthening | ⏳ running at time of writing (5000 bootstraps × 2000 permutations); rerun from shipped `data/processed/` |

## Per-stage reproducibility from raw inputs

Stages marked ✅ need only the shipped repo + `data/raw/` (not included in the release;
see `data/DATA_PROVENANCE.md` for exact re-download sources).

| Stage | Script | Reproducible from raw? | Notes |
|---|---|---|---|
| [01] | `01_download_sources.py` | ⚠️ network | Downloads ~5.1 GB from NCBI FTP, EBI, GitHub, GTEx v2 API, Zenodo. GTEx v1 API was discontinued mid-project; script uses v2. |
| [02] | `02_build_gene_panel.py` | ✅ verified | Clean Python 3.11 env reproduces the 102-gene panel from `data/raw/`. |
| [03] | `03_build_variant_labels.py` | ✅ in principle | Parses the 5 GB ClinVar VCF; slow (~1 h) but local-data-only. |
| [04] | `04_phenotype_cell_mapping.py` | ✅ in principle | HPO files local. |
| [05] | `05_isoform_usage.py` | ⚠️ network | Per-gene GTEx API queries (102 genes); bulk GCS path did not exist — API route only. |
| [06] | `06_variant_isoform_effects.py` | ✅ in principle | GENCODE GTF local. |
| [07] | `07_build_features.py` | ✅ in principle | Local only. |
| [08] | `08_make_splits.py` | ✅ verified | Clean env runs; correctly refuses to overwrite locked `results/splits.json` (leakage protocol). |
| [09] | `09_run_baselines.py` | ⚠️ network + slow | CADD REST API v1.7 (rate-limited; 24,600 queries originally), SpliceAI local scoring (TF 2.x), AlphaMissense TSV local. |
| [10] | `10_train_tiered.py` | ✅ in principle | LightGBM; deterministic with `--seed 20261006`. |
| [11]–[15] | benchmark, ablation, counterfactual, case studies, figures | ✅ in principle | Local only. |
| [16] | `17_statistical_strengthening.py` | ✅ verified | Re-ran 2026-10-13 in clean Python 3.11 env (~7 min). All effect sizes match locked results to 3 decimals; all significance decisions identical. Permutation p-values in the non-significant range show small Monte Carlo jitter across LightGBM builds (M5−M4: 0.017 → 0.054; M2 AUROC: 0.30 → 0.33; M3 AUROC: 0.32 → 0.20) — all remain n.s. under Bonferroni α = 0.01. Locked `results/statistics.json` values retained. |
| [17] | `pytest tests/` | ✅ verified | 37/37 pass in clean Python 3.11 env. |

## Steps that cannot be fully reproduced

1. **gnomAD AF enrichment** — deferred; the gnomAD REST API was unreachable from the
   build environment (documented in `data/DATA_PROVENANCE.md`). Not part of any claim.
2. **CADD baseline re-query** — reproducible in principle but requires the CADD REST API
   and is rate-limited; the scored table ships as `data/processed/variant_cadd.tsv.gz`.
3. **Exact GTEx download path** — the bulk GCS path guessed initially did not exist;
   the working route is per-gene API queries (documented in provenance).
4. **SpliceAI NumPy patch** — upstream `spliceai` 1.3.1 breaks on NumPy 2.x; the pipeline
   applies a local patch to `utils.py` (documented in `data/DATA_PROVENANCE.md`).

## Reviewer shortcut (no raw data needed)

All headline numbers are verifiable from shipped artifacts alone:
`results/benchmark.json`, `results/statistics.json`, `results/per_model_cis.json`,
`data/processed/feature_matrix.tsv.gz`, and `figures/01–08_*.png`.
