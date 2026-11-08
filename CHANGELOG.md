# Changelog

All notable changes to IsoformVUS. Versions follow semantic versioning.
Current version: 0.1.0.

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

## 2026-11-08 - Sync from Drive folder and panel-review fixes

- Synced the repository to the Drive-folder project state. This added `scripts/18_exon_architecture_nmd.py` (4,593 reconstructed MANE Select exons; 4x2x2 exon length x backup skipping x phase matrix; donor/acceptor boundary analysis; mediation; Tavtigian likelihood ratios), `results/exon_architecture_nmd.json`, `results/exon_skipping_nmd_candidates.tsv`, `results/vus_bayesian_stratification.tsv`, and figures 09–12. The version stays 0.1.0; the version heading from the synced changelog was removed.
- A journal review panel rejected the manuscript built from this work (REJECT, confidence 4/5). In response: the M6 model, VUS likelihood ratios, ASO candidate list, donor-edge finding, and "co-transcriptional NMD" mechanism were withdrawn; the test set was disclosed as accessed 13 times as logged (17 reads in all); the reference list was rebuilt (6 wrong, 3 partly wrong, 2 missing in the reviewed version; the panel's proposed SpliceAI and Pozo fixes were both wrong on re-check); the stub pipeline stages were deleted or implemented and the runner rewired; and the docs were rewritten (`docs/methods.md`, `docs/data.md`, `docs/limitations.md`). The review report itself was not copied into the repository.
- Landed 2026-11-08: feature matrix rebuilt with real CADD values (19,014 of 60,521 scored) and benchmark re-run from the shipped matrix. The old M1-M4 gain (+0.0509) is gone (M1 AUPRC 0.9465, M4 0.9395, M2 best at 0.9533); `statistics.json` tie flags corrected and tie-corrected values added; Fig 05 reads values from `statistics.json`; `results/splits.json` pointer fixed; M6 donor-edge offsets fit on train genes only; CMH odds ratios shipped with per-gene tables; sensitivity check 4 saved; script 18 writes its outputs and logs test accesses; `breakthrough/` deleted and its five script-18 inputs moved to `data/interim/script18/`; Makefile targets (`test`, `benchmark`, `figures`, `pipeline`, `lint`, `format`) and CI (pytest + ruff/black) added; the manuscript was rewritten around the rebuilt numbers and the reference list replaced with 11 verified entries.
- New checksums (2026-11-08): `data/processed/feature_matrix.tsv.gz` a67f15d71eba05e31a84ef60d36d8f14177cf04959aa8033360d8a09b847b6a6; `results/benchmark.json` a5b8c7d3d8d19063400ae754f4f0f5949a4ef833075f0beaa8b61739d7e0cad3; `results/statistics.json` 0272ee0f2203ca5894af89fd71adfee44f202f47f3bab378ed13337e20ea1393; `results/exon_architecture_nmd.json` 532f0e877ee2f2fc083cd0521e8430dd68b6189674d7c2052b35dce5aecfca5d.
