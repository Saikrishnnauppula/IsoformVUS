# Changelog

All notable changes to IsoformVUS. Versions follow semantic versioning.

## [0.1.0] — 2026-10-13

First public release.

### Results (locked held-out benchmark, 24 unseen genes / 14,400 test variants)

- Isoform-structure features (M4 vs M1): test AUPRC 0.802 → 0.852 (+0.0509, 95% CI [+0.0473, +0.0546], permutation p < 0.001).
- Phenotype (M2) and cell-type (M3) layers: no statistically significant AUROC gains (permutation p = 0.30, 0.32) but statistically significant AUPRC gains (+0.031 each, p < 0.001).
- Tissue-level isoform usage (M5 vs M4): no gain on either metric (AUPRC Δ ≈ −0.0001, permutation p = 0.017, n.s. under Bonferroni α = 0.01).
- CADD v1.7 baseline (test AUROC 0.9554) outperforms the 26-feature model outright — IsoformVUS is a *complement*, not a competitor. Inside the held-out CADD-weak subgroup (n = 8,601), M5 reaches AUROC 0.861 [0.848, 0.874] with 4.36x top-decile pathogenic enrichment.
- SpliceAI v1.3.1 baseline (test AUROC 0.6961) is beaten by M5 (+0.106).
- AlphaMissense (missense-only, AUROC 0.9435; ClinVar-training caveat disclosed).

### Pre-publication fixes

- CITATION.cff finalized with author and repository URL.
- Manuscript abstract refined: AUROC/AUPRC split stated precisely for M2/M3.
- README dead links removed; figure captions regenerated to match the refined finding.
- `.gitignore` corrected to include `results/`, `figures/`, processed matrices.
- Significance rule made consistent: Bonferroni α = 0.01 applied in the statistics script itself.
- README stage count corrected to 17.

### Reproducibility notes

- Full raw-data pipeline is **not** independently reproducible from this release: `data/raw/` (~5.1 GB) and `data/interim/` (~3.1 GB) are excluded; see `data/DATA_PROVENANCE.md` for exact source versions and re-download instructions.
- Reviewers can verify all headline results without raw data via `results/benchmark.json`, `results/statistics.json`, `results/per_model_cis.json`, and `data/processed/feature_matrix.tsv.gz`.
- Pinned dependencies tested in a clean Python 3.11 environment (`environment.yml`); 37/37 tests passing.

### Data checksums (v0.1.0, sha256)

- `data/processed/feature_matrix.tsv.gz`: 42926422601cd7c808e1e4d8551c469355b82bfda6779eae50f788cb65210bad
- `data/processed/variant_cadd.tsv.gz`: 12c99ed606711a02241787298138187f90a9796a461bee22115280ecba9e750b
- `data/processed/variant_spliceai.tsv.gz`: f9dbaacf1af62c870d8a2f2f2e4395b9c5cbf912bed26869db3b997eda78e0b4
- `results/benchmark.json`: 14eabd7b7bae16bc54a8ca5a52e5a6194b029bfbfae67bff3ca8866fe49f956e
- `results/statistics.json`: 8d3c6a085b106631ac602cd9fce8955215f70bc21ff1700a8890713b96cc591b
