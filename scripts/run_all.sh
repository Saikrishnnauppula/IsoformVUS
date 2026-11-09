#!/usr/bin/env bash
# IsoformVUS end-to-end pipeline. Deterministic stage order; each stage is
# idempotent and records provenance. Exit on first failure.
# Every stage calls a real script; no stage is a stub.
#
# Stage count: 19 numbered stages (01-19) in 23 script invocations.
# Stage 04 runs two scripts (04a, 04b), stage 09 runs five (09b-09f),
# stages 10-11 share one script (scripts/10_11_benchmark.py), and
# stage 19 is the pytest suite.
#
# Network and data requirements (stages that need more than the repo):
# - Stages 01-03 need data/raw (~5.1 GB). The ClinVar weekly the locked
#   results were built from (clinvar_20261004.vcf.gz) now returns HTTP 404
#   from NCBI FTP (see docs/data.md). Stage 01 pins the live weekly
#   clinvar_20261108.vcf.gz (verified 2026-11-09) with a fallback to the
#   unversioned latest file, so a from-scratch run downloads a different
#   ClinVar version and will not reproduce the locked results exactly;
#   data/raw/ and data/interim/ are not in the repository.
# - Stage 04a and stage 05 need live GTEx Portal v2 API calls.
# - Stage 09b needs the live CADD REST API: ~12 h for 60k sequential
#   calls, so the default scope is validation+test variants only.
# - Stages 09d-09f need a SpliceAI venv plus the GRCh38 FASTA and GTF.
# - Stage 03b (gnomAD AF enrichment) is NOT wired in: the gnomAD REST API
#   was unreachable when the pipeline was built, so no AF feature exists
#   in the shipped feature matrix. Run scripts/03b_gnomad_af.py by hand
#   once the API is reachable, then rebuild stage 07.
# - Stage 14 runs the SpliceAI+CADD baseline enrichment (scripts/14_add_baselines.py).
#   The old case-studies stage had no implementing script and was removed;
#   results/case_study_candidates.md exists but nothing generates it.
set -euo pipefail
REPO="$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)"
cd "$REPO"

echo "=== [01] download sources ==="
python3 scripts/01_download_sources.py --source all --out data/raw

echo "=== [02] build gene panel ==="
python3 scripts/02_build_gene_panel.py --config config/base.yaml

echo "=== [03] ingest ClinVar + gnomAD labels ==="
python3 scripts/03_build_variant_labels.py --config config/base.yaml

echo "=== [04] GTEx expression + phenotype -> cell-type relevance ==="
python3 scripts/04a_gtex_expression.py
python3 scripts/04b_cell_relevance.py

echo "=== [05] cell-type -> isoform usage ==="
python3 scripts/05_isoform_usage.py

echo "=== [06] isoform-specific variant effects ==="
python3 scripts/06_variant_isoform_effects.py --config config/base.yaml

echo "=== [07] feature matrix ==="
python3 scripts/07_build_features.py

echo "=== [08] leakage-safe splits ==="
python3 scripts/08_make_splits.py --config config/base.yaml

echo "=== [09] baselines (CADD / SpliceAI / AlphaMissense) ==="
python3 scripts/09b_cadd_api.py
python3 scripts/09c_alphamissense.py
python3 scripts/09d_spliceai_vcf.py
python3 scripts/09e_spliceai_run.py
python3 scripts/09f_spliceai_merge.py

echo "=== [10-11] tiered training + held-out benchmark ==="
# One script covers both stages (scripts/10_11_benchmark.py). It refuses
# to run twice without --allow-rerun; every access is appended to
# results/test_access.log, which records 13 logged accesses (17 reads of
# the test split in all, counting 4 unlogged script-18 evaluations).
python3 scripts/10_11_benchmark.py --allow-rerun

echo "=== [12] ablation (paired analysis) ==="
python3 scripts/12b_paired_analysis.py

echo "=== [13] counterfactual isoform experiment ==="
python3 scripts/13_counterfactual.py

echo "=== [14] SpliceAI + CADD baselines on the locked test set ==="
python3 scripts/14_add_baselines.py

echo "=== [15] figures ==="
python3 scripts/15_make_figures.py

echo "=== [16] remaining manuscript figures ==="
python3 scripts/16_make_missing_figures.py

echo "=== [17] statistical strengthening (AUPRC CIs, permutation p-values, CADD held-out subgroup) ==="
python3 scripts/17_statistical_strengthening.py

echo "=== [18] exon architecture, backup skipping & NMD boundary analysis ==="
python3 scripts/18_exon_architecture_nmd.py

echo "=== [19] tests ==="
PYTHONPATH=src python3 -m pytest tests/ -q

echo "ALL STAGES COMPLETE"
