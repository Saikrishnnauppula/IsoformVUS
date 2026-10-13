#!/usr/bin/env bash
# IsoformVUS end-to-end pipeline. Deterministic stage order; each stage is
# idempotent and records provenance. Exit on first failure.
set -euo pipefail
REPO="$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)"
cd "$REPO"

echo "=== [01] download sources ==="
python3 scripts/01_download_sources.py --source all --out data/raw

echo "=== [02] build gene panel ==="
python3 scripts/02_build_gene_panel.py --config config/base.yaml

echo "=== [03] ingest ClinVar + gnomAD labels ==="
python3 scripts/03_build_variant_labels.py --config config/base.yaml

echo "=== [04] phenotype -> cell-type relevance ==="
python3 scripts/04_phenotype_cell_mapping.py --config config/base.yaml

echo "=== [05] cell-type -> isoform usage ==="
python3 scripts/05_isoform_usage.py --config config/base.yaml

echo "=== [06] isoform-specific variant effects ==="
python3 scripts/06_variant_isoform_effects.py --config config/base.yaml

echo "=== [07] feature matrix ==="
python3 scripts/07_build_features.py --config config/base.yaml

echo "=== [08] leakage-safe splits ==="
python3 scripts/08_make_splits.py --config config/base.yaml

echo "=== [09] baselines (CADD / SpliceAI / AlphaMissense / consequence+AF) ==="
python3 scripts/09_run_baselines.py --config config/base.yaml

echo "=== [10] tiered IsoformVUS models ==="
python3 scripts/10_train_tiered.py --config config/base.yaml

echo "=== [11] held-out benchmark ==="
python3 scripts/11_benchmark.py --config config/base.yaml

echo "=== [12] ablation ==="
python3 scripts/12_ablation.py --config config/base.yaml

echo "=== [13] counterfactual isoform experiment ==="
python3 scripts/13_counterfactual.py --config config/base.yaml

echo "=== [14] case studies ==="
python3 scripts/14_case_studies.py --config config/base.yaml

echo "=== [15] figures ==="
python3 scripts/15_make_figures.py --config config/base.yaml

echo "=== [16] statistical strengthening (AUPRC CIs, permutation p-values, CADD held-out subgroup) ==="
python3 scripts/17_statistical_strengthening.py

echo "=== [17] tests ==="
PYTHONPATH=src python3 -m pytest tests/ -q

echo "ALL STAGES COMPLETE"
