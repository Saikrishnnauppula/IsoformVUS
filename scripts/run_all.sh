#!/usr/bin/env bash
# IsoformVUS pipeline runner. Two modes:
#
#   ./scripts/run_all.sh           VERIFY mode (default). Runs only the
#                                  stages that execute from data shipped in
#                                  this repository, in dependency order, and
#                                  finishes with the pytest suite. This mode
#                                  runs to completion from a fresh clone.
#   ./scripts/run_all.sh --full    FULL mode. The original 19-stage
#                                  from-scratch rebuild (stages 01-19 below),
#                                  preceded by a preflight that lists what
#                                  is missing and fails fast.
#
# Honest limits, stated up front:
# - Bit-exact reproduction of the locked results from raw inputs is NOT
#   possible. The ClinVar weekly the locked results were built from
#   (clinvar_20261004.vcf.gz) returns HTTP 404 from NCBI FTP and data/raw/
#   (~5.1 GB) is not in the repository. Stage 01 pins the live weekly
#   clinvar_20261108.vcf.gz, so a full rebuild uses a different ClinVar
#   version and will not bit-match the locked results.
# - Verify mode excludes two numbered stages, for stated reasons:
#   stage 13 (counterfactual) needs the GENCODE v46 GTF and MANE summary
#   in data/raw/, which are not shipped; stage 17 is superseded by
#   scripts/19_refit_statistics.py, which recomputes its three sections
#   on the no-CADD refit (running 17 after 19 would overwrite
#   results/statistics.json with the pre-refit schema).
# - Verify mode appends real entries to results/test_access.log: stages
#   10/11, 12, 14, and 18 each score the locked test set and log the
#   access. That is the log working as intended.
# - A stage-18 rerun reproduces its committed TSV outputs byte-identically
#   and its JSON values to within 2e-6. The committed
#   results/exon_architecture_nmd.json additionally carries hand-applied
#   note edits (CMH p-values removed) that the script does not re-emit.
# - Stage 03b (gnomAD AF enrichment) is NOT wired into either mode: no
#   gnomAD allele-frequency feature exists in the shipped feature matrix.
#
# Full-mode network and data requirements:
# - Stages 01-03 need data/raw (~5.1 GB) from the sources in docs/data.md.
# - Stage 04a and stage 05 need live GTEx Portal v2 API calls.
# - Stage 09b needs the live CADD REST API: ~12 h for 60k sequential
#   calls, so the default scope is validation+test variants only.
# - Stages 09d-09f need a SpliceAI venv plus the GRCh38 FASTA and GTF.
set -euo pipefail
REPO="$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)"
cd "$REPO"
PY="${PYTHON:-python3}"
MODE="${1:-verify}"

if [ "$MODE" = "verify" ]; then
    echo "=== VERIFY MODE: stages runnable from shipped data ==="
    echo "Excluded: stage 13 (needs data/raw GTF + MANE, not shipped);"
    echo "stage 17 (superseded by scripts/19_refit_statistics.py)."

    echo "=== [10-11] tiered training + held-out benchmark ==="
    $PY scripts/10_11_benchmark.py --allow-rerun

    echo "=== [12] ablation (paired analysis) ==="
    $PY scripts/12b_paired_analysis.py

    echo "=== [14] SpliceAI + CADD baselines on the locked test set ==="
    $PY scripts/14_add_baselines.py

    echo "=== [18] exon architecture, backup skipping & NMD boundary analysis ==="
    $PY scripts/18_exon_architecture_nmd.py

    echo "=== [19s] refit statistics (supersedes stage 17) ==="
    $PY scripts/19_refit_statistics.py

    echo "=== [15] figures ==="
    $PY scripts/15_make_figures.py

    echo "=== [16] remaining manuscript figures (incl. replotted 09-12) ==="
    $PY scripts/16_make_missing_figures.py

    echo "=== [tests] pytest suite ==="
    PYTHONPATH=src $PY -m pytest tests/ -q

    echo "VERIFY MODE COMPLETE"
    exit 0
fi

if [ "$MODE" != "--full" ]; then
    echo "usage: $0 [--full]" >&2
    exit 2
fi

echo "=== FULL MODE PREFLIGHT ==="
missing=0
if [ ! -d data/raw ] || [ -z "$(ls -A data/raw 2>/dev/null)" ]; then
    echo "MISSING: data/raw/ is absent or empty (~5.1 GB of downloads;"
    echo "         see docs/data.md for the source list)."
    missing=1
fi
if [ -z "$(ls data/raw/*.fa* 2>/dev/null)" ]; then
    echo "MISSING: no GRCh38 FASTA in data/raw/ (stages 09d-09f, SpliceAI)."
    missing=1
fi
echo "NOTE: clinvar_20261004.vcf.gz (the locked weekly) returns HTTP 404"
echo "      from NCBI FTP and cannot be restored. Stage 01 will download"
echo "      clinvar_20261108.vcf.gz instead; results will not bit-match"
echo "      the locked results."
if command -v curl >/dev/null 2>&1; then
    code=$(curl -s -o /dev/null -w "%{http_code}" -I \
        https://ftp.ncbi.nlm.nih.gov/pub/clinvar/vcf_GRCh38/clinvar_20261108.vcf.gz \
        || echo "000")
    echo "CHECK: clinvar_20261108.vcf.gz HTTP status: $code"
    if [ "$code" != "200" ]; then
        echo "MISSING: pinned ClinVar weekly is not reachable."
        missing=1
    fi
fi
echo "NOTE: stages 04a/05 need the live GTEx Portal v2 API and stage 09b"
echo "      the live CADD REST API (~12 h); stages 09d-09f need a SpliceAI"
echo "      install. These are not probed by this preflight."
if [ "$missing" -ne 0 ]; then
    echo "PREFLIGHT FAILED: the items above are missing. Failing fast"
    echo "instead of dying mid-run. For a run that works from the shipped"
    echo "repository alone, use verify mode: ./scripts/run_all.sh"
    exit 2
fi

echo "=== [01] download sources ==="
$PY scripts/01_download_sources.py --source all --out data/raw

echo "=== [02] build gene panel ==="
$PY scripts/02_build_gene_panel.py --config config/base.yaml

echo "=== [03] ingest ClinVar + gnomAD labels ==="
$PY scripts/03_build_variant_labels.py --config config/base.yaml

echo "=== [04] GTEx expression + phenotype -> cell-type relevance ==="
$PY scripts/04a_gtex_expression.py
$PY scripts/04b_cell_relevance.py

echo "=== [05] cell-type -> isoform usage ==="
$PY scripts/05_isoform_usage.py

echo "=== [06] isoform-specific variant effects ==="
$PY scripts/06_variant_isoform_effects.py --config config/base.yaml

echo "=== [07] feature matrix ==="
$PY scripts/07_build_features.py

echo "=== [08] leakage-safe splits ==="
$PY scripts/08_make_splits.py --config config/base.yaml

echo "=== [09] baselines (CADD / SpliceAI / AlphaMissense) ==="
$PY scripts/09b_cadd_api.py
$PY scripts/09c_alphamissense.py
$PY scripts/09d_spliceai_vcf.py
$PY scripts/09e_spliceai_run.py
$PY scripts/09f_spliceai_merge.py

echo "=== [10-11] tiered training + held-out benchmark ==="
$PY scripts/10_11_benchmark.py --allow-rerun

echo "=== [12] ablation (paired analysis) ==="
$PY scripts/12b_paired_analysis.py

echo "=== [13] counterfactual isoform experiment ==="
$PY scripts/13_counterfactual.py

echo "=== [14] SpliceAI + CADD baselines on the locked test set ==="
$PY scripts/14_add_baselines.py

echo "=== [15] figures ==="
$PY scripts/15_make_figures.py

echo "=== [16] remaining manuscript figures ==="
$PY scripts/16_make_missing_figures.py

echo "=== [17] statistical strengthening (AUPRC CIs, permutation p-values, CADD held-out subgroup) ==="
$PY scripts/17_statistical_strengthening.py

echo "=== [17b] refit statistics on the no-CADD tiers (supersedes 17's schema) ==="
$PY scripts/19_refit_statistics.py

echo "=== [18] exon architecture, backup skipping & NMD boundary analysis ==="
$PY scripts/18_exon_architecture_nmd.py

echo "=== [19] tests ==="
PYTHONPATH=src $PY -m pytest tests/ -q

echo "ALL STAGES COMPLETE"
