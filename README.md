# IsoformVUS

**Phenotype- and Cell-Type-Aware Isoform-Level Variant Prioritization for Rare-Disease VUS Interpretation**

IsoformVUS is a research project investigating a single question:

> **Does incorporating disease-relevant cell-type information and isoform usage improve
> prioritization of pathogenic rare variants compared with conventional variant-effect predictors?**

The hypothesis: a variant may have different functional consequences depending on which
transcript/isoform the disease-relevant cell type uses. IsoformVUS tests whether phenotype,
cell type, and isoform context recover pathogenic variants that conventional gene-level or
tissue-agnostic approaches rank poorly.

**This is a research project, not a clinical diagnostic tool.** It does not diagnose
patients and its scores must never be auto-converted into ACMG classifications.
See [docs/clinical_safety.md](docs/clinical_safety.md).

## Status

**Complete.** All 17 pipeline stages have been executed; the held-out
benchmark is locked and all results are in `results/`. Nothing in this
repository is fabricated: every number in `paper/manuscript.md` traces to a
file in `results/`, and the full audit trail is in `docs/final_audit.md`.
The honest headline: isoform-structural features add significant
pathogenic-variant ranking signal (AUPRC +0.051, p < 0.001); phenotype and
cell-type layers add modest ranking gains (+0.031 each); CADD v1.7 wins on
raw discrimination but our model complements it in CADD-weak regions.

## Repository layout

```
IsoformVUS/
├── README.md            # this file
├── LICENSE              # MIT (see file)
├── CITATION.cff          # citation metadata
├── requirements.txt      # pinned Python dependencies
├── environment.yml       # conda environment
├── pyproject.toml        # package metadata + build config
├── config/               # YAML configuration files
├── data/                 # raw / interim / processed (gitignored at scale)
│   └── DATA_PROVENANCE.md# source, URL, version, date, build, license per dataset
├── metadata/             # gene lists, HPO mappings, split definitions
├── src/isoformvus/       # the Python package
├── tests/                # automated tests
├── notebooks/           # exploratory notebooks (numbered)
├── scripts/             # runnable pipeline stages
├── results/             # model outputs, metrics, split assignments
├── figures/             # publication figures (built by scripts, never hand-drawn)
├── docs/                # design docs, biology explainers, limitations, PhD materials
├── research/            # literature review + novelty statement
└── paper/               # manuscript
```

## Development order

1. Literature review (`research/literature_review.md`)
2. Novelty assessment (`research/novelty_statement.md`)
3. Feasibility assessment (`docs/feasibility.md`)
4. Dataset selection (`metadata/gene_selection.md`)
5. Small prototype (`notebooks/01_*`, `scripts/01_*`)
6. Phenotype → cell mapping (`src/isoformvus/data/phenotype_cell.py`)
7. Isoform analysis (`src/isoformvus/data/isoform.py`)
8. Variant feature generation (`src/isoformvus/features/`)
9. Baselines (`scripts/09_*`)
10. IsoformVUS model (`scripts/10_*`)
11. Held-out benchmark (`results/benchmark.md`)
12. Ablation
13. Counterfactual isoform experiment
14. Case studies
15. Figures (`figures/`)
16. Manuscript (`paper/manuscript.md`)
17. GitHub cleanup
18. Final scientific audit (`docs/final_audit.md`)

## Reproducibility

- Genome build: **GRCh38** throughout (all sources mapped to GRCh38; conversions documented).
- Transcript annotation: **GENCODE** release pinned in `data/DATA_PROVENANCE.md`.
- Random seeds: every stochastic step takes `--seed` (default `20261006`), recorded in configs.
- Run the pipeline: `bash scripts/run_all.sh` (stages 1–17, deterministic order).

### Verify without re-downloading raw data

Reviewers can verify every reported number without the ~8 GB raw downloads
(`data/raw/`, `data/interim/` — regenerable via `scripts/01_*`, documented
in `data/DATA_PROVENANCE.md`):

> **Note (2026-10-13):** the exact pinned ClinVar weekly file
> (`clinvar_20261004.vcf.gz`) has rotated off NCBI's FTP — only newer weeklies
> remain. A full raw-data rerun today would use a different ClinVar snapshot,
> so exact reproduction of the locked results from raw inputs is no longer
> possible. The shipped processed tables and results JSONs remain the
> verifiable record.

- **Benchmark metrics:** `results/benchmark.json` (tier AUROC/AUPRC, baselines).
  Structure: `tiers.<M1..M5>.test.{auroc, auprc, ...}`, `tiers.<M>.test.auroc_ci95`,
  `baselines.{cadd,spliceai,alphamissense}`, `recovery_M5_vs_M1`.
  Example: `python -c "import json; b=json.load(open('results/benchmark.json')); print(b['tiers']['M4']['test']['auprc'])"`
- **Statistical strengthening:** `results/statistics.json` (paired AUPRC-difference
  CIs, permutation p-values, CADD held-out subgroup analysis) and
  `results/per_model_cis.json` (per-model 95% CIs)
- **Processed feature tables:** `data/processed/feature_matrix.tsv.gz` (26 features),
  `data/processed/variant_cadd.tsv.gz`, `data/processed/variant_spliceai.tsv.gz`
- **Figures:** `figures/01–08_*.png` (all 8, reproducible via `scripts/15_make_figures.py`)
- **Manuscript:** `paper/manuscript.md` — every number traces to a file above
  (12-number spot-check recorded in `docs/final_audit.md`)

## Installation

Requires **Python 3.10+** (tested on 3.11 and 3.12).

```bash
python -m venv .venv && source .venv/bin/activate
pip install -r requirements.txt
pip install -e .
```

Or with conda: `conda env create -f environment.yml` (pins Python 3.11).

Run the test suite:

```bash
pytest tests/ -q        # 37 tests, ~5 s
```

## Citation

See [CITATION.cff](CITATION.cff).
