# IsoformVUS

[![CI](https://github.com/Saikrishnnauppula/IsoformVUS/actions/workflows/ci.yml/badge.svg)](https://github.com/Saikrishnnauppula/IsoformVUS/actions/workflows/ci.yml)
[![License: MIT](https://img.shields.io/badge/License-MIT-yellow.svg)](LICENSE)

IsoformVUS tests whether transcript and exon structure add signal to
rare-disease variant prioritization, on top of CADD, SpliceAI, and
AlphaMissense. The benchmark holds out whole genes: 60,521 ClinVar variants
across 102 Mendelian disease genes, split into 61 train, 17 validation, and
24 test genes (35,921 / 10,200 / 14,400 variants).

**This is research code, not a clinical tool.** Nothing in this repository
diagnoses a patient, and no model output may be turned into an ACMG/AMP
pathogenicity classification. See `docs/limitations.md`.

## Status after peer review

A journal review panel rejected the manuscript built from this work in
November 2026 (REJECT, confidence 4/5). The core M1–M5 benchmark reproduces
on the shipped matrix, but after rebuilding the matrix with real CADD values
the M1–M4 gain is gone: M1 test AUPRC 0.9465, M4 0.9395 (-0.0070), and M2 is
the best tier at 0.9533. Every train row still carries the -1.0 CADD
sentinel, so the rebuilt tier comparison is provisional until CADD is
queried for the train variants.

The rest of the manuscript's headline claims were withdrawn after that
review. The M6 model (test AUROC 0.9923) is consequence-class separation, not
exon biology: a consequence-only model reaches 0.9848, and within missense
M6 scores 0.9321 against 0.9511 for the M1 baseline. The VUS likelihood
ratios, the ASO candidate list, the donor-edge finding, and the
"co-transcriptional NMD" mechanism were all withdrawn. The test set was
accessed 13 times as recorded in results/test_access.log (17 reads in all,
counting 4 unlogged script-18 evaluations): 1 benchmark run on 2026-10-06,
8 baseline-enrichment evaluations on 2026-10-07, and 4 more during this
revision. Full detail, with numbers, is in `docs/limitations.md`.

## Install

```bash
python3.11 -m venv .venv && source .venv/bin/activate
pip install -r requirements.txt
pip install -e .
```

`requirements.txt` pins the tested versions (numpy 1.26.4, pandas 2.2.2,
scikit-learn 1.5.1, pytest 8.3.3). `pyproject.toml` lists the same packages
unpinned, so install from `requirements.txt` for the tested environment.

## Reproduce the results

The intended interface is the Makefile and CI:

```bash
make test       # pytest suite on Python 3.11
make benchmark  # retrain M1–M5 on train genes, evaluate on the 24 test genes
make figures    # regenerate figures from results/*.json
```

The working commands are `make test PYTHON=.venv/bin/python` (or
`PYTHONPATH=src pytest tests/`) and
`python3 scripts/10_11_benchmark.py --allow-rerun` (the benchmark refuses to
rerun while `results/test_access.log` exists unless `--allow-rerun` is
passed, and every run appends a line to that log). CI runs `pytest` on
Python 3.11 plus ruff and black checks on every push.

Bit-exact reproduction from raw inputs is no longer possible: the pinned
ClinVar weekly (`clinvar_20261004.vcf.gz`) now returns HTTP 404 from NCBI
FTP, and `data/raw/` (~5.1 GB) is not in the repository. The shipped
`data/processed/` tables and `results/*.json` files are the verifiable
record. See `docs/data.md` and `docs/methods.md`.

## Repo layout

```text
IsoformVUS/
├── README.md, CHANGELOG.md, CITATION.cff, LICENSE
├── requirements.txt, environment.yml, pyproject.toml
├── config/          # YAML configuration (seed 20261006, GENCODE v46)
├── data/            # processed feature tables (provenance in docs/data.md)
├── metadata/        # 102-gene panel and selection criteria
├── src/isoformvus/  # core package (splits, features)
├── tests/           # pytest suite
├── scripts/         # pipeline stages and run_all.sh (19 stages)
├── results/         # benchmark.json, statistics.json, splits.json, test_access.log
├── figures/         # publication figures
├── docs/            # methods.md, data.md, limitations.md
└── paper/           # manuscript.md (under revision after peer review)
```

## Data and licenses

Code is MIT. The data is not: the AlphaMissense-derived table is
CC BY-NC 4.0, CADD is academic/non-commercial, and the processed tables mix
sources with different terms. Do not redistribute the shipped processed
tables as MIT-licensed. Every input, with version, source, and licence, is
listed in `docs/data.md`; `data/DATA_PROVENANCE.md` records download dates
and checksums.

## Citation

```bibtex
@software{isoformvus,
  author  = {Sai Krishna Uppula},
  title   = {IsoformVUS: Exon Architecture, Alternative Isoform Expression,
             and NMD in Rare-Disease Variant Interpretation},
  version = {0.1.0},
  url     = {https://github.com/Saikrishnnauppula/IsoformVUS}
}
```

## Status and limitations

Version 0.1.0. The M1–M5 benchmark was rebuilt with real CADD values and
re-run; the old +0.0509 gain is gone (M1 AUPRC 0.9465, M4 0.9395, M2 best at
0.9533). The M6, VUS, ASO, donor-edge, and NMD-mechanism claims from the
earlier manuscript were withdrawn after peer review. Every train row still
carries the -1.0 CADD sentinel, so the rebuilt tier comparison is
provisional. All of it, with the numbers behind each withdrawal, is in
`docs/limitations.md`.
