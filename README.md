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
on the shipped matrix. On 2026-11-09 the benchmark was refit with CADD
removed from every fitted tier: the 2026-11-08 rebuilt run had let
test-time CADD signal in through a coefficient fit on a constant in train
(0 of 35,921 train rows scored), inflating every tier to about 0.94, and
its "M2 is the best tier" result is withdrawn. On the refit, tie-corrected
test AUPRC is M1 0.8336 and M4 0.8526 (+0.0189; gene-clustered 95% CI
[+0.0107, +0.0278] over the 24 test genes). M2 and M3 add nothing
(gene-clustered intervals include zero) and M5 adds nothing over M4. Step
AUPRC (M1 0.8015, M4 0.8524) overstates the gain, because step evaluation
penalizes M1's tie block. CADD is an external baseline only; on the 11,240
test variants that carry a CADD score it beats every tier on AUPRC
(0.9590). Refit statistics are recomputed by
`scripts/19_refit_statistics.py`.

The rest of the manuscript's headline claims were withdrawn after that
review. The M6 model (test AUROC 0.9923) is consequence-class separation, not
exon biology: a consequence-only model reaches 0.9848, and within missense
M6 scores 0.9321 against 0.9511 for the M1 baseline. The VUS likelihood
ratios, the ASO candidate list, the donor-edge finding, and the
"co-transcriptional NMD" mechanism were all withdrawn. The test set was
accessed 13 times as recorded in results/test_access.log (17 reads in all,
counting 4 unlogged script-18 evaluations): 1 benchmark run on 2026-10-06,
8 baseline-enrichment evaluations on 2026-10-07, 2 benchmark reruns and
1 baseline run on 2026-11-08, and 1 script-18 M6 evaluation on 2026-11-08.
The log is the original timestamped record. It sat untracked in the
working tree until 2026-11-09 because `.gitignore` carried a bare `*.log`
rule; a negation (`!results/test_access.log`) now keeps it tracked. Full
detail, with numbers, is in `docs/limitations.md`.

## Install

```bash
python3.11 -m venv .venv && source .venv/bin/activate
pip install -r requirements.txt
pip install -e .
```

`requirements.txt` pins the tested versions (numpy 1.26.4, pandas 2.2.2,
scikit-learn 1.5.1, pytest 8.3.3). `pyproject.toml` pins the same packages
with `==`, so both files describe the same tested environment.

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

Bit-exact reproduction from raw inputs is no longer possible: the ClinVar
weekly the locked results were built from (`clinvar_20261004.vcf.gz`) now
returns HTTP 404 from NCBI FTP, and `data/raw/` (~5.1 GB) is not in the
repository. Stage 01 now pins the live weekly `clinvar_20261108.vcf.gz`
(verified 2026-11-09), with a fallback to the unversioned latest file; a
rebuild from any newer release uses a different ClinVar version and will
not reproduce the locked results. The shipped `data/processed/` tables and
`results/*.json` files are the verifiable record. See `docs/data.md` and
`docs/methods.md`.

One pipeline script is deliberately not part of the rebuild:
`scripts/03b_gnomad_af.py` is not wired into `scripts/run_all.sh`, and no
gnomAD allele-frequency feature exists in the shipped feature matrix.

## Verify without re-downloading raw data

The shipped `data/processed/` tables and `results/` files are enough to
check every reported number without the ~5.1 GB raw download:

```bash
make test                                    # pytest suite (includes data-grounded checks on the shipped files)
python3 scripts/10_11_benchmark.py --allow-rerun   # retrain M1–M5 on the shipped matrix, re-evaluate the test genes
make figures                                 # regenerate figures from results/*.json
```

`make test` runs on Python 3.11 with `PYTHONPATH=src`; the data-grounded
tests in `tests/test_shipped_evidence.py` read the shipped matrix, splits,
benchmark, statistics, and test-access log directly. The benchmark rerun
appends a line to `results/test_access.log` (that is what the guard is
for); its outputs should match `results/benchmark.json`.

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
├── scripts/         # pipeline stages and run_all.sh (19 numbered stages, 23 script invocations: stage 04 runs two scripts, stage 09 runs five, stages 10–11 share one)
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
