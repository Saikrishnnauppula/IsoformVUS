# Prototype findings (2026-10-06)

From `notebooks/01_clinvar_prototype.py` on `clinvar_20261004.vcf.gz`
(first 2,000,000 records; full parse in stage 03).

## Label availability (first 2M records)

| Bucket | Count |
|---|---|
| pathogenic (P/LP, no benign/uncertain mixed) | 91,813 |
| benign (B/LB) | 118,157 |
| other (VUS, conflicting, unclassified) | 1,790,030 |

Labels are plentiful: even the 2M-record sample yields ~92k pathogenic and
~118k benign candidates before review-status filtering.

## Review-status distribution

Dominant: `criteria_provided,_single_submitter` (1.45M) — the 1-star tier.
`multiple_submitters,_no_conflicts` (326k) is the ≥2-star pool. Conflicting
classifications: 78k (excluded from labels by the ≥1-star + no-conflict rule
in stage 03).

## Top genes by pathogenic count (first 2M)

APC (2,327), MSH6 (2,056), MSH2 (1,824), MLH1 (1,564), SCN1A (1,547),
USH2A (1,443), CFTR (1,098), ABCA4 (1,085), NEB (898), PMS2 (786),
TTN (760), ALMS1 (726), COL7A1 (722), CHD7 (694), DNAH5 (655).

**Design note for stage 02:** hereditary cancer-predisposition genes (APC,
MSH2/6, MLH1, PMS2) dominate raw pathogenic counts. The gene-selection
criteria exclude genes whose *only* disease association is somatic/cancer;
germline cancer-predisposition genes need an explicit include/exclude
decision in `scripts/02_build_gene_panel.py`, and the 40%-per-category
diversity cap must be enforced or the panel will be cancer-heavy.

## Panel feasibility verdict

FEASIBLE: ≥10 P/LP per gene is easily met for hundreds of genes; the binding
constraints will be (a) multi-isoform requirement, (b) cell-type data
availability, and (c) the diversity cap — not label counts.
