# Methods

We built IsoformVUS as a gene-held-out benchmark: 60,521 ClinVar variants
across 102 Mendelian disease genes, with models trained on train genes and
evaluated on 24 genes the models never saw. A second, descriptive analysis
covered 587,227 unselected ClinVar variants in the same 102 genes. This
document records what we did, in the order we did it, with the defects we
found afterwards. Limits and withdrawn claims live in `docs/limitations.md`;
input versions and licences live in `docs/data.md`.

## Environment

We ran the pipeline in Python 3.11 with the pinned versions in
`requirements.txt` (numpy 1.26.4, pandas 2.2.2, scikit-learn 1.5.1,
pytest 8.3.3). We verified a clean install and 37/37 passing tests in a
fresh Python 3.11.17 virtual environment on 2026-10-13; the suite is
62/62 after the 2026-11-09 revision, which added data-grounded tests
for the points raised in review (tests/test_review_regressions.py). `pyproject.toml`
pins the same packages with `==`, so `pip install -e .` and `pip install
-r requirements.txt` describe the same environment; we used
`requirements.txt` for every locked result.

## Seed and splits

We used seed 20261006 everywhere: `config/base.yaml`, `scripts/08_make_splits.py`,
`scripts/10_11_benchmark.py`, and the stage-18 analysis script.

We split by gene, never by variant. Genes in the same split share constraint,
expression, and disease mechanism, so a variant-level split would let a model
memorize gene identity and score it as variant biology. We assigned genes to
train, validation, and test with a deterministic SHA-256 keyed shuffle
(`src/isoformvus/data/splits.py`, hashing `seed:gene`), stratified by disease
category. The split counts are 61 train genes (35,921 variants), 17
validation genes (10,200 variants), and 24 test genes (14,400 variants),
with zero gene overlap. We wrote the assignment once to `results/splits.json`
before any model training. One pointer in that file was wrong and is now
fixed: its `seed_note` says "see src/isoformvus/data/splits.py", where the
implementation lives.

We also ran a disease-level sensitivity split as a secondary check, holding
out whole disease families instead of genes. We did not use VUS as a training
class at any point; VUS was the application target only.

## Models M1–M5

We fit logistic regression only. Each tier is a `LogisticRegression` with
C=1.0, `class_weight="balanced"`, `max_iter=5000`, fit on train genes. Older
docs said the benchmark used LightGBM; that was wrong: the code fits
`LogisticRegression` and nothing else, and LightGBM is no longer a
dependency of the repository.

We added one feature group per tier, in this order:

- M1, variant scores: AlphaMissense score and missingness indicator, CADD
  PHRED and missingness indicator, and ClinVar review stars. Two of these
  need a warning. The CADD columns in the shipped matrix were constants
  (`cadd_phred = -1.0`, `cadd_missing = 1` for all 60,521 rows), because the
  matrix was built on 2026-10-06 and the CADD table arrived on 2026-10-07.
  The matrix was rebuilt with the real scores on 2026-11-08 (19,014 of
  60,521 rows scored) and the benchmark re-run from it; that run is
  superseded, because every train row still carried the -1.0 sentinel, so
  the fitted CADD coefficient was never fit against varying CADD values.
  On the 2026-11-09 refit the CADD columns were removed from every fitted
  tier, so no tier in the reported results uses CADD; it appears only as
  an external baseline. Review stars come from the same ClinVar records as
  the labels (`CLNREVSTAT` beside `CLNSIG`), so M1 is not a clean
  variant-only baseline; on the shipped matrix, removing review stars and
  AlphaMissense dropped M1 AUPRC to 0.6992, while M4 fell only to 0.8121,
  so the M4 gain did not depend on them.
- M2, phenotype: HPO-derived gene features. M3, cell type: GTEx tissue
  expression relevance per gene. Every feature added by M2 and M3 is a
  gene-level constant: all variants in a gene carry the same value. These
  tiers tested gene-level context offsets, not variant-level phenotype or
  cell-type effects.
- M4, isoform structure: variant-varying flags for MANE Select coding
  membership, alternative-isoform consequence, exon position, and phase.
  This tier carried the result that did not hold up at its first size: on
  the shipped matrix, test step AUPRC 0.8016 to 0.8524, a gain of +0.0509
  with 95% CI [+0.0473, +0.0546], positive in 24 of 24 test genes. After
  the 2026-11-09 refit (CADD removed from all tiers), M4 is +0.0189 above
  M1 in tie-corrected AUPRC (0.8526 against 0.8336; gene-clustered 95% CI
  [+0.0107, +0.0279]). The intervening rebuilt run, which put M4 at -0.0069
  below M1, was inflated by test-time CADD signal through a constant-fit
  coefficient and is superseded.
- M5, isoform usage: the dominant transcript's GTEx usage in the gene's top
  tissue, assigned to every variant of the gene. It is also a gene-level
  constant. On the shipped matrix it added nothing (AUPRC change about
  -0.0001); on the 2026-11-09 refit M5 is -0.0001 from M4 in tie-corrected
  AUPRC (step AUPRC change -0.0002).

We scored average precision two ways. Standard step AUPRC (right-endpoint)
penalizes tied score blocks, and M1 holds a 6,077-SNV tie block of variants
with missing AlphaMissense scores. Under analytical tie-corrected average
precision (McSherry and Najork 2008), on the 2026-11-09 refit M1 scores
0.8336, M2 0.8344, M3 0.8344, M4 0.8526, and M5 0.8525. The +0.031
step-AUPRC jumps in M2 and M3 were a tie artifact; tie-corrected, M2
gains +0.0008 and M3 +0.0008. The rebuilt run's tie-corrected values
(M1 0.9465, M4 0.9395) are superseded: on that run, real CADD values in
validation and test broke the tie block by leaking signal into the tiers
through a coefficient fit on a constant in train.

We set the benchmark prevalence by construction: each test gene contributes
exactly 300 pathogenic and 300 benign variants, so test prevalence is 50%.
All AUPRC and precision values are benchmark quantities at that prevalence,
not estimates for a clinical VUS population. Test prevalence also differs
sharply by variant type (SNV 0.352, deletion 0.915), so global metrics do
not carry across variant mixes.

We reported variant-level bootstrap and permutation intervals (5,000
bootstraps, 2,000 permutations) with Bonferroni alpha 0.01. Gene-clustered
intervals (cluster bootstrap over the 24 test genes, 2,000 resamples, plus
a leave-one-gene-out jackknife) are computed by
scripts/19_refit_statistics.py on the 2026-11-09 refit: M4-M1 +0.0189
[+0.0107, +0.0279]; M2-M1 +0.0008 [-0.0029, +0.0047]. An earlier check put
the M4 minus M1 gene-cluster interval at [+0.0425, +0.0582] on step
AUPRC, which penalizes M1's tie block; the tie-corrected interval above
supersedes it. Tier claims in the manuscript rest on the clustered
intervals, not the variant-level ones.

## Run order

We ran the pipeline in stage order 01 to 19: 19 numbered stages in 23
script invocations (stage 04 runs two scripts, stage 09 runs five, stages
10-11 share one script, and stage 19 is the pytest suite). Every stage
calls a real script: the six stub stages were deleted or implemented, and
the runner wires in `scripts/10_11_benchmark.py` (stages 10-11),
`scripts/12b_paired_analysis.py` (stage 12),
`scripts/14_add_baselines.py` (stage 14), and
`scripts/18_exon_architecture_nmd.py` (stage 18). Since 2026-11-09
`scripts/run_all.sh` has two modes: the default verify mode runs the
stages that execute from shipped data (10-11, 12, 14, 18, the refit
statistics script in place of 17, figures, tests) and runs to completion
from a fresh clone; `--full` keeps the stage order above behind a
preflight that fails fast on missing inputs. Stage 13 is verify-excluded
(it needs the GENCODE GTF and MANE summary in `data/raw/`). Two further
scripts sit outside the runner: `scripts/20_build_exon_atlas.py` and
`scripts/21_validate_exon_reconstruction.py` build the annotation exon
atlas and the annotated tensor reported as primary in manuscript
Section 4.1; they need the same `data/raw/` downloads (provenance in
`data/DATA_PROVENANCE.md`). The intended
entry points are `make test`, `make benchmark`, and `make figures`.
Stage 18 writes `results/exon_architecture_nmd.json`,
`results/exon_skipping_nmd_candidates.tsv`,
`results/vus_bayesian_stratification.tsv`, and
`results/cmh_per_gene_tables.tsv`. A full from-scratch run does not
reproduce the locked results: the ClinVar weekly they were built from
returns HTTP 404 (stage 01 now pins the live `clinvar_20261108.vcf.gz`,
verified 2026-11-09, which is a different ClinVar version), `data/raw/`
is not in the repository, stages 04a and 05 need the live GTEx API, and
stage 09b needs the live CADD API.

## Test-access guard

We locked the test set before training. `scripts/10_11_benchmark.py`
refuses to run when `results/test_access.log` exists unless `--allow-rerun`
is passed, and it appends one line per access. Scripts 12b and 18 log their
test accesses the same way. The log records 16 entries: 1 benchmark run on
2026-10-06, 8 baseline-enrichment evaluations (SpliceAI plus CADD) on
2026-10-07, 2 benchmark reruns, 1 baseline run, and 1 script-18 M6
evaluation on 2026-11-08, and on 2026-11-09 a script-18 verification
rerun, a benchmark rerun for the pinned-stack refit, and 1 baseline run.
Four earlier script-18 M6 evaluations on 2026-11-08 were not logged by
that script, so the test split was read 20 times in all. The committed
log is the original timestamped record; it sat untracked in the working tree until 2026-11-09 because
`.gitignore` carried a bare `*.log` rule, and a negation
(`!results/test_access.log`) now keeps it tracked. Any statement that the
test set was accessed once is false.

## Cutoffs and tuning sets

We set every cutoff in the analysis from data covering all 102 genes,
including the 24 test genes: exon-length bins at 180, 400, and 1,000 bp,
backup-skip usage at 10% with 1.0 median TPM, the 3 bp exon-edge window,
and the AlphaMissense benign cutoff at 0.34. We wrote no pre-registration
document, and we ran no threshold sweep in the locked analysis. The
"held-out" label applies to the final logistic fit only, not to the feature
engineering or the cutoff choices.

## Data inputs and version mismatches

We used GRCh38 throughout, GENCODE v46 for structure, MANE Select v1.5,
GTEx v8, ClinVar weekly `clinvar_20261004.vcf.gz`, CADD v1.7, SpliceAI
v1.3.1, and AlphaMissense v1. Full versions, sources, and licences are in
`docs/data.md`. Four mismatches shaped the results:

- The ClinVar weekly the locked results were built from
  (`clinvar_20261004.vcf.gz`) now returns HTTP 404 from NCBI FTP, so
  bit-exact reproduction from raw inputs is no longer possible. Stage 01
  now pins the live weekly `clinvar_20261108.vcf.gz` (verified 2026-11-09)
  with a fallback to the unversioned latest file; a rebuild uses that
  newer ClinVar version instead. `data/raw/` (~5.1 GB) is not in the
  repository and must be re-downloaded per `docs/data.md`.
- GTEx v8 was quantified on GENCODE v26. We joined it to the v46 structural
  classification on the ENST base (the identifier without its version
  number). We did not quantify the mapping loss: the unmatched-transcript
  counts live in `data/interim/`, which is not in the shipped repository.
  "Backup skipping" in this project is a transcript-abundance fraction from
  GTEx median transcript TPM, not a junction measurement. We computed no
  junction-level PSI anywhere.
- SpliceAI features were scored against the bundled GENCODE V24-era
  canonical annotation while everything else used v46. That gap was noted
  in provenance and nowhere in the manuscript.
- We never built a gnomAD allele-frequency column. The script
  `scripts/03b_gnomad_af.py` targets gnomAD v4.1.1 (`gnomad_r4`), but the
  REST API was unreachable from the build environment, and the enrichment
  was deferred. No allele-frequency feature exists in the benchmark, and
  all 60,521 benchmark rows are ClinVar B/LB controls; no gnomAD population
  controls were used.

AlphaMissense data is CC BY-NC 4.0. It cannot be redistributed as
MIT-licensed, and the AlphaMissense-derived processed table cannot be
shared under the code licence.

## The stage-18 analysis and M6

We reconstructed MANE Select exons and stratified internal coding variants
by exon length, GTEx backup skipping, and reading-frame phase, then built
M6, a logistic model that added consequence-class indicators (`is_ptv`,
`is_syn`, `mane_sev`, `alt_max_sev` and related flags) and exon-architecture
interactions to the M1 logit. For internal missense, M6 fits
exon-architecture offsets on train genes only, with the M1 logit held
fixed. The fitted donor-edge coefficient is +0.31 (95% CI [-0.86, 1.48]).
The SpliceAI interaction terms are not identified: no train missense
variant carries a SpliceAI score. The earlier +3.60 and +1.80 constants
were hard-coded, not fitted. Script 18 refits M1 with a standardized
logistic regression; its pooled M1 (AUROC 0.7981) is a different estimator
from the benchmark's M1 (AUROC 0.7985 on the 2026-11-09 no-CADD refit),
and the within-missense comparison in the manuscript is internal to
script 18. Script 18 reads its five large inputs from
`data/interim/script18/` (see `docs/data.md`). We record this here because
the numbers in `results/` came from that code; the model itself was
withdrawn after peer review, and `docs/limitations.md` gives the reasons.

## Prior work

We placed the project against the closest published methods before building
it. IMPPROVE (Chapman and Lassmann 2025) was the nearest comparator: per-HPO
random forests on CADD plus gene-level expression, with no isoform layer.
Expression-weighted annotation (pext, Cummings et al. 2020) had already
shown that low-expression exons are depleted of pathogenic loss-of-function
variants; our backup-skipping variable is a cruder version of the same
idea, built from annotated-transcript TPM fractions instead of junction
measurements, applied to missense rather than loss of function. The biology
behind the exon gradient was established before we started: tissue-specific
exons often encode disordered protein segments (Buljan et al. 2012),
isoform-selective alleles were demonstrated in patients (Ly et al. 2025),
and disease-associated variants map to tissue-specific alternative isoforms
(Weykopf et al. 2026). Full references for the works cited here are in
paper/manuscript.md, Section 7. Our October follow-up work then falsified the
phenotype and cell-type framing we had started with: the M2 and M3 features
are gene-level constants, their AUPRC gains are a tie artifact. An earlier
version of this paragraph also claimed that a validation-trained
CADD-plus-context stack beat CADD alone by +0.0191 AUPRC, of which the
isoform-specific slice was +0.0015 and the isoform block alone on top of
CADD added exactly 0.0000. Those three numbers appear in no results file,
no shipped script fits a validation-trained stack (every fit in the repo
trains on the train split), and they could not be reproduced from shipped
files; they were removed on 2026-11-09 as unreproducible (audit trail:
results/statistics.json, key methods_209_check: the closest computable
train-fit quantity, M5 minus CADD AUPRC on the CADD-scored test subset, is
-0.2367, not +0.0191). What survived that check was the
tie correction itself, and the M4 structural signal at a smaller size: on
the 2026-11-09 refit M4 is +0.0189 above M1 in tie-corrected AUPRC
(gene-clustered 95% CI [+0.0107, +0.0279]). The rebuilt run's -0.0069
reading was a CADD artifact (manuscript Section 3.1).
