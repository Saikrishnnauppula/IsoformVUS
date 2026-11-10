# Limitations

This document lists what is wrong or unproven in IsoformVUS, with the
numbers behind each item. An internal review of the manuscript in
November 2026 recommended against its original form; the items here are
the reasons.

## Research code only

IsoformVUS is a research project, not a clinical diagnostic tool.
Nothing here diagnoses a patient or determines clinical care. Scores are
research prioritization signals. Model output must never be converted
into ACMG/AMP pathogenicity classifications (Pathogenic, Likely
Pathogenic, VUS, Likely Benign, Benign) without independent expert
review and orthogonal evidence. Keep three things apart: computational
evidence (what the model computed and why), research prioritization (a
ranked shortlist for follow-up experiments), and clinical interpretation
(a judgment by qualified professionals with patient context, family
data, and functional validation). Do not deploy IsoformVUS scores in a
clinical pipeline, do not present them to patients or families as
results, and do not cite them as standalone evidence in a diagnostic
report. The manuscript's earlier clinical assignments were withdrawn,
and this document and the manuscript now agree on that point.

## The tie artifact and the phenotype and cell-type layers

The +0.031 step-AUPRC gains reported for M2 (phenotype) and M3 (cell type)
were a tie artifact. M1 holds a 6,077-SNV tie block of variants with
missing AlphaMissense scores, and standard right-endpoint step AUPRC
penalizes that block by -0.0318. Under analytical tie-corrected average
precision (McSherry and Najork 2008), on the 2026-11-09 refit M1 scores
0.8336, M2 0.8344, and M3 0.8344, so the tie-corrected gains are +0.0008
and +0.0008. On that refit (CADD removed from all fitted tiers), the
tie-corrected values are M1 0.8336, M2 0.8344 (+0.0008), M3 0.8344
(+0.0008), M4 0.8526 (+0.0189), M5 0.8525 (+0.0188); gene-clustered 95%
CIs over the 24 test genes include zero for M2-M1 ([-0.0029, +0.0047])
and exclude zero for M4-M1 ([+0.0107, +0.0279]). The old step-AUPRC gain
(+0.0509) overstates the M4 effect because step evaluation penalizes the
tied M1 baseline; the tie-corrected gain is +0.0189. The rebuilt run's
values (M1 0.9465, M4 0.9395, M2 0.9533 as "best tier") are superseded:
real CADD in validation and test entered through a coefficient fit on a
constant in train (0 of 35,921 train rows scored). Both layers are
gene-level constants: every variant in a gene carries the same value, so
they tested gene-level context offsets, not variant-level phenotype or
cell-type effects. The claim that 94.5% of M4's gain reflects canonical
coding-sequence membership was withdrawn as well: three derivations from
the shipped numbers give 94.0%, 94.9%, and 96.8%, and none gives 94.5%.

## Test-set access and tuning on test genes

The test set was accessed 16 times as recorded in `results/test_access.log`,
not once: 1 benchmark run on 2026-10-06, 8 baseline-enrichment evaluations
on 2026-10-07, 2 benchmark reruns, 1 baseline run, and 1 script-18 M6
evaluation on 2026-11-08, and on 2026-11-09 a script-18 verification rerun,
a benchmark rerun for the pinned-stack refit, and 1 baseline run. Four
earlier script-18 M6 evaluations on 2026-11-08 were not logged by that
script, so the test split was read 20 times in all during this revision.
The committed log is the original timestamped record (16 entries,
2026-10-06 to 2026-11-09). It
sat untracked in the working tree until 2026-11-09 because `.gitignore`
carried a bare `*.log` rule; a negation (`!results/test_access.log`) now
keeps it tracked. The "held-out" label applies
to the final logistic fit only. Every threshold and constant in the analysis
was set from work covering all 102 genes, including the 24 test genes: the
exon bins at 180, 400, and 1,000 bp, backup-skip usage at 10% with 1.0 median
TPM, the 3 bp edge window, and the AlphaMissense cutoff at 0.34. The
donor-edge offsets in M6 are now fit on train genes only (donor-edge +0.31,
95% CI [-0.86, 1.48]); the SpliceAI interaction terms are not identified,
because no train missense variant carries a SpliceAI score. The earlier
+3.60/+1.80 constants were hard-coded, not fitted. No cutoff was
pre-registered and no threshold sweep was run.

## M6 and the consequence-class confound

M6's test AUROC of 0.9923 is consequence-class separation, not exon
biology. A consequence-class-only logistic regression reaches 0.9848 on
the same test set. The test set is 9.5% missense and is dominated by
near-deterministic classes: frameshift 99.97% pathogenic, stop-gained
99.94%, synonymous 0.20%, intronic 3.47%. Within missense, the only class
where exon architecture could matter, M6 scores 0.9321 against 0.9511 for
M1 and 0.9493 for M4: M6 is worse than the baseline. Script 18 refits M1
with a standardized logistic regression; its pooled M1 (AUROC 0.7981) is
a different estimator from the benchmark's M1 (AUROC 0.7985 on the
2026-11-09 refit); the within-missense comparison is internal to
script 18. The donor-edge recall of 42.9% is 3 of 7 variants,
all in validation genes; the test set contributes zero pathogenic edge
variants, so the recall on held-out test genes is 0 of 0. The fitted M6
recalls 0 of 7. No consequence-aware baseline (consequence class plus
gnomAD AF, CADD, SpliceAI, REVEL, and AlphaMissense) was ever run against
M6, and no gnomAD AF column exists in the benchmark at all.

## The 4x2x2 gradient and cell 4D

As ClinVar arithmetic, the gradient holds on both exon geometries. On
the reconstructed frame (script 18), resolved missense pathogenicity
falls from 73.00% (3,741 P/LP of 5,125) in compact, constitutive,
phase-symmetric exons to 1.17% (3 P/LP of 257) in giant, backup-skipped,
symmetric exons, an odds ratio of 0.0044; its decisive cell, 4D, is a
single gene: TTN, 4 exons (3947, 3971, 3981, 4274), 3 pathogenic
variants, 94.0% unresolved (4,041 of 4,298 missense variants are VUS or
conflicting). On the annotated frame (scripts 20–21, GENCODE v46 / MANE
v1.5 atlas), the same cells read 78.11% (5,625 of 7,201) and 1.38%
(8 of 581), and 4D is 7 exons in 3 genes, 91.8% unresolved. The
manuscript reports the annotated tensor as primary. The reconstruction
itself does not survive the atlas: only 3,761 strict-internal CDS exons
exist in these transcripts, against 4,509 reconstructed "internal"
exons; 14.65% of matched reconstructed exons differ in length from
annotation by more than 10%, and 625 have no annotated counterpart
(593 of them single-base artefacts of deletion-POS clustering). Script
18's other geometry consumers — d_junc edge calls, the M6 features, the
Tier 1 nomination — still use reconstructed spans; the atlas measures
their input error but does not correct those outputs.

The within-gene CMH odds ratios for the gradient (0.2005, 52 informative
genes; 0.3461, 42 informative genes, reconstructed frame) are computed
by scripts/18_exon_architecture_nmd.py with per-gene tables in
results/cmh_per_gene_tables.tsv. The p-values that shipped with them
(4.29e-187; 3.53e-42) treated variants as independent within genes when
the effective sample size is genes and exons; script 18 now withholds
the p-values and z statistics at source. Gene-cluster bootstrap
intervals (2,000 resamples of genes, seed 0) are [0.1119, 0.3384] for
length and [0.1280, 1.0299] for backup skipping on the reconstructed
frame, and [0.0946, 0.3247] (OR 0.1860) and [0.1090, 0.9218] (OR 0.2720)
on the annotated frame. The length association survives clustering on
both frames; the backup-skip interval includes 1 on the reconstructed
frame, and the backup odds ratio moves 0.6063 / 0.3461 / 0.3650 across
usage thresholds 0.05 / 0.10 / 0.20, so backup skipping is reported as
descriptive. A bounding analysis over the unresolved variants
(vus_bounding_analysis in the stage-18 JSON; annotated bounds in
results/exon_reconstruction_validation.json) puts the true 1B rate in
[11.07%, 95.91%] and the true 4D rate in [0.07%, 94.09%] on the
reconstructed frame: the intervals overlap completely, and the gradient
is a property of the resolved subset. A jackknife excluding BRCA1,
BRCA2, APC, TTN, and NF1 leaves the pooled gradient at 73.39% versus
7.05%, but that check empties the reconstructed 4D cell entirely, since
that cell is TTN alone. Exon length was never adjusted for disorder,
domain content, or within-gene position, and the median AlphaMissense
score falls 0.943 to 0.104 across the same cells, so length tracks
structural constraint rather than acting as a separate effect.

## The donor-edge finding

The section 4.1 donor-edge finding was withdrawn. The quoted counts
(donor 15 of 23 pathogenic versus acceptor 2 of 15, OR 12.19, Fisher p =
0.0024) do not reproduce from the shipped data under any filter tried:
all splits give donor 54 (21 pathogenic) and acceptor 50 (13), and test
genes alone give donor 10 (0) and acceptor 5 (0). The generating code
does not exist in the repository, the counts arrive through a missing
intermediate file, and the shipped JSON contradicts itself on the
synonymous counts for the same window (2,832 versus 5,390). The Fisher p
is also mis-rounded: recomputation gives two-sided p = 0.0025. No
threshold sweep, no multiplicity correction, and no SpliceAI or
MaxEntScan comparison was run for these variants.

## The three isoform regimes

All three isoform-regime rates are withdrawn in full (2026-11-09). The
regime classification sat on the 587k variant table in
`data/interim/script18/` (shipped since 2026-11-09), and the code that
built it was deleted with the experiment archive, so no regime count — neither
the original rates (Regime I 99.93% benign, Regime II 89.20%, Regime III
28.57%) nor the corrected rates published in the previous revision
(98.38%, 730 of 742; 35.02%, 145 of 414) — can be regenerated from
shipped materials. A correction to an untraceable number is still
untraceable, so the manuscript reports none of them (Section 6.1,
item 8).

## The NMD mechanism

The "co-transcriptional NMD" claim and the 60-90% figure were withdrawn.
EJC-dependent NMD is translation-dependent and cytoplasmic; nothing in
this project shows NMD occurring co-transcriptionally, and no RNA data of
any kind was presented. The sentence claiming nuclear splicing
incorporates the exon into 60-90% of nascent transcripts has no source,
citation, or derivation anywhere in the repository. The long-exon NMD-escape rule is now cited in the manuscript (Section
4.2): Lindeboom et al. 2016 established that PTCs in exons over about
400 nt escape NMD (61% versus 98% efficiency), and the start-proximal
rule (PTCs within 150 nt of the start codon) is cited there too. Neither
rule is modeled: the flagship giant exons are exactly where NMD is
predicted to be weakest, yet the manuscript claimed PTVs there stay
99.85-100% pathogenic because NMD fires uniformly. The PTV
pattern itself cannot test NMD: across all 16 tensor cells, resolved PTVs
are 99.93% pathogenic (78,284 P/LP of 78,337), which is what the PVS1
curation rule predicts for a panel enriched for haploinsufficient
loss-of-function genes, with no variance left for any exon feature to
explain. "Backup skipping" is a transcript-abundance fraction from GTEx
median transcript TPM, not a junction measurement; we computed no
junction-level PSI anywhere.

## VUS likelihood ratios

The VUS likelihood ratios were withdrawn as clinical assignments. They
were derived from the same 587,227-variant ClinVar table they were then
applied to, so optimism is built in rather than possible. The derivation
is now stated (LR = P(cell|P/LP) / P(cell|B/LB), resolved missense only,
95% CI on the log scale), and the files agree: cell 4D is 0.0088 (95% CI
[0.0028, 0.0274]) in both the exon JSON and the per-variant VUS file. The
earlier version listed three inconsistent values (0.0102, 0.0085, 0.0088)
with no stated formula and no confidence interval. The neighbouring 4C
cell shows how fragile the assignments were: its interval runs 0.0275 to
0.0511, and its upper bound sits 0.002 below the Strong Benign threshold
of 0.053, so one more pathogenic variant breaks the assignment. The ratios also
double-count the evidence curators already use: the withdrawn mediation
analysis claimed 77.9% of the exon-length effect and 94.7% of the
backup-skipping effect runs through AlphaMissense (those percentages are
themselves withdrawn — no generating code is in the repository), so a
variant in a giant backup-skipped exon would attract
BP4 from its low AlphaMissense score and a further -4 points from an exon
ratio that is mostly the same signal. The ratio was never computed within
AlphaMissense bands, so no residual value was ever shown. Applying the
resolved-cell ratio to the unresolved pool assumes resolved variants are
exchangeable with VUS; the data contradict that, since resolvability
itself tracks evidence strength and exon context. The 3,279 cell-4D VUS
relabeled Likely Benign are all TTN, and 3 pathogenic variants already
sit inside cell 4D. Counted on 2026-11-09: none of the 23,485 rows in
results/vus_bayesian_stratification.tsv sits in a stratum whose
likelihood-ratio interval crosses 1.0, so that particular failure mode
does not occur in the shipped file; the file stays research-only for the
reasons above.

## ASO candidates

The 302 Tier 1 nominated exons (Exon_Nomination_Score in
results/exon_skipping_nmd_candidates.tsv) are a computational nomination
with no functional evidence. Missense tolerance does not equal tolerance of
whole-exon deletion: the list includes BRCA2 exon 956 (4,932 bp), BRCA1
exon 11 (3,423 bp), and TTN exon 3981 (17,097 bp). No candidate has
protein-level evidence of a functional skipped isoform, and in
haploinsufficient genes skipping yields a hypomorph at best. The
"10,826 pathogenic patient PTVs" counts ClinVar records, not patients.
The manuscript's disease-category labels for the flagship genes
contradict its own source table, which lists BRCA1 as digestive_system
and BRCA2 as nervous_system. The tolerance score behind the tiers is
computed from ClinVar missense labels over all genes including test
genes, so it is a label-derived proxy, not a functional measurement.

## Mediation

The mediation analysis is associational, not causal. AlphaMissense is a
trained predictor, not a measured biological mediator: it is a function of
the same sequence features that define the exposure, so sequential
ignorability cannot hold, and every unmeasured confounder of the
exon-pathogenicity relationship is also a confounder of the
AlphaMissense-pathogenicity relationship. The reported overlap
percentages (77.9% and 94.7%) are withdrawn: the generating code for
those numbers is not in the repository, so the mediation result is not
verifiable from shipped materials, and no percentage is reported in the
manuscript. No swap test with ESM1b or EVE was run to check whether the
overlap is specific to 3D structure.

## Labels and prevalence

Two thirds of the benchmark labels are single-submitter: 40,404 of 60,521
rows are 1-star (66.8%), with 18,531 at 2-star, 1,584 at 3-star, and 2 at
4-star. Single-submitter classifications are the most likely to be
PVS1-only calls, and the 587,227-row table carries no review-star column
at all, so the manuscript's 2-star sensitivity check cannot be verified
from the shipped ClinVar table. Test prevalence is constructed at 50% by
design (300 pathogenic and 300 benign per test gene), so AUPRC and
precision values do not transfer to a clinical VUS population. Review
stars are target-adjacent: they come from the same ClinVar records as the
labels and sit in M1 as a feature. The top 6 genes (TTN, BRCA2, ATM, NF1,
APC, BRCA1) hold about 22% of all 587,227 variants, and ClinVar counts
have no testing-volume denominator, so pooled per-kilobase densities are
submission rates as much as biology.

## Feature matrix and benchmark provenance

The feature matrix was rebuilt with real CADD values on 2026-11-08 and
the benchmark re-run from it. The old matrix carried `cadd_phred = -1.0`
and `cadd_missing = 1` for all 60,521 rows while `variant_cadd.tsv.gz`
held real scores for 24,600 rows with 100% key overlap, and the old
benchmark could not have been produced from it. The rebuilt matrix
scores 19,014 of 60,521 rows; every train row still carries the -1.0
sentinel, so the fitted CADD coefficient was never fit against varying
CADD values, and the rebuilt tier comparison was superseded the next day
by the 2026-11-09 refit, which removed CADD from every fitted tier.
`statistics.json` marks the rebuilt values as superseded and the
M2-minus-M1 (+0.0310) and M3-minus-M1 (+0.0317) pre-rebuild differences
not significant, with a tie-artifact note, and `results/splits.json`
points to `src/isoformvus/data/splits.py` (the generator string in
`src/isoformvus/data/splits.py` now matches it).

## References and counts

The reference list in the reviewed manuscript was rebuilt from primary
sources: of 15 entries, 6 were wrong, 3 partly wrong, and 2 cited works
were missing from it entirely (CADD, Kircher et al. 2014; CEP290, Taylor
et al. 2022). The review's proposed SpliceAI and Pozo fixes were both wrong
on re-check: the SpliceAI DOI 10.1016/j.cell.2019.01.017 is not that paper's
DOI, and the Pozo DOI 10.1186/s12864-018-5013-2 belongs to Bhuiyan et al.
2018. The rebuilt list has 11 entries, only works the manuscript text
cites. Two entries pointed at the wrong PMC article. Three cited entries
were multi-field chimeras with wrong authors, journal, pages, DOI, and
title. Smaller count defects from the same review: the section 3.1
component counts sum to 311,249, not the stated 301,264, and the difference
is exactly the 9,985 in-frame indels; the shipped exon counts are now
4,593 reconstructed and 4,509 internal, while the earlier 4,132/3,955
matched no reconstruction run and the counts move with sort tie-ordering
(a stable sort gives 4,578/4,494); and the manuscript header's GENCODE
v44 use is now stated in Methods (E8 annotation-fragility experiment
only).
