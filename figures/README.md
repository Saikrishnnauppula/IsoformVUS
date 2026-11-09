# Figures — provenance and status

## Regenerable (scripts/15_make_figures.py, scripts/16_make_missing_figures.py)

| Figure | Content | Source |
| --- | --- | --- |
| 01_architecture.png | Schematic of the M1–M5 tiers and feature counts (7/9/12/21/24 on the 2026-11-09 refit) | data/processed/feature_tiers.json, results/splits.json |
| 02_workflow.png | Schematic of the data and analysis workflow | results/*.json, data/DATA_PROVENANCE.md |
| 03_isoform_usage.png | Per-variant effect classification and GTEx isoform-usage evidence | data/processed/variant_isoform_effects.tsv.gz, data/processed/isoform_usage.tsv.gz |
| 04_benchmark.png | Held-out test AUROC/AUPRC by tier | results/benchmark.json |
| 05_ablation.png | Paired AUROC differences, 2026-11-09 no-CADD refit | results/statistics.json (paired_auroc_diffs_refit_no_cadd) |
| 06_topk_recovery.png | M5 vs M1 top-k recovery | results/benchmark.json (recovery_M5_vs_M1) |
| 07_counterfactual.png | LMNA chr1:156115000:C>T per-transcript effects and M5 score delta | results/counterfactual.json |
| 08_case_study.png | LMNA case study (same counterfactual, plus interpretation chain) | results/counterfactual.json |

Note on Figure 08, panel C: the panel now reads its counts from
results/counterfactual.json in code ("LMNA → 45 transcripts; coding-specific
on 18"). The 18 is the file's total coding_specific count across biotypes;
the narrower count of protein-coding transcripts with a coding-specific
effect is 10, a different quantity. The disease list and cell types in the
panel remain descriptive text. The manuscript (Section 5.3) cites only the
verified parts of this figure: the +0.0417 score delta, the 45 GENCODE v46
transcripts, and the M5 rank of 8 in results/recovered_variants.tsv.

Note on Figure 02: the top banner is generated in code from
results/benchmark.json and results/statistics.json (refit values: test
AUPRC 0.8015 → 0.8524 step; tie-corrected 0.8336 → 0.8526, +0.0189,
gene-clustered 95% CI [+0.0107, +0.0279]; 23/50 top M5 pathogenic variants
missed by M1). No banner number is hard-coded.

## Regenerable since 2026-11-09 (scripts/16_make_missing_figures.py)

Figures 09–12 were archival and unregenerable until 2026-11-09: no script
in the repository generated them and their stage-18 inputs were not
shipped. The inputs now ship in data/interim/script18/, and
scripts/16_make_missing_figures.py replots all four from the shipped
stage-18 outputs. The replots are new renderings of the same result
files, not the original images (layout and palette differ), and a few
stale values in the originals are corrected by the shipped data: figure
11 uses the shipped donor-edge LR of 1.817, not 4.52, and drops the
ACMG-point threshold lines; figure 12 uses the shipped M6
internal-missense AUROC of 0.9321, not 95.1%. No manuscript claim rests
on any of them. Do not cite them as evidence.

| Figure | Content | Source |
| --- | --- | --- |
| 09_exon_architecture_nmd.png | Exon tensor (missense resolved P/LP by size x skip x symmetry), NMD zones, in-frame indel cells | results/exon_architecture_nmd.json (exon_phase_symmetry, inframe_indels_and_nmd_zones). This is the reconstructed-geometry tensor; the manuscript's primary tensor is the annotated one in results/exon_reconstruction_validation.json (Section 4.1). |
| 10_structural_mediation_controls.png | Label-free biophysical controls; mediation and sensitivity panel | results/exon_architecture_nmd.json (biophysical_and_density_controls, sensitivity_checks) and results/withdrawn_carried_over_values.json (causal_mediation). The mediation values are carried over, not recomputed, and the mediation percentages are withdrawn (manuscript Section 6.1); the figure says so on its face. |
| 11_donor_acceptor_boundary_asymmetry.png | Donor–acceptor boundary asymmetry; per-cell likelihood ratios | results/withdrawn_carried_over_values.json (splice_junction_boundary_asymmetry) and results/exon_architecture_nmd.json (tiered_model_and_bayesian_lr). Depicts the withdrawn donor-edge finding (manuscript Section 6.1); the distance profile is carried over, not recomputed. Retained only as a record of what was withdrawn, not as evidence. |
| 12_exon_skipping_candidates.png | Top-15 nominated exons (Tier 1/2 PTV records); stage-18 model comparison | results/exon_skipping_nmd_candidates.tsv (Exon_Nomination_Score column); results/exon_architecture_nmd.json (tiered_model_and_bayesian_lr) |
