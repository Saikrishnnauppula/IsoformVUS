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

Note on Figure 08, panel C: the disease list, cell types, and "5 affected
protein-coding isoforms" text are hard-coded in scripts/16_make_missing_figures.py
and are not backed by a shipped results file. The manuscript (Section 5.3) cites
only the verified parts of this figure: the +0.0417 score delta, the 45 GENCODE v46
transcripts, and the M5 rank of 8 in results/recovered_variants.tsv. The
counterfactual.json file itself records 10 protein-coding transcripts with a
coding-specific effect, not 5.

## Archival — unregenerable, no manuscript claim rests on them

| Figure | Content | Status |
| --- | --- | --- |
| 09_exon_architecture_nmd.png | Exon architecture and NMD | Archival. No script in the repository generates it. |
| 10_structural_mediation_controls.png | Structural mediation controls | Archival. No script in the repository generates it. |
| 11_donor_acceptor_boundary_asymmetry.png | Donor–acceptor boundary asymmetry | Archival. Depicts the withdrawn donor-edge finding (manuscript Section 6.1). Retained only as a record of what was withdrawn, not as evidence. |
| 12_exon_skipping_candidates.png | Exon-skipping candidates | Archival. No script in the repository generates it. |

Figures 09–12 came from the stage-18 run. The five input files stage 18 needs
sit in data/interim/script18/, which is not in the shipped tree, and the code
that built them was removed with the experiment archive, so these figures cannot
currently be regenerated from the repository. Do not cite them as evidence.
