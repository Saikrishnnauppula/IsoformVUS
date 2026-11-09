# Held-out benchmark (test set)

Test: 14400 variants (7200 pathogenic), 24 genes. Gene-level splits; test set accessed per results/test_access.log.

Refit 2026-11-09: cadd_phred and cadd_missing REMOVED from all fitted tiers. CADD was queried only for validation+test variants, so every train row (n=35,921) carried the -1 sentinel and the fitted cadd_phred coefficient was fit on a constant; test-time CADD signal entered through that coefficient. CADD is now a pure external baseline (see Baselines). The 2026-11-08 rebuilt-matrix values (M1 AUPRC 0.9465, M4 0.9395) are SUPERSEDED: they were inflated by test-time CADD signal flowing through the constant-fit coefficient. The refit values agree with the pre-rebuild run (CADD constant in all splits; M1 0.8016, M4 0.8524) to within 0.0001, which cross-validates the refit. Environment: pinned requirements stack (scikit-learn 1.5.1 / pandas 2.2.2 / numpy 1.26.4).

| Model | Feats | Val AUROC | Test AUROC [95% CI] | Test AUPRC | Top-10 recall [95% CI] |
|---|---|---|---|---|---|
| M1 | 7 | 0.8253 | 0.7985 [0.7916,0.8060] | 0.8015 | 0.0014 [0.0014,0.0014] |
| M2 | 9 | 0.8242 | 0.7969 [0.7893,0.8053] | 0.8328 | 0.0014 [0.0014,0.0014] |
| M3 | 12 | 0.8303 | 0.7962 [0.7889,0.8044] | 0.8328 | 0.0014 [0.0014,0.0014] |
| M4 | 21 | 0.8328 | 0.8029 [0.7955,0.8106] | 0.8524 | 0.0014 [0.0014,0.0014] |
| M5 | 24 | 0.8325 | 0.8028 [0.7955,0.8106] | 0.8522 | 0.0014 [0.0014,0.0014] |

## Baselines (test)

- alphamissense_missense_only (n=1441): AUROC=0.9435 AUPRC=0.9750
- cadd_phred_partial (n=11240): AUROC=0.9554 AUPRC=0.9590

## Recovery: M5 vs M1 (test)

- top10: M5 top-k pathogenic=10, M1=10, recovered missed by M1=7
- top20: M5 top-k pathogenic=20, M1=20, recovered missed by M1=12
- top50: M5 top-k pathogenic=50, M1=50, recovered missed by M1=23

## Baselines added post-hoc (enrichment completed)

- SpliceAI max delta (n=13687): AUROC=0.6961 AUPRC=0.7192 | M5 minus SpliceAI AUROC=+0.1057
- CADD PHRED (n=11240): AUROC=0.9554 AUPRC=0.9590 | M5 minus CADD AUROC=-0.2406 (SAME-SUBSET head-to-head only, n=11240; M5 AUROC on this subset=0.7148). On the FULL test set (n=14400) M5 AUROC=0.8028, i.e. -0.1526 vs CADD's subset AUROC -- different populations, not a head-to-head.
