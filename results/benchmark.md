# Held-out benchmark (test set)

Test: 14400 variants (7200 pathogenic), 24 genes. Gene-level splits; test set accessed per results/test_access.log.

Rerun 2026-11-08 from the rebuilt feature matrix with real CADD values (19,014 scored of 60,521; 11,240 of 14,400 test). CADD was queried only for validation+test variants, so every train row carries the -1 sentinel and the fitted cadd_phred coefficient is fit on a constant; test-time CADD signal enters through that coefficient. The pre-rebuild values (M1 AUPRC 0.8016, M4 0.8524) were computed with CADD constant at -1 in all splits.

| Model | Feats | Val AUROC | Test AUROC [95% CI] | Test AUPRC | Top-10 recall [95% CI] |
|---|---|---|---|---|---|
| M1 | 9 | 0.9410 | 0.9484 [0.9444,0.9525] | 0.9465 | 0.0014 [0.0014,0.0014] |
| M2 | 11 | 0.9504 | 0.9584 [0.9545,0.9622] | 0.9533 | 0.0014 [0.0014,0.0014] |
| M3 | 14 | 0.9387 | 0.9454 [0.9410,0.9500] | 0.9312 | 0.0014 [0.0014,0.0014] |
| M4 | 23 | 0.9380 | 0.9396 [0.9355,0.9436] | 0.9395 | 0.0014 [0.0014,0.0014] |
| M5 | 26 | 0.9431 | 0.9461 [0.9420,0.9498] | 0.9454 | 0.0014 [0.0014,0.0014] |

## Baselines (test)

- alphamissense_missense_only (n=1441): AUROC=0.9435 AUPRC=0.9750

## Recovery: M5 vs M1 (test)

- top10: M5 top-k pathogenic=10, M1=10, recovered missed by M1=5
- top20: M5 top-k pathogenic=20, M1=20, recovered missed by M1=9
- top50: M5 top-k pathogenic=50, M1=50, recovered missed by M1=19

## Baselines added post-hoc (enrichment completed)

- SpliceAI max delta (n=13687): AUROC=0.6961 AUPRC=0.7192 | M5 minus SpliceAI AUROC=+0.2499
- CADD PHRED (n=11240): AUROC=0.9554 AUPRC=0.9590 | M5 minus CADD AUROC=+0.0322
