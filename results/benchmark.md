# Held-out benchmark (test set)

Test: 14400 variants (7200 pathogenic), 24 genes. Gene-level splits; test set accessed per results/test_access.log.

| Model | Feats | Val AUROC | Test AUROC [95% CI] | Test AUPRC | Top-10 recall [95% CI] |
|---|---|---|---|---|---|
| M1 | 9 | 0.8253 | 0.7985 [0.7917,0.8061] | 0.8016 | 0.0014 [0.0014,0.0014] |
| M2 | 11 | 0.8244 | 0.7967 [0.7892,0.8051] | 0.8327 | 0.0014 [0.0014,0.0014] |
| M3 | 14 | 0.8297 | 0.7964 [0.7891,0.8045] | 0.8330 | 0.0014 [0.0014,0.0014] |
| M4 | 23 | 0.8327 | 0.8029 [0.7955,0.8105] | 0.8524 | 0.0014 [0.0014,0.0014] |
| M5 | 26 | 0.8325 | 0.8029 [0.7955,0.8106] | 0.8523 | 0.0014 [0.0014,0.0014] |

## Baselines (test)

- alphamissense_missense_only (n=1441): AUROC=0.9435 AUPRC=0.9750

## Recovery: M5 vs M1 (test)

- top10: M5 top-k pathogenic=10, M1=10, recovered missed by M1=7
- top20: M5 top-k pathogenic=20, M1=20, recovered missed by M1=10
- top50: M5 top-k pathogenic=50, M1=50, recovered missed by M1=23

## Baselines added post-hoc (enrichment completed)

- SpliceAI max delta (n=14310): AUROC=0.6844 AUPRC=0.7097 | M5 minus SpliceAI AUROC=+0.1162
- CADD PHRED (n=11240): AUROC=0.9554 AUPRC=0.9590 | M5 minus CADD AUROC=-0.2408

## Baselines added post-hoc (enrichment completed)

- SpliceAI max delta (n=13687): AUROC=0.6961 AUPRC=0.7192 | M5 minus SpliceAI AUROC=+0.1057
- CADD PHRED (n=11240): AUROC=0.9554 AUPRC=0.9590 | M5 minus CADD AUROC=-0.2408
