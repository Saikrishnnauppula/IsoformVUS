# Case-study candidacy decisions (2026-10-06)

Rule from the brief: a gene is a case study only if public evidence and
available data support it — never forced.

## TRNT1 — INCLUDED (conditional on isoform annotation check)

- **ClinVar evidence (real data, clinvar_20261004.vcf.gz):** 75 P/LP records,
  265 B/LB, 310 VUS, 18 conflicting. Label pool is sufficient.
- **Published mechanism:** Ly et al., *Molecular Cell* 2025 — nonsense/
  frameshift/missense variants between two translation-initiation sites
  selectively eliminate the longer TRNT1 protein isoform; isoform-specific
  TRNT1 mutations explain heterogeneous SIFD (sideroblastic anemia,
  immunodeficiency, fevers, developmental delay) phenotypes. Even the
  structurally benign M30L is pathogenic via alternative-start-codon
  disruption.
- **Why it fits IsoformVUS:** textbook isoform-selective alleles with a
  phenotype (SIFD) → cell-type story to build.
- **Remaining check (stage 06):** GENCODE v46 must annotate the alternative
  start-codon isoforms; if the annotation lacks them, the case study is
  dropped or reframed as an annotation-gap finding (which is itself a
  result, cf. Weykopf et al. 2026).
- **Annotation check result (2026-10-06, real GENCODE v46 primary assembly):**
  TRNT1 has 13 protein-coding transcripts with **13 distinct exon
  structures** — alternative isoforms ARE annotated. The Ly et al. mechanism
  (alternative translation-initiation sites) is a CDS-level difference;
  stage 06 will verify CDS-level distinctness per transcript. PROCEED.

## DNM1 — INCLUDED (literature-anchored)

- **MANE evidence (real data, MANE v1.5):** DNM1 has both a MANE Select
  (ENST00000372923.8) and a **MANE Plus Clinical** transcript
  (ENST00000486160.4) — the annotation itself acknowledges a clinically
  relevant alternative isoform.
- **Published mechanism:** AlphaGenome Atlas/AVI (2026 preprint) rescued a
  real epileptic-encephalopathy VUS in DNM1 intron 10 affecting a
  brain-specific isoform (exon 10a, absent from MANE Select).
- **IsoformVUS angle:** the counterfactual experiment (stage 13) can show
  the same variant scored under MANE Select vs MANE Plus Clinical —
  exactly the hypothesis test.

## Other candidates — PENDING benchmark

Additional case studies will be selected from stage-11 top-k recovered
variants (pathogenic variants IsoformVUS ranks highly that baselines miss),
not hand-picked.
