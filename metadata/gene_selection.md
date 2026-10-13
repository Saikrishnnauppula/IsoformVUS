# Gene selection criteria

Target: **~100 Mendelian disease genes** (acceptable range 50–200; 100 chosen to
balance statistical power against SpliceAI/CADD compute and single-cell data
wrangling). Final list: `metadata/gene_panel.tsv` [PENDING — built by
`scripts/04_select_genes.py`].

## Inclusion criteria (all must hold)

1. **Mendelian disease association.** Gene appears in OMIM with a monogenic
   phenotype AND has HPO gene→phenotype annotations (`phenotype.hpoa`), so the
   phenotype→cell-type chain has measurable input.
2. **Sufficient ClinVar pathogenic variants.** ≥10 Pathogenic/Likely-pathogenic
   submissions with review status ≥ 1 star (criteria provided), GRCh38-mappable.
   Guarantees the gene contributes meaningful positive labels.
3. **Benign controls available.** ≥10 Benign/Likely-benign ClinVar variants OR
   ≥20 gnomAD variants with AF > 1% in the gene region (documented per gene;
   common ≠ benign is respected — these are *controls for the classifier*, and
   the AF>1% set is labeled "population-common", analyzed separately from
   ClinVar B/LB).
4. **Multiple isoforms.** ≥2 protein-coding GENCODE v46 transcripts with distinct
   exon structures (not just UTR differences), so isoform-selectivity is possible.
5. **Expression data available.** Gene detected (median TPM > 1) in ≥1 GTEx v8
   tissue relevant to its HPO phenotypes.
6. **Cell-type data available.** Gene detected in ≥1 CELLxGENE/HCA dataset covering
   a phenotype-relevant tissue, with Cell Ontology cell-type labels.
7. **Biologically meaningful alternative isoforms where possible.** Preference
   (not a hard filter) for genes where alternative exons change coding sequence
   (cassette exons in CDS), recorded per gene in the panel table.

## Exclusion criteria

- Genes with a single annotated protein-coding transcript (no isoform contrast).
- Genes whose only disease association is somatic/cancer (germline rare-disease
  focus; different biology and ascertainment).
- Mitochondrial genes (different inheritance, annotation, and frequency logic).
- Genes failing GRCh38 liftover consistency (ClinVar GRCh38 VCF is primary; any
  gene with >5% unmappable variants is dropped and logged).

## Diversity guardrails

- Disease categories are derived from the HPO ontology itself
  (`src/isoformvus/data/hpo.py`): each gene's HPO terms are walked up to the
  direct children of HP:0000118 (Phenotypic abnormality); the most frequent
  child names the category (e.g. `nervous_system`, `cardiovascular_system`).
  This replaced an earlier hand-written term→category hint map that
  mislabeled genes (e.g. APC as cardiovascular). The ontology vote is still
  a heuristic for pleiotropic genes and is documented as such.
- Cap any single disease category at 40% of the panel so one biology does
  not dominate the benchmark.
- **Per-gene variant caps** (enforced in stage 03): TTN contributes
  5,946 P / 17,567 B variants — without caps the benchmark would be a
  TTN benchmark. Cap at 300 variants per class per gene (seeded); the cap
  value and seed are recorded in the labels provenance.
- Include ≥10 genes with a *known* isoform-selective pathogenic mechanism from
  literature (positive-control genes for the counterfactual experiment), flagged
  in the panel table — these also seed the case studies.

## What "sufficient" means statistically

With ~100 genes × ~20–60 labeled variants ≈ 2,000–6,000 variants, gene-level
splits leave ~20 genes (400–1,200 variants) in held-out test — enough for AUROC
with bootstrap CIs but not for rare-subgroup claims beyond the isoform-selective
subset. Subgroup analyses are reported with CIs and labeled exploratory.
