# Feasibility assessment

Honest accounting of what IsoformVUS can and cannot do with public data,
decided *before* large-scale compute. Updated as evidence arrives.

## Data availability verdicts

| Dependency | Verdict | Basis |
|---|---|---|
| ClinVar labels (P/LP vs B/LB) | FEASIBLE | Public VCF + variant_summary; GRCh38; well-trodden |
| gnomAD allele frequencies | PARTIAL (2026-10-06: the gnomAD GraphQL API times out from this
  environment; AF enrichment deferred — see below) | Public browser/downloads v4.1; AF fields documented |

## 2026-10-06 update: gnomAD API unreachable

Direct tests (Python urllib + curl, 60s timeouts) against
`https://gnomad.broadinstitute.org/api` timed out with zero bytes received.
The AF enrichment script (`scripts/03b_gnomad_af.py`) is implemented and
documented but cannot run from here. Consequences:

- `af_gnomad` / `af_popmax` features are marked PENDING in the feature
  matrix; the missing-data policy (explicit missingness indicators, no
  zero-imputation) covers them.
- gnomAD-common controls are deferred; the benchmark uses ClinVar B/LB
  controls only, and this is disclosed in results.
- The AF feature can be filled in later from any network where the API or
  bulk VCFs are reachable; the script is ready.

This is an infrastructure limitation, not a data-availability one, and it
is recorded rather than worked around with fabricated frequencies.
| GENCODE v46 transcripts | FEASIBLE | Public GTF; exon structures parseable |
| MANE Select mapping | FEASIBLE | Public TSV; one canonical transcript per gene |
| HPO ontology + phenotype annotations | FEASIBLE | Public OBO/JSON; phenotype.hpoa links genes↔phenotypes |
| GTEx bulk tissue expression | FEASIBLE (summary) | Open summary matrices; gene-level TPM |
| Single-cell expression (CELLxGENE/HCA) | FEASIBLE with caveats | Public H5ADs; cell-type labels vary by dataset; 3′ bias limits isoform resolution |
| Long-read full-length isoform evidence | PARTIAL | Public PacBio/ONT datasets exist (e.g. ENCODE LR-RNA-seq, GTEx LR pilot) but coverage is tissue-limited, not cell-type-resolved at scale |
| SpliceAI scores | FEASIBLE | Open model v1.3.1; runnable on GRCh38 locally |
| CADD scores | FEASIBLE | Precomputed v1.7 GRCh38 TSVs downloadable |
| AlphaMissense scores | FEASIBLE (missense only) | Public per-protein TSVs; comparisons restricted to missense |

## Key risks and mitigations
1. **Isoform usage from short-read scRNA-seq is unreliable.** 3′/5′ bias and
   ambiguous multi-mapping mean per-cell-type isoform *quantification* is noisy.
   Mitigation: use isoform usage as *categorical evidence* (dominant / supported /
   low / uncertain) aggregated across datasets, never as precise PSI; record
   uncertainty; long-read datasets used as supporting evidence where available.
2. **Cell-type labels are dataset-specific.** "Excitatory neuron" in one atlas ≠
   another's. Mitigation: map to Cell Ontology terms; keep dataset-of-origin in
   provenance; cell-relevance score aggregates across datasets.
3. **ClinVar label bias.** Well-studied genes/diseases dominate; VUSs are not a
   random sample. Mitigation: gene-level splits, disease-stratified sampling,
   explicit bias discussion in limitations; never train on VUS as a label.
4. **Phenotype→cell-type mapping is the weakest causal link.** HPO annotations
   are gene→phenotype, not cell→phenotype. Mitigation: measurable evidence chain
   (HPO → affected system → tissue expression (GTEx) → cell-type expression
   (CELLxGENE)); score, don't hard-code; ablate this layer (Model 3 vs Model 2)
   to test whether it adds anything.
5. **Compute.** SpliceAI on ~10k variants is fine on CPU/GPU-single; full-genome
   scoring is out of scope. Keep the panel at ~100 genes.

## Scope decision (recorded 2026-10-06)

Proceed with the ~100-gene pilot. If the ablation shows isoform features add no
signal on the pilot, the honest outcome is a *negative result with a diagnosed
cause* (data resolution? biology? modeling?) — not scope creep into more genes.

## CADD v1.7 REST API coverage (verified 2026-10-07)
- val_test scope (24,600 variants): 19,014 scored (77.3%), 5,586 null (22.7%).
- Nulls are SYSTEMATIC, not transient (retry of 1,031 non-deletion nulls recovered 0):
  - All 4,555 deletions return `[]` from the API regardless of query format.
  - Remaining nulls are complex/MNVs (e.g. AG>GA) and longer insertions.
- Conclusion: the GRCh38-v1.7 REST endpoint reliably scores SNVs + short simple
  indels only. CADD encoded as missing for deletions/complex (vartype features
  already capture frameshift signal; ablation will quantify CADD's contribution).
- Cache: data/interim/cadd_cache.jsonl (24,600 unique keys, 0 dups);
  merged table: data/processed/variant_cadd.tsv.gz.
