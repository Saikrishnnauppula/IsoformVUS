# Feature specification

Every feature below is computed per (variant, phenotype-set, cell-type) context
unless marked variant-only. Each entry documents *why* it is included — a
feature without a biological justification does not enter the model.

## Tier 1 — variant (variant-only; Baseline A / Model 1)

| Feature | Why included |
|---|---|
| `consequence` (VEP-style: missense, stop_gained, splice_*, intronic, …) | Mechanism prior: protein-truncating and splice variants have higher base rates of pathogenicity |
| `gene_symbol`, `gene_id` | Identity for grouping; NOT one-hot encoded (leakage — see docs/leakage_prevention.md) |
| `transcript_id` (MANE Select + overlapping) | Which annotation the consequence is called against |
| `exon_number`, `exon_total` | Positional context; terminal-exon truncations may escape NMD |
| `af_gnomad` (max across populations) | Population frequency: common variants are unlikely causal for rare severe disease |
| `af_popmax` | Ancestry-aware frequency; rare in all populations matters |
| `conservation` (phyloP / phastCons at position) | Conserved positions tolerate change poorly |
| `gene_constraint` (LOEUF / pLI of gene) | Constrained genes are less tolerant of damaging variation |
| `spliceai_max_delta` | Splice-disruption potential from sequence (variant-only score) |
| `cadd_phred` | Baseline deleteriousness integrating many annotations |
| `alphamissense_class` (missense only; else missing) | Structure-aware missense pathogenicity prior |

## Tier 2 — phenotype (+phenotype; Baseline B / Model 2)

| Feature | Why included |
|---|---|
| `phenotype_gene_score` | Does this gene have curated HPO association with the patient's phenotypes? (phenotype.hpoa) |
| `n_matching_hpo_terms` | Breadth of phenotype match |
| `phenotype_specificity` | Rare/specific HPO terms carry more signal than "seizure" alone (information content) |

## Tier 3 — cell type (+cell type; Baseline C / Model 3)

| Feature | Why included |
|---|---|
| `cell_relevance_score` | Measurable phenotype→cell-type relevance (see `data/phenotype_cell.py`) |
| `gene_expressed_in_cell` (binary + TPM) | A variant cannot act through a cell type that never expresses the gene |
| `cell_type_specificity_tau` | Housekeeping-expressed genes vs cell-type-restricted genes behave differently |

## Tier 4 — isoform (+isoform; IsoformVUS / Model 4)

| Feature | Why included |
|---|---|
| `isoform_usage_in_cell` (0..1 + width) | The core hypothesis variable: is the affected isoform actually used where it matters? |
| `isoform_usage_category` (one-hot) | Dominant/supported/low/uncertain — the model sees our uncertainty |
| `n_isoforms_gene` | More isoforms → more room for isoform-selective effects |
| `is_mane_select` | Whether the affected transcript is the canonical one (annotation-quality prior) |

## Tier 5 — isoform-specific effect (+functional; IsoformVUS+ / Model 5)

| Feature | Why included |
|---|---|
| `gene_effect_summary` (all/some/isoform_specific/none) | Structural classification from `data/transcripts.py` |
| `affects_dominant_isoform` (binary) | Pathogenicity requires the *relevant* isoform to be hit |
| `isoform_specific_coding` (binary) | Coding-region hit unique to a subset of isoforms |
| `isoform_specific_junction` (binary) | Splice-junction hit unique to a subset of isoforms |
| `spliceai_in_isoform_context` | SpliceAI delta re-weighted by usage of transcripts containing the affected junction |
| `nmd_escape_prediction` | Truncation in last exon / <50bp from last junction may escape NMD (weaker effect) |

## Missing-data policy

- AlphaMissense: missing for non-missense → explicit `is_missing` indicator + 0 imputation is **forbidden**; the model gets a missingness flag so "no score" is not confused with "benign score".
- Same policy for any feature unavailable for a variant class.
- Features missing >30% in train are dropped and logged; the drop is reported in results.

## The critical comparison

Model 4 vs Model 3 (and Model 5 vs Model 4) use **identical** variants, splits,
and learners — the only difference is the isoform feature tiers. Any performance
delta is attributable to isoform information, not to data or modeling choices.
