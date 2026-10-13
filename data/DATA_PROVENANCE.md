# Data provenance

Every dataset used in IsoformVUS. No analysis result is valid unless its inputs are
listed here with exact versions. Genome build is **GRCh38** unless a row says otherwise.

| # | Dataset | Source URL | Version | Download date | Genome build | Transcript annotation | License | Preprocessing |
|---|---------|-----------|---------|---------------|--------------|----------------------|---------|---------------|
| 1 | ClinVar VCF (weekly) | https://ftp.ncbi.nlm.nih.gov/pub/clinvar/vcf_GRCh38/clinvar_20261004.vcf.gz | clinvar_20261004 (released 2026-10-05) | 2026-10-06 | GRCh38 | n/a | NCBI public domain | none yet; sha256=4bff32d1…81aa12 |
| 2 | gnomAD | https://gnomad.broadinstitute.org/downloads | v4.1.0 (planned) | DEFERRED — REST API unreachable from this environment throughout 2026-10-06–07; AF enrichment not performed | GRCh38 | n/a | ODC-ODbL / CC-BY per file | not performed; disclosed in docs/feasibility.md |
| 3 | GENCODE human annotation (primary assembly) | https://ftp.ebi.ac.uk/pub/databases/gencode/Gencode_human/release_46/gencode.v46.primary_assembly.annotation.gtf.gz | v46 (2024-05-13) | 2026-10-06 | GRCh38 | GENCODE v46 | free academic/commercial | none yet |
| 4 | MANE Select + Plus Clinical | https://ftp.ncbi.nlm.nih.gov/refseq/MANE/MANE_human/current/MANE.GRCh38.v1.5.summary.txt.gz | v1.5 | 2026-10-06 | GRCh38 | MANE v1.5 | NCBI public domain | none yet; v1.5 supersedes the v1.4 plan (current/ serves v1.5 as of 2026-10-06) |
| 5 | HPO ontology + annotations | https://github.com/obophenotype/human-phenotype-ontology/releases/download/v2026-02-16/ (hp.obo, phenotype.hpoa, phenotype_to_genes.txt) | v2026-02-16 | 2026-10-06 | n/a | n/a | HPO license | column order of phenotype_to_genes.txt verified: hpo_id, hpo_name, ncbi_gene_id, gene_symbol, disease_id; legacy hpo.jax.org download page 404s — use GitHub releases |
| 6 | GTEx v8 median gene expression | https://gtexportal.org/api/v2 (medianGeneExpression; v1 API discontinued, v2 verified 2026-10-06) | v8 (54 tissues, 948 donors; GENCODE v26) | 2026-10-06 | GRCh38 | GENCODE v26 | GTEx summary data open | per-gene API queries; 102/102 panel genes retrieved; bulk GCS path guessed initially did not exist — corrected to API route |
| 7 | CELLxGENE / HCA single-cell | https://cellxgene.cziscience.com/ | DEFERRED — stage 05 uses GTEx v8 tissue-level transcript data as the measured isoform-usage source; single-cell isoform refinement identified as the key future experiment (manuscript §9) | GRCh38 | per-dataset | CC-BY per dataset | not performed |
| 8 | SpliceAI | https://github.com/Illumina/SpliceAI | v1.3.1 (pip; TF 2.21, local scoring) | 2026-10-07 | GRCh38 | **bundled GENCODE V24-era canonical annotation** (not v46 — disclosed limitation) | MIT (code) | 13,687/14,310 test variants scored; 623 unscorable (no annotated transcript window); NumPy 2.x patch to utils.py |
| 9 | CADD | https://cadd.gs.washington.edu/ (REST API v1.7) | v1.7 | 2026-10-07 | GRCh38 | n/a | academic/non-commercial | 24,600 val+test variants queried; 19,014 scored (indels/complex often unscorable) |
| 10 | AlphaMissense | https://zenodo.org/records/8208688/files/AlphaMissense_hg38.tsv.gz | v1 (2023-09-19), md5 9fd167735f16a1b87da6eb3e4c25fcb5 verified | 2026-10-06 | GRCh38 | n/a | CC-BY-NC 4.0 | missense-only; non-missense get explicit missing indicator |

## Provenance rules

1. Every `scripts/` stage that downloads data appends a row-completion entry here
   (version string, download date, md5/sha256, exact URL/file).
2. Derived files in `data/processed/` record their inputs in a sidecar `.provenance.json`.
3. If a source version changes mid-project, the change is logged here with date and reason.
4. Never silently substitute one build or annotation for another; all coordinates are
   GRCh38 and all transcripts GENCODE v46 unless this file documents an exception.
