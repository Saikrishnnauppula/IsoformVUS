# Data

Every input to IsoformVUS, with version, source, and licence. Genome build
is GRCh38 for all rows. `data/DATA_PROVENANCE.md` records download dates
and checksums; the table below is the readable copy of that record, with
the corrections the review found.

## Inputs

| Resource | Version | Source | Licence |
|---|---|---|---|
| ClinVar VCF, weekly | `clinvar_20261004.vcf.gz`, released 2026-10-05, downloaded 2026-10-06 | https://ftp.ncbi.nlm.nih.gov/pub/clinvar/vcf_GRCh38/clinvar_20261004.vcf.gz | NCBI public domain |
| ClinVar variant_summary | Monthly tab-delimited release, first Thursday | https://ftp.ncbi.nlm.nih.gov/pub/clinvar/tab_delimited/variant_summary.txt.gz | NCBI public domain |
| GENCODE human annotation | v46, released 2024-05-13, primary assembly GTF | https://ftp.ebi.ac.uk/pub/databases/gencode/Gencode_human/release_46/gencode.v46.primary_assembly.annotation.gtf.gz | Free for academic and commercial use |
| MANE Select + Plus Clinical | v1.5 summary, downloaded 2026-10-06 | https://ftp.ncbi.nlm.nih.gov/refseq/MANE/MANE_human/current/MANE.GRCh38.v1.5.summary.txt.gz | NCBI public domain |
| HPO ontology and annotations | v2026-02-16 (`hp.obo`, `phenotype.hpoa`, `phenotype_to_genes.txt`) | https://github.com/obophenotype/human-phenotype-ontology/releases/download/v2026-02-16/ | HPO licence |
| GTEx | v8, 54 tissues, 948 donors, quantified on GENCODE v26 | https://gtexportal.org/api/v2 (median transcript expression) | GTEx summary data open |
| SpliceAI | v1.3.1, local scoring, TensorFlow 2.21 | https://github.com/Illumina/SpliceAI | MIT (code); models CC BY-NC 4.0 |
| CADD | v1.7, REST API | https://cadd.gs.washington.edu/ | Academic / non-commercial |
| AlphaMissense | v1, 2023-09-19, Zenodo record 8208688, md5 `9fd167735f16a1b87da6eb3e4c25fcb5` | https://zenodo.org/records/8208688/files/AlphaMissense_hg38.tsv.gz | CC BY-NC 4.0 |
| gnomAD | v4.1.1 planned (`gnomad_r4`), never run | https://gnomad.broadinstitute.org/api | ODC-ODbL / CC-BY per file |

Three notes on the table. First, GENCODE v44 appears in the manuscript
header as "v46 / v44" and in no other analysis: a v44 GTF was used only for
an annotation-fragility experiment, and the benchmark itself is v46. Second,
SpliceAI ran against its bundled GENCODE V24-era canonical annotation, not
v46; the scores in the shipped tables carry that mismatch. Third, gnomAD
was deferred because the REST API was unreachable from the build
environment throughout 2026-10-06 and 2026-10-07. The script
`scripts/03b_gnomad_af.py` exists and targets v4.1.1, but no
allele-frequency column exists in the feature matrix.

## What ships in the repository

`data/processed/` holds the small derived tables: `feature_matrix.tsv.gz`
(60,521 rows, 29 columns), `variant_labels.tsv.gz`, `variant_cadd.tsv.gz`,
`variant_spliceai.tsv.gz`, `variant_alphamissense.tsv.gz`,
`variant_isoform_effects.tsv.gz`, `isoform_usage.tsv.gz`,
`gtex_median_tpm.tsv.gz`, `cell_relevance.tsv.gz`, and
`feature_tiers.json`, plus `THIRD_PARTY_DATA_LICENSES.md`. `results/`
holds `benchmark.json`, `statistics.json`, `per_model_cis.json`,
`splits.json`, `test_access.log`, and the stage-18 outputs.

## What must be re-downloaded

`data/raw/` (~5.1 GB) and `data/interim/` (~3.1 GB) are not in the
repository. A rebuild from raw inputs needs the ClinVar VCF, the GENCODE
v46 GTF, the MANE v1.5 summary, HPO v2026-02-16, the GTEx v8 API, the
SpliceAI package with a GRCh38 FASTA and GTF, the CADD v1.7 REST API, and
the AlphaMissense Zenodo file, at the URLs in the table above.

Two rebuilds from raw inputs will not reproduce the locked results exactly.
The pinned ClinVar weekly `clinvar_20261004.vcf.gz` now returns HTTP 404
from NCBI FTP; weeklies rotate and only newer snapshots remain, so a
rebuild today uses a different ClinVar version. And the GTEx v8 to GENCODE
v46 join loses an unquantified number of transcripts: GTEx v8 used GENCODE
v26 identifiers, the join used the ENST base only, and the
unmatched-transcript counts sat in `data/interim/`, which is not shipped.

## The CADD situation

The CADD scores exist and the feature matrix now carries them.
`variant_cadd.tsv.gz` holds real v1.7 scores for 24,600 validation and
test variants, with 100% variant-key overlap against the matrix. The
rebuilt matrix merges those scores in (19,014 of 60,521 rows scored), and
the benchmark was re-run from it. Every train row still carries
`cadd_phred = -1.0` with `cadd_missing = 1`: CADD was queried only for
validation and test variants, so the fitted CADD coefficient was never fit
against varying CADD values. The benchmark's CADD baseline (test AUROC
0.9554, n = 11,240) was scored from the CADD table, not from the matrix.

CADD scoring itself was partial: 24,600 val+test variants were queried and
19,014 scored. Indels and complex variants were often unscorable, and
scoring is skewed by label: CADD scores 92.96% of benign but 63.15% of
pathogenic test variants, and 0% of deletions. SpliceAI scoring was also
partial: 13,687 of 14,310 test variants scored, 623 unscorable for lack of
an annotated transcript window.

## Script-18 inputs (not shipped)

Stage 18 (`scripts/18_exon_architecture_nmd.py`) needs five large files
that live in `data/interim/script18/` and are not committed:

| File | Size | Provenance |
|---|---|---|
| `supp_clinvar_all_table.tsv` | 67 MB | The 587,227-row ClinVar table with dual GENCODE v44/v46 consequences, from the project archive |
| `variant_molecular_consequences.tsv.gz` | 512 KB | Per-variant molecular consequences from the experiment archive |
| `ptc_features.tsv.gz` | 300 KB | Premature-termination-codon features from the experiment archive |
| `variant_transcript_sev.tsv.gz` | 2.6 MB | Per-transcript severity scores from the experiment archive |
| `labels_2024_vs_2026.tsv.gz` | 1.5 MB | 2024 vs 2026 ClinVar labels for the temporal holdout, from the experiment archive |

The code that built these five files was removed with the experiment
archive (`breakthrough/`, deleted 2026-11-08), so they cannot currently
be regenerated from the repository. A fresh clone cannot run stage 18
until they are restored to `data/interim/script18/`.

## Licences

The MIT licence in `LICENSE` covers the code only. The processed tables
mix sources with different terms, and `data/processed/
THIRD_PARTY_DATA_LICENSES.md` says so file by file. The binding cases are
AlphaMissense at CC BY-NC 4.0 (`variant_alphamissense.tsv.gz` must not be
represented as MIT-licensed or used commercially), CADD at academic /
non-commercial, and SpliceAI models at CC BY-NC 4.0. ClinVar, MANE, and
GENCODE carry no such restriction. Check each source's current terms
before redistributing or reusing any processed table.
