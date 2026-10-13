# Biology in plain English

Written for a Computer Science + ML reader with limited formal biology coursework.
Everything here is prerequisite vocabulary for the IsoformVUS manuscript.

## The core vocabulary

**Gene.** A stretch of DNA that carries instructions to build a functional product,
usually a protein. Think of it as one function in a codebase. Humans have ~20,000
protein-coding genes.

**Variant.** A difference in DNA at a specific position relative to a reference genome
— like a diff against `main`. A single-letter change (SNV), insertion, or deletion.
Everyone carries millions of variants; almost all are harmless.

**Pathogenic variant.** A variant with credible evidence that it *causes* or strongly
contributes to disease. "Pathogenic" is an evidence-based label, not a synonym for "rare".

**VUS (Variant of Uncertain Significance).** A variant found in a patient where the
evidence is currently insufficient to call it pathogenic *or* benign. VUSs are the
central clinical pain point: sequencing finds them by the hundreds, and each one is
a question mark. IsoformVUS exists to *prioritize* which VUSs deserve follow-up — it
does not diagnose.

**Transcript.** The RNA copy made from a gene. DNA → (transcription) → RNA →
(translation) → protein. A gene's DNA is the source file; a transcript is one
compiled artifact built from it.

**Isoform.** One of several *alternative* transcripts a single gene can produce.
Through **splicing** (below), the cell can include or skip different exons, so one
gene yields multiple distinct proteins. ~95% of human multi-exon genes are
alternatively spliced. If a gene is a function, isoforms are its overloaded
signatures — same name, different behavior.

**Cell type.** Cells with a shared identity and job: neurons, cardiomyocytes,
hepatocytes, T cells, etc. The key fact for this project: **different cell types
prefer different isoforms of the same gene**. A variant sitting in an exon that the
brain's isoform includes but the liver's isoform skips will matter in brain cells
and be silent in liver cells — same DNA, different consequence.

**Splicing.** After transcription, the raw RNA (pre-mRNA) contains coding chunks
(**exons**) and filler (**introns**). The spliceosome cuts out introns and joins
exons at **splice junctions**. Alternative splicing = choosing different exon
combinations. A variant near a splice junction can reroute this choice — include
an exon that should be skipped, or skip one that should be kept.

**ClinVar.** NCBI's public archive of variant ↔ disease interpretations submitted by
labs worldwide. Each record: variant, condition, classification
(Pathogenic / Likely pathogenic / VUS / Likely benign / Benign), evidence summary.
It is our *label source* — and a biased one (well-studied genes are overrepresented;
see docs/limitations.md).

**gnomAD.** The Genome Aggregation Database: allele frequencies from ~800k
exomes/genomes. Tells you how common a variant is in the general population. Rare
in gnomAD ≠ pathogenic (most rare variants are benign); common in gnomAD ≈ almost
certainly not causal for a rare severe disease. We use it for frequency features
and for benign controls — never as a label by itself.

**GTEx.** Genotype-Tissue Expression: bulk RNA-seq across ~50 human tissues from
donors. Answers "is this gene expressed in heart?" at *tissue* resolution — useful
as a prior, but tissues are mixtures of cell types, so it cannot resolve isoform
usage per cell type.

**CELLxGENE.** A public atlas of single-cell RNA-seq datasets (Human Cell Atlas
data and more). Answers "which *cell types* express this gene?" — the resolution
IsoformVUS needs. Limitation: standard short-read scRNA-seq is 3′-biased and
cannot reliably distinguish full-length isoforms; we treat isoform assignments
from it as *supported/uncertain*, not ground truth.

**HPO (Human Phenotype Ontology).** A controlled vocabulary of human phenotypic
abnormalities ("Seizure" = HP:0001250, "Cardiomyopathy" = HP:0001638), organized
as a hierarchy. It lets us turn a patient's clinical description into computable
terms, then link phenotypes → affected tissues → cell types through curated
annotations rather than hand-waving.

**SpliceAI.** A deep-learning model that predicts, from DNA sequence alone, whether
a variant disrupts splicing (creates/breaks donor/acceptor sites). Outputs delta
scores 0–1 per variant. Strong tool — but it is *transcript-agnostic* in practice:
it scores the sequence change, not which isoform the relevant cell actually uses.
IsoformVUS puts its score in isoform context.

**CADD.** Combined Annotation Dependent Depletion: a single deleteriousness score
per variant integrating conservation, functional annotations, etc. The field's
default "how bad is this variant" number — gene-level, tissue-agnostic. Our
baseline to beat.

---

## The IsoformVUS concept, by example

Meet **a hypothetical patient**: a child with an undiagnosed neurodevelopmental
disorder — seizures (HP:0001250) and intellectual disability (HP:0001249).
Exome sequencing finds a rare heterozygous variant in gene **GENE-X**, currently a VUS.

**What conventional tools say.** CADD gives the variant a middling score (say 18 —
below the usual "look here" thresholds). SpliceAI predicts a moderate splice effect
(0.4). AlphaMissense calls the amino-acid change "ambiguous". Gene-level reasoning
says: GENE-X is expressed in many tissues; nothing stands out. The variant stays a
VUS, and the diagnostic odyssey continues.

**What IsoformVUS adds, step by step.**

1. **Phenotype → cell type.** HPO terms (seizure, intellectual disability) point to
   the brain. Measurable evidence — HPO tissue annotations + single-cell atlases —
   narrows this to specific cell types, e.g. cortical excitatory neurons, with a
   *cell-relevance score* rather than a hard-coded guess.

2. **Cell type → isoform.** In cortical neurons, CELLxGENE/long-read data show
   GENE-X predominantly uses **isoform B**, which includes exon 7. In liver and
   blood, the dominant isoform **A** skips exon 7. This is a *measured* usage
   difference, with uncertainty recorded.

3. **Isoform-specific variant effect.** The patient's variant sits inside exon 7
   and creates a premature stop codon *in isoform B*. In isoform A, exon 7 doesn't
   exist — the variant is intronic there and harmless. So the variant is
   **isoform-selective**: devastating in the brain's isoform, silent elsewhere.

4. **Prioritization.** IsoformVUS combines the variant features (consequence,
   frequency, conservation, SpliceAI/CADD/AlphaMissense) with phenotype relevance,
   cell-type relevance, isoform usage in the relevant cell type, and the
   isoform-specific consequence. The same variant that looked mediocre in a
   tissue-agnostic view now ranks at the top — because the model knows *which
   cells matter* and *which isoform those cells use*.

**The counterfactual that makes the point.** Ask the model: "score this variant
assuming isoform A is the relevant one" vs "assuming isoform B". Same DNA, two
different biological interpretations. When the answer flips with the isoform, the
isoform information was load-bearing — that's the experiment in §15 of the brief,
and the honest test of the hypothesis.

**What IsoformVUS does NOT do.** It does not prove the variant causes the disease,
does not replace functional validation, and does not output a clinical diagnosis.
It produces a ranked, explainable shortlist: *"start here, and here's why."*
