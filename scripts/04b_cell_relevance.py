#!/usr/bin/env python3
"""Stage 04b: phenotype -> cell-type relevance scores (real implementation).

For each panel gene:
    HPO terms -> top-level systems (ontology, measured from annotations)
      -> GTEx tissues (curated anatomical map, tier='curated')
      -> median TPM per tissue (measured via GTEx v8 API, tier='expression')
      -> representative cell types (tier='inferred')

cell_relevance(gene, cell_type) =
    w1 * phenotype_tissue_match + w2 * expression_support

- phenotype_tissue_match: fraction of the gene's HPO top-level systems that
  map to the cell type's tissue (0..1).
- expression_support: log1p(median TPM) of the gene in the tissue, scaled to
  0..1 across tissues for that gene (relative expression profile).

Every edge carries its evidence tier. Cell-type assignments are 'inferred'
until single-cell data (stage 05) refines them — this is disclosed, not
hidden.

Writes data/processed/cell_relevance.tsv.gz:
    gene_symbol, cell_ontology_id, cell_type_name, tissue,
    score, phenotype_match, expression_support, evidence_tier
"""

from __future__ import annotations

import argparse
import gzip
import json
import math
import sys
from collections import defaultdict
from pathlib import Path

REPO = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(REPO / "src"))

from isoformvus.data.hpo import categorize_gene, parse_obo  # noqa: E402
from isoformvus.data.phenotype_cell import combine_evidence  # noqa: E402

# Curated anatomical map: HPO top-level system -> GTEx tissueSiteDetailIds.
# Tier: 'curated' (anatomical knowledge, documented here — NOT hard-coded
# disease->tissue guesses; the *gene-specific* part comes from measured data).
SYSTEM_TO_TISSUES: dict[str, list[str]] = {
    "nervous_system": ["Brain_Amygdala", "Brain_Anterior_cingulate_cortex_BA24",
                       "Brain_Caudate_basal_ganglia", "Brain_Cerebellar_Hemisphere",
                       "Brain_Cerebellum", "Brain_Cortex",
                       "Brain_Frontal_Cortex_BA9", "Brain_Hippocampus",
                       "Brain_Hypothalamus", "Brain_Nucleus_accumbens_basal_ganglia",
                       "Brain_Putamen_basal_ganglia", "Brain_Spinal_cord_cervical_c-1",
                       "Brain_Substantia_nigra", "Nerve_Tibial", "Pituitary"],
    "cardiovascular_system": ["Heart_Atrial_Appendage", "Heart_Left_Ventricle",
                              "Artery_Aorta", "Artery_Coronary", "Artery_Tibial"],
    "musculoskeletal_system": ["Muscle_Skeletal"],
    "digestive_system": ["Liver", "Pancreas", "Stomach", "Colon_Sigmoid",
                         "Colon_Transverse", "Small_Intestine_Terminal_Ileum",
                         "Esophagus_Gastroesophageal_Junction",
                         "Esophagus_Mucosa", "Esophagus_Muscularis",
                         "Adrenal_Gland"],
    "genitourinary_system": ["Kidney_Cortex", "Kidney_Medulla", "Bladder",
                             "Prostate", "Testis", "Ovary", "Uterus", "Vagina",
                             "Breast_Mammary_Tissue"],
    "respiratory_system": ["Lung"],
    "integument": ["Skin_Not_Sun_Exposed_Suprapubic", "Skin_Sun_Exposed_Lower_leg"],
    "blood_and_blood-forming_tissues": ["Whole_Blood", "Spleen"],
    "immune_system": ["Spleen", "Whole_Blood"],
    "head_or_neck": ["Thyroid", "Pituitary", "Minor_Salivary_Gland"],
    "eye": [],          # no GTEx eye tissue: expression_support unavailable
    "ear": [],          # no GTEx ear tissue
    "limbs": ["Muscle_Skeletal", "Nerve_Tibial", "Artery_Tibial"],
    "neoplasm": [],     # pan-tissue: handled as max-expression below
}

# Representative cell types per tissue (Cell Ontology ids where standard).
# Tier: 'inferred' — single-cell refinement is stage 05.
TISSUE_TO_CELLS: dict[str, list[tuple[str, str]]] = {
    "Brain_Cortex": [("CL:0000540", "neuron"), ("CL:0000127", "astrocyte"),
                     ("CL:0000236", "oligodendrocyte"), ("CL:0000129", "microglial cell")],
    "Brain_Frontal_Cortex_BA9": [("CL:0000540", "neuron"), ("CL:0000127", "astrocyte")],
    "Brain_Cerebellum": [("CL:0000540", "neuron"), ("CL:0000127", "astrocyte")],
    "Heart_Left_Ventricle": [("CL:0000746", "cardiac muscle cell"),
                             ("CL:0000115", "endothelial cell")],
    "Heart_Atrial_Appendage": [("CL:0000746", "cardiac muscle cell")],
    "Muscle_Skeletal": [("CL:0000188", "skeletal muscle cell")],
    "Liver": [("CL:0000182", "hepatocyte")],
    "Kidney_Cortex": [("CL:0000653", "renal tubular epithelial cell")],
    "Lung": [("CL:0002063", "pneumocyte"), ("CL:0000236", "ciliated cell")],
    "Whole_Blood": [("CL:0000084", "T cell"), ("CL:0000236", "B cell"),
                    ("CL:0000235", "monocyte")],
    "Skin_Sun_Exposed_Lower_leg": [("CL:0000312", "keratinocyte")],
    "Pancreas": [("CL:0000171", "pancreatic acinar cell")],
    "Testis": [("CL:0000020", "spermatocyte")],
    "Thyroid": [("CL:0000402", "thyroid follicular cell")],
}
DEFAULT_CELLS = [("CL:0000000", "cell (tissue-level, unspecified)")]


def load_gene_hpo(path: Path) -> dict[str, set[str]]:
    gene_hpo: dict[str, set[str]] = defaultdict(set)
    with open(path) as fh:
        for line in fh:
            if line.startswith("#") or not line.strip():
                continue
            parts = line.rstrip("\n").split("\t")
            if len(parts) < 4 or parts[0] == "hpo_id":
                continue
            gene_hpo[parts[3]].add(parts[0])
    return dict(gene_hpo)


def load_tpm(path: Path) -> dict[str, dict[str, float]]:
    tpm: dict[str, dict[str, float]] = defaultdict(dict)
    with gzip.open(path, "rt") as fh:
        fh.readline()
        for line in fh:
            gene, _ensg, tissue, val = line.rstrip("\n").split("\t")
            tpm[gene][tissue] = float(val)
    return dict(tpm)


def main() -> int:
    ap = argparse.ArgumentParser()
    ap.add_argument("--panel", default=str(REPO / "metadata/gene_panel.tsv"))
    ap.add_argument("--tpm", default=str(REPO / "data/processed/gtex_median_tpm.tsv.gz"))
    ap.add_argument("--out", default=str(REPO / "data/processed/cell_relevance.tsv.gz"))
    ap.add_argument("--w1", type=float, default=0.6)
    ap.add_argument("--w2", type=float, default=0.4)
    args = ap.parse_args()

    raw = REPO / "data" / "raw"
    if not Path(args.tpm).exists():
        print(f"[04b] missing {args.tpm}; run scripts/04a_gtex_expression.py first",
              file=sys.stderr)
        return 1

    names, parents = parse_obo(raw / "hp.obo")
    gene_hpo = load_gene_hpo(raw / "phenotype_to_genes.txt")
    tpm = load_tpm(Path(args.tpm))
    panel = [l.split("\t")[0] for l in open(args.panel).read().splitlines()[1:]]

    out = Path(args.out)
    n_rows = 0
    with gzip.open(out, "wt") as fh:
        fh.write("\t".join(["gene_symbol", "cell_ontology_id", "cell_type_name",
                             "tissue", "score", "phenotype_match",
                             "expression_support", "evidence_tier"]) + "\n")
        for gene in panel:
            hpos = gene_hpo.get(gene, set())
            # Gene's HPO top-level systems (multi-system genes keep all).
            systems: dict[str, int] = defaultdict(int)
            for hpo in hpos:
                cat = categorize_gene({hpo}, parents, names)
                systems[cat] += 1
            gene_tpm = tpm.get(gene, {})
            # Relative expression profile for this gene (0..1 across tissues).
            logtpm = {t: math.log1p(v) for t, v in gene_tpm.items()}
            mx = max(logtpm.values()) if logtpm else 1.0

            # Candidate tissues: union of mapped tissues, or all tissues for
            # neoplasm/unmapped genes (max-expression fallback).
            mapped_tissues: set[str] = set()
            for s in systems:
                mapped_tissues.update(SYSTEM_TO_TISSUES.get(s, []))
            use_fallback = not mapped_tissues
            tissues = (set(gene_tpm) if use_fallback else mapped_tissues)

            for tissue in sorted(tissues):
                expr_sup = logtpm.get(tissue, 0.0) / mx if mx else 0.0
                # Phenotype match: fraction of gene's systems mapping to tissue.
                n_sys_match = sum(1 for s in systems
                                  if tissue in SYSTEM_TO_TISSUES.get(s, []))
                ph_match = (n_sys_match / len(systems)) if systems else 0.0
                if use_fallback:
                    ph_match = 0.0  # no phenotype link: expression only
                score = combine_evidence(ph_match, expr_sup, w1=args.w1, w2=args.w2)
                tier = "curated+expression" if not use_fallback else "expression_only"
                cells = TISSUE_TO_CELLS.get(tissue, DEFAULT_CELLS)
                for cl_id, cl_name in cells:
                    cell_tier = tier + "+inferred_celltype"
                    fh.write("\t".join([gene, cl_id, cl_name, tissue,
                                         f"{score:.4f}", f"{ph_match:.4f}",
                                         f"{expr_sup:.4f}", cell_tier]) + "\n")
                    n_rows += 1
    print(f"[04b] wrote {out} ({n_rows} gene-cell rows)")
    print("[04b] evidence tiers: curated (HPO->tissue map) + expression (GTEx TPM); "
          "cell-type assignments are 'inferred' pending single-cell data (stage 05)")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
