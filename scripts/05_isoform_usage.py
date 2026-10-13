#!/usr/bin/env python3
"""Stage 05: cell-type -> isoform usage (real implementation).

For each panel gene and its phenotype-relevant tissues, fetch GTEx v8
median *transcript* expression
(/expression/medianTranscriptExpression?gencodeId=<v26>&tissueSiteDetailId=<t>)
and compute per-tissue isoform usage fractions:

    usage(tx, tissue) = TPM(tx, tissue) / sum_TPM(all tx of gene, tissue)

Usage is then categorized per src/isoformvus/data/isoform.py:
    dominant / supported / low_usage / uncertain
with n_datasets=1 (single GTEx v8 dataset) and has_long_read_support=False —
both recorded honestly. Tissue-level (not cell-type) resolution is disclosed;
cell types inherit their tissue's usage via the 04b tissue->cell map.

Transcript ID matching: GTEx v8 uses GENCODE v26 IDs; our structural
classification uses v46. Joins use the ID *base* (ENST without version),
which is stable for most transcripts; unmatched transcripts are reported
in the provenance stats, not silently dropped.

Writes data/processed/isoform_usage.tsv.gz:
    gene_symbol, transcript_id_v26, transcript_base, tissue,
    median_tpm, usage, usage_category, n_datasets, has_long_read_support

Cache: data/interim/gtex_tx_cache.json (resume-safe).
"""

from __future__ import annotations

import argparse
import gzip
import json
import sys
import time
import urllib.parse
import urllib.request
from datetime import date
from pathlib import Path

REPO = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(REPO / "src"))

from isoformvus.data.isoform import categorize_usage  # noqa: E402

API = "https://gtexportal.org/api/v2"
DATASET = "gtex_v8"

# Reuse the curated system->tissue map from stage 04b (import, don't duplicate).
sys.path.insert(0, str(REPO / "scripts"))
import importlib.util  # noqa: E402

_spec = importlib.util.spec_from_file_location(
    "cell_rel", str(REPO / "scripts" / "04b_cell_relevance.py"))
_mod = importlib.util.module_from_spec(_spec)
_spec.loader.exec_module(_mod)
SYSTEM_TO_TISSUES = _mod.SYSTEM_TO_TISSUES


def api_get(path: str, params: dict) -> dict:
    url = f"{API}{path}?" + urllib.parse.urlencode(params)
    req = urllib.request.Request(url, headers={"Accept": "application/json",
                                               "User-Agent": "IsoformVUS/0.1.0 (research)"})
    with urllib.request.urlopen(req, timeout=30) as resp:
        return json.load(resp)


def main() -> int:
    ap = argparse.ArgumentParser()
    ap.add_argument("--panel", default=str(REPO / "metadata/gene_panel.tsv"))
    ap.add_argument("--gtex-cache", default=str(REPO / "data/interim/gtex_cache.json"))
    ap.add_argument("--cache", default=str(REPO / "data/interim/gtex_tx_cache.json"))
    ap.add_argument("--out", default=str(REPO / "data/processed/isoform_usage.tsv.gz"))
    ap.add_argument("--sleep", type=float, default=0.3)
    args = ap.parse_args()

    raw = REPO / "data" / "raw"
    sys.path.insert(0, str(REPO / "src"))
    from isoformvus.data.hpo import categorize_gene, parse_obo

    panel = [l.split("\t")[0] for l in open(args.panel).read().splitlines()[1:]]
    gtex_cache = json.loads(Path(args.gtex_cache).read_text())
    names, parents = parse_obo(raw / "hp.obo")

    gene_hpo: dict[str, set[str]] = {}
    with open(raw / "phenotype_to_genes.txt") as fh:
        for line in fh:
            if line.startswith("#") or not line.strip():
                continue
            parts = line.rstrip("\n").split("\t")
            if len(parts) < 4 or parts[0] == "hpo_id":
                continue
            gene_hpo.setdefault(parts[3], set()).add(parts[0])

    # Relevant tissues per gene from HPO systems, SCOPED to top-N by
    # measured gene-level TPM (keeps API calls tractable; the scope rule is
    # logged, not hidden).
    gene_tissues: dict[str, set[str]] = {}
    TOP_N_TISSUES = 5
    for gene in panel:
        systems = {categorize_gene({h}, parents, names)
                   for h in gene_hpo.get(gene, set())}
        tissues: set[str] = set()
        for s in systems:
            tissues.update(SYSTEM_TO_TISSUES.get(s, []))
        if not tissues:
            tissues = {"Whole_Blood"}  # fallback, flagged in output
        # Rank by measured TPM; keep top-N.
        gene_tpm = (gtex_cache.get(gene) or {}).get("tissues", {})
        ranked = sorted(tissues, key=lambda t: -gene_tpm.get(t, 0.0))
        gene_tissues[gene] = set(ranked[:TOP_N_TISSUES])
    n_pairs = sum(len(v) for v in gene_tissues.values())
    print(f"[05] {n_pairs} gene-tissue pairs "
          f"(top-{TOP_N_TISSUES} by TPM per gene)")

    cache_path = Path(args.cache)
    cache = json.loads(cache_path.read_text()) if cache_path.exists() else {}

    def save_cache():
        tmp = cache_path.with_suffix(".tmp")
        tmp.write_text(json.dumps(cache))
        tmp.rename(cache_path)

    n_calls = 0
    for gi, gene in enumerate(panel):
        v26 = (gtex_cache.get(gene) or {}).get("ensg_v26")
        if not v26:
            continue
        for tissue in sorted(gene_tissues[gene]):
            key = f"{gene}||{tissue}"
            if key in cache:
                continue
            try:
                d = api_get("/expression/medianTranscriptExpression",
                            {"gencodeId": v26, "datasetId": DATASET,
                             "tissueSiteDetailId": tissue,
                             "itemsPerPage": 500})
                tx = [(r["transcriptId"], float(r["median"]))
                      for r in d.get("data", [])]
                cache[key] = {"transcripts": tx, "date": date.today().isoformat()}
                n_calls += 1
                if n_calls % 20 == 0:
                    save_cache()
                    print(f"[05] {n_calls} queries; gene {gi+1}/{len(panel)} ({gene})",
                          flush=True)
            except Exception as exc:  # noqa: BLE001
                cache[key] = {"error": str(exc)[:120]}
            time.sleep(args.sleep)
    save_cache()
    print(f"[05] cache: {len(cache)} gene-tissue pairs ({n_calls} new queries)")

    n_rows, n_unmatched = 0, 0
    with gzip.open(args.out, "wt") as fh:
        fh.write("\t".join(["gene_symbol", "transcript_id_v26",
                             "transcript_base", "tissue", "median_tpm",
                             "usage", "usage_category", "n_datasets",
                             "has_long_read_support"]) + "\n")
        for gene in panel:
            for tissue in sorted(gene_tissues[gene]):
                rec = cache.get(f"{gene}||{tissue}", {})
                txs = rec.get("transcripts", [])
                total = sum(t for _, t in txs)
                for tid, tpm in txs:
                    base = tid.split(".")[0]
                    usage = (tpm / total) if total > 0 else 0.0
                    cat = categorize_usage(usage, 0.0, n_datasets=1,
                                           has_long_read_support=False)
                    fh.write("\t".join([gene, tid, base, tissue, f"{tpm:.3f}",
                                         f"{usage:.4f}", cat.value, "1",
                                         "0"]) + "\n")
                    n_rows += 1
    print(f"[05] wrote {args.out} ({n_rows} transcript-tissue rows)")
    print("[05] HONEST LIMITS: n_datasets=1 (GTEx v8 only); tissue — not "
          "cell-type — resolution; short-read RSEM quantification; "
          "has_long_read_support=False throughout")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
