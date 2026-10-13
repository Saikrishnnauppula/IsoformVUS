#!/usr/bin/env python3
"""Stage 04a: GTEx v8 median gene expression for panel genes (real measurement).

Two-step route via the GTEx Portal v2 API (https://gtexportal.org/api/v2):
    1. /reference/gene?geneId=<ENSG base> -> v26 versioned gencodeId
    2. /expression/medianGeneExpression?gencodeId=<v26 id>&datasetId=gtex_v8
       -> median TPM per tissue (54 tissues)

ENSG base IDs come from the GENCODE v46 GTF (stable across GENCODE releases;
only the version suffix changes, resolved in step 1).

Writes data/processed/gtex_median_tpm.tsv.gz:
    gene_symbol, ensg_v26, tissue, median_tpm
Cache: data/interim/gtex_cache.json (resume-safe).

Provenance: GTEx v8 (54 tissues, 948 donors), GENCODE v26 annotation,
API query date recorded.
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

API = "https://gtexportal.org/api/v2"
DATASET = "gtex_v8"


def api_get(path: str, params: dict) -> dict:
    url = f"{API}{path}?" + urllib.parse.urlencode(params)
    req = urllib.request.Request(url, headers={"Accept": "application/json",
                                               "User-Agent": "IsoformVUS/0.1.0 (research)"})
    with urllib.request.urlopen(req, timeout=30) as resp:
        return json.load(resp)


def ensg_for_symbols(gtf_path: Path, symbols: set[str]) -> dict[str, str]:
    """gene symbol -> unversioned ENSG id from the GENCODE GTF."""
    import re
    attr_re = re.compile(r'(\S+)\s+"([^"]+)"')
    mapping: dict[str, str] = {}
    with gzip.open(gtf_path, "rt") as fh:
        for line in fh:
            if line.startswith("#"):
                continue
            cols = line.split("\t")
            if len(cols) != 9 or cols[2] != "gene":
                continue
            attrs = dict(attr_re.findall(cols[8]))
            sym = attrs.get("gene_name", "")
            if sym in symbols and sym not in mapping:
                mapping[sym] = attrs["gene_id"].split(".")[0]
            if len(mapping) == len(symbols):
                break
    return mapping


def main() -> int:
    ap = argparse.ArgumentParser()
    ap.add_argument("--panel", default=str(REPO / "metadata/gene_panel.tsv"))
    ap.add_argument("--cache", default=str(REPO / "data/interim/gtex_cache.json"))
    ap.add_argument("--out", default=str(REPO / "data/processed/gtex_median_tpm.tsv.gz"))
    ap.add_argument("--sleep", type=float, default=0.3)
    args = ap.parse_args()

    symbols = [l.split("\t")[0] for l in
               open(args.panel).read().splitlines()[1:]]
    print(f"[04a] {len(symbols)} panel genes")

    cache_path = Path(args.cache)
    cache_path.parent.mkdir(parents=True, exist_ok=True)
    cache = json.loads(cache_path.read_text()) if cache_path.exists() else {}

    gtf = REPO / "data/raw/gencode.v46.primary_assembly.annotation.gtf.gz"
    ensg = ensg_for_symbols(gtf, set(symbols))
    print(f"[04a] resolved {len(ensg)}/{len(symbols)} ENSG ids from GENCODE v46")

    n_new = 0
    for sym in symbols:
        if sym in cache:
            continue
        base = ensg.get(sym)
        if not base:
            cache[sym] = {"error": "no ENSG in GTF", "date": date.today().isoformat()}
            continue
        try:
            ref = api_get("/reference/gene", {"geneId": base})
            data = ref.get("data", [])
            if not data:
                cache[sym] = {"error": "no v26 record", "date": date.today().isoformat()}
                continue
            v26 = data[0]["gencodeId"]
            expr = api_get("/expression/medianGeneExpression",
                           {"gencodeId": v26, "datasetId": DATASET,
                            "itemsPerPage": 100})
            tissues = {r["tissueSiteDetailId"]: r["median"]
                       for r in expr.get("data", [])}
            cache[sym] = {"ensg_v26": v26, "tissues": tissues,
                          "date": date.today().isoformat(),
                          "dataset": DATASET}
            n_new += 1
            if n_new % 10 == 0:
                print(f"[04a] {n_new} genes fetched")
        except Exception as exc:  # noqa: BLE001 - network; cached as error
            cache[sym] = {"error": str(exc)[:120], "date": date.today().isoformat()}
        time.sleep(args.sleep)

    cache_path.write_text(json.dumps(cache))
    print(f"[04a] cache: {len(cache)} genes ({n_new} newly fetched)")

    with gzip.open(args.out, "wt") as fh:
        fh.write("gene_symbol\tensg_v26\ttissue\tmedian_tpm\n")
        n_rows = 0
        for sym in symbols:
            rec = cache.get(sym, {})
            for tissue, tpm in sorted(rec.get("tissues", {}).items()):
                fh.write(f"{sym}\t{rec.get('ensg_v26','')}\t{tissue}\t{tpm}\n")
                n_rows += 1
    print(f"[04a] wrote {args.out} ({n_rows} rows)")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
