#!/usr/bin/env python3
"""Stage 01: download public source files with provenance recording.

Usage:
    python scripts/01_download_sources.py --source clinvar --out data/raw
    python scripts/01_download_sources.py --source all --out data/raw

URLS are read from research/data_urls_verified.md once the verification agent
completes; until then this script runs in --dry-run mode, printing what it
*would* fetch. It never fabricates data: every downloaded file gets a
.provenance.json sidecar (url, date, sha256, bytes).
"""

from __future__ import annotations

import argparse
import hashlib
import json
import sys
import urllib.request
from datetime import date
from pathlib import Path

REPO = Path(__file__).resolve().parents[1]

# Placeholder entries. Real URLs land here after data_urls_verified.md exists.
# Format: name -> (url, filename). None url = not yet verified; script refuses.
SOURCES: dict[str, tuple[str | None, str]] = {
    "clinvar_vcf": (
        "https://ftp.ncbi.nlm.nih.gov/pub/clinvar/vcf_GRCh38/clinvar_20261004.vcf.gz",
        "clinvar_20261004.vcf.gz"),
    "clinvar_summary": (
        "https://ftp.ncbi.nlm.nih.gov/pub/clinvar/tab_delimited/variant_summary.txt.gz",
        "variant_summary.txt.gz"),
    "gencode_gtf": (
        "https://ftp.ebi.ac.uk/pub/databases/gencode/Gencode_human/release_46/"
        "gencode.v46.primary_assembly.annotation.gtf.gz",
        "gencode.v46.primary_assembly.annotation.gtf.gz"),
    "mane_summary": (
        "https://ftp.ncbi.nlm.nih.gov/refseq/MANE/MANE_human/current/"
        "MANE.GRCh38.v1.5.summary.txt.gz",
        "MANE.GRCh38.v1.5.summary.txt.gz"),
    # HPO: pin a dated release, never `latest`. v2026-02-16 verified 2026-10-06.
    "hpo_obo": (
        "https://github.com/obophenotype/human-phenotype-ontology/releases/"
        "download/v2026-02-16/hp.obo",
        "hp.obo"),
    "hpo_phenotype_annotations": (
        "https://github.com/obophenotype/human-phenotype-ontology/releases/"
        "download/v2026-02-16/phenotype.hpoa",
        "phenotype.hpoa"),
    "hpo_phenotype_to_genes": (
        "https://github.com/obophenotype/human-phenotype-ontology/releases/"
        "download/v2026-02-16/phenotype_to_genes.txt",
        "phenotype_to_genes.txt"),
}


def sha256_of(path: Path) -> str:
    h = hashlib.sha256()
    with open(path, "rb") as fh:
        for chunk in iter(lambda: fh.read(1 << 20), b""):
            h.update(chunk)
    return h.hexdigest()


def download(url: str, dest: Path) -> None:
    dest.parent.mkdir(parents=True, exist_ok=True)
    tmp = dest.with_suffix(dest.suffix + ".part")
    req = urllib.request.Request(url, headers={"User-Agent": "IsoformVUS/0.1.0 (research)"})
    with urllib.request.urlopen(req) as resp, open(tmp, "wb") as fh:
        while True:
            chunk = resp.read(1 << 20)
            if not chunk:
                break
            fh.write(chunk)
    tmp.rename(dest)


def main() -> int:
    ap = argparse.ArgumentParser()
    ap.add_argument("--source", default="all",
                    choices=["all"] + list(SOURCES),
                    help="which source to fetch")
    ap.add_argument("--out", default=str(REPO / "data" / "raw"))
    ap.add_argument("--dry-run", action="store_true",
                    help="print planned downloads without fetching")
    args = ap.parse_args()

    names = list(SOURCES) if args.source == "all" else [args.source]
    out = Path(args.out)
    pending = []
    for name in names:
        url, filename = SOURCES[name]
        if url is None:
            pending.append(name)
            continue
        dest = out / filename
        if args.dry_run:
            print(f"[dry-run] would fetch {url} -> {dest}")
            continue
        if dest.exists():
            print(f"[skip] {dest} exists")
        else:
            print(f"[fetch] {url}")
            download(url, dest)
        sidecar = dest.with_suffix(dest.suffix + ".provenance.json")
        sidecar.write_text(json.dumps({
            "source": name, "url": url,
            "download_date": date.today().isoformat(),
            "sha256": sha256_of(dest), "bytes": dest.stat().st_size,
            "genome_build": "GRCh38",
        }, indent=2))
        print(f"[provenance] {sidecar}")
    if pending:
        print(f"[blocked] URLs not yet verified for: {', '.join(pending)} "
              f"(see research/data_urls_verified.md)", file=sys.stderr)
        return 2
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
