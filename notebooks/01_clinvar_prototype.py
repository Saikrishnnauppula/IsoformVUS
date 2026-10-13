# %% [markdown]
# 01 — ClinVar prototype exploration
#
# Smallest meaningful experiment: can we extract a clean, labeled,
# gene-grouped ClinVar variant set suitable for the IsoformVUS panel?
# Run after `scripts/01_download_sources.py` completes.

# %% imports
import gzip
import re
from collections import Counter
from pathlib import Path

REPO = Path(__file__).resolve().parents[1]
RAW = REPO / "data" / "raw"
VCF = RAW / "clinvar_20261004.vcf.gz"

# %% helpers
SIG_MAP = {
    "Pathogenic": 1, "Likely_pathogenic": 1,
    "Pathogenic/Likely_pathogenic": 1,
    "Benign": 0, "Likely_benign": 0, "Benign/Likely_benign": 0,
}

def parse_info(info: str) -> dict:
    d = {}
    for kv in info.split(";"):
        if "=" in kv:
            k, v = kv.split("=", 1)
            d[k] = v
    return d

def significance_bucket(clnsig: str) -> str | None:
    """Map CLNSIG to pathogenic / benign / other (VUS, conflicting, ...)."""
    s = clnsig.replace(" ", "")
    if "Pathogenic" in s and "Benign" not in s and "Uncertain" not in s:
        return "pathogenic"
    if "Benign" in s and "Pathogenic" not in s and "Uncertain" not in s:
        return "benign"
    return "other"

# %% stream the VCF (no full load)
assert VCF.exists(), f"missing {VCF}; run scripts/01_download_sources.py first"
buckets = Counter()
gene_pathogenic = Counter()
stars = Counter()
n = 0
with gzip.open(VCF, "rt") as fh:
    for line in fh:
        if line.startswith("#"):
            continue
        n += 1
        cols = line.rstrip("\n").split("\t")
        info = parse_info(cols[7])
        b = significance_bucket(info.get("CLNSIG", ""))
        buckets[b] += 1
        stars[info.get("CLNREVSTAT", "?")] += 1
        if b == "pathogenic":
            geneinfo = info.get("GENEINFO", "")
            gene = geneinfo.split(":")[0] if geneinfo else "?"
            gene_pathogenic[gene] += 1
        if n >= 2_000_000:  # prototype cap; full parse in stage 03
            break

print(f"records scanned: {n}")
print("buckets:", dict(buckets))
print("top review statuses:", stars.most_common(5))
print("top-15 genes by pathogenic count:", gene_pathogenic.most_common(15))

# %% prototype questions (answered once full parse runs in stage 03)
# Q1. How many genes have >=10 P/LP with >=1 star?  -> panel feasibility
# Q2. What fraction of records are VUS/conflicting? -> label noise estimate
# Q3. Do top pathogenic genes have multiple protein-coding isoforms in GENCODE v46?
print("\n[prototype] full answers pending stage 02/03 implementation")
