"""Cell type -> isoform usage.

For each (gene, cell type), estimate which isoforms are used and how confident
we are. Because short-read scRNA-seq cannot reliably resolve full-length
isoforms (3' bias, ambiguous multi-mapping), usage is reported as *evidence
categories*, not precise point estimates:

    dominant   - consistently the highest-used isoform across supporting datasets
    supported  - detected above background in the cell type, not dominant
    low_usage  - detected at low level; NOT assumed biologically irrelevant
                 (low-abundance isoforms can be functional, e.g. poison exons,
                 NMD targets, or tightly regulated switches)
    uncertain  - conflicting evidence across datasets, or below detection with
                 inadequate depth to call absence

Each assignment carries:
    - n_datasets: number of independent datasets behind the call
    - has_long_read_support: bool (full-length evidence exists)
    - notes: free-text caveats (e.g. "3' bias; isoforms differ only in 5' exons")

The quantitative substrate is a per-(gene, cell type, transcript) usage matrix
with values in [0, 1] (fraction of gene expression attributed to the isoform),
plus a width (uncertainty interval). Downstream features consume the category
(one-hot) AND the usage value with its width, so the model sees both the
estimate and our uncertainty about it.
"""

from __future__ import annotations

from dataclasses import dataclass, field
from enum import Enum


class UsageCategory(str, Enum):
    DOMINANT = "dominant"
    SUPPORTED = "supported"
    LOW_USAGE = "low_usage"
    UNCERTAIN = "uncertain"


@dataclass(frozen=True)
class IsoformUsage:
    """Usage of one transcript in one cell type."""

    gene_symbol: str
    transcript_id: str
    cell_ontology_id: str
    cell_type_name: str
    usage: float                  # 0..1 fraction of gene expression
    usage_width: float = 0.0      # half-width of uncertainty interval
    category: UsageCategory = UsageCategory.UNCERTAIN
    n_datasets: int = 0
    has_long_read_support: bool = False
    notes: str = ""

    def __post_init__(self) -> None:
        if not 0.0 <= self.usage <= 1.0:
            raise ValueError("usage must be in [0, 1]")
        if self.usage_width < 0:
            raise ValueError("usage_width must be >= 0")


def categorize_usage(usage: float, usage_width: float, n_datasets: int,
                     has_long_read_support: bool,
                     dominant_threshold: float = 0.5,
                     detection_threshold: float = 0.05) -> UsageCategory:
    """Assign a usage category from the estimate and its support.

    Rules (thresholds in config; defaults documented here):
    - < detection_threshold and a single noisy dataset -> UNCERTAIN, not "absent".
      Absence of evidence is not evidence of absence at scRNA-seq depth.
    - >= dominant_threshold with >= 2 datasets (or 1 + long-read) -> DOMINANT.
    - >= detection_threshold otherwise -> SUPPORTED.
    - < detection_threshold but detected in >= 2 datasets -> LOW_USAGE
      (reproducible low signal; kept, not discarded).
    """
    if n_datasets < 1:
        return UsageCategory.UNCERTAIN
    if usage >= dominant_threshold and (n_datasets >= 2 or has_long_read_support):
        return UsageCategory.DOMINANT
    if usage >= detection_threshold:
        return UsageCategory.SUPPORTED
    if n_datasets >= 2:
        return UsageCategory.LOW_USAGE
    return UsageCategory.UNCERTAIN


def check_usage_sums(usages: list[IsoformUsage], tol: float = 0.05) -> list[str]:
    """Sanity check: per (gene, cell type), isoform usages should ~sum to 1.

    Returns warning strings for groups violating the tolerance. A systematic
    violation indicates unannotated isoforms or quantification leakage — a
    finding, not a bug to silently normalize away.
    """
    from collections import defaultdict
    groups: dict[tuple[str, str], float] = defaultdict(float)
    for u in usages:
        groups[(u.gene_symbol, u.cell_ontology_id)] += u.usage
    warnings = []
    for (gene, ct), total in sorted(groups.items()):
        if abs(total - 1.0) > tol:
            warnings.append(
                f"{gene} @ {ct}: isoform usages sum to {total:.3f} (tolerance {tol})")
    return warnings
