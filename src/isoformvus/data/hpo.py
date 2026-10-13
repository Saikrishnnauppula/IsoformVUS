"""HPO ontology utilities: parse hp.obo, walk is_a hierarchy.

Used to assign each gene a *coarse disease category*: for the gene's HPO
term set, walk each term up to the direct children of HP:0000118
(Phenotypic abnormality) and take the most frequent child. This replaces
hand-written term->category hints with the ontology itself.
"""

from __future__ import annotations

from collections import Counter
from pathlib import Path

PHENOTYPIC_ABNORMALITY = "HP:0000118"


def parse_obo(path: Path) -> tuple[dict[str, str], dict[str, list[str]]]:
    """Parse hp.obo. Returns (id->name, id->[parent ids])."""
    names: dict[str, str] = {}
    parents: dict[str, list[str]] = {}
    cur: str | None = None
    with open(path) as fh:
        for line in fh:
            line = line.rstrip("\n")
            if line == "[Term]":
                cur = None
            elif line.startswith("id: HP:"):
                cur = line.split("id: ")[1].strip()
                names[cur] = ""
                parents[cur] = []
            elif cur and line.startswith("name: "):
                names[cur] = line.split("name: ", 1)[1].strip()
            elif cur and line.startswith("is_a: HP:"):
                parents[cur].append(line.split("is_a: ")[1].split(" ")[0].strip())
    return names, parents


def top_level_categories(hpo_ids: set[str],
                         parents: dict[str, list[str]]) -> Counter:
    """Map each HPO term to its ancestor that is a direct child of
    HP:0000118. Terms that never reach it are counted as 'unmapped'."""
    cats: Counter = Counter()
    for hpo in hpo_ids:
        seen: set[str] = set()
        frontier = [hpo]
        found = None
        while frontier and found is None:
            node = frontier.pop()
            if node in seen:
                continue
            seen.add(node)
            for p in parents.get(node, []):
                if p == PHENOTYPIC_ABNORMALITY:
                    found = node
                    break
                frontier.append(p)
        cats[found if found else "unmapped"] += 1
    return cats


def categorize_gene(hpo_ids: set[str],
                    parents: dict[str, list[str]],
                    names: dict[str, str]) -> str:
    """Coarse category = most frequent top-level HPO category, as its name
    slugified. Ties broken alphabetically (deterministic)."""
    cats = top_level_categories(hpo_ids, parents)
    cats.pop("unmapped", None)
    if not cats:
        return "unmapped"
    top = sorted(cats.items(), key=lambda kv: (-kv[1], kv[0]))[0][0]
    label = names.get(top, top)
    return label.lower().replace("abnormality of the ", "").replace("abnormality of ", "") \
        .replace(" ", "_").replace(",", "")
