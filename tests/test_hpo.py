"""Tests for HPO ontology utilities."""

from isoformvus.data.hpo import categorize_gene, parse_obo, top_level_categories

TINY_OBO = """\
format-version: 1.2

[Term]
id: HP:0000118
name: Phenotypic abnormality

[Term]
id: HP:0000707
name: Abnormality of the nervous system
is_a: HP:0000118 ! Phenotypic abnormality

[Term]
id: HP:0001250
name: Seizure
is_a: HP:0000707 ! Abnormality of the nervous system

[Term]
id: HP:0001638
name: Cardiomyopathy
is_a: HP:0001626 ! Abnormality of the cardiovascular system

[Term]
id: HP:0001626
name: Abnormality of the cardiovascular system
is_a: HP:0000118 ! Phenotypic abnormality

[Term]
id: HP:9999999
name: Orphan term
"""


def test_parse_obo(tmp_path):
    p = tmp_path / "hp.obo"
    p.write_text(TINY_OBO)
    names, parents = parse_obo(p)
    assert names["HP:0001250"] == "Seizure"
    assert parents["HP:0001250"] == ["HP:0000707"]


def test_categorize_gene(tmp_path):
    p = tmp_path / "hp.obo"
    p.write_text(TINY_OBO)
    names, parents = parse_obo(p)
    assert categorize_gene({"HP:0001250"}, parents, names) == "nervous_system"
    assert categorize_gene({"HP:0001638"}, parents, names) == "cardiovascular_system"
    # majority vote
    assert categorize_gene({"HP:0001250", "HP:0001638", "HP:0001626"},
                            parents, names) == "cardiovascular_system"
    # orphan term -> unmapped
    assert categorize_gene({"HP:9999999"}, parents, names) == "unmapped"
    assert categorize_gene(set(), parents, names) == "unmapped"


def test_top_level_categories_counts(tmp_path):
    p = tmp_path / "hp.obo"
    p.write_text(TINY_OBO)
    _, parents = parse_obo(p)
    cats = top_level_categories({"HP:0001250", "HP:0001638"}, parents)
    assert cats["HP:0000707"] == 1
    assert cats["HP:0001626"] == 1
