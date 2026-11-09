"""Data-grounded tests for the 2026-11-09 journal-panel fix round.

These pin the numbers the fix round shipped: the corrected tensor sums,
the gene-cluster bootstrap intervals for the CMH odds ratios, the
bounding analysis over unresolved variants, the recomputed likelihood
ratio table, the renamed statistics block, the CLNSIG composite
sensitivity, the quarantined carried-over values, and the annotated
tensor that manuscript Section 4.1 now reports as primary. Expected
literals are the values produced by scripts 18, 19, and 21 on the
pinned stack on 2026-11-09.
"""

import json
from pathlib import Path

import pandas as pd
import pytest

REPO = Path(__file__).resolve().parents[1]
RESULTS = REPO / "results"


@pytest.fixture(scope="module")
def stage18():
    return json.loads((RESULTS / "exon_architecture_nmd.json").read_text())


@pytest.fixture(scope="module")
def statistics():
    return json.loads((RESULTS / "statistics.json").read_text())


@pytest.fixture(scope="module")
def validation():
    return json.loads((RESULTS / "exon_reconstruction_validation.json").read_text())


def test_tensor_sums(stage18):
    tensor = stage18["exon_phase_symmetry"]["tensor_4x2x2"]
    assert len(tensor) == 16
    assert sum(c["n_exons"] for c in tensor) == 4509
    assert sum(c["missense_n"] for c in tensor) == 220209
    assert sum(c["missense_plp"] + c["missense_blb"] for c in tensor) == 27965


def test_cmh_ors_and_clustered_cis(stage18):
    cmh = stage18["cmh_within_gene_missense"]
    length = cmh["giant_large_exon_gt400bp"]
    backup = cmh["expr_backup_skip"]
    assert length["cmh_or"] == pytest.approx(0.2005, abs=1e-4)
    assert backup["cmh_or"] == pytest.approx(0.3461, abs=1e-4)
    # p-values and z statistics are withheld at source.
    for block in (length, backup):
        assert "p" not in block and "p_value" not in block
        assert block["z"] is None
        boot = block["gene_cluster_bootstrap"]
        assert boot["n_boot"] == 2000
    lo, hi = length["gene_cluster_bootstrap"]["ci95"]
    assert lo == pytest.approx(0.1119, abs=2e-3)
    assert hi == pytest.approx(0.3384, abs=2e-3)
    assert hi < 1.0  # length association survives clustering
    blo, bhi = backup["gene_cluster_bootstrap"]["ci95"]
    assert blo == pytest.approx(0.1280, abs=2e-3)
    assert bhi == pytest.approx(1.0299, abs=2e-3)
    assert blo < 1.0 < bhi  # backup CI includes 1, as shipped


def test_bounding_block(stage18):
    b = stage18["vus_bounding_analysis"]
    cells = {
        (c["size_bin"], c["backup_skip"], c["symmetric_phase"]): c for c in b["cells"]
    }
    c1b = cells[("1_Compact_<=180", 0, 1)]
    assert c1b["true_rate_min"] == pytest.approx(0.1107, abs=1e-3)
    assert c1b["true_rate_max"] == pytest.approx(0.9591, abs=1e-3)
    c4d = cells[("4_Giant_>1000", 1, 1)]
    assert c4d["true_rate_min"] == pytest.approx(0.0007, abs=1e-3)
    assert c4d["true_rate_max"] == pytest.approx(0.9409, abs=1e-3)
    # The pooled intervals overlap completely: the gradient is a
    # property of the resolved subset.
    cc = b["pooled"]["compact_constitutive"]
    gs = b["pooled"]["giant_backup_skipped"]
    assert cc["true_rate_min"] < gs["true_rate_max"]
    assert gs["true_rate_min"] < cc["true_rate_max"]


def test_usage_threshold_sweep(stage18):
    sweep = stage18["sensitivity_checks"]["usage_threshold_sweep"]
    by_thr = {row["usage_threshold"]: row for row in sweep}
    ors = [by_thr[t]["backup_skip_cmh_or"] for t in (0.05, 0.10, 0.20)]
    assert ors == pytest.approx([0.6063, 0.3461, 0.3650], abs=1e-3)


def test_lr_table_matches_tensor(stage18):
    tensor = stage18["exon_phase_symmetry"]["tensor_4x2x2"]
    table = stage18["tiered_model_and_bayesian_lr"]["bayesian_lr_table"]
    assert sum(r["n_exons"] for r in table) == sum(c["n_exons"] for c in tensor)
    tmap = {(c["size_bin"], c["backup_skip"], c["symmetric_phase"]): c for c in tensor}
    for row in table:
        cell = tmap[(row["size_bin"], row["backup_skip"], row["symmetric_phase"])]
        assert row["n_exons"] == cell["n_exons"]
        assert row["n_genes"] == cell["n_genes"]
    cell4d = [
        r
        for r in table
        if r["size_bin"] == "4_Giant_>1000"
        and r["backup_skip"] == 1
        and r["symmetric_phase"] == 1
    ][0]
    assert cell4d["likelihood_ratio"] == pytest.approx(0.0088, abs=1e-4)


def test_carried_over_values_quarantined(stage18):
    assert "causal_mediation" not in stage18
    sj = stage18["splice_junction_boundary_asymmetry"]
    assert "distance_profile" not in sj
    assert "donor_vs_acceptor_edge_0_3bp" not in sj
    withdrawn = json.loads((RESULTS / "withdrawn_carried_over_values.json").read_text())
    assert withdrawn["status"].startswith("withdrawn 2026-11-09")
    assert "causal_mediation" in withdrawn
    assert "distance_profile" in withdrawn["splice_junction_boundary_asymmetry"]


def test_candidate_tsv_renamed():
    df = pd.read_csv(RESULTS / "exon_skipping_nmd_candidates.tsv", sep="\t")
    assert "Exon_Nomination_Score" in df.columns
    assert not any("ASO" in c for c in df.columns)
    tiers = df["aso_tier"].value_counts().to_dict()
    tier1 = [k for k in tiers if k.startswith("Tier 1")][0]
    tier2 = [k for k in tiers if k.startswith("Tier 2")][0]
    assert tiers[tier1] == 302
    assert tiers[tier2] == 626
    assert "computational nomination only" in tier1
    assert "Gold-Standard" not in tier1


def test_test7_rates(stage18):
    t7 = stage18["sensitivity_checks"]["test7_prospective_2024_to_2026"]
    assert t7["compact_constit_p_rate"] == pytest.approx(0.7861, abs=1e-3)
    assert t7["giant_skip_p_rate"] == pytest.approx(0.0645, abs=1e-3)
    assert t7["n_new_2026_internal_missense"] == 2408


def test_statistics_demo_block_gated(statistics):
    assert "m2_m3_clarification" not in statistics
    demo = statistics["step_metric_tie_artifact_demo"]
    assert "not evidence" in demo["note"]


def test_clnsig_sensitivity(statistics):
    sens = statistics["clnsig_composite_sensitivity"]
    assert sens["n_excluded_total"] == 19
    assert sens["n_excluded_test"] == 17
    assert sens["max_abs_diff"] == pytest.approx(0.000862, abs=1e-5)


def test_topk_ci_suppressed():
    bench = json.loads((RESULTS / "benchmark.json").read_text())
    for tier in ("M1", "M2", "M3", "M4", "M5"):
        assert bench["tiers"][tier]["test"]["top10_recall_ci95"] is None
    assert "saturated" in bench["topk_recall_ci_note"]


def test_figure_scripts_have_no_dead_literals():
    for name in ("15_make_figures.py", "16_make_missing_figures.py"):
        text = (REPO / "scripts" / name).read_text()
        for dead in ("0.947", "19/50", "5 affected"):
            assert dead not in text, f"{dead} still in {name}"


def test_annotated_tensor_primary(validation):
    atc = validation["annotated_tensor_comparison"]
    ann = atc["annotated_tensor_16_strict_internal"]
    assert sum(c["n_exons"] for c in ann) == 3758
    assert sum(c["missense_n"] for c in ann) == 209508
    assert sum(c["plp"] + c["blb"] for c in ann) == 26801
    cells = {(c["size_bin"], c["backup_skip"], c["symmetric_phase"]): c for c in ann}
    c1b = cells[("1_Compact_<=180", 0, 1)]
    assert c1b["resolved_plp_rate"] == pytest.approx(0.7811, abs=1e-3)
    c4d = cells[("4_Giant_>1000", 1, 1)]
    assert c4d["resolved_plp_rate"] == pytest.approx(0.0138, abs=1e-3)
    assert c4d["n_exons"] == 7
    assert c4d["n_genes"] == 3
    cmh = atc["annotated_cmh"]
    assert cmh["giant_large_exon_gt400bp"]["cmh_or"] == pytest.approx(0.1860, abs=1e-3)
    assert cmh["expr_backup_skip"]["cmh_or"] == pytest.approx(0.2720, abs=1e-3)
    blo, bhi = cmh["expr_backup_skip"]["gene_cluster_bootstrap"]["ci95"]
    assert bhi < 1.0  # annotated-frame backup CI excludes 1, as shipped


def test_atlas_validation_counts(validation):
    assert validation["atlas_summary"]["genes_covered"] == 102
    assert validation["atlas_summary"]["cds_exons"] == 3965
    assert (
        validation["atlas_summary"][
            "internal_coding_exons_strict_excl_first_and_last_cds"
        ]
        == 3761
    )
    strand = validation["strand_validation"]
    assert strand["n_discordant"] == 0
    fidelity = validation["length_fidelity"]["coordinate_matching"]
    assert fidelity["pct_unique_with_rel_diff_gt_10pct"] == pytest.approx(
        14.65, abs=0.1
    )
    assert fidelity["n_no_overlap"] == 625
