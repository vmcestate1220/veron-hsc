"""
Pipeline smoke tests for veron-hsc.

Verifies pipeline components work end-to-end using real AF3 data
(if present in results/af3_outputs/) or skips gracefully.
"""
import json
import os
import tempfile

import numpy as np
import pandas as pd
import pytest

PROJECT_ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
REAL_CIF = os.path.join(
    PROJECT_ROOT, "results", "af3_outputs", "cd34_vs_ch02_v10a",
    "fold_cd34_vs_ch02_v10a_model_0.cif",
)
SCREENING_CSV = os.path.join(PROJECT_ROOT, "results", "screening_results.csv")
LEADS_CSV = os.path.join(PROJECT_ROOT, "results", "veron_prioritized_leads.csv")

has_real_cif = os.path.exists(REAL_CIF)
has_screening = os.path.exists(SCREENING_CSV)
has_leads = os.path.exists(LEADS_CSV)


# ── CIF Parser ──────────────────────────────────


class TestCIFParser:
    @pytest.mark.skipif(not has_real_cif, reason="Real AF3 data not available")
    def test_parse_real_cif_has_two_chains(self):
        from scripts.utils import parse_cif_atoms
        chains = parse_cif_atoms(REAL_CIF, ca_only=True)
        assert "A" in chains, "Missing receptor chain A"
        assert "B" in chains, "Missing ligand chain B"

    @pytest.mark.skipif(not has_real_cif, reason="Real AF3 data not available")
    def test_parse_real_cif_receptor_length(self):
        from scripts.utils import parse_cif_atoms
        chains = parse_cif_atoms(REAL_CIF, ca_only=True)
        assert len(chains["A"]) == 141, f"Expected 141 receptor CAs, got {len(chains['A'])}"

    @pytest.mark.skipif(not has_real_cif, reason="Real AF3 data not available")
    def test_parse_real_cif_ligand_length(self):
        from scripts.utils import parse_cif_atoms
        chains = parse_cif_atoms(REAL_CIF, ca_only=True)
        assert len(chains["B"]) == 12, f"Expected 12 ligand CAs, got {len(chains['B'])}"

    @pytest.mark.skipif(not has_real_cif, reason="Real AF3 data not available")
    def test_cif_column_count_validation(self):
        """Real CIF ATOM lines should have exactly 18 columns."""
        with open(REAL_CIF) as f:
            for line in f:
                if line.startswith("ATOM"):
                    parts = line.split()
                    assert len(parts) == 18, (
                        f"Expected 18 columns, got {len(parts)}: {line.strip()}"
                    )
                    break


# ── Stealth Scoring ──────────────────────────────


class TestStealthScoring:
    def test_perfect_stealth(self):
        from scripts.run_screening import compute_stealth_scores
        df = pd.DataFrame([{"candidate_id": "X", "n_9mers": 4, "n_hits": 0}])
        result = compute_stealth_scores(df)
        assert result.iloc[0]["stealth_score"] == 1.0

    def test_partial_stealth(self):
        from scripts.run_screening import compute_stealth_scores
        df = pd.DataFrame([{"candidate_id": "X", "n_9mers": 10, "n_hits": 5}])
        result = compute_stealth_scores(df)
        assert result.iloc[0]["stealth_score"] == pytest.approx(0.5)

    def test_no_9mers_defaults_to_1(self):
        from scripts.run_screening import compute_stealth_scores
        df = pd.DataFrame([{"candidate_id": "X", "n_9mers": 0, "n_hits": 0}])
        result = compute_stealth_scores(df)
        assert result.iloc[0]["stealth_score"] == 1.0

    def test_ranking_order(self):
        from scripts.run_screening import compute_stealth_scores
        df = pd.DataFrame([
            {"candidate_id": "BAD", "n_9mers": 10, "n_hits": 8},
            {"candidate_id": "GOOD", "n_9mers": 10, "n_hits": 0},
            {"candidate_id": "MID", "n_9mers": 10, "n_hits": 3},
        ])
        result = compute_stealth_scores(df)
        assert list(result["candidate_id"]) == ["GOOD", "MID", "BAD"]


# ── Lead Score ───────────────────────────────────


class TestLeadScore:
    def test_lead_score_formula(self):
        from scripts.postprocess import merge_and_rank
        structural = pd.DataFrame([{
            "candidate_id": "TEST",
            "iptm": 0.80,
            "ligand_plddt": 90.0,
        }])
        with tempfile.NamedTemporaryFile(mode="w", suffix=".csv", delete=False) as f:
            f.write("candidate_id,stealth_score\n")
            f.write("TEST,1.0\n")
            screening_path = f.name
        try:
            result = merge_and_rank(structural, screening_path)
            # 0.4*1.0 + 0.4*0.8 + 0.2*(90/100) = 0.4+0.32+0.18 = 0.90
            assert result.iloc[0]["lead_score"] == pytest.approx(0.90, abs=0.01)
        finally:
            os.unlink(screening_path)

    def test_lead_score_ranking_order(self):
        from scripts.postprocess import merge_and_rank
        structural = pd.DataFrame([
            {"candidate_id": "LOW", "iptm": 0.20, "ligand_plddt": 30.0},
            {"candidate_id": "HIGH", "iptm": 0.90, "ligand_plddt": 95.0},
        ])
        with tempfile.NamedTemporaryFile(mode="w", suffix=".csv", delete=False) as f:
            f.write("candidate_id,stealth_score\n")
            f.write("LOW,1.0\n")
            f.write("HIGH,1.0\n")
            screening_path = f.name
        try:
            result = merge_and_rank(structural, screening_path)
            assert result.iloc[0]["candidate_id"] == "HIGH"
            assert result.iloc[1]["candidate_id"] == "LOW"
        finally:
            os.unlink(screening_path)


# ── MWSP Motif Distance ─────────────────────────


class TestMWSPMotifDistance:
    def test_motif_not_found_returns_nan(self):
        from scripts.postprocess import compute_motif_distance
        result = compute_motif_distance("/nonexistent.cif", "AAAAAAAA")
        assert np.isnan(result)

    @pytest.mark.skipif(not has_real_cif, reason="Real AF3 data not available")
    def test_motif_distance_v10a_is_positive(self):
        from scripts.postprocess import compute_motif_distance
        dist = compute_motif_distance(REAL_CIF, "THRPPMWSPAWP")
        assert not np.isnan(dist)
        assert 0 < dist < 20


# ── End-to-End Data Validation ───────────────────


class TestEndToEnd:
    @pytest.mark.skipif(not has_screening, reason="screening_results.csv not available")
    def test_screening_csv_has_all_candidates(self):
        df = pd.read_csv(SCREENING_CSV)
        assert len(df) == 12

    @pytest.mark.skipif(not has_leads, reason="veron_prioritized_leads.csv not available")
    def test_prioritized_leads_has_lead_score(self):
        df = pd.read_csv(LEADS_CSV)
        assert "lead_score" in df.columns
        assert "candidate_id" in df.columns
        assert all(df["lead_score"].notna())

    @pytest.mark.skipif(not has_leads, reason="veron_prioritized_leads.csv not available")
    def test_8g12_control_penalized_for_immunogenicity(self):
        """8G12 scFv should not be #1 despite high iPTM — stealth penalty."""
        df = pd.read_csv(LEADS_CSV)
        top = df.iloc[0]
        assert top["candidate_id"] != "8G12_scFv", (
            "8G12 control should not be top lead"
        )
