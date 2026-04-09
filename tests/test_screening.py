"""Tests for MHC screening with multi-allele support."""
import pandas as pd
import pytest

from scripts.screening_utils import DEFAULT_ALLELES, tile_9mers


def test_tile_9mers_basic():
    result = tile_9mers("ABCDEFGHIJK")  # 11-mer -> 3 9-mers
    assert result == ["ABCDEFGHI", "BCDEFGHIJ", "CDEFGHIJK"]


def test_tile_9mers_exact_9():
    result = tile_9mers("ABCDEFGHI")
    assert result == ["ABCDEFGHI"]


def test_tile_9mers_too_short():
    result = tile_9mers("ABCDEFGH")
    assert result == []


def test_default_alleles_includes_common():
    assert "HLA-A*02:01" in DEFAULT_ALLELES
    assert "HLA-B*07:02" in DEFAULT_ALLELES
    assert "HLA-A*01:01" in DEFAULT_ALLELES


def test_screen_candidates_short_sequence():
    """Sequences too short for 9-mer tiling should not crash."""
    from scripts.screening_utils import screen_candidates
    seqs = {"TEST_SHORT": "AAAA"}
    summary, _ = screen_candidates(seqs, alleles=["HLA-A*01:01"])
    assert "n_hits" in summary.columns
    assert len(summary) == 1
    assert summary.iloc[0]["n_hits"] == 0
    assert not summary.iloc[0]["flagged"]


def test_stealth_score_computation():
    from scripts.run_screening import compute_stealth_scores
    df = pd.DataFrame([
        {"candidate_id": "A", "n_9mers": 10, "n_hits": 0},
        {"candidate_id": "B", "n_9mers": 10, "n_hits": 3},
        {"candidate_id": "C", "n_9mers": 0, "n_hits": 0},
    ])
    result = compute_stealth_scores(df)
    assert result.loc[result["candidate_id"] == "A", "stealth_score"].iloc[0] == 1.0
    assert result.loc[result["candidate_id"] == "B", "stealth_score"].iloc[0] == pytest.approx(0.7)
    assert result.loc[result["candidate_id"] == "C", "stealth_score"].iloc[0] == 1.0
