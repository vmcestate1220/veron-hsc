"""Tests for shared utility functions."""
import numpy as np
import pytest

from scripts.utils import AA3TO1, extract_receptor_sequence, parse_cif_atoms


def test_aa3to1_has_20_standard():
    assert len(AA3TO1) == 20
    assert AA3TO1["ALA"] == "A"
    assert AA3TO1["TRP"] == "W"


def test_aa3to1_missing_returns_none():
    assert AA3TO1.get("UNK") is None


def test_extract_receptor_sequence(tmp_path):
    """Minimal PDB with 3 residues should produce a 3-letter sequence."""
    pdb = tmp_path / "test.pdb"
    pdb.write_text(
        "ATOM      1  CA  ALA A   1       0.000   0.000   0.000  1.00  0.00\n"
        "ATOM      2  CA  GLY A   2       3.800   0.000   0.000  1.00  0.00\n"
        "ATOM      3  CA  TRP A   3       7.600   0.000   0.000  1.00  0.00\n"
        "END\n"
    )
    seq = extract_receptor_sequence(str(pdb))
    assert seq == "AGW"


def test_parse_cif_atoms_ca_only(tmp_path):
    """parse_cif_atoms with ca_only=True returns only CA atoms."""
    cif = tmp_path / "test.cif"
    cif.write_text(
        "data_test\n"
        "#\n"
        "loop_\n"
        "_atom_site.group_PDB\n"
        "_atom_site.id\n"
        "_atom_site.type_symbol\n"
        "_atom_site.label_atom_id\n"
        "_atom_site.label_alt_id\n"
        "_atom_site.label_comp_id\n"
        "_atom_site.label_asym_id\n"
        "_atom_site.label_entity_id\n"
        "_atom_site.label_seq_id\n"
        "_atom_site.pdbx_PDB_ins_code\n"
        "_atom_site.Cartn_x\n"
        "_atom_site.Cartn_y\n"
        "_atom_site.Cartn_z\n"
        "_atom_site.occupancy\n"
        "_atom_site.B_iso_or_equiv\n"
        "_atom_site.auth_seq_id\n"
        "_atom_site.auth_asym_id\n"
        "_atom_site.pdbx_PDB_model_num\n"
        "ATOM 1  N  N   . ALA A 1 1 ? 1.0 2.0 3.0 1.00 30.0 1 A 1\n"
        "ATOM 2  C  CA  . ALA A 1 1 ? 4.0 5.0 6.0 1.00 45.0 1 A 1\n"
        "ATOM 3  C  C   . ALA A 1 1 ? 7.0 8.0 9.0 1.00 40.0 1 A 1\n"
        "#\n"
    )
    chains = parse_cif_atoms(str(cif), ca_only=True)
    assert "A" in chains
    assert len(chains["A"]) == 1
    _, atom_name, _, seq_id, coord = chains["A"][0]
    assert atom_name == "CA"
    assert seq_id == 1
    np.testing.assert_allclose(coord, [4.0, 5.0, 6.0])


def test_parse_cif_atoms_all(tmp_path):
    """parse_cif_atoms with ca_only=False returns all heavy atoms."""
    cif = tmp_path / "test.cif"
    cif.write_text(
        "data_test\n"
        "#\n"
        "loop_\n"
        "_atom_site.group_PDB\n"
        "_atom_site.id\n"
        "_atom_site.type_symbol\n"
        "_atom_site.label_atom_id\n"
        "_atom_site.label_alt_id\n"
        "_atom_site.label_comp_id\n"
        "_atom_site.label_asym_id\n"
        "_atom_site.label_entity_id\n"
        "_atom_site.label_seq_id\n"
        "_atom_site.pdbx_PDB_ins_code\n"
        "_atom_site.Cartn_x\n"
        "_atom_site.Cartn_y\n"
        "_atom_site.Cartn_z\n"
        "_atom_site.occupancy\n"
        "_atom_site.B_iso_or_equiv\n"
        "_atom_site.auth_seq_id\n"
        "_atom_site.auth_asym_id\n"
        "_atom_site.pdbx_PDB_model_num\n"
        "ATOM 1  N  N   . ALA A 1 1 ? 1.0 2.0 3.0 1.00 30.0 1 A 1\n"
        "ATOM 2  C  CA  . ALA A 1 1 ? 4.0 5.0 6.0 1.00 45.0 1 A 1\n"
        "ATOM 3  C  C   . ALA A 1 1 ? 7.0 8.0 9.0 1.00 40.0 1 A 1\n"
        "#\n"
    )
    chains = parse_cif_atoms(str(cif), ca_only=False)
    assert len(chains["A"]) == 3


def test_parse_cif_skips_short_lines(tmp_path):
    """Lines with fewer than 18 columns should be silently skipped."""
    cif = tmp_path / "test.cif"
    cif.write_text(
        "data_test\n"
        "ATOM 1  C  CA  . ALA A 1 1 ? 4.0 5.0 6.0 1.00 45.0 1 A 1\n"
        "ATOM short line\n"
        "#\n"
    )
    chains = parse_cif_atoms(str(cif), ca_only=True)
    assert len(chains["A"]) == 1
