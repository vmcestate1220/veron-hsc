"""Shared pytest fixtures for veron-hsc tests."""
import os
import pytest

PROJECT_ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))


@pytest.fixture
def project_root():
    return PROJECT_ROOT


@pytest.fixture
def af3_outputs_dir(project_root):
    return os.path.join(project_root, "results", "af3_outputs")


@pytest.fixture
def candidates_fasta(project_root):
    return os.path.join(project_root, "data", "ligands", "candidates.fasta")
