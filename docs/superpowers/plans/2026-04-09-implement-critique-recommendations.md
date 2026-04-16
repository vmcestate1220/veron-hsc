# Implement Critique Recommendations

> **For agentic workers:** REQUIRED SUB-SKILL: Use superpowers:subagent-driven-development (recommended) or superpowers:executing-plans to implement this plan task-by-task. Steps use checkbox (`- [ ]`) syntax for tracking.

**Goal:** Address all 6 prioritized recommendations from the project critique: expand MHC screening to multi-allele, clarify data provenance, fix environment.yml, extract shared utilities, add smoke tests, and reconcile the dual AF3 JSON generation paths.

**Architecture:** Extract duplicated code (AA3TO1, sequence extraction, CIF parsing) into a shared `scripts/utils.py`. Expand `screening_utils.py` to support a configurable allele panel (default: HLA-A\*02:01 + HLA-B\*07:02 + HLA-A\*01:01). Add a `data_source` metadata field to postprocessor output. Delete the standalone `af3_generator.py` (redundant with `run_screening.py`). Add `tests/` with pytest smoke tests. Fill `environment.yml` with actual dependencies.

**Tech Stack:** Python 3.11, pytest, pandas, numpy, mhcflurry, biopython, openmm

---

### Task 1: Extract shared utilities into `scripts/utils.py`

**Files:**
- Create: `scripts/utils.py`
- Modify: `scripts/af3_generator.py` (will be deleted in Task 6, but update imports first)
- Modify: `scripts/run_screening.py:42-47` (remove AA3TO1 + extract_receptor_sequence)
- Modify: `scripts/postprocess.py:162-198` (remove parse_cif_ca_atoms)
- Modify: `scripts/generate_test_af3_data.py:27-33` (remove AA3TO1)
- Modify: `scripts/lead_profile_v10a.py:42-47,62-89` (remove AA3TO1 + parse_cif_all_atoms)
- Test: `tests/test_utils.py`

- [ ] **Step 1: Write failing tests for utils module**

Create `tests/__init__.py` (empty) and `tests/test_utils.py`:

```python
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
```

- [ ] **Step 2: Run tests to verify they fail**

Run: `./veron-hsc/bin/python -m pytest tests/test_utils.py -v`
Expected: FAIL with `ModuleNotFoundError: No module named 'scripts.utils'`

- [ ] **Step 3: Create `scripts/utils.py` with shared code**

```python
"""
Shared utilities for the veron-hsc pipeline.

Consolidates code that was previously duplicated across multiple scripts:
  - AA3TO1 mapping
  - Receptor sequence extraction from PDB
  - mmCIF atom parsing (AF3 Server format)
"""

import numpy as np
from Bio.PDB import PDBParser

AA3TO1 = {
    "ALA": "A", "ARG": "R", "ASN": "N", "ASP": "D", "CYS": "C",
    "GLN": "Q", "GLU": "E", "GLY": "G", "HIS": "H", "ILE": "I",
    "LEU": "L", "LYS": "K", "MET": "M", "PHE": "F", "PRO": "P",
    "SER": "S", "THR": "T", "TRP": "W", "TYR": "Y", "VAL": "V",
}


def extract_receptor_sequence(pdb_path):
    """Extract 1-letter amino acid sequence from a PDB file (standard residues only)."""
    parser = PDBParser(QUIET=True)
    structure = parser.get_structure("receptor", pdb_path)
    return "".join(
        AA3TO1.get(r.get_resname(), "X")
        for r in structure.get_residues()
        if r.get_id()[0] == " "
    )


def parse_cif_atoms(cif_path, ca_only=True):
    """
    Parse atoms from AF3 Server mmCIF files (18-column format).

    AF3 CIF columns (0-indexed):
      [0] group_PDB  [1] id  [2] type_symbol  [3] label_atom_id
      [4] label_alt_id  [5] label_comp_id  [6] label_asym_id
      [7] label_entity_id  [8] label_seq_id  [9] pdbx_PDB_ins_code
      [10] Cartn_x  [11] Cartn_y  [12] Cartn_z
      [13] occupancy  [14] B_iso_or_equiv  [15] auth_seq_id
      [16] auth_asym_id  [17] pdbx_PDB_model_num

    Args:
        cif_path: Path to mmCIF file.
        ca_only: If True, return only CA atoms. If False, return all ATOM records.

    Returns:
        dict: {chain_id: [(atom_serial, atom_name, resname, seq_id, np.array([x,y,z])), ...]}
    """
    chains = {}
    with open(cif_path) as f:
        for line in f:
            if not line.startswith("ATOM"):
                continue
            parts = line.split()
            if len(parts) < 18:
                continue

            atom_name = parts[3]
            if ca_only and atom_name != "CA":
                continue

            atom_serial = int(parts[1])
            resname = parts[5]
            chain = parts[6]
            seq_id = int(parts[8])
            x = float(parts[10])
            y = float(parts[11])
            z = float(parts[12])

            chains.setdefault(chain, []).append(
                (atom_serial, atom_name, resname, seq_id, np.array([x, y, z]))
            )
    return chains
```

- [ ] **Step 4: Run tests to verify they pass**

Run: `./veron-hsc/bin/python -m pytest tests/test_utils.py -v`
Expected: 5 tests PASS

- [ ] **Step 5: Update `run_screening.py` to use shared utils**

Remove lines 42-47 (AA3TO1 dict) and lines 55-62 (extract_receptor_sequence function). Add import and use it:

Replace the `AA3TO1` dict and `extract_receptor_sequence` function with:
```python
from scripts.utils import AA3TO1, extract_receptor_sequence
```

Remove the `AA3TO1 = { ... }` block (lines 42-47) and the `def extract_receptor_sequence(pdb_path):` function (lines 55-62).

- [ ] **Step 6: Update `postprocess.py` to use shared utils**

Replace the `parse_cif_ca_atoms` function (lines 162-198) with an import and thin wrapper.

At the top of the file, add:
```python
from scripts.utils import parse_cif_atoms
```

Replace the `parse_cif_ca_atoms` function body. In `compute_motif_distance`, change:
```python
chains = parse_cif_ca_atoms(cif_path)
```
to:
```python
chains = parse_cif_atoms(cif_path, ca_only=True)
```

And update the unpacking since parse_cif_atoms returns 5-tuples `(serial, name, resname, seq_id, coord)`:
```python
receptor_coords = np.array([c for _, _, _, _, c in chains["A"]])
motif_coords = np.array(
    [c for _, _, _, sid, c in chains["B"] if sid in motif_ids
])
```

- [ ] **Step 7: Update `generate_test_af3_data.py` to use shared utils**

Remove the `AA3TO1` dict (lines 27-33). Add import:
```python
from scripts.utils import AA3TO1
```

- [ ] **Step 8: Update `lead_profile_v10a.py` to use shared utils**

Remove the `AA3TO1` dict (lines 42-47) and the `parse_cif_all_atoms` function (lines 62-89). Add imports:
```python
from scripts.utils import AA3TO1, parse_cif_atoms
```

Replace all calls to `parse_cif_all_atoms(path)` with `parse_cif_atoms(path, ca_only=False)`.

In `find_binding_pocket`, update unpacking:
```python
motif_atoms = [(a[3], a[4]) for a in chains["B"] if a[3] in motif_ids]
```
stays the same (index 3 = seq_id, index 4 = coord).

In `per_residue_contacts`, update:
```python
receptor_coords = np.array([a[4] for a in chains["A"]])
```
stays the same (index 4 = coord).

Ligand residue grouping update:
```python
for atom in chains["B"]:
    _, _, resname, seq_id, xyz = atom
```
stays the same (5-tuple matches).

- [ ] **Step 9: Verify all existing scripts still work**

Run: `./veron-hsc/bin/python -c "from scripts.utils import AA3TO1, extract_receptor_sequence, parse_cif_atoms; print('utils OK')"` 
Run: `./veron-hsc/bin/python -m pytest tests/test_utils.py -v`
Expected: All pass

- [ ] **Step 10: Commit**

```bash
git add scripts/utils.py tests/__init__.py tests/test_utils.py scripts/run_screening.py scripts/postprocess.py scripts/generate_test_af3_data.py scripts/lead_profile_v10a.py
git commit -m "refactor: extract shared utils (AA3TO1, seq extraction, CIF parsing)"
```

---

### Task 2: Expand MHC screening to multi-allele panel

**Files:**
- Modify: `scripts/screening_utils.py:17` (add default allele panel)
- Modify: `scripts/run_screening.py:97-108` (update stealth scoring for multi-allele)
- Test: `tests/test_screening.py`

- [ ] **Step 1: Write failing tests for multi-allele screening**

Create `tests/test_screening.py`:

```python
"""Tests for MHC screening with multi-allele support."""
import pytest


def test_tile_9mers_basic():
    from scripts.screening_utils import tile_9mers
    result = tile_9mers("ABCDEFGHIJK")  # 11-mer → 3 9-mers
    assert result == ["ABCDEFGHI", "BCDEFGHIJ", "CDEFGHIJK"]


def test_tile_9mers_exact_9():
    from scripts.screening_utils import tile_9mers
    result = tile_9mers("ABCDEFGHI")
    assert result == ["ABCDEFGHI"]


def test_tile_9mers_too_short():
    from scripts.screening_utils import tile_9mers
    result = tile_9mers("ABCDEFGH")
    assert result == []


def test_default_alleles_includes_common():
    from scripts.screening_utils import DEFAULT_ALLELES
    assert "HLA-A*02:01" in DEFAULT_ALLELES
    assert "HLA-B*07:02" in DEFAULT_ALLELES
    assert "HLA-A*01:01" in DEFAULT_ALLELES


def test_screen_candidates_returns_per_allele_hits():
    """screen_candidates summary should have allele-level detail columns."""
    from scripts.screening_utils import screen_candidates
    seqs = {"TEST_SHORT": "AAAA"}  # too short for 9-mer
    summary, _ = screen_candidates(seqs, alleles=["HLA-A*01:01"])
    assert "n_hits" in summary.columns
    assert len(summary) == 1
    assert summary.iloc[0]["n_hits"] == 0


def test_stealth_score_multi_allele():
    """Stealth score should reflect worst allele, not just one."""
    from scripts.run_screening import compute_stealth_scores
    import pandas as pd
    df = pd.DataFrame([
        {"candidate_id": "A", "n_9mers": 10, "n_hits": 0},
        {"candidate_id": "B", "n_9mers": 10, "n_hits": 3},
        {"candidate_id": "C", "n_9mers": 0, "n_hits": 0},
    ])
    result = compute_stealth_scores(df)
    assert result.loc[result["candidate_id"] == "A", "stealth_score"].iloc[0] == 1.0
    assert result.loc[result["candidate_id"] == "B", "stealth_score"].iloc[0] == 0.7
    assert result.loc[result["candidate_id"] == "C", "stealth_score"].iloc[0] == 1.0
```

- [ ] **Step 2: Run tests to verify they fail**

Run: `./veron-hsc/bin/python -m pytest tests/test_screening.py -v`
Expected: `test_default_alleles_includes_common` fails (DEFAULT_ALLELES doesn't exist yet). Others may pass or fail.

- [ ] **Step 3: Update `screening_utils.py` to support multi-allele panels**

Replace the module-level constants:
```python
ALLELE = "HLA-A*01:01"
```
with:
```python
DEFAULT_ALLELES = [
    "HLA-A*01:01",
    "HLA-A*02:01",  # most common worldwide (~25-50% of populations)
    "HLA-B*07:02",  # common Caucasian allele (~10-15%)
]
```

Update `screen_candidates` signature from:
```python
def screen_candidates(sequences):
```
to:
```python
def screen_candidates(sequences, alleles=None):
```

At the top of the function body, add:
```python
if alleles is None:
    alleles = DEFAULT_ALLELES
```

Replace the prediction loop body. For each candidate, predict against ALL alleles, and track the worst (lowest IC50) across alleles:

```python
    for cand_id, seq in sequences.items():
        if len(seq) < KMER_LEN:
            summary_rows.append(
                {
                    "candidate_id": cand_id,
                    "sequence_length": len(seq),
                    "n_9mers": 0,
                    "n_hits": 0,
                    "worst_ic50_nM": float("nan"),
                    "worst_peptide": "",
                    "worst_allele": "",
                    "flagged": False,
                }
            )
            continue

        ninemers = tile_9mers(seq)

        # Predict across all alleles
        all_allele_preds = []
        for allele in alleles:
            preds = predictor.predict(
                peptides=ninemers,
                alleles=[allele],
                verbose=0,
            )
            preds["candidate_id"] = cand_id
            preds["position"] = list(range(len(ninemers)))
            all_allele_preds.append(preds)

        combined = pd.concat(all_allele_preds, ignore_index=True)
        all_rows.append(combined)

        # Summarize across all alleles
        hits = combined[combined["affinity"] < IC50_THRESHOLD_NM]
        best_binder = combined.loc[combined["affinity"].idxmin()]
        summary_rows.append(
            {
                "candidate_id": cand_id,
                "sequence_length": len(seq),
                "n_9mers": len(ninemers),
                "n_hits": len(hits),
                "worst_ic50_nM": best_binder["affinity"],
                "worst_peptide": best_binder["peptide"],
                "worst_allele": best_binder["allele"],
                "flagged": len(hits) > 0,
            }
        )
```

- [ ] **Step 4: Update `run_screening.py` to pass through alleles**

In `run_mhc_screening`, update:
```python
def run_mhc_screening(candidates):
    from screening_utils import screen_candidates
    return screen_candidates(candidates)
```
No change needed — `screen_candidates` now defaults to multi-allele.

Update the print statement in `main()` to reflect multi-allele:
```python
    print(f"\n[2/3] Running MHC-I screening ({len(summary_df)} candidates, multi-allele panel) ...")
```

- [ ] **Step 5: Run tests to verify they pass**

Run: `./veron-hsc/bin/python -m pytest tests/test_screening.py -v`
Expected: All 6 tests PASS

- [ ] **Step 6: Commit**

```bash
git add scripts/screening_utils.py scripts/run_screening.py tests/test_screening.py
git commit -m "feat: expand MHC screening to multi-allele panel (HLA-A*01:01, A*02:01, B*07:02)"
```

---

### Task 3: Add data provenance tracking to postprocessor

**Files:**
- Modify: `scripts/postprocess.py:306-381` (add data_source detection + metadata column)
- Modify: `scripts/generate_test_af3_data.py` (add marker file to synthetic zips)
- Test: `tests/test_provenance.py`

- [ ] **Step 1: Write failing tests for provenance detection**

Create `tests/test_provenance.py`:

```python
"""Tests for data provenance detection in postprocessor."""
import json
import os
import pytest


def test_detect_synthetic_data(tmp_path):
    """Directories with .synthetic marker should be detected as synthetic."""
    from scripts.postprocess import detect_data_source
    job_dir = tmp_path / "cd34_vs_test"
    job_dir.mkdir()
    (job_dir / ".synthetic").write_text("generated by generate_test_af3_data.py")
    assert detect_data_source(str(job_dir)) == "synthetic"


def test_detect_real_data(tmp_path):
    """Directories without .synthetic marker should be detected as af3_server."""
    from scripts.postprocess import detect_data_source
    job_dir = tmp_path / "cd34_vs_test"
    job_dir.mkdir()
    assert detect_data_source(str(job_dir)) == "af3_server"


def test_detect_real_data_with_terms(tmp_path):
    """Directories with terms_of_use.md are real AF3 Server results."""
    from scripts.postprocess import detect_data_source
    job_dir = tmp_path / "cd34_vs_test"
    job_dir.mkdir()
    (job_dir / "terms_of_use.md").write_text("AF3 terms")
    assert detect_data_source(str(job_dir)) == "af3_server"
```

- [ ] **Step 2: Run tests to verify they fail**

Run: `./veron-hsc/bin/python -m pytest tests/test_provenance.py -v`
Expected: FAIL with `ImportError: cannot import name 'detect_data_source'`

- [ ] **Step 3: Add `detect_data_source` to `postprocess.py`**

Add this function after the imports section (around line 52):

```python
def detect_data_source(job_dir):
    """
    Detect whether AF3 results are from real server or synthetic generation.

    Returns 'synthetic' if .synthetic marker exists, 'af3_server' otherwise.
    """
    if os.path.exists(os.path.join(job_dir, ".synthetic")):
        return "synthetic"
    return "af3_server"
```

- [ ] **Step 4: Add data_source column to postprocessor output**

In `main()`, in the loop that builds rows (around line 326), add after `cand_id = id_map.get(...)`:
```python
        data_source = detect_data_source(paths["dir"])
```

And add to the row dict:
```python
            "data_source": data_source,
```

Update the CSV header print in the final ranking table to include data_source.

- [ ] **Step 5: Update `generate_test_af3_data.py` to write `.synthetic` marker**

In the `main()` function, after extracting each zip (or inside the zip write block), add a `.synthetic` marker. Since the test generator creates zips that `ingest_results.py` extracts, add the marker inside the zip:

In the `main()` loop, after `zf.writestr(f"fold_{job_name}_model_0.cif", cif_content)`, add:
```python
            zf.writestr(".synthetic", "generated by generate_test_af3_data.py\n")
```

- [ ] **Step 6: Run tests to verify they pass**

Run: `./veron-hsc/bin/python -m pytest tests/test_provenance.py -v`
Expected: 3 tests PASS

- [ ] **Step 7: Commit**

```bash
git add scripts/postprocess.py scripts/generate_test_af3_data.py tests/test_provenance.py
git commit -m "feat: add data provenance tracking (synthetic vs af3_server)"
```

---

### Task 4: Fix `environment.yml`

**Files:**
- Modify: `environment.yml`

- [ ] **Step 1: Write the environment.yml with actual dependencies**

```yaml
name: veron-hsc
channels:
  - conda-forge
  - defaults
dependencies:
  - python=3.11
  - openmm>=8.5
  - pdbfixer
  - biopython>=1.87
  - pandas
  - numpy
  - matplotlib
  - seaborn
  - pip
  - pip:
    - mhcflurry>=2.2
    - torch
```

- [ ] **Step 2: Verify it parses correctly**

Run: `./veron-hsc/bin/python -c "import yaml; yaml.safe_load(open('environment.yml')); print('OK')"` or if pyyaml not available: `./veron-hsc/bin/python -c "import re; text=open('environment.yml').read(); assert 'python=3.11' in text; assert 'openmm' in text; assert 'mhcflurry' in text; print('OK')"`

- [ ] **Step 3: Commit**

```bash
git add environment.yml
git commit -m "fix: populate environment.yml with actual conda/pip dependencies"
```

---

### Task 5: Add pipeline smoke tests

**Files:**
- Create: `tests/test_pipeline_smoke.py`
- Create: `tests/conftest.py`

This task creates end-to-end smoke tests that: generate synthetic data, ingest it, run postprocessing, and verify the output ranking makes sense. Also tests the screening stealth score computation and CIF parser against real AF3 data.

- [ ] **Step 1: Create `tests/conftest.py` with shared fixtures**

```python
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
```

- [ ] **Step 2: Write smoke tests**

Create `tests/test_pipeline_smoke.py`:

```python
"""
Pipeline smoke tests for veron-hsc.

These tests verify the pipeline components work end-to-end using
either real AF3 data (if present) or synthetic test data.
"""
import json
import os

import numpy as np
import pandas as pd
import pytest

PROJECT_ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))


class TestCIFParser:
    """Test CIF parsing against real AF3 Server output (if available)."""

    REAL_CIF = os.path.join(
        PROJECT_ROOT, "results", "af3_outputs", "cd34_vs_ch02_v10a",
        "fold_cd34_vs_ch02_v10a_model_0.cif",
    )

    @pytest.mark.skipif(
        not os.path.exists(REAL_CIF),
        reason="Real AF3 data not available",
    )
    def test_parse_real_cif_has_two_chains(self):
        from scripts.utils import parse_cif_atoms
        chains = parse_cif_atoms(self.REAL_CIF, ca_only=True)
        assert "A" in chains, "Missing receptor chain A"
        assert "B" in chains, "Missing ligand chain B"

    @pytest.mark.skipif(
        not os.path.exists(REAL_CIF),
        reason="Real AF3 data not available",
    )
    def test_parse_real_cif_receptor_length(self):
        from scripts.utils import parse_cif_atoms
        chains = parse_cif_atoms(self.REAL_CIF, ca_only=True)
        # Receptor is 141 residues (original 150-290, renumbered 1-141)
        assert len(chains["A"]) == 141, f"Expected 141 receptor CAs, got {len(chains['A'])}"

    @pytest.mark.skipif(
        not os.path.exists(REAL_CIF),
        reason="Real AF3 data not available",
    )
    def test_parse_real_cif_ligand_length(self):
        from scripts.utils import parse_cif_atoms
        chains = parse_cif_atoms(self.REAL_CIF, ca_only=True)
        # CH02_V10A ligand is 12 residues
        assert len(chains["B"]) == 12, f"Expected 12 ligand CAs, got {len(chains['B'])}"

    @pytest.mark.skipif(
        not os.path.exists(REAL_CIF),
        reason="Real AF3 data not available",
    )
    def test_cif_column_count_validation(self):
        """Verify real CIF ATOM lines have exactly 18 columns."""
        with open(self.REAL_CIF) as f:
            for line in f:
                if line.startswith("ATOM"):
                    parts = line.split()
                    assert len(parts) == 18, (
                        f"Expected 18 columns, got {len(parts)}: {line.strip()}"
                    )
                    break  # Only need to check one line


class TestStealthScoring:
    """Test stealth score computation."""

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
        # Should be sorted descending by stealth_score
        assert list(result["candidate_id"]) == ["GOOD", "MID", "BAD"]


class TestLeadScore:
    """Test the weighted lead score computation."""

    def test_lead_score_formula(self):
        from scripts.postprocess import merge_and_rank
        structural = pd.DataFrame([{
            "candidate_id": "TEST",
            "iptm": 0.80,
            "ligand_plddt": 90.0,
        }])
        # Create a minimal screening CSV
        import tempfile
        with tempfile.NamedTemporaryFile(mode="w", suffix=".csv", delete=False) as f:
            f.write("candidate_id,stealth_score\n")
            f.write("TEST,1.0\n")
            screening_path = f.name

        try:
            result = merge_and_rank(structural, screening_path)
            # lead_score = 0.4 * 1.0 + 0.4 * 0.8 + 0.2 * 0.9 = 0.4 + 0.32 + 0.18 = 0.90
            assert result.iloc[0]["lead_score"] == pytest.approx(0.90, abs=0.01)
        finally:
            os.unlink(screening_path)

    def test_lead_score_ranking_order(self):
        from scripts.postprocess import merge_and_rank
        structural = pd.DataFrame([
            {"candidate_id": "LOW", "iptm": 0.20, "ligand_plddt": 30.0},
            {"candidate_id": "HIGH", "iptm": 0.90, "ligand_plddt": 95.0},
        ])
        import tempfile
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


class TestMWSPMotifDistance:
    """Test MWSP motif distance computation."""

    def test_motif_not_found_returns_nan(self):
        from scripts.postprocess import compute_motif_distance
        # Sequence without MWSP → NaN
        result = compute_motif_distance("/nonexistent.cif", "AAAAAAAA")
        assert np.isnan(result)

    @pytest.mark.skipif(
        not os.path.exists(os.path.join(
            PROJECT_ROOT, "results", "af3_outputs", "cd34_vs_ch02_v10a",
            "fold_cd34_vs_ch02_v10a_model_0.cif",
        )),
        reason="Real AF3 data not available",
    )
    def test_motif_distance_v10a_is_positive(self):
        from scripts.postprocess import compute_motif_distance
        cif = os.path.join(
            PROJECT_ROOT, "results", "af3_outputs", "cd34_vs_ch02_v10a",
            "fold_cd34_vs_ch02_v10a_model_0.cif",
        )
        dist = compute_motif_distance(cif, "THRPPMWSPAWP")
        assert not np.isnan(dist)
        assert dist > 0
        assert dist < 20  # sanity: should be within 20 Angstroms


class TestEndToEnd:
    """End-to-end: generate synthetic → ingest → postprocess → verify ranking."""

    SCREENING_CSV = os.path.join(PROJECT_ROOT, "results", "screening_results.csv")

    @pytest.mark.skipif(
        not os.path.exists(os.path.join(PROJECT_ROOT, "results", "screening_results.csv")),
        reason="screening_results.csv not available",
    )
    def test_screening_csv_has_all_candidates(self):
        df = pd.read_csv(self.SCREENING_CSV)
        assert len(df) == 12  # 11 CH02 variants + 8G12 control

    @pytest.mark.skipif(
        not os.path.exists(os.path.join(PROJECT_ROOT, "results", "veron_prioritized_leads.csv")),
        reason="veron_prioritized_leads.csv not available",
    )
    def test_prioritized_leads_has_lead_score(self):
        csv_path = os.path.join(PROJECT_ROOT, "results", "veron_prioritized_leads.csv")
        df = pd.read_csv(csv_path)
        assert "lead_score" in df.columns
        assert "candidate_id" in df.columns
        assert all(df["lead_score"].notna())

    @pytest.mark.skipif(
        not os.path.exists(os.path.join(PROJECT_ROOT, "results", "veron_prioritized_leads.csv")),
        reason="veron_prioritized_leads.csv not available",
    )
    def test_8g12_control_penalized_for_immunogenicity(self):
        """8G12 scFv should not be #1 despite high iPTM — stealth penalty."""
        csv_path = os.path.join(PROJECT_ROOT, "results", "veron_prioritized_leads.csv")
        df = pd.read_csv(csv_path)
        top = df.iloc[0]
        assert top["candidate_id"] != "8G12_scFv", (
            "8G12 control should not be top lead — stealth score should penalize it"
        )
```

- [ ] **Step 3: Run tests**

Run: `./veron-hsc/bin/python -m pytest tests/ -v`
Expected: All tests pass (some may skip if AF3 data not present)

- [ ] **Step 4: Commit**

```bash
git add tests/conftest.py tests/test_pipeline_smoke.py
git commit -m "test: add pipeline smoke tests (CIF parser, stealth scoring, lead score, e2e)"
```

---

### Task 6: Reconcile AF3 JSON generation paths — delete `af3_generator.py`

**Files:**
- Delete: `scripts/af3_generator.py`
- Delete: `scripts/fix_jsons.py`
- Modify: `scripts/run_screening.py:49,67-83` (ensure it's the canonical path, uses `modelSeeds: []`)
- Modify: `CLAUDE.md` (remove af3_generator.py references)

The standalone `af3_generator.py` is redundant with `run_screening.py`. The `fix_jsons.py` script was a one-time patch for a format mismatch that no longer applies. `consolidate_jsons.py` already forces `modelSeeds: []`, making the `MODEL_SEEDS = [42, ...]` in `run_screening.py` dead code.

- [ ] **Step 1: Update `run_screening.py` to use empty seeds (matching consolidator)**

Change line 49:
```python
MODEL_SEEDS = [42, 123, 256, 789, 1024]
```
to:
```python
MODEL_SEEDS = []  # AF3 Server picks random seeds; consolidate_jsons.py also enforces this
```

Also fix the JSON wrapping to be consistent — `generate_af3_jsons` should wrap in a list (AF3 Server expects `[{...}]`):

In `generate_af3_jsons`, change:
```python
        with open(out_path, "w") as f:
            json.dump(payload, f, indent=2)
```
to:
```python
        with open(out_path, "w") as f:
            json.dump([payload], f, indent=2)
```

- [ ] **Step 2: Delete redundant scripts**

```bash
git rm scripts/af3_generator.py scripts/fix_jsons.py
```

- [ ] **Step 3: Update `CLAUDE.md` pipeline description**

Remove the `af3_generator.py` entry from the "Pipeline Architecture" section. The pipeline list should go:
1. `preprocess.py`
2. `run_screening.py` (handles AF3 JSON generation + MHC screening)  
3. `screening_utils.py` (module, not standalone)
4. `consolidate_jsons.py`
5. `ingest_results.py`
6. `postprocess.py`
7. `generate_test_af3_data.py`

- [ ] **Step 4: Verify nothing imports from deleted files**

Run: `grep -r "af3_generator\|fix_jsons" scripts/ tests/`
Expected: No matches

- [ ] **Step 5: Run all tests to confirm nothing broke**

Run: `./veron-hsc/bin/python -m pytest tests/ -v`
Expected: All tests PASS

- [ ] **Step 6: Commit**

```bash
git add scripts/run_screening.py CLAUDE.md
git commit -m "refactor: reconcile AF3 JSON paths — delete redundant af3_generator.py and fix_jsons.py"
```

---

### Task 7: Remove empty `notebooks/` directory

**Files:**
- Delete: `notebooks/` directory

- [ ] **Step 1: Remove the empty directory**

```bash
rmdir notebooks/
```

- [ ] **Step 2: Remove from `.gitignore` if referenced**

The `.gitignore` references `ipynb_checkpoints` which is fine to keep (generic). No change needed.

- [ ] **Step 3: Update CLAUDE.md if it references notebooks/**

Check `CLAUDE.md` — it mentions `notebooks/ — Jupyter notebooks for analysis`. Remove that line.

- [ ] **Step 4: Commit**

```bash
git rm -r notebooks/ 2>/dev/null; git add CLAUDE.md
git commit -m "chore: remove empty notebooks/ directory"
```

---

### Task 8: Fix README claims about WT HLA-C*07:01 IC50

**Files:**
- Modify: `README.md:15-16`

The README claims WT has HLA-C\*07:01 IC50 = 117 nM, but this value doesn't appear in any data file. The `lead_profile_v10a.py` only profiles V10A. Remove the unsubstantiated comparison.

- [ ] **Step 1: Fix the README table**

Replace the HLA-C\*07:01 row:
```
| HLA-C\*07:01 IC50 | **449 nM** | 117 nM | 3.8x less immunogenic |
```
with:
```
| HLA-C\*07:01 IC50 | **449 nM** | — | Borderline (threshold 500 nM) |
```

Also add a note that the data comes from real AF3 Server results (not synthetic):
After the table, update the paragraph to clarify that structural metrics are from real AF3 Server runs.

- [ ] **Step 2: Commit**

```bash
git add README.md
git commit -m "fix: correct README — remove unsubstantiated WT HLA-C*07:01 IC50 claim"
```

---

### Execution order and dependencies

```
Task 1 (shared utils) ← no deps, do first
Task 2 (multi-allele) ← depends on Task 1 (uses utils imports)
Task 3 (provenance)   ← depends on Task 1 (uses utils imports)
Task 4 (env.yml)      ← independent
Task 5 (smoke tests)  ← depends on Tasks 1-3
Task 6 (reconcile)    ← depends on Task 1
Task 7 (notebooks)    ← independent
Task 8 (README)       ← independent
```

Recommended order: 1 → 2 → 3 → 4 → 5 → 6 → 7 → 8

Tasks 4, 7, 8 are independent and can be parallelized with anything.
