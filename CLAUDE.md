# CLAUDE.md

This file provides guidance to Claude Code (claude.ai/code) when working with code in this repository.

## Project Overview

LNP-ligand binding pipeline ("veron-hsc") targeting the CD34+ receptor (UniProt P28906) on hematopoietic stem cells. Uses AlphaFold for structure prediction and OpenMM for molecular dynamics.

## Environment Setup

- **Conda environment** at `veron-hsc/` (Python 3.11, ARM64/Apple Silicon)
- Activate: `conda activate ./veron-hsc`
- Key packages: `openmm` 8.5 (conda-forge `_apple` build), `pdbfixer`, `biopython`, `mhcflurry` 2.2 (pip), `pandas`, `torch`, `matplotlib`, `seaborn`
- GPU platform: **OpenCL** on Apple M3 (Metal plugin not yet shipped in packages; script auto-selects Metal > OpenCL > CPU)
- Run scripts with: `./veron-hsc/bin/python scripts/<name>.py`

## Pipeline Architecture

Pipeline scripts in `scripts/`:

1. **`preprocess.py`** — Fetches CD34 AlphaFold structure (AF-P28906-F1, v6), truncates to globular/stalk ectodomain (residues 150-290, removing disordered mucin domain), repairs with PDBFixer, energy-minimizes with AMBER14-all + OBC2 implicit solvent. Output: `data/processed/cd34_relaxed.pdb`
2. **`run_screening.py`** — Master orchestrator: AF3 JSON generation (receptor × each candidate ligand) → multi-allele MHC-I screening → stealth scoring → prioritized CSV at `results/screening_results.csv`
3. **`screening_utils.py`** — MHC-I immunogenicity module. Tiles each candidate into 9-mers, predicts binding across HLA-A\*01:01, A\*02:01, and B\*07:02 via MHCflurry, flags candidates with any IC50 < 500 nM
4. **`consolidate_jsons.py`** — Merges per-candidate AF3 JSONs into a single bulk upload manifest (`results/veron_master_manifest.json`)
5. **`ingest_results.py`** — Monitors `results/af3_outputs/` for new `.zip` files, extracts into named sub-directories, validates presence of `summary_confidences*.json` + `.cif`. Idempotent (skips already-extracted).
6. **`postprocess.py`** — Core quality gate: extracts iPTM/pLDDT/fraction_disordered from AF3 results, computes MWSP motif distance to receptor surface, merges with MHC screening, computes weighted lead score (40% stealth + 40% iPTM + 20% pLDDT), tracks data provenance (synthetic vs af3_server), outputs `results/veron_prioritized_leads.csv` and a Stealth-vs-Affinity scatter plot.
7. **`lead_profile_v10a.py`** — Deep structural + stealth analysis for the V10A lead candidate.
8. **`generate_test_af3_data.py`** — Generates synthetic AF3 Server output zips (summary_confidences + mmCIF) for pipeline validation. Includes `.synthetic` marker for provenance tracking.
9. **`utils.py`** — Shared utilities: AA3TO1 mapping, receptor sequence extraction, mmCIF atom parsing.

## Directory Layout

- `data/raw/` — Downloaded structures (AlphaFold PDBs)
- `data/processed/` — Pipeline outputs (truncated + relaxed PDBs)
- `data/ligands/` — Candidate FASTA (`candidates.fasta`: 11 CH02 variants + 8G12 scFv control)
- `results/af3_inputs/` — Generated AF3 Server JSON inputs (one per candidate)
- `results/af3_outputs/` — AlphaFold 3 prediction outputs
- `results/screening_results.csv` — Prioritized candidates ranked by stealth score
- `results/veron_prioritized_leads.csv` — Final master report with weighted lead scores
- `results/figures/stealth_vs_affinity.png` — Stealth-vs-Affinity scatter plot
- `results/figures/` — Generated visualizations
- `tests/` — pytest smoke tests

## Full pipeline execution

```bash
./veron-hsc/bin/python scripts/preprocess.py           # 1. Fetch + minimize CD34
./veron-hsc/bin/python scripts/run_screening.py         # 2. AF3 JSONs + MHC screening
# ... submit JSONs to AF3 Server, download result zips to results/af3_outputs/ ...
./veron-hsc/bin/python scripts/ingest_results.py        # 3. Extract + validate AF3 results
./veron-hsc/bin/python scripts/postprocess.py           # 4. Quality gate + final ranking
```

## Notes

- OpenMM residue numbering: the relaxed PDB is renumbered 1-141 by OpenMM's writer (original residues 150-290). All downstream scripts (AF3 generator, postprocessor) use this renumbered sequence.
- Glycan removal is not yet implemented; current minimization is protein-backbone only.
- MWSP motif distance is NaN for candidates where alanine-scanning destroyed the motif (M7A, W8A) and for the 8G12 scFv control. This is expected.
- Lead score formula: `0.4 × stealth + 0.4 × iPTM + 0.2 × (pLDDT/100)`
