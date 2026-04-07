#!/usr/bin/env python
"""
AF3 JSON input generator for the veron-hsc pipeline.

Reads the receptor sequence from the preprocessed CD34 PDB and pairs it
with every candidate ligand in the FASTA library.  Produces one JSON per
candidate, each with 5 model seeds for structural convergence.

Usage:
    python scripts/af3_generator.py
"""

import json
import os

from Bio import SeqIO
from Bio.PDB import PDBParser

# ──────────────────────────────────────────────
# Paths
# ──────────────────────────────────────────────
PROJECT_ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
RECEPTOR_PDB = os.path.join(PROJECT_ROOT, "data", "processed", "cd34_relaxed.pdb")
CANDIDATES_FASTA = os.path.join(PROJECT_ROOT, "data", "ligands", "candidates.fasta")
OUTPUT_DIR = os.path.join(PROJECT_ROOT, "results", "af3_inputs")

MODEL_SEEDS = []  # empty = server picks a random seed (reliable for bulk uploads)

AA3TO1 = {
    "ALA": "A", "ARG": "R", "ASN": "N", "ASP": "D", "CYS": "C",
    "GLN": "Q", "GLU": "E", "GLY": "G", "HIS": "H", "ILE": "I",
    "LEU": "L", "LYS": "K", "MET": "M", "PHE": "F", "PRO": "P",
    "SER": "S", "THR": "T", "TRP": "W", "TYR": "Y", "VAL": "V",
}


def extract_receptor_sequence(pdb_path):
    """Extract 1-letter amino acid sequence from a PDB file."""
    parser = PDBParser(QUIET=True)
    structure = parser.get_structure("receptor", pdb_path)
    seq = ""
    for residue in structure.get_residues():
        if residue.get_id()[0] == " ":
            seq += AA3TO1.get(residue.get_resname(), "X")
    return seq


def build_af3_json(job_name, receptor_seq, ligand_seq, seeds):
    """
    Build an AlphaFold 3 Server-compatible JSON input.

    Format follows the AF3 API specification:
      - Two protein chains (receptor + ligand)
      - modelSeeds for reproducibility / convergence assessment
    """
    return {
        "name": job_name,
        "modelSeeds": seeds,
        "sequences": [
            {"proteinChain": {"sequence": receptor_seq, "count": 1}},
            {"proteinChain": {"sequence": ligand_seq, "count": 1}},
        ],
    }


def main():
    print("=" * 56)
    print("  veron-hsc  |  AF3 input generation")
    print("=" * 56)

    # 1. Extract receptor sequence
    receptor_seq = extract_receptor_sequence(RECEPTOR_PDB)
    print(f"[receptor] CD34 ectodomain: {len(receptor_seq)} residues")

    # 2. Load candidates
    candidates = list(SeqIO.parse(CANDIDATES_FASTA, "fasta"))
    print(f"[candidates] Loaded {len(candidates)} sequences from {CANDIDATES_FASTA}")

    # 3. Generate one JSON per candidate
    os.makedirs(OUTPUT_DIR, exist_ok=True)
    generated = []
    for record in candidates:
        candidate_id = record.id
        ligand_seq = str(record.seq)
        job_name = f"cd34_vs_{candidate_id}"

        payload = build_af3_json(job_name, receptor_seq, ligand_seq, MODEL_SEEDS)

        out_path = os.path.join(OUTPUT_DIR, f"{job_name}.json")
        with open(out_path, "w") as f:
            json.dump([payload], f, indent=2)

        generated.append(out_path)
        print(f"  -> {os.path.basename(out_path)}  (ligand: {len(ligand_seq)} aa, seeds: {len(MODEL_SEEDS)})")

    print(f"\n[af3_gen] Wrote {len(generated)} JSON files to {OUTPUT_DIR}")


if __name__ == "__main__":
    main()
