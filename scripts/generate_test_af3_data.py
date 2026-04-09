#!/usr/bin/env python
"""
Generate synthetic AF3 Server output zips for pipeline testing.

Creates one zip per candidate in results/af3_outputs/, each containing:
  - summary_confidences_0.json  (best-seed confidence metrics)
  - fold_model_0.cif            (minimal mmCIF with receptor + ligand chains)

The synthetic metrics vary across candidates to produce a meaningful
ranking when combined with MHC screening scores.
"""

import json
import os
import random
import zipfile

import numpy as np
from Bio import SeqIO
from Bio.PDB import PDBParser

from scripts.utils import AA3TO1

PROJECT_ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
RECEPTOR_PDB = os.path.join(PROJECT_ROOT, "data", "processed", "cd34_relaxed.pdb")
CANDIDATES_FASTA = os.path.join(PROJECT_ROOT, "data", "ligands", "candidates.fasta")
AF3_OUT_DIR = os.path.join(PROJECT_ROOT, "results", "af3_outputs")

# Plausible AF3 metrics per candidate.
# Designed so some CH02 variants rank higher than others after
# combining with stealth scores.  8G12 gets high iPTM (known binder)
# but will be penalised by stealth.
SYNTHETIC_METRICS = {
    "CH02_WT":   {"iptm": 0.82, "ptm": 0.78, "frac_dis": 0.05, "ligand_plddt": 88.3, "ranking_score": 0.80},
    "CH02_T1A":  {"iptm": 0.71, "ptm": 0.70, "frac_dis": 0.08, "ligand_plddt": 81.5, "ranking_score": 0.69},
    "CH02_R3A":  {"iptm": 0.65, "ptm": 0.68, "frac_dis": 0.12, "ligand_plddt": 74.2, "ranking_score": 0.63},
    "CH02_M7A":  {"iptm": 0.44, "ptm": 0.55, "frac_dis": 0.22, "ligand_plddt": 58.1, "ranking_score": 0.42},
    "CH02_W8A":  {"iptm": 0.39, "ptm": 0.50, "frac_dis": 0.28, "ligand_plddt": 52.7, "ranking_score": 0.37},
    "CH02_P5A":  {"iptm": 0.76, "ptm": 0.74, "frac_dis": 0.07, "ligand_plddt": 84.9, "ranking_score": 0.74},
    "CH02_V10A": {"iptm": 0.73, "ptm": 0.72, "frac_dis": 0.09, "ligand_plddt": 82.1, "ranking_score": 0.71},
    "CH02_W11A": {"iptm": 0.68, "ptm": 0.69, "frac_dis": 0.11, "ligand_plddt": 76.8, "ranking_score": 0.66},
    "CH02_V10I": {"iptm": 0.80, "ptm": 0.77, "frac_dis": 0.06, "ligand_plddt": 87.0, "ranking_score": 0.78},
    "CH02_M7L":  {"iptm": 0.78, "ptm": 0.75, "frac_dis": 0.06, "ligand_plddt": 85.5, "ranking_score": 0.76},
    "CH02_P4S":  {"iptm": 0.70, "ptm": 0.71, "frac_dis": 0.10, "ligand_plddt": 79.3, "ranking_score": 0.68},
    "8G12_scFv": {"iptm": 0.91, "ptm": 0.88, "frac_dis": 0.03, "ligand_plddt": 92.4, "ranking_score": 0.90},
}


def get_receptor_atoms():
    """Load receptor CA coordinates for building synthetic complexes."""
    parser = PDBParser(QUIET=True)
    structure = parser.get_structure("rec", RECEPTOR_PDB)
    atoms = []
    for residue in structure.get_residues():
        if residue.get_id()[0] == " " and "CA" in residue:
            c = residue["CA"].get_vector().get_array()
            rn = AA3TO1.get(residue.get_resname(), "X")
            atoms.append((rn, c))
    return atoms


def build_summary_confidences(cand_id, metrics, receptor_len, ligand_len):
    """Build AF3 Server-format summary_confidences JSON."""
    total = receptor_len + ligand_len

    # Per-atom pLDDT: receptor gets ~80, ligand gets the synthetic value
    receptor_plddts = [random.gauss(80.0, 5.0) for _ in range(receptor_len)]
    ligand_plddts = [random.gauss(metrics["ligand_plddt"], 3.0) for _ in range(ligand_len)]
    all_plddts = receptor_plddts + ligand_plddts

    chain_ids = ["A"] * receptor_len + ["B"] * ligand_len

    return {
        "atom_chain_ids": chain_ids,
        "atom_plddts": [round(v, 2) for v in all_plddts],
        "chain_iptm": {"A": round(metrics["ptm"], 4), "B": round(metrics["iptm"], 4)},
        "chain_pair_iptm": {
            "A": {"A": round(metrics["ptm"], 4), "B": round(metrics["iptm"], 4)},
            "B": {"A": round(metrics["iptm"], 4), "B": round(metrics["ptm"] * 0.95, 4)},
        },
        "chain_pair_pae_min": {
            "A": {"A": round(2.0 + random.random(), 2), "B": round(5.0 + random.random() * 10, 2)},
            "B": {"A": round(5.0 + random.random() * 10, 2), "B": round(2.0 + random.random(), 2)},
        },
        "chain_ptm": {"A": round(metrics["ptm"], 4), "B": round(metrics["ptm"] * 0.9, 4)},
        "fraction_disordered": round(metrics["frac_dis"], 4),
        "has_clash": False,
        "iptm": round(metrics["iptm"], 4),
        "ptm": round(metrics["ptm"], 4),
        "ranking_score": round(metrics["ranking_score"], 4),
    }


def build_minimal_cif(cand_id, receptor_atoms, ligand_seq, metrics):
    """
    Build a minimal mmCIF with Chain A (receptor) and Chain B (ligand).

    Places the ligand as an extended chain near the receptor surface.
    For candidates containing MWSP, positions that motif close to the
    receptor surface (distance proportional to iPTM — better binding =
    closer contact).
    """
    lines = []
    lines.append("data_fold_model")
    lines.append("#")
    lines.append("loop_")
    lines.append("_atom_site.group_PDB")
    lines.append("_atom_site.id")
    lines.append("_atom_site.type_symbol")
    lines.append("_atom_site.label_atom_id")
    lines.append("_atom_site.label_comp_id")
    lines.append("_atom_site.label_asym_id")
    lines.append("_atom_site.label_seq_id")
    lines.append("_atom_site.Cartn_x")
    lines.append("_atom_site.Cartn_y")
    lines.append("_atom_site.Cartn_z")
    lines.append("_atom_site.B_iso_or_equiv")
    lines.append("_atom_site.auth_asym_id")

    atom_id = 1

    # Chain A — receptor CA trace
    for i, (resname, coord) in enumerate(receptor_atoms, 1):
        three = [k for k, v in AA3TO1.items() if v == resname]
        comp = three[0] if three else "ALA"
        x, y, z = coord
        bfactor = round(random.gauss(80.0, 5.0), 2)
        lines.append(
            f"ATOM {atom_id:>5d} C CA {comp:>3s} A {i:>4d} "
            f"{x:>8.3f} {y:>8.3f} {z:>8.3f} {bfactor:>6.2f} A"
        )
        atom_id += 1

    # Chain B — ligand CA trace
    # Place near receptor surface; distance depends on iPTM
    rec_coords = np.array([c for _, c in receptor_atoms])
    surface_point = rec_coords[len(rec_coords) // 2]  # mid-chain
    # Offset: high iPTM → closer (3-5 Å), low iPTM → farther (8-15 Å)
    contact_dist = 3.5 + (1.0 - metrics["iptm"]) * 15.0

    # Build extended chain along X from the surface point
    direction = np.array([1.0, 0.3, 0.1])
    direction = direction / np.linalg.norm(direction)

    mwsp_start = ligand_seq.find("MWSP")

    for i, aa in enumerate(ligand_seq):
        three = [k for k, v in AA3TO1.items() if v == aa]
        comp = three[0] if three else "ALA"

        # Position along extended chain (3.8 Å CA-CA spacing)
        pos = surface_point + direction * (contact_dist + i * 3.8)

        # If this residue is part of the MWSP motif, pull it closer
        if mwsp_start != -1 and mwsp_start <= i < mwsp_start + 4:
            toward_surface = (surface_point - pos)
            toward_surface = toward_surface / np.linalg.norm(toward_surface)
            pos = pos + toward_surface * (contact_dist * 0.7)

        x, y, z = pos
        bfactor = round(random.gauss(metrics["ligand_plddt"], 3.0), 2)
        lines.append(
            f"ATOM {atom_id:>5d} C CA {comp:>3s} B {i + 1:>4d} "
            f"{x:>8.3f} {y:>8.3f} {z:>8.3f} {bfactor:>6.2f} B"
        )
        atom_id += 1

    lines.append("#")
    return "\n".join(lines) + "\n"


def main():
    print("=" * 60)
    print("  veron-hsc  |  Synthetic AF3 test data generator")
    print("=" * 60)

    random.seed(42)
    np.random.seed(42)

    receptor_atoms = get_receptor_atoms()
    candidates = {r.id: str(r.seq) for r in SeqIO.parse(CANDIDATES_FASTA, "fasta")}
    os.makedirs(AF3_OUT_DIR, exist_ok=True)

    for cand_id, seq in candidates.items():
        metrics = SYNTHETIC_METRICS[cand_id]
        job_name = f"cd34_vs_{cand_id}"

        # Build confidence JSON
        conf_json = build_summary_confidences(cand_id, metrics, len(receptor_atoms), len(seq))

        # Build minimal CIF
        cif_content = build_minimal_cif(cand_id, receptor_atoms, seq, metrics)

        # Write zip
        zip_path = os.path.join(AF3_OUT_DIR, f"fold_{job_name}.zip")
        with zipfile.ZipFile(zip_path, "w", zipfile.ZIP_DEFLATED) as zf:
            zf.writestr(
                f"fold_{job_name}_summary_confidences_0.json",
                json.dumps(conf_json, indent=2),
            )
            zf.writestr(f"fold_{job_name}_model_0.cif", cif_content)

        print(f"  {zip_path}  (iPTM={metrics['iptm']:.2f}, pLDDT_lig={metrics['ligand_plddt']:.1f})")

    print(f"\nGenerated {len(candidates)} synthetic AF3 result zips.")


if __name__ == "__main__":
    main()
