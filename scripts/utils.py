"""
Shared utilities for the veron-hsc pipeline.

Consolidates code previously duplicated across multiple scripts:
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
