#!/usr/bin/env python
"""
Structural preprocessing for the veron-hsc pipeline.

Fetches the CD34 ectodomain (UniProt P28906), truncates to the
globular/stalk region (residues 150-290), repairs missing atoms,
adds hydrogens, and performs energy minimization with AMBER14-all.

Usage:
    python scripts/preprocess.py
"""

import os
import sys
import urllib.request

from Bio.PDB import PDBParser, PDBIO, Select
from pdbfixer import PDBFixer
from openmm import Platform, LangevinMiddleIntegrator
from openmm import unit as openmm_unit
from openmm.app import (
    PDBFile,
    ForceField,
    Modeller,
    HBonds,
    NoCutoff,
    Simulation,
)

# ──────────────────────────────────────────────
# Paths
# ──────────────────────────────────────────────
PROJECT_ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
RAW_DIR = os.path.join(PROJECT_ROOT, "data", "raw")
PROCESSED_DIR = os.path.join(PROJECT_ROOT, "data", "processed")
RAW_PDB = os.path.join(RAW_DIR, "AF-P28906-F1-model_v6.pdb")
TRUNCATED_PDB = os.path.join(PROCESSED_DIR, "cd34_truncated.pdb")
OUTPUT_PDB = os.path.join(PROCESSED_DIR, "cd34_relaxed.pdb")

ALPHAFOLD_URL = (
    "https://alphafold.ebi.ac.uk/files/AF-P28906-F1-model_v6.pdb"
)

# Globular/stalk ectodomain boundaries (1-indexed, inclusive)
RES_START = 150
RES_END = 290


# ──────────────────────────────────────────────
# 1. Fetch the structure
# ──────────────────────────────────────────────
def fetch_structure():
    """Download AlphaFold-predicted CD34 structure if not cached."""
    if os.path.exists(RAW_PDB):
        print(f"[fetch] Using cached structure: {RAW_PDB}")
        return
    os.makedirs(RAW_DIR, exist_ok=True)
    print(f"[fetch] Downloading CD34 from AlphaFold DB ...")
    urllib.request.urlretrieve(ALPHAFOLD_URL, RAW_PDB)
    print(f"[fetch] Saved to {RAW_PDB}")


# ──────────────────────────────────────────────
# 2. Truncate to globular/stalk region
# ──────────────────────────────────────────────
class ResidueRangeSelect(Select):
    """Accept only residues within [RES_START, RES_END]."""

    def accept_residue(self, residue):
        resseq = residue.get_id()[1]
        return RES_START <= resseq <= RES_END


def truncate_structure():
    """Keep only residues 150-290 of chain A."""
    parser = PDBParser(QUIET=True)
    structure = parser.get_structure("cd34", RAW_PDB)

    os.makedirs(PROCESSED_DIR, exist_ok=True)

    io = PDBIO()
    io.set_structure(structure)
    io.save(TRUNCATED_PDB, ResidueRangeSelect())

    # Verify residue count
    trunc = parser.get_structure("cd34_trunc", TRUNCATED_PDB)
    n_res = sum(1 for _ in trunc.get_residues())
    print(f"[truncate] Kept {n_res} residues (target range {RES_START}-{RES_END})")
    return TRUNCATED_PDB


# ──────────────────────────────────────────────
# 3. PDBFixer: repair missing atoms + hydrogens
# ──────────────────────────────────────────────
def fix_structure(pdb_path):
    """Add missing heavy atoms and hydrogens at pH 7.0."""
    fixer = PDBFixer(filename=pdb_path)
    fixer.findMissingResidues()
    fixer.findMissingAtoms()
    fixer.addMissingAtoms()
    fixer.addMissingHydrogens(pH=7.0)

    n_atoms = fixer.topology.getNumAtoms()
    print(f"[pdbfixer] Structure now has {n_atoms} atoms (with hydrogens)")
    return fixer


# ──────────────────────────────────────────────
# 4. Select best available GPU platform
# ──────────────────────────────────────────────
def get_best_platform():
    """Prefer Metal > OpenCL > CPU > Reference."""
    preference = ["Metal", "OpenCL", "CPU"]
    available = {
        Platform.getPlatform(i).getName(): Platform.getPlatform(i)
        for i in range(Platform.getNumPlatforms())
    }
    for name in preference:
        if name in available:
            plat = available[name]
            print(f"[platform] Selected: {name} (speed factor {plat.getSpeed()})")
            return plat
    plat = available["Reference"]
    print(f"[platform] Falling back to Reference (speed factor {plat.getSpeed()})")
    return plat


# ──────────────────────────────────────────────
# 5. Energy minimization
# ──────────────────────────────────────────────
def minimize(fixer):
    """AMBER14-all minimization with implicit solvent (OBC2)."""
    forcefield = ForceField("amber14-all.xml", "implicit/obc2.xml")

    modeller = Modeller(fixer.topology, fixer.positions)
    system = forcefield.createSystem(
        modeller.topology,
        nonbondedMethod=NoCutoff,
        constraints=HBonds,
    )

    integrator = LangevinMiddleIntegrator(
        300 * openmm_unit.kelvin,
        1.0 / openmm_unit.picosecond,
        0.002 * openmm_unit.picoseconds,
    )

    platform = get_best_platform()
    simulation = Simulation(modeller.topology, system, integrator, platform)
    simulation.context.setPositions(modeller.positions)

    # Report initial energy
    state = simulation.context.getState(getEnergy=True)
    e_before = state.getPotentialEnergy().value_in_unit(
        openmm_unit.kilojoules_per_mole
    )
    print(f"[minimize] Initial PE: {e_before:.1f} kJ/mol")

    # Minimize
    simulation.minimizeEnergy(
        tolerance=10 * openmm_unit.kilojoules_per_mole / openmm_unit.nanometer
    )

    state = simulation.context.getState(getEnergy=True, getPositions=True)
    e_after = state.getPotentialEnergy().value_in_unit(
        openmm_unit.kilojoules_per_mole
    )
    print(
        f"[minimize] Final PE:   {e_after:.1f} kJ/mol"
        f"  (delta = {e_after - e_before:.1f})"
    )

    # Write output
    os.makedirs(PROCESSED_DIR, exist_ok=True)
    with open(OUTPUT_PDB, "w") as f:
        PDBFile.writeFile(simulation.topology, state.getPositions(), f)
    print(f"[minimize] Relaxed structure saved to {OUTPUT_PDB}")


# ──────────────────────────────────────────────
# Main
# ──────────────────────────────────────────────
def main():
    print("=" * 56)
    print("  veron-hsc  |  CD34 structural preprocessing")
    print("=" * 56)

    fetch_structure()
    truncated_path = truncate_structure()
    fixer = fix_structure(truncated_path)
    minimize(fixer)

    print("\nDone.")


if __name__ == "__main__":
    main()
