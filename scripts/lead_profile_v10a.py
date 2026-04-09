#!/usr/bin/env python
"""
Lead Candidate Profile Generator for CH02_V10A.

Performs:
  1. Structural deep-dive — CD34 residues within 5Å of MWSP motif
  2. Binding pose comparison — V10A vs WT (C-terminal approach, contacts)
  3. Global stealth expansion — multi-allele MHCflurry screening
  4. Outputs structured JSON for report generation

Usage:
    python scripts/lead_profile_v10a.py
"""

import json
import os
import sys
from collections import defaultdict

import numpy as np
import pandas as pd
from mhcflurry import Class1PresentationPredictor

from scripts.utils import AA3TO1, parse_cif_atoms

PROJECT_ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
AF3_OUT_DIR = os.path.join(PROJECT_ROOT, "results", "af3_outputs")

V10A_CIF = os.path.join(
    AF3_OUT_DIR, "cd34_vs_ch02_v10a", "fold_cd34_vs_ch02_v10a_model_0.cif"
)
WT_CIF = os.path.join(
    AF3_OUT_DIR, "cd34_vs_ch02_wt", "fold_cd34_vs_ch02_wt_model_0.cif"
)

V10A_SEQ = "THRPPMWSPAWP"   # V10A: Val→Ala at position 10
WT_SEQ   = "THRPPMWSPVWP"   # Wild type

BINDING_MOTIF = "MWSP"
CONTACT_CUTOFF = 5.0  # Angstroms

# Global HLA panel for stealth expansion
GLOBAL_HLA_PANEL = [
    "HLA-A*01:01",  # original screening allele
    "HLA-A*02:01",  # most common worldwide
    "HLA-B*07:02",  # common Caucasian allele
    "HLA-C*07:01",  # common Class I allele
]
IC50_THRESHOLD = 500  # nM


# ──────────────────────────────────────────────
# 1. Binding pocket: receptor residues within 5Å of MWSP
# ──────────────────────────────────────────────
def find_binding_pocket(cif_path, ligand_seq, motif=BINDING_MOTIF, cutoff=CONTACT_CUTOFF):
    """
    Identify all Chain A (receptor) residues with any heavy atom
    within `cutoff` Å of any heavy atom in the MWSP motif on Chain B.

    Returns:
      - pocket_residues: list of (seq_id, resname, one_letter, min_dist)
      - motif_atoms: number of motif atoms used
      - receptor_atoms: number of receptor atoms used
    """
    motif_pos = ligand_seq.find(motif)
    if motif_pos == -1:
        return [], 0, 0

    # 1-indexed residue IDs for the motif residues
    motif_ids = set(range(motif_pos + 1, motif_pos + 1 + len(motif)))

    chains = parse_cif_atoms(cif_path, ca_only=False)
    if "A" not in chains or "B" not in chains:
        return [], 0, 0

    # Gather motif atom coordinates (Chain B, within motif residues)
    motif_atoms = [(a[3], a[4]) for a in chains["B"] if a[3] in motif_ids]
    if not motif_atoms:
        return [], 0, 0

    motif_coords = np.array([c for _, c in motif_atoms])

    # For each receptor residue, find minimum distance to any motif atom
    receptor_residues = defaultdict(list)
    for atom in chains["A"]:
        _, _, resname, seq_id, xyz = atom
        receptor_residues[seq_id].append((resname, xyz))

    pocket = []
    for seq_id in sorted(receptor_residues.keys()):
        resname = receptor_residues[seq_id][0][0]
        res_coords = np.array([xyz for _, xyz in receptor_residues[seq_id]])

        # Pairwise distances: res_atoms x motif_atoms
        diffs = res_coords[:, np.newaxis, :] - motif_coords[np.newaxis, :, :]
        dists = np.linalg.norm(diffs, axis=2)
        min_dist = float(dists.min())

        if min_dist <= cutoff:
            one_letter = AA3TO1.get(resname, "X")
            pocket.append((seq_id, resname, one_letter, round(min_dist, 2)))

    return pocket, len(motif_atoms), len(chains["A"])


# ──────────────────────────────────────────────
# 2. Binding pose comparison — V10A vs WT
# ──────────────────────────────────────────────
def per_residue_contacts(cif_path, ligand_seq, cutoff=CONTACT_CUTOFF):
    """
    For each ligand residue (Chain B), compute minimum distance to
    the receptor surface (any Chain A heavy atom).

    Returns list of (seq_id, resname, one_letter, min_dist_to_receptor).
    """
    chains = parse_cif_atoms(cif_path, ca_only=False)
    if "A" not in chains or "B" not in chains:
        return []

    receptor_coords = np.array([a[4] for a in chains["A"]])

    # Group ligand atoms by residue
    lig_residues = defaultdict(list)
    for atom in chains["B"]:
        _, _, resname, seq_id, xyz = atom
        lig_residues[seq_id].append((resname, xyz))

    contacts = []
    for seq_id in sorted(lig_residues.keys()):
        resname = lig_residues[seq_id][0][0]
        res_coords = np.array([xyz for _, xyz in lig_residues[seq_id]])

        diffs = res_coords[:, np.newaxis, :] - receptor_coords[np.newaxis, :, :]
        dists = np.linalg.norm(diffs, axis=2)
        min_dist = float(dists.min())

        one_letter = AA3TO1.get(resname, "X")
        contacts.append((seq_id, resname, one_letter, round(min_dist, 2)))

    return contacts


def compare_poses(v10a_cif, wt_cif, v10a_seq, wt_seq):
    """Compare per-residue receptor contact distances for V10A vs WT."""
    v10a_contacts = per_residue_contacts(v10a_cif, v10a_seq)
    wt_contacts = per_residue_contacts(wt_cif, wt_seq)
    return v10a_contacts, wt_contacts


# ──────────────────────────────────────────────
# 3. Global stealth — multi-allele MHCflurry
# ──────────────────────────────────────────────
def global_stealth_screen(sequence, alleles, threshold=IC50_THRESHOLD):
    """
    Run MHCflurry Class1PresentationPredictor against multiple alleles.

    Returns per-allele summary and overall stealth verdict.
    """
    if len(sequence) < 9:
        return [], True

    ninemers = [sequence[i:i+9] for i in range(len(sequence) - 9 + 1)]
    predictor = Class1PresentationPredictor.load()

    allele_results = []
    any_hit = False

    for allele in alleles:
        preds = predictor.predict(
            peptides=ninemers,
            alleles=[allele],
            verbose=0,
        )
        hits = preds[preds["affinity"] < threshold]
        best = preds.loc[preds["affinity"].idxmin()]

        allele_results.append({
            "allele": allele,
            "n_9mers": len(ninemers),
            "n_hits": len(hits),
            "worst_ic50_nM": round(float(best["affinity"]), 2),
            "worst_peptide": best["peptide"],
            "flagged": len(hits) > 0,
        })

        if len(hits) > 0:
            any_hit = True

    return allele_results, not any_hit


# ──────────────────────────────────────────────
# Main
# ──────────────────────────────────────────────
def main():
    print("=" * 60)
    print("  CH02_V10A Lead Candidate Profile")
    print("=" * 60)

    results = {}

    # ── 1. Binding Pocket ──
    print("\n[1/3] Identifying CD34 binding pocket (Chain A residues within 5Å of MWSP) ...")
    pocket, n_motif_atoms, n_receptor_atoms = find_binding_pocket(V10A_CIF, V10A_SEQ)
    print(f"      Motif atoms: {n_motif_atoms}, Receptor atoms: {n_receptor_atoms}")
    print(f"      Found {len(pocket)} receptor residues in binding pocket:\n")
    print(f"      {'Res#':<6s} {'AA3':<5s} {'AA1':<4s} {'Min Dist (Å)':<12s}")
    print(f"      {'-'*30}")
    for seq_id, resname, aa1, dist in pocket:
        print(f"      {seq_id:<6d} {resname:<5s} {aa1:<4s} {dist:<12.2f}")

    results["binding_pocket"] = [
        {"seq_id": s, "resname": r, "one_letter": o, "min_dist_A": d}
        for s, r, o, d in pocket
    ]

    # ── 2. Pose Comparison ──
    print("\n[2/3] Comparing V10A vs WT per-residue contact distances ...")
    v10a_contacts, wt_contacts = compare_poses(V10A_CIF, WT_CIF, V10A_SEQ, WT_SEQ)

    print(f"\n      {'Pos':<5s} {'V10A':>10s}  {'WT':>10s}  {'Δ(Å)':>8s}  {'Residue'}")
    print(f"      {'-'*50}")
    deltas = []
    for v, w in zip(v10a_contacts, wt_contacts):
        v_id, _, v_aa, v_dist = v
        w_id, _, w_aa, w_dist = w
        delta = v_dist - w_dist
        deltas.append(delta)
        marker = " ←closer" if delta < -0.5 else (" →farther" if delta > 0.5 else "")
        label = f"{v_aa}{v_id}" + (f" (V10A: A, WT: V)" if v_id == 10 else "")
        print(f"      {v_id:<5d} {v_dist:>10.2f}  {w_dist:>10.2f}  {delta:>+8.2f}  {label}{marker}")

    results["pose_comparison"] = {
        "v10a": [{"pos": s, "aa": o, "min_dist_A": d} for s, _, o, d in v10a_contacts],
        "wt":   [{"pos": s, "aa": o, "min_dist_A": d} for s, _, o, d in wt_contacts],
    }

    # C-terminus analysis (last 3 residues)
    print(f"\n      C-terminal approach (last 3 residues):")
    for i in range(-3, 0):
        v = v10a_contacts[i]
        w = wt_contacts[i]
        delta = v[3] - w[3]
        print(f"        {v[2]}{v[0]}: V10A={v[3]:.2f}Å  WT={w[3]:.2f}Å  Δ={delta:+.2f}Å")

    # ── 3. Global Stealth ──
    print(f"\n[3/3] Running global HLA panel stealth screen ...")
    print(f"      Panel: {', '.join(GLOBAL_HLA_PANEL)}")
    allele_results, stealth_pass = global_stealth_screen(V10A_SEQ, GLOBAL_HLA_PANEL)

    print(f"\n      {'Allele':<16s} {'Hits':>5s} {'Worst IC50 (nM)':>16s} {'Peptide':<12s} {'Flag'}")
    print(f"      {'-'*60}")
    for ar in allele_results:
        flag = "FLAG" if ar["flagged"] else "CLEAN"
        print(f"      {ar['allele']:<16s} {ar['n_hits']:>5d} {ar['worst_ic50_nM']:>16.2f} {ar['worst_peptide']:<12s} {flag}")

    verdict = "GLOBAL STEALTH CONFIRMED" if stealth_pass else "STEALTH BREACH DETECTED"
    print(f"\n      Verdict: {verdict}")

    results["global_stealth"] = {
        "allele_results": allele_results,
        "overall_stealth": stealth_pass,
    }

    # ── Save structured output ──
    out_json = os.path.join(PROJECT_ROOT, "results", "leads", "CH02_V10A_analysis.json")
    os.makedirs(os.path.dirname(out_json), exist_ok=True)
    with open(out_json, "w") as f:
        json.dump(results, f, indent=2)
    print(f"\n      Analysis JSON saved to {out_json}")
    print("\nDone.")


if __name__ == "__main__":
    main()
