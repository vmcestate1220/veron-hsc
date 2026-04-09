#!/usr/bin/env python
"""
Post-processing quality gate for the veron-hsc pipeline.

For each ingested AF3 result:
  1. Best-seed selection — pick seed with highest ranking_score
  2. Confidence extraction — iptm, fraction_disordered, ranking_score
  3. Ligand pLDDT — mean pLDDT for Chain B atoms (from full_data JSON)
  4. Spatial validation — minimum distance (Å) between the ligand's
     MWSP motif and the CD34 receptor surface (Chain A)
  5. Merge with MHC screening results
  6. Compute weighted "Actionable Lead" score:
       40% stealth + 40% iPTM + 20% normalised pLDDT
  7. Generate Stealth-vs-Affinity scatter plot

Handles:
  - Lowercase AF3 Server directory names → uppercase candidate IDs
  - Renumbered receptor (residues 1-141, original 150-290)
  - Real mmCIF format (18 columns) from AF3 Server
  - pLDDT stored in full_data JSONs (not summary_confidences)

Usage:
    python scripts/postprocess.py
"""

import glob
import json
import os

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import seaborn as sns
from Bio import SeqIO
from scripts.utils import parse_cif_atoms

# ──────────────────────────────────────────────
# Paths
# ──────────────────────────────────────────────
PROJECT_ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
AF3_OUT_DIR = os.path.join(PROJECT_ROOT, "results", "af3_outputs")
SCREENING_CSV = os.path.join(PROJECT_ROOT, "results", "screening_results.csv")
CANDIDATES_FASTA = os.path.join(PROJECT_ROOT, "data", "ligands", "candidates.fasta")
OUTPUT_CSV = os.path.join(PROJECT_ROOT, "results", "veron_prioritized_leads.csv")
FIGURES_DIR = os.path.join(PROJECT_ROOT, "results", "figures")

BINDING_MOTIF = "MWSP"
W_STEALTH = 0.40
W_IPTM = 0.40
W_PLDDT = 0.20


# ──────────────────────────────────────────────
# Candidate ID mapping (lowercase dir → FASTA ID)
# ──────────────────────────────────────────────
def build_id_map(fasta_path):
    """
    Build a mapping from lowercase dir names (as AF3 Server produces)
    to the original FASTA candidate IDs.
      cd34_vs_ch02_wt → CH02_WT
      cd34_vs_8g12_scfv → 8G12_scFv
    """
    id_map = {}
    for r in SeqIO.parse(fasta_path, "fasta"):
        dir_name = f"cd34_vs_{r.id}".lower()
        id_map[dir_name] = r.id
    return id_map


# ──────────────────────────────────────────────
# 1. Discover results + select best seed
# ──────────────────────────────────────────────
def find_result_dirs(af3_dir):
    """
    Find extracted AF3 result directories.  For each, identify the
    best seed (highest ranking_score) and return paths to its files.
    """
    results = {}
    for entry in sorted(os.listdir(af3_dir)):
        full = os.path.join(af3_dir, entry)
        if not os.path.isdir(full):
            continue
        conf_files = sorted(glob.glob(os.path.join(full, "*summary_confidences*.json")))
        if not conf_files:
            continue

        # Find best seed by ranking_score
        best_seed = 0
        best_ranking = -1
        for cf in conf_files:
            with open(cf) as f:
                data = json.load(f)
            rs = data.get("ranking_score", 0)
            if rs > best_ranking:
                best_ranking = rs
                best_seed_file = cf
                # Extract seed number from filename
                basename = os.path.basename(cf)
                # ..._summary_confidences_N.json → N
                seed_str = basename.rsplit("_", 1)[1].replace(".json", "")
                best_seed = int(seed_str)
                best_ranking = rs

        # Find corresponding model CIF and full_data JSON for best seed
        cif_pattern = os.path.join(full, f"*model_{best_seed}.cif")
        fd_pattern = os.path.join(full, f"*full_data_{best_seed}.json")
        cif_files = glob.glob(cif_pattern)
        fd_files = glob.glob(fd_pattern)

        if not cif_files:
            continue

        results[entry] = {
            "dir": full,
            "conf_json": best_seed_file,
            "full_data_json": fd_files[0] if fd_files else None,
            "cif_file": cif_files[0],
            "best_seed": best_seed,
        }
    return results


# ──────────────────────────────────────────────
# 2. Confidence extraction
# ──────────────────────────────────────────────
def extract_confidences(conf_json_path):
    with open(conf_json_path) as f:
        data = json.load(f)
    return {
        "iptm": data.get("iptm", float("nan")),
        "ptm": data.get("ptm", float("nan")),
        "fraction_disordered": data.get("fraction_disordered", float("nan")),
        "ranking_score": data.get("ranking_score", float("nan")),
        "has_clash": data.get("has_clash", 0.0),
    }


# ──────────────────────────────────────────────
# 3. Ligand pLDDT (Chain B) — from full_data JSON
# ──────────────────────────────────────────────
def extract_ligand_plddt(full_data_path):
    """Compute mean pLDDT for Chain B atoms from full_data JSON."""
    if not full_data_path or not os.path.exists(full_data_path):
        return float("nan")

    with open(full_data_path) as f:
        data = json.load(f)

    chain_ids = data.get("atom_chain_ids", [])
    plddts = data.get("atom_plddts", [])

    ligand_plddts = [p for c, p in zip(chain_ids, plddts) if c == "B"]
    if not ligand_plddts:
        return float("nan")
    return float(np.mean(ligand_plddts))


# ──────────────────────────────────────────────
# 4. Spatial validation — MWSP motif distance
# ──────────────────────────────────────────────
def compute_motif_distance(cif_path, ligand_seq):
    """
    Compute minimum CA-CA distance (Å) between the MWSP motif on the
    ligand (Chain B) and the CD34 receptor (Chain A).
    """
    motif_pos = ligand_seq.find(BINDING_MOTIF)
    if motif_pos == -1:
        return float("nan")

    # 1-indexed residue IDs for the motif
    motif_ids = set(range(motif_pos + 1, motif_pos + 1 + len(BINDING_MOTIF)))

    chains = parse_cif_atoms(cif_path, ca_only=True)
    if "A" not in chains or "B" not in chains:
        return float("nan")

    receptor_coords = np.array([c for _, _, _, _, c in chains["A"]])
    motif_coords = np.array(
        [c for _, _, _, sid, c in chains["B"] if sid in motif_ids
    ])

    if len(motif_coords) == 0:
        return float("nan")

    diffs = motif_coords[:, np.newaxis, :] - receptor_coords[np.newaxis, :, :]
    dists = np.linalg.norm(diffs, axis=2)
    return float(dists.min())


# ──────────────────────────────────────────────
# 5. Merge with screening and compute lead score
# ──────────────────────────────────────────────
def merge_and_rank(structural_df, screening_csv):
    screening = pd.read_csv(screening_csv)
    merged = pd.merge(structural_df, screening, on="candidate_id", how="left")
    merged["plddt_norm"] = merged["ligand_plddt"].clip(0, 100) / 100.0
    merged["lead_score"] = (
        W_STEALTH * merged["stealth_score"]
        + W_IPTM * merged["iptm"]
        + W_PLDDT * merged["plddt_norm"]
    )
    return merged.sort_values("lead_score", ascending=False)


# ──────────────────────────────────────────────
# 6. Visualization
# ──────────────────────────────────────────────
def plot_stealth_vs_affinity(df, output_dir):
    os.makedirs(output_dir, exist_ok=True)
    fig, ax = plt.subplots(figsize=(10, 7))
    sns.set_style("whitegrid")

    ch02 = df[df["candidate_id"].str.startswith("CH02")]
    ctrl = df[~df["candidate_id"].str.startswith("CH02")]

    # Clamp point sizes so small pLDDT values are still visible
    size_scale = ch02["ligand_plddt"].clip(20, 100) * 2

    scatter = ax.scatter(
        ch02["stealth_score"], ch02["iptm"],
        s=size_scale,
        c=ch02["lead_score"], cmap="RdYlGn", vmin=0.2, vmax=0.8,
        edgecolors="black", linewidth=0.8, zorder=3,
        label="CH02 variants",
    )
    if not ctrl.empty:
        ax.scatter(
            ctrl["stealth_score"], ctrl["iptm"],
            s=ctrl["ligand_plddt"].clip(20, 100) * 2,
            c="royalblue", marker="D", edgecolors="black", linewidth=0.8,
            zorder=3, label="8G12 scFv (control)",
        )

    for _, row in df.iterrows():
        label = row["candidate_id"].replace("CH02_", "")
        ax.annotate(
            label,
            (row["stealth_score"], row["iptm"]),
            textcoords="offset points", xytext=(6, 6),
            fontsize=8, fontweight="bold",
        )

    ax.set_xlabel("Stealth Score (1 = no immunogenic 9-mers)", fontsize=12)
    ax.set_ylabel("iPTM (interface predicted TM-score)", fontsize=12)
    ax.set_title(
        "veron-hsc: Stealth vs Affinity Trade-off (Real AF3 Data)",
        fontsize=14, fontweight="bold",
    )
    ax.legend(loc="best", fontsize=10)

    sm = plt.cm.ScalarMappable(cmap="RdYlGn", norm=plt.Normalize(0.2, 0.8))
    sm.set_array([])
    cbar = plt.colorbar(sm, ax=ax, shrink=0.6, pad=0.02)
    cbar.set_label("Lead Score", fontsize=10)

    fig.tight_layout()
    out_path = os.path.join(output_dir, "stealth_vs_affinity.png")
    fig.savefig(out_path, dpi=150)
    plt.close(fig)
    return out_path


# ──────────────────────────────────────────────
# Main
# ──────────────────────────────────────────────
def main():
    print("=" * 60)
    print("  veron-hsc  |  Post-processing quality gate (real AF3)")
    print("=" * 60)

    # Build lowercase dir → FASTA ID mapping
    id_map = build_id_map(CANDIDATES_FASTA)
    cand_seqs = {r.id: str(r.seq) for r in SeqIO.parse(CANDIDATES_FASTA, "fasta")}

    # 1. Discover results
    result_dirs = find_result_dirs(AF3_OUT_DIR)
    print(f"\n[1/5] Found {len(result_dirs)} ingested AF3 results")

    if not result_dirs:
        print("      No results to process. Run ingest_results.py first.")
        return

    # 2. Extract structural metrics
    print("[2/5] Extracting confidence metrics (best seed) + ligand pLDDT ...")
    rows = []
    for dir_name, paths in result_dirs.items():
        cand_id = id_map.get(dir_name, dir_name)
        conf = extract_confidences(paths["conf_json"])
        lig_plddt = extract_ligand_plddt(paths["full_data_json"])
        lig_seq = cand_seqs.get(cand_id, "")
        motif_dist = compute_motif_distance(paths["cif_file"], lig_seq)

        row = {
            "candidate_id": cand_id,
            "iptm": conf["iptm"],
            "ptm": conf["ptm"],
            "fraction_disordered": conf["fraction_disordered"],
            "ranking_score": conf["ranking_score"],
            "has_clash": conf["has_clash"],
            "ligand_plddt": round(lig_plddt, 2),
            "mwsp_min_dist_A": round(motif_dist, 2) if not np.isnan(motif_dist) else float("nan"),
            "best_seed": paths["best_seed"],
        }
        rows.append(row)

        dist_str = f"{motif_dist:.1f}Å" if not np.isnan(motif_dist) else "N/A"
        print(
            f"      {cand_id:18s}  seed={paths['best_seed']}  "
            f"iPTM={conf['iptm']:.3f}  pLDDT_lig={lig_plddt:.1f}  "
            f"MWSP={dist_str}"
        )

    structural_df = pd.DataFrame(rows)

    # 3. Merge with MHC screening
    print(f"\n[3/5] Merging with MHC screening results ...")
    ranked = merge_and_rank(structural_df, SCREENING_CSV)
    ranked.to_csv(OUTPUT_CSV, index=False, float_format="%.4f")
    print(f"      Master report saved to {OUTPUT_CSV}")

    # 4. Visualization
    print(f"\n[4/5] Generating Stealth vs Affinity plot ...")
    fig_path = plot_stealth_vs_affinity(ranked, FIGURES_DIR)
    print(f"      Plot saved to {fig_path}")

    # 5. Print top leads
    print(f"\n[5/5] Final Actionable Lead Ranking")
    print("=" * 70)
    print(f"  {'Rank':<5s} {'Candidate':<18s} {'Lead':>6s} {'Stealth':>8s} "
          f"{'iPTM':>6s} {'pLDDT':>6s} {'MWSP Å':>7s} {'MHC':>6s}")
    print("-" * 70)
    for rank, (_, row) in enumerate(ranked.iterrows(), 1):
        mhc_flag = "CLEAN" if not row.get("flagged", False) else "FLAG"
        dist_str = f"{row['mwsp_min_dist_A']:.1f}" if not pd.isna(row["mwsp_min_dist_A"]) else "N/A"
        print(
            f"  {rank:<5d} {row['candidate_id']:<18s} {row['lead_score']:>6.3f} "
            f"{row['stealth_score']:>8.3f} {row['iptm']:>6.3f} "
            f"{row['ligand_plddt']:>6.1f} {dist_str:>7s} {mhc_flag:>6s}"
        )

    print("\nDone.")


if __name__ == "__main__":
    main()
