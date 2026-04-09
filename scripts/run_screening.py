#!/usr/bin/env python
"""
Master screening orchestrator for the veron-hsc pipeline.

Flow:
  1. Generate AF3 input JSONs (receptor × each candidate ligand)
  2. Run MHC-I immunogenicity screening (HLA-A*01:01, 9-mer tiling)
  3. Produce a prioritized CSV with stealth scores

Stealth Score rationale:
  Candidates that bind CD34 (structural docking) without triggering
  MHC-I presentation are "stealthier".  The score penalises candidates
  whose 9-mers bind HLA strongly (low IC50 = high immunogenicity risk).

  stealth_score = 1 - (n_hits / n_9mers)        for peptides with ≥1 9-mer
                = 1.0                             for sequences < 9 residues

  A score of 1.0 means no 9-mer sub-peptide was predicted to bind
  HLA-A*01:01 below the 500 nM threshold.

Usage:
    python scripts/run_screening.py
"""

import os
import sys
import json

from Bio import SeqIO

from scripts.utils import AA3TO1, extract_receptor_sequence

# ── project paths ──────────────────────────────
PROJECT_ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
sys.path.insert(0, os.path.join(PROJECT_ROOT, "scripts"))

RECEPTOR_PDB = os.path.join(PROJECT_ROOT, "data", "processed", "cd34_relaxed.pdb")
CANDIDATES_FASTA = os.path.join(PROJECT_ROOT, "data", "ligands", "candidates.fasta")
AF3_OUT_DIR = os.path.join(PROJECT_ROOT, "results", "af3_inputs")
RESULTS_DIR = os.path.join(PROJECT_ROOT, "results")
CSV_OUTPUT = os.path.join(RESULTS_DIR, "screening_results.csv")

MODEL_SEEDS = []  # AF3 Server picks random seeds; consolidate_jsons.py also enforces this


def generate_af3_jsons(receptor_seq, candidates):
    """Write one AF3 JSON per candidate. Returns list of output paths."""
    os.makedirs(AF3_OUT_DIR, exist_ok=True)
    paths = []
    for cand_id, ligand_seq in candidates.items():
        job_name = f"cd34_vs_{cand_id}"
        payload = {
            "name": job_name,
            "modelSeeds": MODEL_SEEDS,
            "sequences": [
                {"proteinChain": {"sequence": receptor_seq, "count": 1}},
                {"proteinChain": {"sequence": ligand_seq, "count": 1}},
            ],
        }
        out_path = os.path.join(AF3_OUT_DIR, f"{job_name}.json")
        with open(out_path, "w") as f:
            json.dump([payload], f, indent=2)
        paths.append(out_path)
    return paths


# ──────────────────────────────────────────────
# Phase 2: MHC screening
# ──────────────────────────────────────────────
def run_mhc_screening(candidates):
    from screening_utils import screen_candidates
    return screen_candidates(candidates)


# ──────────────────────────────────────────────
# Phase 3: Stealth scoring & CSV output
# ──────────────────────────────────────────────
def compute_stealth_scores(summary_df):
    """
    stealth_score = 1 - (n_hits / n_9mers)
    Higher = better (fewer immunogenic 9-mers).
    """
    df = summary_df.copy()
    df["stealth_score"] = df.apply(
        lambda row: 1.0 if row["n_9mers"] == 0
        else 1.0 - (row["n_hits"] / row["n_9mers"]),
        axis=1,
    )
    return df.sort_values("stealth_score", ascending=False)


# ──────────────────────────────────────────────
# Main
# ──────────────────────────────────────────────
def main():
    print("=" * 60)
    print("  veron-hsc  |  Ligand screening orchestrator")
    print("=" * 60)

    # Load candidates
    records = list(SeqIO.parse(CANDIDATES_FASTA, "fasta"))
    candidates = {r.id: str(r.seq) for r in records}
    print(f"\n[1/3] Loaded {len(candidates)} candidates from FASTA")

    # Phase 1 — AF3 JSON generation
    receptor_seq = extract_receptor_sequence(RECEPTOR_PDB)
    print(f"      Receptor: {len(receptor_seq)} residues")
    json_paths = generate_af3_jsons(receptor_seq, candidates)
    print(f"      Generated {len(json_paths)} AF3 JSONs "
          f"({len(MODEL_SEEDS)} seeds each) → {AF3_OUT_DIR}/")

    # Phase 2 — MHC screening
    print(f"\n[2/3] Running MHC-I screening (HLA-A*01:01, 9-mer tiling) ...")
    summary_df, detail_df = run_mhc_screening(candidates)
    n_flagged = summary_df["flagged"].sum()
    print(f"      Screened {len(summary_df)} candidates, "
          f"{n_flagged} flagged for immunogenicity")

    # Phase 3 — Stealth scoring
    print(f"\n[3/3] Computing stealth scores ...")
    ranked = compute_stealth_scores(summary_df)
    os.makedirs(RESULTS_DIR, exist_ok=True)
    ranked.to_csv(CSV_OUTPUT, index=False, float_format="%.4f")
    print(f"      Results saved to {CSV_OUTPUT}")

    # Print top-3 leads summary
    print("\n" + "=" * 60)
    print("  TOP 3 LEADS — Stealth Score (Binding vs Immunity)")
    print("=" * 60)
    top3 = ranked.head(3)
    for _, row in top3.iterrows():
        flag_str = "CLEAN" if not row["flagged"] else "FLAGGED"
        print(
            f"  {row['candidate_id']:18s}  "
            f"stealth={row['stealth_score']:.3f}  "
            f"hits={int(row['n_hits'])}/{int(row['n_9mers'])} 9-mers  "
            f"worst_IC50={row['worst_ic50_nM']:>8.1f} nM  "
            f"[{flag_str}]"
        )

    print("\nDone.")


if __name__ == "__main__":
    main()
