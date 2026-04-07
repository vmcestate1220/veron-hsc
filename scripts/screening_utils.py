#!/usr/bin/env python
"""
MHC-I immunogenicity screening for the veron-hsc pipeline.

Predicts HLA-A*01:01 binding affinity for all 9-mer sub-peptides of
candidate sequences using MHCflurry.  Flags candidates with any 9-mer
predicted IC50 < 500 nM as potential immunogenicity risks.

Usage (as module):
    from screening_utils import screen_candidates
    results = screen_candidates(sequences_dict)
"""

import pandas as pd
from mhcflurry import Class1PresentationPredictor

ALLELE = "HLA-A*01:01"
IC50_THRESHOLD_NM = 500  # nM; below this → immunogenicity flag
KMER_LEN = 9


def tile_9mers(sequence):
    """Generate all 9-mer sub-peptides from a sequence."""
    return [sequence[i : i + KMER_LEN] for i in range(len(sequence) - KMER_LEN + 1)]


def screen_candidates(sequences):
    """
    Screen a dict of {candidate_id: amino_acid_sequence} for MHC-I binding.

    Returns a DataFrame with per-candidate immunogenicity summary:
      - candidate_id
      - sequence_length
      - n_9mers: total 9-mer windows
      - n_hits: 9-mers with IC50 < threshold
      - worst_ic50_nM: lowest (strongest-binding) IC50 among 9-mers
      - worst_peptide: the 9-mer with lowest IC50
      - flagged: True if any hit
      - all_9mer_details: DataFrame of all per-9mer predictions
    """
    predictor = Class1PresentationPredictor.load()

    all_rows = []
    summary_rows = []

    for cand_id, seq in sequences.items():
        if len(seq) < KMER_LEN:
            # Sequence too short for 9-mer tiling — flag as non-screenable
            summary_rows.append(
                {
                    "candidate_id": cand_id,
                    "sequence_length": len(seq),
                    "n_9mers": 0,
                    "n_hits": 0,
                    "worst_ic50_nM": float("nan"),
                    "worst_peptide": "",
                    "flagged": False,
                }
            )
            continue

        ninemers = tile_9mers(seq)
        preds = predictor.predict(
            peptides=ninemers,
            alleles=[ALLELE],
            verbose=0,
        )

        preds["candidate_id"] = cand_id
        preds["position"] = list(range(len(ninemers)))
        all_rows.append(preds)

        # Summarize
        hits = preds[preds["affinity"] < IC50_THRESHOLD_NM]
        best_binder = preds.loc[preds["affinity"].idxmin()]
        summary_rows.append(
            {
                "candidate_id": cand_id,
                "sequence_length": len(seq),
                "n_9mers": len(ninemers),
                "n_hits": len(hits),
                "worst_ic50_nM": best_binder["affinity"],
                "worst_peptide": best_binder["peptide"],
                "flagged": len(hits) > 0,
            }
        )

    summary_df = pd.DataFrame(summary_rows)
    detail_df = pd.concat(all_rows, ignore_index=True) if all_rows else pd.DataFrame()

    return summary_df, detail_df
