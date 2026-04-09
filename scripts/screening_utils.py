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

DEFAULT_ALLELES = [
    "HLA-A*01:01",
    "HLA-A*02:01",  # most common worldwide (~25-50% of populations)
    "HLA-B*07:02",  # common Caucasian allele (~10-15%)
]
IC50_THRESHOLD_NM = 500  # nM; below this → immunogenicity flag
KMER_LEN = 9


def tile_9mers(sequence):
    """Generate all 9-mer sub-peptides from a sequence."""
    return [sequence[i : i + KMER_LEN] for i in range(len(sequence) - KMER_LEN + 1)]


def screen_candidates(sequences, alleles=None):
    """
    Screen a dict of {candidate_id: amino_acid_sequence} for MHC-I binding.

    Predicts binding against all alleles in the panel (default: DEFAULT_ALLELES).
    Hits and worst-binder stats are aggregated across all alleles.

    Returns a DataFrame with per-candidate immunogenicity summary:
      - candidate_id
      - sequence_length
      - n_9mers: total 9-mer windows
      - n_hits: 9-mers×alleles with IC50 < threshold
      - worst_ic50_nM: lowest (strongest-binding) IC50 across all alleles
      - worst_peptide: the 9-mer with lowest IC50
      - worst_allele: the allele with lowest IC50
      - flagged: True if any hit
    """
    if alleles is None:
        alleles = DEFAULT_ALLELES

    predictor = Class1PresentationPredictor.load()

    all_rows = []
    summary_rows = []

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
        allele_preds = []
        for allele in alleles:
            preds = predictor.predict(
                peptides=ninemers,
                alleles=[allele],
                verbose=0,
            )
            preds["candidate_id"] = cand_id
            preds["position"] = list(range(len(ninemers)))
            allele_preds.append(preds)

        combined = pd.concat(allele_preds, ignore_index=True)
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

    summary_df = pd.DataFrame(summary_rows)
    detail_df = pd.concat(all_rows, ignore_index=True) if all_rows else pd.DataFrame()

    return summary_df, detail_df
