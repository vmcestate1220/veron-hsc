#!/usr/bin/env python
"""
Consolidate individual AF3 input JSONs into a single master manifest
for bulk upload to the AlphaFold Server.

Reads every .json in results/af3_inputs/, extracts the job object
(handling both bare-dict and list-wrapped formats), validates that
each job has exactly two protein entities (id A + B), enforces
modelSeeds=[], and writes the combined list to
results/veron_master_manifest.json.

Jobs are ordered to match the FASTA candidate order (CH02_WT first).

Usage:
    python scripts/consolidate_jsons.py
"""

import glob
import json
import os

from Bio import SeqIO

PROJECT_ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
AF3_INPUTS_DIR = os.path.join(PROJECT_ROOT, "results", "af3_inputs")
CANDIDATES_FASTA = os.path.join(PROJECT_ROOT, "data", "ligands", "candidates.fasta")
OUTPUT_FILE = os.path.join(PROJECT_ROOT, "results", "veron_master_manifest.json")


def load_job(path):
    """Load a job dict from a JSON file, handling both list and dict formats."""
    with open(path) as f:
        data = json.load(f)
    return data[0] if isinstance(data, list) else data


def validate_job(job):
    """Ensure the job has exactly 2 proteinChain entities with sequence + count."""
    seqs = job.get("sequences", [])
    if len(seqs) != 2:
        return False, f"expected 2 sequences, got {len(seqs)}"
    for i, s in enumerate(seqs):
        if "proteinChain" not in s:
            return False, f"entity {i}: missing 'proteinChain' key"
        p = s["proteinChain"]
        if "sequence" not in p or "count" not in p:
            return False, f"entity {i}: missing 'sequence' or 'count'"
    return True, "ok"


def main():
    # Build ordered list of expected job names from FASTA
    candidates = [r.id for r in SeqIO.parse(CANDIDATES_FASTA, "fasta")]
    expected_order = [f"cd34_vs_{cid}" for cid in candidates]

    # Index all available JSONs by job name
    json_files = glob.glob(os.path.join(AF3_INPUTS_DIR, "*.json"))
    available = {}
    for path in json_files:
        job = load_job(path)
        available[job["name"]] = job

    print(f"Found {len(available)} JSON files in {AF3_INPUTS_DIR}")

    # Assemble in FASTA order
    jobs = []
    for job_name in expected_order:
        if job_name not in available:
            print(f"  MISSING  {job_name}")
            continue

        job = available[job_name]
        job["modelSeeds"] = []

        ok, reason = validate_job(job)
        status = "OK" if ok else f"WARN ({reason})"
        print(f"  + {job['name']:40s} [{status}]")

        jobs.append(job)

    with open(OUTPUT_FILE, "w") as f:
        json.dump(jobs, f, indent=2)

    print(f"\nWrote {len(jobs)} jobs to {OUTPUT_FILE}")


if __name__ == "__main__":
    main()
