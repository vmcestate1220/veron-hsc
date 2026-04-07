#!/usr/bin/env python
"""
AF3 result ingestion for the veron-hsc pipeline.

Handles two zip formats:
  1. Bulk zip — a single zip containing multiple job subdirectories
     (e.g., folds_2026_04_07.zip → cd34_vs_ch02_wt/, cd34_vs_ch02_m7a/, ...)
  2. Per-job zips — one zip per job (fold_cd34_vs_CH02_WT.zip)

Extracts into results/af3_outputs/<job_dir>/ and validates each job
directory contains summary_confidences JSON + .cif structure files.

Usage:
    python scripts/ingest_results.py
"""

import glob
import logging
import os
import zipfile

PROJECT_ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
AF3_OUT_DIR = os.path.join(PROJECT_ROOT, "results", "af3_outputs")

logging.basicConfig(
    level=logging.INFO,
    format="%(asctime)s [%(levelname)s] %(message)s",
    datefmt="%H:%M:%S",
)
log = logging.getLogger("ingest")


def validate_job_dir(job_dir):
    """
    Verify a job directory contains at least one summary_confidences
    JSON and at least one .cif structure file (searching recursively).
    """
    has_conf = False
    has_cif = False
    for root, _, files in os.walk(job_dir):
        for f in files:
            if "summary_confidences" in f and f.endswith(".json"):
                has_conf = True
            if f.endswith(".cif") and "template" not in root:
                has_cif = True
    missing = []
    if not has_conf:
        missing.append("summary_confidences*.json")
    if not has_cif:
        missing.append("*.cif")
    return len(missing) == 0, missing


def is_bulk_zip(zf):
    """Detect if a zip is a bulk archive (contains job subdirectories)."""
    top_dirs = set()
    for name in zf.namelist():
        parts = name.split("/")
        if len(parts) > 1 and parts[0]:
            top_dirs.add(parts[0])
    # Bulk zip has multiple top-level directories matching job names
    job_dirs = [d for d in top_dirs if d.startswith("cd34_vs_")]
    return len(job_dirs) > 1, job_dirs


def ingest_zip(zip_path):
    """Extract a zip and validate. Handles both bulk and per-job formats."""
    log.info("Opening %s ...", os.path.basename(zip_path))

    with zipfile.ZipFile(zip_path, "r") as zf:
        bulk, job_dirs = is_bulk_zip(zf)

        if bulk:
            log.info("  Bulk zip detected — %d job directories", len(job_dirs))
            zf.extractall(AF3_OUT_DIR)

            success = 0
            for jd in sorted(job_dirs):
                full_path = os.path.join(AF3_OUT_DIR, jd)
                ok, missing = validate_job_dir(full_path)
                if ok:
                    log.info("    %s — OK", jd)
                    success += 1
                else:
                    log.warning("    %s — MISSING: %s", jd, ", ".join(missing))
            return success, len(job_dirs)
        else:
            # Per-job zip
            extract_dir = zip_path.replace(".zip", "")
            os.makedirs(extract_dir, exist_ok=True)
            zf.extractall(extract_dir)
            ok, missing = validate_job_dir(extract_dir)
            name = os.path.basename(extract_dir)
            if ok:
                log.info("  %s — OK", name)
            else:
                log.warning("  %s — MISSING: %s", name, ", ".join(missing))
            return (1 if ok else 0), 1


def main():
    print("=" * 60)
    print("  veron-hsc  |  AF3 result ingestion")
    print("=" * 60)

    zips = sorted(glob.glob(os.path.join(AF3_OUT_DIR, "*.zip")))
    if not zips:
        log.info("No zip files found in %s", AF3_OUT_DIR)
        return

    total_success = 0
    total_jobs = 0
    for zp in zips:
        s, t = ingest_zip(zp)
        total_success += s
        total_jobs += t

    log.info("Ingestion complete: %d/%d jobs verified", total_success, total_jobs)


if __name__ == "__main__":
    main()
