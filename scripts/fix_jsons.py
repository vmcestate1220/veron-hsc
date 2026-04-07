#!/usr/bin/env python
"""
One-time fix: wrap existing AF3 input JSONs in a list for Server compatibility.

The AlphaFold Server expects the top-level JSON to be a list containing the
job dictionary.  This script reads every .json in results/af3_inputs/ and
wraps any bare dict in a list, saving back in place.

Usage:
    python scripts/fix_jsons.py
"""

import glob
import json
import os

PROJECT_ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
AF3_INPUTS_DIR = os.path.join(PROJECT_ROOT, "results", "af3_inputs")


def main():
    json_files = sorted(glob.glob(os.path.join(AF3_INPUTS_DIR, "*.json")))
    print(f"Found {len(json_files)} JSON files in {AF3_INPUTS_DIR}")

    fixed = 0
    skipped = 0
    for path in json_files:
        with open(path) as f:
            data = json.load(f)

        if isinstance(data, dict):
            with open(path, "w") as f:
                json.dump([data], f, indent=2)
            fixed += 1
            print(f"  FIXED  {os.path.basename(path)}")
        else:
            skipped += 1
            print(f"  OK     {os.path.basename(path)}")

    print(f"\nDone: {fixed} fixed, {skipped} already correct.")


if __name__ == "__main__":
    main()
