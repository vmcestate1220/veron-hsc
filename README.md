# veron-hsc: AI-Driven Engineering of Stealth LNP Ligands

Computational optimization of the CH02 peptide motif for CD34+ hematopoietic stem cell targeting, combining **AlphaFold 3** structure prediction with **MHC-I immunogenicity screening** to identify ligands that bind with high affinity while evading immune detection.

## Key Result

**CH02_V10A** (Val→Ala at position 10) emerged as the top lead from a panel of 12 candidates:

| Metric | V10A | Wild Type | Improvement |
|--------|:----:|:---------:|:-----------:|
| Lead Score | **0.695** | 0.626 | +11% |
| iPTM | **0.490** | 0.350 | +40% |
| MWSP Motif Distance | **4.1 Å** | 5.3 Å | 1.2 Å closer |
| Met6 Burial Depth | **1.79 Å** | 5.70 Å | 3.91 Å deeper |
| HLA-C\*07:01 IC50 | **449 nM** | — | Borderline (threshold 500 nM) |
| HLA-A Stealth | Clean | Clean | — |

The single Val→Ala mutation eliminates a C-terminal steric clash, triggering a backbone pivot that buries the MWSP pharmacophore 3.91 Å deeper into the CD34 receptor surface — converting a shallow, uniform binding mode into a polarized, motif-anchored pose.

## Pipeline Architecture

```
┌─────────────────┐    ┌──────────────────┐    ┌─────────────────┐
│   preprocess.py  │    │  af3_generator.py │    │ run_screening.py│
│                  │    │                   │    │                 │
│ AlphaFold DB     │───▶│ 12 receptor-ligand│───▶│ MHCflurry       │
│ CD34 ectodomain  │    │ AF3 Server JSONs  │    │ HLA-A*01:01     │
│ OpenMM minimize  │    │ (5 seeds each)    │    │ Stealth scoring │
└─────────────────┘    └──────────────────┘    └────────┬────────┘
                                                        │
┌─────────────────┐    ┌──────────────────┐             │
│  postprocess.py  │◀───│ ingest_results.py│◀── AF3 Server results
│                  │    │                   │             │
│ iPTM, pLDDT      │    │ Bulk zip extract  │    ┌───────▼────────┐
│ MWSP distance    │    │ Validation        │    │screening_results│
│ Lead scoring     │    └──────────────────┘    │     .csv        │
│ Visualization    │                            └────────────────┘
└────────┬────────┘
         │
         ▼
  veron_prioritized_leads.csv
  stealth_vs_affinity.png
```

## Stealth vs. Affinity Landscape

![Stealth vs Affinity](results/figures/stealth_vs_affinity.png)

## Quick Start

### Prerequisites

- Conda (Miniforge recommended for Apple Silicon)
- macOS with Apple M-series GPU (OpenCL) or Linux with CUDA

### Setup

```bash
# Create conda environment
conda env create -f environment.yml -p ./veron-hsc
conda activate ./veron-hsc

# Install pip-only dependencies
./veron-hsc/bin/pip install mhcflurry
mhcflurry-downloads fetch models_class1_presentation
```

### Run the Pipeline

```bash
# Stage 1: Fetch and minimize CD34 receptor
./veron-hsc/bin/python scripts/preprocess.py

# Stage 2: Generate AF3 inputs + MHC screening
./veron-hsc/bin/python scripts/run_screening.py

# Stage 3: Consolidate for AF3 Server bulk upload
./veron-hsc/bin/python scripts/consolidate_jsons.py
# → Upload results/veron_master_manifest.json to alphafoldserver.com

# Stage 4: Ingest AF3 results (download zips to results/af3_outputs/)
./veron-hsc/bin/python scripts/ingest_results.py

# Stage 5: Post-processing quality gate + final ranking
./veron-hsc/bin/python scripts/postprocess.py

# Stage 6: Lead candidate profiling
./veron-hsc/bin/python scripts/lead_profile_v10a.py
```

## Project Structure

```
veron-hsc/
├── scripts/
│   ├── preprocess.py           # CD34 structure preparation (AlphaFold DB + OpenMM)
│   ├── af3_generator.py        # AF3 Server JSON input generation
│   ├── screening_utils.py      # MHCflurry immunogenicity screening module
│   ├── run_screening.py        # Master orchestrator (AF3 gen + MHC screening)
│   ├── consolidate_jsons.py    # Merge JSONs into bulk upload manifest
│   ├── ingest_results.py       # AF3 result zip extraction + validation
│   ├── postprocess.py          # Quality gate: iPTM, pLDDT, MWSP, lead scoring
│   └── lead_profile_v10a.py    # V10A deep structural + stealth analysis
├── data/
│   ├── raw/                    # Downloaded AlphaFold structures
│   ├── processed/              # Truncated + minimized receptor PDBs
│   └── ligands/                # Candidate FASTA library (12 sequences)
├── results/
│   ├── af3_inputs/             # Generated AF3 Server JSON inputs
│   ├── screening_results.csv   # MHC-I screening results
│   ├── veron_prioritized_leads.csv  # Final ranked lead table
│   ├── veron_master_manifest.json   # Bulk AF3 Server upload manifest
│   ├── figures/                # Visualizations
│   ├── leads/                  # Lead candidate profiles
│   └── presentation/           # Marp slide deck
├── environment.yml             # Conda environment specification
├── CLAUDE.md                   # Claude Code project instructions
└── README.md
```

## Scoring Methodology

Candidates are ranked by a weighted **Lead Score**:

```
Lead Score = 0.4 × Stealth + 0.4 × iPTM + 0.2 × (pLDDT / 100)
```

| Component | Weight | Source | Rationale |
|-----------|:------:|--------|-----------|
| Stealth Score | 40% | MHCflurry | 1 − (immunogenic 9-mers / total 9-mers) |
| iPTM | 40% | AlphaFold 3 | Interface predicted TM-score (binding confidence) |
| pLDDT | 20% | AlphaFold 3 | Per-residue confidence of the ligand chain |

## Technology Stack

| Tool | Version | Purpose |
|------|---------|---------|
| [AlphaFold 3 Server](https://alphafoldserver.com) | 2026 | Protein complex structure prediction |
| [OpenMM](https://openmm.org) | 8.5 | Molecular dynamics / energy minimization |
| [MHCflurry](https://github.com/openvax/mhcflurry) | 2.2 | MHC Class I binding prediction |
| [BioPython](https://biopython.org) | 1.87 | Sequence and structure parsing |
| Python | 3.11 | Pipeline runtime |

## Candidate Library

The CH02 peptide (`THRPPMWSPVWP`) targets the CD34 ectodomain via the MWSP binding motif. The library includes:

- **7 alanine-scanning variants** (T1A, R3A, P4S, P5A, M7A, W8A, V10A, W11A) — identify dispensable positions
- **3 conservative substitutions** (M7L, P4S, V10I) — probe steric tolerance
- **1 positive control** (8G12 scFv) — anti-CD34 monoclonal antibody single-chain variable fragment

## License

This project is for research purposes only.
