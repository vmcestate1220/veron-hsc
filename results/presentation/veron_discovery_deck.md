---
marp: true
theme: default
paginate: true
backgroundColor: #fdfdfd
style: |
  section {
    font-family: 'Helvetica Neue', Arial, sans-serif;
  }
  section.lead h1 {
    font-size: 2.4em;
    color: #1a1a2e;
  }
  section.lead h2 {
    color: #16213e;
    font-weight: 400;
  }
  h1 {
    color: #1a1a2e;
  }
  h2 {
    color: #0f3460;
  }
  h3 {
    color: #e94560;
  }
  table {
    font-size: 0.72em;
    margin: 0 auto;
  }
  th {
    background-color: #1a1a2e;
    color: #ffffff;
  }
  blockquote {
    border-left: 4px solid #e94560;
    background: #f8f8ff;
    padding: 0.5em 1em;
    font-size: 0.9em;
  }
  em {
    color: #e94560;
  }
  code {
    background: #eef;
    font-size: 0.85em;
  }
  .columns {
    display: grid;
    grid-template-columns: 1fr 1fr;
    gap: 1.5em;
  }
  footer {
    font-size: 0.55em;
    color: #888;
  }
footer: "veron-hsc  |  Confidential  |  2026-04-07"
---

<!-- _class: lead -->
<!-- _paginate: false -->
<!-- _footer: "" -->

# The Stealth LNP Strategy

## Engineering Immune-Invisible Ligands for CD34+ HSC Targeting

<br>

**veron-hsc Discovery Pipeline**
AlphaFold 3 Structure Prediction | MHCflurry Immunogenicity Screening

<br>

April 2026

---

# The Problem

## Immunogenicity Is the Bottleneck in HSC-Targeted Gene Therapy

<div class="columns">
<div>

### The clinical challenge

- LNP-based gene therapies require **surface ligands** to target CD34+ hematopoietic stem cells
- Peptide ligands are potent but trigger **anti-drug antibodies** (ADAs) in 30-60% of patients
- Even a single immunogenic 9-mer presented on MHC-I can prime a cytotoxic T-cell response
- Current approach: empirical screening in animal models -- slow, expensive, species-dependent

</div>
<div>

### What we need

- A ligand that **binds CD34 with high affinity** (interface confidence iPTM > 0.4)
- **Zero immunogenic epitopes** across major HLA alleles (stealth score = 1.0)
- A **computational pipeline** that can screen thousands of variants before touching a pipette

> *"The best antibody response is the one that never happens."*

</div>
</div>

---

# The Methodology

## A Three-Stage Computational Pipeline

<div class="columns">
<div>

### Stage 1: Structure & Screening

1. **CD34 receptor prep** -- AlphaFold DB structure (P28906), truncated to ectodomain (res 150-290), energy-minimized with OpenMM AMBER14-all
2. **Ligand library** -- 11 CH02 peptide variants (alanine scan + conservative substitutions) + 8G12 scFv control
3. **MHC-I screening** -- MHCflurry Class I presentation predictor tiles each 12-mer into 9-mer windows, predicts IC50 against HLA-A\*01:01

</div>
<div>

### Stage 2: Structure Prediction

4. **AlphaFold 3 Server** -- 12 receptor-ligand complexes, 5 seeds each (60 models total)
5. **Best-seed selection** -- Highest ranking_score per candidate

### Stage 3: Quality Gate

6. **Metric extraction** -- iPTM, pLDDT, MWSP motif distance
7. **Lead scoring** -- Weighted composite:

$$\text{Lead} = 0.4 \times \text{Stealth} + 0.4 \times \text{iPTM} + 0.2 \times \frac{\text{pLDDT}}{100}$$

</div>
</div>

---

# The Screening Landscape

## Stealth vs. Affinity Trade-off Across 12 Candidates

![w:820 center](../figures/stealth_vs_affinity.png)

CH02_V10A achieves the highest iPTM (0.49) while maintaining perfect HLA-A\*01:01 stealth. 8G12 scFv control (diamond) shows high pLDDT but poor interface score (iPTM = 0.15), confirming it does not engage CD34 as a short peptide.

---

# The Lead Ranking

## Top 5 Actionable Candidates

| Rank | Candidate | Lead Score | Stealth | iPTM | pLDDT | MWSP Dist | MHC |
|:----:|-----------|:----------:|:-------:|:----:|:-----:|:---------:|:---:|
| **1** | **CH02_V10A** | **0.695** | **1.000** | **0.490** | **49.6** | **4.1 A** | **CLEAN** |
| 2 | CH02_T1A | 0.674 | 1.000 | 0.440 | 49.1 | 4.7 A | CLEAN |
| 3 | CH02_M7L | 0.660 | 1.000 | 0.410 | 47.8 | -- | CLEAN |
| 4 | CH02_W8A | 0.639 | 1.000 | 0.380 | 43.3 | -- | CLEAN |
| 5 | CH02_M7A | 0.634 | 1.000 | 0.360 | 45.1 | -- | CLEAN |

<br>

- All top 5 candidates are **MHC-CLEAN** on HLA-A\*01:01 with perfect stealth scores
- CH02_WT (wild type) ranks **7th** -- multiple engineered variants outperform the parent
- The 8G12 scFv antibody control ranks 8th, flagged for 7 immunogenic 9-mers

---

# The Discovery: CH02_V10A

## A Single Mutation Unlocks a Motif-Anchored Binding Pose

<div class="columns">
<div>

### The mutation: Val &rarr; Ala at position 10

**Sequence:** `THRPPMWSP`*A*`WP`

Removing valine's branched side chain (~27 A^3 volume reduction) triggers a **global conformational rearrangement** of the peptide on the CD34 surface.

### Two critical structural shifts

| Metric | V10A | WT | Delta |
|--------|:----:|:--:|:-----:|
| **Met6 burial** | 1.79 A | 5.70 A | **-3.91 A** |
| **Ala10 approach** | 2.80 A | 3.41 A | **-0.61 A** |
| MWSP CA-CA dist | 4.13 A | 5.26 A | -1.13 A |
| iPTM | 0.490 | 0.350 | +0.140 |

</div>
<div>

### What happened at the atomic level

1. **C-terminal relief** -- The smaller Ala10 methyl eliminates steric clash with the receptor, pulling the C-terminus 0.61 A closer

2. **Backbone pivot** -- The freed volume allows the entire peptide to rotate, swinging Met6 **3.91 A deeper** into a hydrophobic pocket on CD34

3. **N-terminal release** -- Thr1 and His2 swing outward (+4.8, +3.5 A), trading dispensable N-terminal contacts for a tighter MWSP grip

> **Net effect:** Shallow, uniform WT pose &rarr; **polarized, motif-anchored** V10A pose with deeply buried pharmacophore

</div>
</div>

---

# The Binding Pocket

## 9 CD34 Residues Form the MWSP Recognition Site

<div class="columns">
<div>

### Primary contact ridge (res 21-26)

| Res # | AA | UniProt | Dist (A) | Contact type |
|:-----:|:--:|:-------:|:--------:|:-------------|
| D21 | Asp | 170 | 3.13 | Electrostatic |
| **I22** | **Ile** | **171** | **1.79** | **Anchor (closest)** |
| K23 | Lys | 172 | 2.75 | Electrostatic |
| A24 | Ala | 173 | 2.79 | Hydrophobic |
| E25 | Glu | 174 | 2.91 | Electrostatic |
| I26 | Ile | 175 | 3.72 | Hydrophobic |

### Secondary contacts

| K54 (203) | 3.59 A | Distal salt bridge |
|:-:|:-:|:-|
| G58 (207) | 4.58 A | Peripheral |
| E59 (208) | 4.26 A | Peripheral charged |

</div>
<div>

### Pocket architecture

The **DIKAEI** sequence (UniProt 170-175) presents an alternating charged/hydrophobic surface:

```
D - I - K - A - E - I
-   h   +   h   -   h
```

This complements the MWSP motif's mixed-character side chains:

- **Met6** &rarr; buries into the I22/I26 hydrophobic cleft
- **Trp7** &rarr; stacks against the K23/A24 ridge
- **Ser8** &rarr; hydrogen bonds to D21/E25
- **Pro9** &rarr; packs against K54 distal contact

**Anchor residue: Ile22 (1.79 A)** -- this intimate contact between receptor Ile22 and ligand Met6 defines the primary binding hotspot.

</div>
</div>

---

# Global Stealth Assessment

## The 4-Fold Improvement on HLA-C*07:01

<div class="columns">
<div>

### Multi-allele MHCflurry panel

| HLA Allele | V10A IC50 | WT IC50 | V10A vs WT |
|:-----------|:---------:|:-------:|:----------:|
| HLA-A\*01:01 | 26,696 nM | 26,749 nM | Equivalent |
| HLA-A\*02:01 | 21,566 nM | 15,058 nM | Improved |
| HLA-B\*07:02 | *463 nM* | 823 nM | Borderline |
| HLA-C\*07:01 | *449 nM* | ***117 nM*** | **3.8x improved** |

### Key finding

The WT peptide carries a **strong HLA-C\*07:01 epitope** (117 nM) from its HRPPMWSP**V** 9-mer. The V10A mutation converts this to HRPPMWSP**A**, raising IC50 to 449 nM -- a **3.8-fold reduction** in binding affinity to MHC.

</div>
<div>

### Risk stratification

**HLA-A alleles: CLEAR**
Both A\*01:01 and A\*02:01 show IC50 > 20,000 nM. No immunogenicity concern for >80% of the global population.

**HLA-C\*07:01: IMPROVED**
WT had a definitive hit at 117 nM. V10A pushed this above the 400 nM gray zone -- an inherited MWSP-core liability that V10A *ameliorates*.

**HLA-B\*07:02: MONITOR**
New borderline hit at 463 nM (WT was clean at 823 nM). Affects ~10% Caucasian carriers. Near the 500 nM threshold where MHCflurry has limited predictive value.

> **Risk classification: LOW**
> Both flagged values are in the 400-500 nM gray zone. Confirmatory T-cell assays recommended for HLA-B\*07:02 carriers.

</div>
</div>

---

<!-- _class: lead -->

# Conclusion

## Three Actionable Leads for HSC-Targeted Gene Therapy

<br>

| | CH02_V10A | CH02_T1A | CH02_M7L |
|:--|:---------:|:--------:|:--------:|
| **Lead Score** | **0.695** | 0.674 | 0.660 |
| **iPTM** | **0.490** | 0.440 | 0.410 |
| **MWSP Distance** | **4.1 A** | 4.7 A | -- |
| **HLA-A Stealth** | Clean | Clean | Clean |
| **Key advantage** | Motif-anchored pose | N-term flexibility | Conservative M&rarr;L |

<br>

### Recommended next steps

1. **T-cell validation** -- ELISPOT against HLA-B\*07:02+ PBMCs for V10A neopeptide RPPMWSPAW
2. **MD simulation** -- 100 ns explicit-solvent dynamics to confirm motif-anchored pose stability
3. **LNP formulation** -- Conjugate top 2 leads to ionizable lipid nanoparticles for in vitro CD34+ binding assay

---

<!-- _class: lead -->
<!-- _paginate: false -->

# Appendix

---

# Appendix A: Full Ranking Table

| Rank | Candidate | Lead Score | Stealth | iPTM | pLDDT | MWSP (A) | MHC |
|:----:|-----------|:----------:|:-------:|:----:|:-----:|:--------:|:---:|
| 1 | CH02_V10A | 0.695 | 1.000 | 0.490 | 49.6 | 4.1 | CLEAN |
| 2 | CH02_T1A | 0.674 | 1.000 | 0.440 | 49.1 | 4.7 | CLEAN |
| 3 | CH02_M7L | 0.660 | 1.000 | 0.410 | 47.8 | -- | CLEAN |
| 4 | CH02_W8A | 0.639 | 1.000 | 0.380 | 43.3 | -- | CLEAN |
| 5 | CH02_M7A | 0.634 | 1.000 | 0.360 | 45.1 | -- | CLEAN |
| 6 | CH02_P5A | 0.626 | 1.000 | 0.350 | 43.1 | 5.5 | CLEAN |
| 7 | CH02_WT | 0.626 | 1.000 | 0.350 | 42.9 | 5.3 | CLEAN |
| 8 | 8G12_scFv | 0.621 | 0.971 | 0.150 | 86.1 | -- | FLAG |
| 9 | CH02_R3A | 0.610 | 1.000 | 0.310 | 42.9 | 5.0 | CLEAN |
| 10 | CH02_P4S | 0.606 | 1.000 | 0.300 | 43.2 | 4.8 | CLEAN |
| 11 | CH02_V10I | 0.584 | 1.000 | 0.250 | 41.8 | 4.2 | CLEAN |
| 12 | CH02_W11A | 0.562 | 1.000 | 0.230 | 34.9 | 6.0 | CLEAN |

---

# Appendix B: Per-Residue Contact Map -- V10A vs WT

| Pos | AA | V10A (A) | WT (A) | Delta (A) | Note |
|:---:|:--:|:--------:|:------:|:---------:|:-----|
| 1 | T | 10.88 | 6.05 | +4.83 | N-term released |
| 2 | H | 6.55 | 3.10 | +3.45 | N-term released |
| 3 | R | 3.99 | 4.86 | -0.87 | Closer in V10A |
| 4 | P | 4.18 | 4.41 | -0.23 | |
| 5 | P | 3.17 | 3.04 | +0.13 | |
| **6** | **M** | **1.79** | **5.70** | **-3.91** | **Pharmacophore burial** |
| 7 | W | 3.18 | 3.32 | -0.14 | |
| 8 | S | 2.79 | 2.74 | +0.05 | |
| 9 | P | 3.59 | 2.89 | +0.70 | |
| **10** | **A/V** | **2.80** | **3.41** | **-0.61** | **Mutation site** |
| 11 | W | 2.86 | 2.21 | +0.65 | |
| 12 | P | 5.04 | 5.17 | -0.13 | |

---

# Appendix C: Pipeline & Methods

<div class="columns">
<div>

### Software stack

| Tool | Version | Purpose |
|------|---------|---------|
| AlphaFold 3 | Server | Structure prediction |
| OpenMM | 8.5 | Energy minimization |
| MHCflurry | 2.2 | MHC-I binding prediction |
| BioPython | 1.87 | Sequence/structure parsing |
| Python | 3.11 | Pipeline runtime |

### Receptor preparation

- CD34 ectodomain (UniProt P28906)
- AlphaFold DB v6 structure
- Truncated to residues 150-290 (globular + stalk)
- AMBER14-all + OBC2 implicit solvent minimization
- Final: 141 residues, renumbered 1-141

</div>
<div>

### AlphaFold 3 parameters

- 12 receptor-ligand pairs
- 5 model seeds per pair (60 total models)
- Best seed selected by ranking_score
- iPTM as primary affinity metric

### Scoring formula

$$\text{Lead Score} = 0.4S + 0.4I + 0.2\frac{P}{100}$$

Where:
- $S$ = Stealth score (1 - hit_fraction)
- $I$ = iPTM (interface predicted TM-score)
- $P$ = Ligand chain pLDDT (mean)

</div>
</div>
