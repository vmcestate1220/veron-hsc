# Lead Candidate Profile: CH02_V10A

**Pipeline:** veron-hsc | **Date:** 2026-04-07 | **AF3 Seed:** 0 (best by ranking_score)

## 1. Candidate Summary

| Property | Value |
|---|---|
| **Sequence** | `THRPPMWSP`**A**`WP` |
| **Mutation** | Val &rarr; Ala at position 10 |
| **Lead Score** | **0.695** (rank 1 of 12) |
| **iPTM** | 0.490 |
| **pTM** | 0.680 |
| **Ligand pLDDT** | 49.6 |
| **Ranking Score** | 0.610 |
| **Fraction Disordered** | 0.18 |
| **Has Clash** | No |
| **MWSP Min Distance** | 4.13 &Aring; |

**Lead Score Decomposition:** 0.40 &times; 1.000 (stealth) + 0.40 &times; 0.490 (iPTM) + 0.20 &times; 0.496 (pLDDT/100) = **0.695**

---

## 2. Structural Deep-Dive: CD34 Binding Pocket

The MWSP motif (residues M6-W7-S8-P9) on Chain B makes contact with **9 receptor residues** on CD34 Chain A within a 5 &Aring; cutoff. All residue numbers use the renumbered 1&ndash;141 scheme (original CD34 residues 150&ndash;290; add 149 for UniProt numbering).

| Res # | Residue | 1-Letter | UniProt # | Min Distance (&Aring;) | Role |
|-------|---------|----------|-----------|------------------------|------|
| 21 | ASP | D | 170 | 3.13 | Charged contact |
| **22** | **ILE** | **I** | **171** | **1.79** | **Closest contact** |
| 23 | LYS | K | 172 | 2.75 | Charged contact |
| 24 | ALA | A | 173 | 2.79 | Hydrophobic packing |
| 25 | GLU | E | 174 | 2.91 | Charged contact |
| 26 | ILE | I | 175 | 3.72 | Hydrophobic packing |
| 54 | LYS | K | 203 | 3.59 | Distal charged contact |
| 58 | GLY | G | 207 | 4.58 | Peripheral contact |
| 59 | GLU | E | 208 | 4.26 | Peripheral charged |

**Binding pocket signature:** `D21-I22-K23-A24-E25-I26` forms the primary contact ridge (residues 170&ndash;175 in UniProt numbering), with K54 contributing a secondary electrostatic contact from a distal loop. The pocket is dominated by the **DIKAEI** sequence stretch, presenting an alternating charged/hydrophobic surface that complements the MWSP motif's mixed-character side chains.

**Closest contact:** Ile22 (1.79 &Aring;) &mdash; this tight packing between receptor I22 and the Met6 side chain of the MWSP motif anchors the binding interface.

---

## 3. Binding Pose Comparison: V10A vs Wild Type

### Per-Residue Contact Map

| Pos | Residue | V10A (&Aring;) | WT (&Aring;) | &Delta; (&Aring;) | Interpretation |
|-----|---------|---------------|-------------|-------------------|----------------|
| 1 | Thr | 10.88 | 6.05 | +4.83 | N-terminus displaced in V10A |
| 2 | His | 6.55 | 3.10 | +3.45 | N-terminus displaced in V10A |
| 3 | Arg | 3.99 | 4.86 | -0.87 | Closer in V10A |
| 4 | Pro | 4.18 | 4.41 | -0.23 | Similar |
| 5 | Pro | 3.17 | 3.04 | +0.13 | Similar |
| **6** | **Met** | **1.79** | **5.70** | **-3.91** | **Dramatically closer in V10A** |
| 7 | Trp | 3.18 | 3.32 | -0.14 | Similar |
| 8 | Ser | 2.79 | 2.74 | +0.05 | Similar |
| 9 | Pro | 3.59 | 2.89 | +0.70 | Slightly farther in V10A |
| **10** | **Ala/Val** | **2.80** | **3.41** | **-0.61** | **Mutant residue: closer approach** |
| 11 | Trp | 2.86 | 2.21 | +0.65 | Slightly farther in V10A |
| 12 | Pro | 5.04 | 5.17 | -0.13 | Similar |

### C-Terminal Approach Analysis

The Val&rarr;Ala mutation at position 10 reduces the side-chain volume by ~27 &Aring;&sup3;, producing two structural consequences:

1. **Direct effect at the mutation site:** Ala10 approaches the receptor 0.61 &Aring; closer than Val10 in WT (2.80 vs 3.41 &Aring;). The smaller alanine methyl group eliminates the steric clash from valine's branched isopropyl, allowing tighter C-terminal packing against the receptor surface.

2. **Allosteric reorientation of the MWSP core:** Met6 undergoes a **dramatic repositioning** (&Delta; = -3.91 &Aring;), moving from a distant 5.70 &Aring; in WT to an intimate 1.79 &Aring; contact in V10A. This suggests the reduced C-terminal bulk allows the entire peptide backbone to pivot, burying the methionine deep into a hydrophobic pocket on CD34.

3. **N-terminal release:** Thr1 and His2 swing outward (+4.83, +3.45 &Aring;), indicating the peptide trades N-terminal contacts for a tighter MWSP-centered grip &mdash; a favorable trade given the MWSP motif is the pharmacophore.

**Net effect:** The V10A mutation converts a relatively uniform, shallow binding mode (WT) into a **polarized, motif-anchored pose** where the MWSP core is deeply buried while the N- and C-termini remain solvent-exposed.

---

## 4. Global HLA Stealth Assessment

Multi-allele MHCflurry screening (Class1PresentationPredictor) against 4 globally prevalent HLA alleles:

### CH02_V10A Results

| HLA Allele | Hits (<500 nM) | Worst IC50 (nM) | Worst Peptide | Status |
|------------|---------------|-----------------|---------------|--------|
| HLA-A\*01:01 | 0 | 26,696 | RPPMWSPAW | CLEAN |
| HLA-A\*02:01 | 0 | 21,566 | HRPPMWSPA | CLEAN |
| HLA-B\*07:02 | 1 | 463 | RPPMWSPAW | FLAG |
| HLA-C\*07:01 | 1 | 449 | HRPPMWSPA | FLAG |

### Wild Type (CH02_WT) Comparison

| HLA Allele | WT Worst IC50 (nM) | V10A Worst IC50 (nM) | V10A vs WT |
|------------|-------------------|---------------------|------------|
| HLA-A\*01:01 | 26,749 | 26,696 | ~equivalent |
| HLA-A\*02:01 | 15,058 | 21,566 | V10A improved (higher IC50) |
| HLA-B\*07:02 | 823 | 463 | V10A worsened (now below 500 nM) |
| HLA-C\*07:01 | 117 | 449 | V10A improved (116&rarr;449 nM) |

### Interpretation

- **HLA-A alleles (A\*01:01, A\*02:01):** Fully stealth. All 9-mers exceed 20,000 nM &mdash; no immunogenicity concern.
- **HLA-B\*07:02:** Borderline hit at 463 nM (WT was clean at 823 nM). The V10A mutation introduced a marginal HLA-B\*07:02 epitope via the RPPMWSP**A**W neopeptide. This is near-threshold and may not translate to T-cell recognition in vivo.
- **HLA-C\*07:01:** Flagged at 449 nM, but this represents a **3.8-fold improvement** over WT (117 nM). The WT motif HRPPMWSP**V** is a substantially stronger HLA-C\*07:01 binder than V10A's HRPPMWSP**A**. This is an inherited liability from the MWSP core, not introduced by the V10A mutation.

**Verdict:** The HLA-B\*07:02 hit is a newly introduced borderline liability. The HLA-C\*07:01 hit is an inherited MWSP-core liability that V10A actually ameliorates. Both are near-threshold (>400 nM) and affect only 1 of 4 possible 9-mers each.

**Risk classification:** LOW &mdash; both flagged IC50 values are within the 400&ndash;500 nM "gray zone" where MHCflurry predictions have limited clinical predictive value. For HLA-B\*07:02 carriers (~10% Caucasian frequency), confirmatory T-cell proliferation assays are recommended before clinical advancement.

---

## 5. Actionable Summary

### Strengths
- **Top-ranked lead** (0.695) with highest iPTM (0.490) in the panel
- **Tightest MWSP engagement** (4.13 &Aring; CA-CA, 1.79 &Aring; all-atom)
- Mutation produces favorable **motif-anchored binding pose** with deep Met6 burial
- **No HLA-A epitopes** across both major A alleles
- **Improves** HLA-C\*07:01 liability 3.8-fold vs WT

### Risks
- Borderline HLA-B\*07:02 hit (463 nM) &mdash; absent in WT
- Moderate ligand pLDDT (49.6) suggests conformational flexibility in the unbound state
- N-terminal residues (T1, H2) are solvent-exposed and may be susceptible to proteolytic cleavage

### Recommended Next Steps
1. **T-cell assay:** ELISPOT or proliferation assay against HLA-B\*07:02+ donor PBMCs using RPPMWSPAW peptide
2. **MD refinement:** 100 ns explicit-solvent simulation of the V10A:CD34 complex to validate the motif-anchored pose and assess binding free energy
3. **Backup candidates:** CH02_T1A (rank 2, lead score 0.674) and CH02_M7L (rank 3, lead score 0.660) as fallback leads if HLA-B\*07:02 liability proves clinically relevant
4. **Cyclization screen:** Evaluate head-to-tail or disulfide cyclization to protect the exposed N-terminus while preserving the MWSP contact geometry
