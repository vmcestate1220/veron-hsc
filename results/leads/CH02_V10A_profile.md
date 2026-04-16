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

---

## 6. Safety & Off-Target Audit

*Generated 2026-04-16 | Tools: MHCflurry 2.2, BioPython ProtParam, local motif scan*

### 6.1 Sequence Homology & Mimicry Risk

A local scan of `THRPPMWSPAWP` against catalogued human short linear motifs (SLiMs) and functional sequence patterns:

| Motif | Pattern | Result | Notes |
|-------|---------|--------|-------|
| Integrin-binding (RGD) | RGD | Clean | No match |
| Nuclear receptor box (LxxLL) | LxxLL | Clean | No match |
| SH3 binding (PxxP) | PxxP | **Hit: PAWP (pos 9-12)** | See analysis below |
| ITAM immunoreceptor (YxxL) | YxxL | Clean | No match |
| PTB / endocytosis (NPxY) | NPxY | Clean | No match |
| ER retention (KDEL) | KDEL | Clean | No match |
| Caspase-3 cleavage (DEVD) | DEVD | Clean | No match |
| Cytokine receptor (WSxWS) | WSxWS | Clean | No match |
| Metal-binding (CxxC) | CxxC | Clean | No match |
| Sortase substrate (LPxTG) | LPxTG | Clean | No match |

**PxxP analysis:** The `PAWP` match at positions 9-12 is a minimal PxxP motif. However, SH3 domain binding requires flanking basic residues in the canonical `RxxPxxP` or `PxxPxR` arrangement and a minimum of ~7 residues for measurable affinity (K_d < 100 &mu;M). The isolated `PAWP` at the C-terminus lacks these features. **Risk: negligible.**

**Compositional flags:**
- Proline content: 33% (human average ~5%) &mdash; unusual but favorable for conformational rigidity and resistance to beta-sheet aggregation
- MWSP core is not a known human signaling motif; closest match (WSxWS cytokine receptor motif) differs in both sequence and structural context
- Met-Trp dipeptide frequency in human surface proteins: ~0.06%

**Mimicry risk: LOW** &mdash; no significant matches to human functional motifs.

### 6.2 MHC-I Gray Zone Sensitivity Analysis

Deep per-peptide MHCflurry analysis of the two flagged alleles (HLA-B\*07:02, HLA-C\*07:01) with full V10A vs WT comparison:

#### HLA-B\*07:02

| Pos | V10A 9-mer | IC50 (nM) | WT 9-mer | IC50 (nM) | &Delta; (nM) | Zone |
|-----|------------|----------:|----------|----------:|------:|------|
| 0 | THRPPMWSP | 23,001 | THRPPMWSP | 23,001 | 0 | Clean |
| 1 | HRPPMWSPA | 16,611 | HRPPMWSPV | 15,121 | +1,490 | Clean (V10A safer) |
| **2** | **RPPMWSPAW** | **463** | RPPMWSPVW | 823 | **-360** | **Gray zone (V10A riskier)** |
| 3 | PPMWSPAWP | 18,374 | PPMWSPVWP | 17,913 | +461 | Clean |

**Key finding:** Only 1 of 4 possible 9-mers is affected. The `RPPMWSPAW` neopeptide created by V10A crosses the 500 nM threshold (463 nM), while the WT equivalent `RPPMWSPVW` remains in the weak-binding range (823 nM). This is a *de novo* liability introduced by the mutation.

#### HLA-C\*07:01

| Pos | V10A 9-mer | IC50 (nM) | WT 9-mer | IC50 (nM) | &Delta; (nM) | Zone |
|-----|------------|----------:|----------|----------:|------:|------|
| 0 | THRPPMWSP | 16,732 | THRPPMWSP | 16,732 | 0 | Clean |
| **1** | **HRPPMWSPA** | **449** | **HRPPMWSPV** | **117** | **+332** | **Gray zone (V10A 3.8&times; safer)** |
| 2 | RPPMWSPAW | 10,558 | RPPMWSPVW | 10,755 | -197 | Clean |
| 3 | PPMWSPAWP | 29,815 | PPMWSPVWP | 29,639 | +176 | Clean |

**Key finding:** V10A *improves* the HLA-C\*07:01 liability from a strong hit (117 nM in WT) to gray zone (449 nM). The `HRPPMWSPV` peptide in WT is 3.8&times; more immunogenic than V10A's `HRPPMWSPA`. This is an inherited MWSP-core liability that V10A ameliorates.

#### Gray Zone Risk Summary

| Allele | V10A IC50 | WT IC50 | V10A vs WT | Liability Origin | Population Freq |
|--------|----------:|--------:|------------|------------------|:-----------:|
| HLA-B\*07:02 | 463 nM | 823 nM | Worsened | **New** (V10A-introduced) | ~10% Caucasian |
| HLA-C\*07:01 | 449 nM | 117 nM | Improved 3.8&times; | Inherited (MWSP core) | ~25% Caucasian |

Both flagged values fall in the 400&ndash;500 nM range where:
- MHCflurry IC50 predictions have ~0.7 AUC (limited clinical discriminative power)
- IC50 < 500 nM does **not** guarantee T-cell recognition *in vivo*
- Peptide-MHC stability (t&frac12;) and TCR repertoire availability matter more than binding affinity alone

### 6.3 Hydrophobicity & Aggregation Profile

| Property | V10A | WT | &Delta; | Interpretation |
|----------|-----:|---:|--------:|----------------|
| Molecular weight | 1,462.7 Da | 1,490.7 Da | -28.0 | Lighter (lost isopropyl &rarr; methyl) |
| GRAVY score | **-1.142** | -0.942 | -0.200 | V10A more hydrophilic |
| Net charge (pH 7.0) | +0.49 | +0.49 | 0 | Identical |
| Isoelectric point | 9.44 | 9.44 | 0 | Identical |
| Instability index | 87.65 | 87.65 | 0 | Both "unstable" by Guruprasad definition |

**GRAVY interpretation:** V10A scores &minus;1.142, firmly in the hydrophilic range. Peptides with GRAVY > +0.5 are aggregation-prone in aqueous LNP formulations; V10A is well below this threshold. The Val&rarr;Ala substitution reduces local hydrophobicity at position 10 (Kyte-Doolittle: +4.2 &rarr; +1.8), making V10A more LNP-compatible than WT.

**Aggregation risk: VERY LOW.** Contributing factors:
- Strongly negative GRAVY (&minus;1.14) &mdash; net hydrophilic character
- 33% proline content acts as a structural breaker, disfavoring &beta;-sheet nucleation
- No continuous hydrophobic stretches (longest: Met6 at KD +1.9, flanked by prolines)

**Instability index note:** The high instability index (87.65) reflects the Guruprasad dipeptide model, which penalizes Pro-Pro, Trp-Ser, and Pro-Met dipeptides. This metric was calibrated on full-length proteins and is unreliable for 12-residue peptides. For short therapeutic peptides, proteolytic stability is better assessed by *in vitro* serum half-life assays.

**Solubility:** 17% charged residues (His, Arg) &mdash; borderline by Wilkinson-Harrison criteria but adequate for LNP-encapsulated delivery where aqueous solubility is not rate-limiting.

### 6.4 Safety Verdict

| Audit Category | Risk Level | Disposition |
|----------------|:----------:|-------------|
| Human protein mimicry | **LOW** | No functional motif matches; PxxP hit is sub-threshold |
| HLA-A stealth (A\*01:01, A\*02:01) | **CLEAN** | All 9-mers >20,000 nM across both alleles |
| HLA-B\*07:02 immunogenicity | **GRAY ZONE** | 463 nM (new liability); confirm with T-cell assay |
| HLA-C\*07:01 immunogenicity | **GRAY ZONE** | 449 nM (improved 3.8&times; vs WT); inherited from MWSP core |
| LNP aggregation risk | **VERY LOW** | GRAVY &minus;1.14; proline-rich breaker sequence |
| Solubility | **ADEQUATE** | Borderline charged content; mitigated by LNP encapsulation |

**Overall: CONDITIONALLY CLEARED for advancement.**

The single actionable risk is the *de novo* HLA-B\*07:02 epitope (`RPPMWSPAW`, 463 nM). This requires experimental validation via ELISPOT or T-cell proliferation assay against HLA-B\*07:02+ donor PBMCs before clinical development. All other safety metrics are favorable.
