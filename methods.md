# Defect Detection Methods in TPE Calculator

## Overview

Four defect detection methods are available in `run_Ed.sh` via the `DEFECT_METHOD` parameter.
Each method determines, for a given recoil simulation, whether the PKA energy was sufficient
to create a stable crystallographic defect (Frenkel pair). The binary DEFECT/NO_DEFECT signal
drives the binary search algorithm that converges to the threshold displacement energy Ed.

---

## Method Comparison

| Method | Detection Principle | Measures | Typical Ed |
|--------|-------------------|----------|------------|
| **displace** | Atomic displacement >1.0 Å from reference | Any large atomic motion | Lower (most sensitive) |
| **CNA** | Common Neighbor Analysis of local structure | Atoms in non-BCC/FCC environments | Intermediate |
| **PTM** | Polyhedral Template Matching of local structure | Atoms in non-BCC/FCC environments | Intermediate |
| **OVITO/WS** | Wigner-Seitz cell occupancy analysis | Vacancies + interstitials (Frenkel pairs) | Standard (most specific) |

---

## Detailed Mechanisms

### 1. displace — Atomic displacement threshold

**Implementation** (`in.ed.recoil.lmp:120-124`):
```lammps
variable        ndisp atom c_disp[4]>1.0
compute         ndef all reduce sum v_ndisp
print           "NDEFECTS: ${NDEF}"
```

**Principle**: The PKA simulation runs; at the end, each atom's displacement from its
initial position is checked. Atoms displaced more than 1.0 Å contribute to the count.
If count > 0 → DEFECT.

**Physical meaning**: This is the most sensitive method. Any large atomic motion —
including **replacement collision sequences** where atoms move to neighboring sites
without creating a Frenkel pair — is flagged.

**Limitation**: 2 atoms can move >1Å during a replacement collision without creating
any stable vacancy or interstitial. The method produces a lower (artifactually
low) Ed compared to methods that specifically detect Frenkel pair formation.

**Example** (Direction 1, PKA type=2, HfNbZrTiTa):
```
E=46 eV: NDEFECTS=0   → NO_DEFECT
E=54 eV: NDEFECTS=2   → DEFECT  (but PTM=0, OVITO=NO_DEFECT — replacement collision)
E=58 eV: NDEFECTS=6   → DEFECT  (genuine Frenkel pair)
```

---

### 2. CNA — Common Neighbor Analysis

**Implementation** (`in.ed.recoil.lmp:130-134`):
```lammps
compute         cna_all all cna/atom 3.5
variable        is_cna_defect atom c_cna_all!=${expected_type}
compute         cna_ndef all reduce sum v_is_cna_defect
print           "CNA_DEFECTS: ${CNA_NDEF}"
```

**Cutoff radius**: 3.5 Å (hardcoded). Covers first-nearest-neighbor shell
for both BCC (~2.94 Å) and FCC (~2.55 Å), plus partial second-neighbor
information for structure identification.

**Expected structure types**:
- BCC → expected_type = 3
- FCC → expected_type = 1
- Custom/non-bcc-fcc → expected_type = 0 (all atoms are "defects")

**Caveat**: In practice, CNA may flag a very large number of atoms as "defects"
even in pristine configurations because the 3.5 Å cutoff and BCC/FCC template
are not perfectly selective for random solid solutions (especially multi-element
alloys). This makes CNA less reliable for TDE calculations.

---

### 3. PTM — Polyhedral Template Matching

**Implementation** (`in.ed.recoil.lmp:137-142`):
```lammps
compute         ptm_all all ptm/atom default 0.1
variable        is_ptm_defect atom c_ptm_all[1]!=${expected_type}
compute         ptm_ndef all reduce sum v_is_ptm_defect
print           "PTM_DEFECTS: ${PTM_NDEF}"
```

**RMSD threshold**: 0.1 (hardcoded). Atoms whose local structure deviates
from the template by more than this RMSD are flagged as defects.

**Principle**: PTM identifies the local crystal structure of each atom by
matching its neighbor shell against ideal templates. `c_ptm_all[1]` returns
the structure type (0=unknown, 1=FCC, 2=HCP, 3=BCC, 4-6=other). Atoms
whose structure type does not match the expected bulk type are counted.

**Advantages** over CNA:
- Better noise immunity against thermal vibrations
- Clearer structure type assignment
- More reliable for FCC identification

**Relation to Frenkel pairs**: A Frenkel pair (vacancy + interstitial)
creates two types of structural disruption:
- The vacancy: no atom at the original lattice site, so the vacancy itself
  is not a PTM "atom" (it's a missing atom). However, the nearest neighbors
  of the vacancy have fewer neighbors in their shell, so their local
  structure assignment changes (typically → unknown = 0).
- The interstitial: an extra atom in a non-lattice position is detected
  as "unknown" structure.

Both effects contribute to PTM_DEFECTS > 0. The count is typically
larger than the number of Frenkel pairs because several neighbors
of the defect are also flagged.

---

### 4. OVITO/WS — Wigner-Seitz Defect Analysis

**Implementation** (`check_defects.py`):
```python
pipeline = import_file(final_file)
ws = WignerSeitzAnalysisModifier()
ws.reference = FileSource()
ws.reference.load(ref_file)
pipeline.modifiers.append(ws)
data = pipeline.compute()
ni = data.attributes['WignerSeitz.interstitial_count']
nv = data.attributes['WignerSeitz.vacancy_count']
if ni > 0 or nv > 0: print("DEFECT")
```

**Principle**: The simulation cell is divided into Wigner-Seitz cells centered
on each reference lattice site. After the cascade:
- If a reference site is empty → **vacancy** (nv++)
- If an extra atom is found in a WS cell → **interstitial** (ni++)
- If an atom moved but still occupies its own WS cell → no defect counted

**This is the most physically rigorous method.** It directly counts the
number of Frenkel pairs (equal to both nv and ni in a perfect crystal,
though surface effects may cause slight differences in practice).

**Requirement**: Needs OVITO Python API (`ovitos`). Falls back to other
methods when OVITO is not available.

---

## Why PTM and OVITO/WS Give Identical Binary Results

### Physical Equivalence

PTM and WS detect the same physical event: **formation of a Frenkel pair**.

| Energy | Event | NDEFECTS (displace) | PTM_DEFECTS | WS (OVITO) |
|--------|-------|---------------------|-------------|------------|
| Below Ed | No defect | 0 | 0 | NO_DEFECT |
| ~Ed | Replacement collision only | 2 (moved >1Å) | 0 | NO_DEFECT |
| >Ed | Frenkel pair created | 6+ | 26+ (>0) | DEFECT (ni>0 or nv>0) |

### Key Observations

1. **At E=54 eV (Direction 1)**: NDEFECTS=2 but PTM=0. Two atoms moved >1Å,
   indicating a replacement collision sequence. No Frenkel pair was created
   because the atoms ended up in valid lattice sites. Both PTM (structure=BCC)
   and WS (no vacancy/interstitial) correctly report NO_DEFECT. The displace
   method falsely reports DEFECT.

2. **At E=58 eV**: NDEFECTS=6, PTM=26. A genuine Frenkel pair was created.
   All three methods report DEFECT. The numerical difference (6 vs 26) is due
   to PTM counting the structural disruption around the defect.

3. **The binary search tests only one energy value at a time**, and the
   DEFECT/NO_DEFECT signal determines the search direction. If both PTM and WS
   agree at every tested energy (which they do, because both detect Frenkel
   pairs), the binary search converges identically and produces the same Ed.

### When PTM and WS Would Differ

In rare edge cases, the two methods may disagree:

- **Very close to threshold**: A single Frenkel pair might be marginally stable.
  WS directly counts one vacancy and one interstitial. PTM requires that the
  structural disruption exceeds the 0.1 RMSD threshold. If the defect is
  "shallow" (barely displaced), PTM might not flag it.

- **High-temperature simulations**: Thermal fluctuations may cause PTM to
  mis-identify a few atoms near expected_type transitions. WS is immune to
  this because it compares to a reference configuration at the same temperature.

- **Multi-component alloys with chemical disorder**: In random solid solutions,
  chemical disorder alone can distort the local structure enough that PTM
  assigns "unknown" to pristine atoms. WS is unaffected because it only
  checks whether a lattice site is occupied vs empty, regardless of the atom type.

---

## Practical Recommendations

| Context | Recommended Method | Rationale |
|---------|-------------------|-----------|
| OVITO available | **OVITO/WS** | Most physically rigorous; directly counts Frenkel pairs |
| OVITO unavailable, BCC or FCC | **PTM** (with `expected_type` correctly set) | Good sensitivity, noise-resistant |
| OVITO unavailable, single-element | **displace** | Simple, no template dependency |
| Multi-element / custom structure | **displace** or **OVITO/WS** | PTM/CNA may confuse chemical disorder for structural defects |
| Quick preview runs | **displace** | Fastest (no OVITO overhead); results are ~lower bound |
| Publication-quality results | **OVITO/WS** | Industry standard; lowest systematic error |

---

## References

- Wigner-Seitz defect analysis: Nordlund, K., et al. (2018). *Nat. Commun.*, 9, 1084.
- CNA: Honeycutt, J.D. & Andersen, H.C. (1987). *J. Phys. Chem.*, 91, 4950.
- PTM: Larsen, P.M., et al. (2016). *Modelling Simul. Mater. Sci. Eng.*, 24, 055007.
