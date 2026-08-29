# Computational Design of the HCMV Pentamer

This repository contains the computational design workflow used to stabilize the human cytomegalovirus (HCMV) pentamer through structure-guided optimization of inter-subunit interfaces.

The design strategy was performed in two sequential stages:

1. Optimization of the **gL–UL128–UL130** ternary interface, resulting in **PC-Int1**.
2. Optimization of the **UL128–UL130–UL131A** ternary interface on the PC-Int1 background, resulting in **PC-Int2**.

The computational workflow combines structural modeling, interface residue selection, ProteinMPNN sequence design, and AlphaFold3-based structural evaluation.

---

## 1. Structural Template Preparation

The initial structural template of the HCMV pentamer was derived from the crystal structure **PDB 5VOB**.

Using **Coot**, the amino acid sequence was modified to match the wild-type (WT) sequence used in this study. Missing loop regions were manually modeled to generate a complete structural template for subsequent computational design.

---

## 2. Design of the gL–UL128–UL130 Interface

### 2.1 Interface definition

The center of the gL–UL128–UL130 ternary interface was defined as the **NH1 atom of residue R66 in chain D (UL130)** of the 5VOB structure.

Designable residues were identified according to three criteria:

1. At least one atom of the residue was located within **10 Å** of the interface center.
2. The residue participated in the gL–UL128, UL128–UL130, or gL–UL130 interface, defined by at least one interatomic contact within **4 Å**.
3. Residues protected by existing GSK patent (CN111094323B and CN113271968B) claims were excluded from redesign.

The final designable residues submitted to ProteinMPNN were:
```text
gL:      141, 144, 146, 152, 154, 159, 161
UL128:   139, 140, 141, 143, 144, 145
UL130:   65, 66, 70, 97, 98, 99
```



### 2.2 Designable residues





### 2.3 ProteinMPNN design

ProteinMPNN was used to generate **100 candidate sequences** while restricting mutations to the predefined designable positions.

Candidate substitutions were subsequently filtered according to two criteria:

1. The mutation occurred in more than **25%** of the generated sequences at that position.
2. Structural inspection indicated the introduction of additional favorable interfacial interactions, including hydrogen bonds, salt bridges, or improved hydrophobic packing.

The following substitutions were retained:

```text
UL130:P65E
UL130:R66F
UL130:Q70W
UL130:T98I
```

The following substitutions were excluded after structural inspection:

```text
UL128:K139P
UL128:I140P
UL128:T141G
```

### 2.4 AlphaFold3 screening

All two-, three-, and four-mutation combinations of the retained substitutions were subsequently modeled using **AlphaFold3**.

Candidate models were ranked according to:

1. The average pLDDT of the candidate mutation sites.
2. For tied candidates, the backbone RMSD of the designed residues relative to the WT structure was used as a secondary criterion.

The highest-ranked design contained:

```text
UL130:P65E
UL130:R66F
UL130:T98I
```

This design was designated **PC-Int1** and selected for subsequent experimental validation.

---

## 3. Design of the UL128–UL130–UL131A Interface

To further improve pentamer stability, a second computational design round was performed at the **UL128–UL130–UL131A ternary interface**.

The three stabilizing mutations identified in PC-Int1 were retained as the structural background:

```text
UL130:P65E
UL130:R66F
UL130:T98I
```

The UL128–UL130–UL131A interface was then independently optimized on the PC-Int1 background.

### 3.1 Interface definition

The center of the UL128–UL130–UL131A interface was defined as the **Cβ atom of residue D210 in chain D (UL130)** of the 5VOB structure.

Designable residues were selected using the same criteria as described above:

1. Spatial proximity to the interface center.
2. Participation in inter-subunit contacts.
3. Exclusion of residues covered by existing GSK patent claims.

### 3.2 Designable residues

The final designable residues submitted to ProteinMPNN were:

```text
UL128:
86, 100, 123, 124, 125, 126, 127, 128, 129, 130

UL130:
161, 163, 164, 165, 196, 206, 207, 208, 209, 210, 211, 212, 213

UL131A:
37, 82, 83, 85, 89, 113, 118
```

### 3.3 ProteinMPNN design

ProteinMPNN was used to generate **100 candidate sequences** while restricting mutations to the predefined designable positions.

Candidate substitutions were filtered using the same criteria described above:

- Mutation frequency >25%
- Favorable structural effects at the interface

The following substitutions were retained:

```text
UL128:V128I
UL128:Y130C
UL130:V163C
UL131A:V82L
UL131A:R118M
```

The following substitutions were excluded:

```text
UL128:A125R
UL128:G126R
UL128:S127E
UL128:P129D
UL130:N164K
UL130:F196L
UL130:F206L
UL131A:L113I
```

These substitutions were excluded because structural analysis suggested that they were unlikely to introduce additional favorable interfacial interactions or could potentially perturb the local packing environment.

### 3.4 AlphaFold3 screening

All combinations of the retained mutations were subsequently modeled using **AlphaFold3** together with the three stabilizing mutations from PC-Int1.

Candidate designs were ranked according to:

1. The average pLDDT of the designed residues.
2. The local backbone RMSD of the designed residues relative to the WT structure.

The highest-ranked design contained four additional substitutions:

```text
UL128:V128I
UL128:Y130C
UL130:V163C
UL131A:V82L
```

on the PC-Int1 background.

This design was designated **PC-Int2** and selected for subsequent experimental characterization.

---

## 4. Final Designs

| Design | Background | Optimized interface | Mutations |
|---|---|---|---|
| **WT** | 5VOB-derived WT | — | — |
| **PC-Int1** | WT | gL–UL128–UL130 | UL130:P65E, R66F, T98I |
| **PC-Int2** | PC-Int1 | UL128–UL130–UL131A | UL128:V128I, Y130C; UL130:V163C; UL131A:V82L |

### PC-Int1

```text
UL130:P65E
UL130:R66F
UL130:T98I
```

### PC-Int2

```text
UL130:P65E
UL130:R66F
UL130:T98I

UL128:V128I
UL128:Y130C
UL130:V163C
UL131A:V82L
```

---

## 5. Software

The computational design workflow used the following software:

| Software | Purpose |
|---|---|
| [Coot](https://www2.mrc-lmb.cam.ac.uk/personal/pemsley/coot/) | Structural model preparation and manual loop modeling |
| [ProteinMPNN](https://github.com/dauparas/ProteinMPNN) | Protein sequence design |
| [AlphaFold3](https://github.com/google-deepmind/alphafold3) | Structural prediction and candidate evaluation |

---

## 6. Repository Structure

```text
HCMV-pentamer-design/
│
├── README.md
│
├── structures/
│   ├── 5VOB/
│   ├── WT/
│   ├── PC-Int1/
│   └── PC-Int2/
│
├── proteinmpnn/
│   ├── PC-Int1/
│   └── PC-Int2/
│
├── alphafold3/
│   ├── PC-Int1/
│   └── PC-Int2/
│
├── scripts/
│   ├── interface_selection.py
│   ├── analyze_mpnn_sequences.py
│   ├── filter_mutations.py
│   └── analyze_af3_models.py
│
└── results/
    ├── PC-Int1/
    └── PC-Int2/
```

---

## 7. Key Parameters

| Parameter | Value |
|---|---:|
| Interface distance cutoff | 10 Å |
| Interatomic contact cutoff | 4 Å |
| ProteinMPNN sequences per round | 100 |
| Mutation frequency threshold | >25% |
| Structural prediction | AlphaFold3 |

---

## 8. Citation

If you use the computational workflow or designs described in this repository, please cite the associated publication and the corresponding software.

### Structural template

**PDB 5VOB**

### ProteinMPNN

Dauparas et al.  
*Robust deep learning–based protein sequence design using ProteinMPNN.*

### AlphaFold3

Abramson et al.  
*Accurate structure prediction of biomolecular interactions with AlphaFold 3.*

### Coot

Emsley et al.  
*Coot: model-building tools for molecular graphics.*

---

## 9. License

Please see the `LICENSE` file for information regarding the use and redistribution of the code and data contained in this repository.
