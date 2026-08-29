# HCMV-Pentamer-interface-stabilization-design
This repository contains the computational design workflow used to improve the stability of the human cytomegalovirus (HCMV) pentameric complex through structure-guided interface optimization.
The design strategy consists of two sequential interface-optimization stages:
1.	Optimization of the gL–UL128–UL130 ternary interface, resulting in PC-Int1.
2.	Optimization of the UL128–UL130–UL131A ternary interface on the PC-Int1 background, resulting in PC-Int2.
The workflow combines structural modeling, interface residue selection, ProteinMPNN sequence design, structural inspection, and AlphaFold3-based structural evaluation.
________________________________________
Overview
The HCMV pentamer is composed of five viral proteins:
•	gL
•	UL128
•	UL130
•	UL131A
•	gH
The computational design described here focuses on two ternary interfaces within the pentamer:
                 UL131A
                    │
                    │
              UL128 ─── UL130
                │          │
                │          │
                └── gL ────┘
The design procedure was performed sequentially.
Stage 1
gL–UL128–UL130 interface
5VOB
  │
  ├── WT sequence reconstruction
  ├── Missing-loop modeling
  ├── Interface residue identification
  ├── Patent-protected residue exclusion
  ├── ProteinMPNN sequence design
  ├── Mutation-frequency filtering
  ├── Structural inspection
  └── AlphaFold3 screening
          │
          ▼
       PC-Int1
PC-Int1 contains:
UL130:P65E
UL130:R66F
UL130:T98I
Stage 2
UL128–UL130–UL131A interface
PC-Int1 was used as the structural background for a second round of interface optimization.
PC-Int1
   │
   ├── Interface residue identification
   ├── Patent-protected residue exclusion
   ├── ProteinMPNN sequence design
   ├── Mutation-frequency filtering
   ├── Structural inspection
   └── AlphaFold3 screening
           │
           ▼
        PC-Int2
PC-Int2 contains the three PC-Int1 substitutions plus:
UL128:V128I
UL128:Y130C
UL130:V163C
UL131A:V82L
________________________________________
1. Structural Template
The initial structural template of the HCMV pentamer was derived from the crystal structure:
PDB: 5VOB
The structure was processed using Coot.
The amino acid sequences in the structural model were first modified to match the WT sequences used in this study. Missing loop regions were subsequently modeled manually to generate a complete structural template for computational design.
The resulting WT model was used as the structural reference for both rounds of interface optimization.
________________________________________
2. Design Strategy I: gL–UL128–UL130 Interface
2.1 Interface definition
The center of the gL–UL128–UL130 ternary interface was defined as the NH1 atom of R66 in chain D (UL130) of the 5VOB structure.
Candidate designable residues were selected according to three criteria.
Criterion 1 — Spatial proximity
At least one atom of the residue had to be located within:
10 Å
of the interface center.
Criterion 2 — Interfacial contact
The residue had to participate in at least one of the following interfaces:
•	gL–UL128
•	UL128–UL130
•	gL–UL130
An inter-subunit contact was defined using an interatomic distance cutoff of:
4 Å
Criterion 3 — Exclusion of protected residues
Residues covered by the relevant GSK patent claims were excluded from redesign.
The excluded residues were:
gL:
140, 145, 150, 160, 166

UL128:
83, 98, 142

UL130:
52, 62, 64, 95, 167, 204

UL131A:
86
Relevant patent references:
•	CN111094323B
•	CN113271968B
________________________________________
2.2 Designable residues
After applying the above criteria, the final designable positions were:
gL:
141, 144, 146, 152, 154, 159, 161

UL128:
139, 140, 141, 143, 144, 145

UL130:
65, 66, 70, 97, 98, 99
These positions were provided to ProteinMPNN as designable residues.
________________________________________
3. ProteinMPNN Sequence Design
ProteinMPNN was used to generate sequence variants while restricting mutations to the predefined designable positions.
For each design round:
Number of generated sequences = 100
Only the predefined interface residues were allowed to change.
The resulting sequences were analyzed position-by-position.
________________________________________
3.1 Mutation-frequency filtering
A candidate substitution was retained for further evaluation when the substitution occurred in more than:
25%
of the generated ProteinMPNN sequences at the corresponding position.
Thus, the mutation-frequency criterion was:
mutation frequency > 0.25
________________________________________
3.2 Structural filtering
Mutation-frequency filtering was followed by structural inspection.
Candidate substitutions were evaluated for their potential to introduce or improve favorable interfacial interactions, including:
•	Hydrogen bonds
•	Salt bridges
•	Hydrophobic packing
•	Local side-chain packing
•	Inter-subunit contacts
Mutations that did not produce an apparent structural stabilization were excluded from subsequent combinatorial modeling.
________________________________________
4. Identification of PC-Int1
For the first interface-design round, the following substitutions were retained:
UL130:P65E
UL130:R66F
UL130:Q70W
UL130:T98I
Three additional substitutions proposed by ProteinMPNN were excluded after structural inspection:
UL128:K139P
UL128:I140P
UL128:T141G
These substitutions did not introduce apparent stabilizing interactions at the interface.
________________________________________
4.1 Combinatorial modeling
All two-, three-, and four-mutation combinations of the four retained substitutions were subsequently modeled using AlphaFold3.
The resulting models were evaluated using two structural criteria.
Primary criterion
The average pLDDT of the candidate mutation sites.
Secondary criterion
For models with comparable pLDDT values, candidates were further evaluated according to the backbone RMSD of the designed residues relative to the WT structure.
The highest-ranked design contained:
UL130:P65E
UL130:R66F
UL130:T98I
This design was designated:
PC-Int1
PC-Int1 was subsequently selected for experimental validation and served as the structural background for the second interface-optimization round.
________________________________________
5. Design Strategy II: UL128–UL130–UL131A Interface
To further improve pentamer stability, a second computational design round was performed at the:
UL128–UL130–UL131A ternary interface
Because the gL–UL128–UL130 interface was identified as an important structural module for pentamer stabilization, the three mutations identified in PC-Int1 were retained during the second round:
UL130:P65E
UL130:R66F
UL130:T98I
Thus, all second-round designs were generated on the PC-Int1 background.
________________________________________
6. Interface Definition
The center of the UL128–UL130–UL131A interface was defined as the:
Cβ atom of D210
in chain D (UL130) of the 5VOB structure.
The same three-step residue-selection strategy used for the first interface was applied:
1.	Spatial proximity to the interface center
2.	Participation in inter-subunit contacts
3.	Exclusion of patent-protected residues
________________________________________
7. Designable Residues
The final designable positions for the second design round were:
UL128:
86, 100, 123, 124, 125, 126, 127, 128, 129, 130

UL130:
161, 163, 164, 165, 196, 206, 207, 208, 209, 210, 211, 212, 213

UL131A:
37, 82, 83, 85, 89, 113, 118
These positions were submitted to ProteinMPNN while maintaining the PC-Int1 mutations as the structural background.
________________________________________
8. ProteinMPNN Design and Filtering
ProteinMPNN was again used to generate:
100 candidate sequences
Candidate substitutions were filtered using the same two-stage strategy:
Step 1
Mutation frequency:
>25%
Step 2
Structural inspection of the predicted substitution.
Particular attention was given to:
•	Hydrogen-bond formation
•	Salt-bridge formation
•	Hydrophobic packing
•	Local steric compatibility
•	Inter-subunit packing
________________________________________
8.1 Retained substitutions
The following substitutions were retained for combinatorial evaluation:
UL128:V128I
UL128:Y130C
UL130:V163C
UL131A:V82L
UL131A:R118M
________________________________________
8.2 Excluded substitutions
The following ProteinMPNN-derived substitutions were excluded:
UL128:A125R
UL128:G126R
UL128:S127E
UL128:P129D

UL130:N164K
UL130:F196L
UL130:F206L

UL131A:L113I
These substitutions were excluded because structural inspection suggested that they were unlikely to introduce additional favorable interfacial interactions or could potentially perturb the local packing environment.
________________________________________
9. Identification of PC-Int2
All combinations of the retained second-round mutations were modeled using AlphaFold3 together with the three stabilizing mutations from PC-Int1.
The candidate models were evaluated according to:
1.	Average pLDDT of the designed residues
2.	Local backbone RMSD of the designed residues relative to the WT structure
The highest-ranked design contained four additional substitutions:
UL128:V128I
UL128:Y130C
UL130:V163C
UL131A:V82L
on the PC-Int1 background.
The resulting design was designated:
PC-Int2
Therefore, the final PC-Int2 mutation set is:
UL130:P65E
UL130:R66F
UL130:T98I

UL128:V128I
UL128:Y130C

UL130:V163C

UL131A:V82L
________________________________________
10. Computational Workflow
The complete design workflow can be summarized as:
                         PDB 5VOB
                            │
                            ▼
                    Structural preparation
                            │
                  ┌─────────┴─────────┐
                  │                   │
                  ▼                   ▼
        gL–UL128–UL130         Structural reference
          interface
                  │
                  ▼
       Interface residue selection
                  │
                  ├── 10 Å proximity
                  ├── 4 Å contact
                  └── patent exclusion
                  │
                  ▼
             ProteinMPNN
             100 sequences
                  │
                  ▼
        Mutation frequency >25%
                  │
                  ▼
        Structural inspection
                  │
                  ▼
        Candidate substitutions
                  │
                  ▼
             AlphaFold3
                  │
                  ▼
              PC-Int1
                  │
                  │
                  ▼
        UL128–UL130–UL131A
             interface
                  │
                  ▼
       Interface residue selection
                  │
                  ▼
             ProteinMPNN
                  │
                  ▼
        Mutation frequency >25%
                  │
                  ▼
        Structural inspection
                  │
                  ▼
             AlphaFold3
                  │
                  ▼
              PC-Int2
________________________________________
11. Design Summary
Design	Structural background	Interface	Additional mutations
WT	5VOB-derived WT model	—	—
PC-Int1	WT	gL–UL128–UL130	UL130:P65E, R66F, T98I
PC-Int2	PC-Int1	UL128–UL130–UL131A	UL128:V128I, Y130C; UL130:V163C; UL131A:V82L
