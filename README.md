# StruS


## Getting Started

These instructions will get you a copy of the project up and running on your local machine for development and testing purposes.

### Prerequisites

- Python 3.11+
- g++ compilator

### Installing

A step by step series of examples that tell you how to get a development env running.

#### [Linux]

Clone the repo and create virtual environment.

```bash
git clone https://github.com/KNBi-Bioinformatyka-Strukturalna-RNA/rna-model-error-detector.git rna-model-error-detector
cd rna-model-error-detector/StruS
python -m venv .venv
source .venv/bin/activate
(.venv) pip install -r requirements.txt
```

#### [Windows]

Clone the repo and create virtual environment.

```pwsh
git clone https://github.com/KNBi-Bioinformatyka-Strukturalna-RNA/rna-model-error-detector.git rna-model-error-detector
cd rna-model-error-detector/StruS
python -m venv .venv
.venv\Scripts\activate
(.venv) pip install -r requirements.txt
```

## Running the tool

### Single prediction RTBS:

```bash
(.venv) python StruS.py RTBS target.pdb prediction.pdb
```

### Multiple predictions structRMSD:

```bash
(.venv) python StruS.py structRMSD target.pdb -p predictions
```

### Run both tools with different output folder:

```bash
(.venv) python StruS.py target.pdb prediction.pdb -o results
```

# structRMSD

StructRMSD compares structural RNA motifs between one reference RNA structure (target) and one or more predicted RNA structures.

The program performs the following steps:

1. identifies or loads structural motifs for the reference RNA,
2. evaluates the global similarity between the target and every prediction using USalign,
3. filters out predictions whose TM-score is below the selected threshold,
4. for every remaining prediction, performs local superposition and computes RMSD for each structural motif,
5. exports detailed and summary statistics as CSV files.

---

## Features

- automatic motif detection using **annotator** (rnapolis)
- motif extraction directly from **Dot-Bracket** or **BPSEQ**
- reuse of previously generated motif lists
- automatic filtering of predictions using global TM-score
- all-atom local RMSD calculation for every motif
- CSV export of detailed and summary results

## Assumptions

RNAMotifFinder assumes that:

- each execution analyses **one reference RNA structure (target)**;
- the target and all predictions contain the **same RNA chain** (default: chain `A`);
- residue numbering is identical between the target and every prediction;
- corresponding residues represent the same nucleotides.

The program matches atoms using

```
(residue_id, atom_name)
```

and therefore does **not** perform sequence alignment or residue renumbering automatically.

If residue numbering differs between the target and predictions, RMSD values may be incorrect or motifs may be skipped because matching atoms cannot be found.

To analyse multiple targets, run the program separately for each target.

---

## Input files

The following arguments are always required:

- reference RNA structure (`--target`)
- one prediction or a directory containing multiple predictions 

The target structure must contain 3D atomic coordinates because they are required for RMSD calculation.

Dot-Bracket and BPSEQ files contain only secondary structure information. They are used exclusively to identify structural motifs and do **not** replace the target PDB structure.

---

## Motif sources

Structural motifs can be obtained from four different sources.

## 1. Annotator (default)

If no additional option is provided, motifs are generated automatically using the **annotator** tool from the **rnapolis** package.

```
target.pdb
      │
      ▼
annotator
      │
      ▼
annotator.json
      │
      ▼
motif list
```

---

## 2. Previously generated motif list

Previously generated motif lists can be reused.

```
--motif-tree target.structure_tree.json
```

No motif detection is performed.

---

## 3. Dot-Bracket file

Motifs can be extracted directly from a Dot-Bracket file.

```
target.dbn
      │
      ▼
DotBracket
      │
      ▼
BpSeq
      │
      ▼
bpseq.elements
      │
      ▼
motif list
```

The program uses the same underlying mechanism (`BpSeq.elements`) as the `motif-extractor` utility from rnapolis, but accesses it directly through the Python API instead of executing the external program.

---

## 4. BPSEQ file

Motifs can also be extracted from a BPSEQ file.

```
target.bpseq
      │
      ▼
BpSeq
      │
      ▼
bpseq.elements
      │
      ▼
motif list
```

---

## Workflow

```
                     target.pdb
                          │
          ┌───────────────┼─────────────────┐
          │               │                 │
          │               │                 │
     annotator         DBN/BPSEQ     structure_tree.json
          │               │                 │
          └───────────────┴─────────────────┘
                          │
                    motif list
                          │
                          ▼
                 prediction PDB(s)
                          │
                 global TM-score
                     (USalign)
                          │
              TM-score >= threshold ?
                   │
            no ────┘
                   │
                  yes
                   │
                   ▼
        local motif superposition
       (Bio.PDB.Superimposer)
                   │
                   ▼
               motif RMSD
                   │
                   ▼
        per_motif_rmsd.csv
        motif_summary.csv
```

---

## RMSD calculation

For every structural motif the program:

1. collects all residues belonging to the motif,
2. extracts all atoms from both structures,
3. matches atoms using

```
(residue_id, atom_name)
```

4. performs optimal local superposition using the Kabsch algorithm (`Bio.PDB.Superimposer`),
5. computes all-atom RMSD.

If the fraction of matched atoms is below the selected coverage threshold (90% by default), RMSD is not calculated for that motif.

---

## Prediction filtering

Before motif analysis, every prediction is compared against the reference structure using **USalign**.

Only predictions satisfying

```
TM-score >= 0.45
```

are analysed by default.

The threshold can be changed using

```
--tm-threshold
```

---

## Command-line options

| Argument | Description |
|----------|-------------|
| `--motif-tree` | previously generated motif list |
| `--dbn` | Dot-Bracket secondary structure |
| `--bpseq` | BPSEQ secondary structure |
| `--tm-threshold` | minimum accepted TM-score (default: 0.45) |
| `--usalign-bin` | path to USalign executable |
| `--out-per-motif` | output CSV containing per-motif RMSD values |
| `--out-summary` | output CSV containing summary statistics |

---

## Output files

## `per_motif_rmsd.csv`

Each row corresponds to one pair:

```
(motif, prediction)
```

Columns:

- `motif_id`
- `motif_type`
- `residue_range`
- `prediction_file`
- `tm_score`
- `rmsd`

---

## `motif_summary.csv`

Each row corresponds to one motif.

Columns:

- `motif_id`
- `motif_type`
- `residue_range`
- `n_predictions`
- `mean_rmsd`
- `std_rmsd`

---

## Motif types

The `motif_type` column classifies motifs into one of the following categories.

| motif_type | Description |
|------------|-------------|
| `stem` | RNA helix (paired region) |
| `hairpin` | Hairpin loop |
| `single_strand` | Unpaired single-stranded region |
| `internal_loop` | Internal loop connecting two helices |
| `bulge` | Bulged loop |
| `junction_3way` | Three-way junction |
| `junction_4way` | Four-way junction |
| `junction_nway` | Multi-way junction with *n* branches |

Loop motifs are classified automatically according to the number and length of their strands.

---

## Notes

- Structural motifs can be obtained from annotator, Dot-Bracket, BPSEQ or a previously generated motif list.
- Dot-Bracket and BPSEQ files are used **only to identify structural motifs**.
- RMSD is always calculated from 3D atomic coordinates stored in the target and prediction structures.
- The program assumes consistent residue numbering between the target and all predictions.
- The program performs local superposition independently for every motif using the Kabsch algorithm implemented in BioPython.
- All atoms belonging to a motif are included in the RMSD calculation (all-atom RMSD).
- Predictions whose global TM-score is below the selected threshold are excluded before motif-level analysis.