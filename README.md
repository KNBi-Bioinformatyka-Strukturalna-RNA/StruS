# StruS


## Getting Started

These instructions will get you a copy of the project up and running on your local machine for development and testing purposes.

### Prerequisites

- Python 3.11+
- g++ compilator

### Installing

A step by step series of examples that tell you how to get a development env running.

#### [Linux]

Clone the repo and run the installer.

```bash
git clone [https://github.com/KNBi-Bioinformatyka-Strukturalna-RNA/rna-model-error-detector.git](https://github.com/KNBi-Bioinformatyka-Strukturalna-RNA/StruS) StruS
cd StruS
chmod +x install_strus.sh StruS
./install_strus.sh
```

#### [Windows]

Clone the repo and create virtual environment.

```pwsh
git clone [https://github.com/KNBi-Bioinformatyka-Strukturalna-RNA/rna-model-error-detector.git](https://github.com/KNBi-Bioinformatyka-Strukturalna-RNA/StruS) StruS
cd StruS
python -m venv .venv
.venv\Scripts\activate
(.venv) pip install -r requirements.txt
```

## Running the tool

### Single prediction RTBS:

```bash
StruS RTBS target.pdb prediction.pdb
```

### Multiple predictions structRMSD:

```bash
StruS structRMSD target.pdb -p predictions
```

### Run both tools with different output folder:

```bash
StruS target.pdb prediction.pdb -o results
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


# RTBS

Compares a target RNA secondary structure tree with a predicted tree, annotates each matched node of the prediction with a penalty score derived from the real MBR (Matrix of BEAR-encoded RNA secondary structures) from Mattei et al. 2014 (doi:10.1093/nar/gku283). Calculates two versions of the RTBS measure, which determines the structural similarity of structures based on a tree structure, and creates a visualization of penalties for structural differences between structures for each node of the tree of structural elements.

---
The program performs the following steps:
1. Finds structural elements in a PDB file (both target and prediction) using the RNApolis annotator tool.
2. Creates a graph (usually a tree structure) based on structural elements, indicating their relative positions.
3. Compares the target and prediction graphs, trying to match as many nodes as possible.
4. Based on the node matching, it calculates penalties for differences between them based on the BEAR matrix.
5. Sums penalties from all nodes and calculates normalized measures in two variants: target-oriented and symmetrical.
6. Returns the results of the measures and a visualization of the node matching and the penalties imposed on them.

If RTBS is run in multi-prediction mode, it also returns a text file with the prediction ranking according to the RTBS target-oriented measure.

---
## Penalty system
A penalty is assigned to a node based on the BEAR matrix (doi:10.1093/nar/gku283). If the prediction node is not matched with any node in the target, the highest possible penalty is assigned, i.e., UNMATCHED PENALTY (approximately 3.68). 

If matched:

Penalty = -MBR_score(BEAR_char_target, BEAR_char_pred)
 - negative MBR scores become positive penalties (bad)
 - positive MBR scores become negative penalties (reward, node turns green)

MBR (Matrix of BEAR-encoded RNA secondary structures) is a substitution matrix that assigns scores to substitutions between RNA secondary-structure elements encoded using BEAR, based on transition rates observed in related RNA families.

BEAR alphabet:
```
  STEM (non-branching):             a b c d e f g h i =   (lengths 1–10)
  STEM (branching/multi-loop):      A B C D E F G H I J   (lengths 1–10)
  LOOP (hairpin):                   j k l m n o p q r s t u v w x y z ^  (lengths 3–18+)
  LEFTINTERNALLOOP (non-branch):    ! " # $ % & ' ( ) +   (lengths 1–11+)
  RIGHTINTERNALLOOP (non-branch):   2 3 4 5 6 7 8 9 0 >   (lengths 1–11+)
  LEFTINTERNALLOOP (branch):        K L M N O P Q R S T U V W   (lengths 2–14+)
  RIGHTINTERNALLOOP (branch):       Y Z ~ _ | / \\ @         (lengths 2–10+)
  BULGELEFT:                        [
  BULGERIGHT:                       ]
  BULGELFETBRANCH:                  {
  BULGERIGTHBRANCH:                 }
  UNSTRUCTURED (singleStrand):                     :
  UNKNOWN:                          ?
```

The exact penalties charged for different nodes are in the mbr_matrix.json file.

---

## RTBS normalisation:
```
RTBS = 1 - (sum_penalty − best_possible) / (worst_possible − best_possible)

sum_penalty = the sum of penalties from all matched nodes and the penalty for mismatches from the target (tgt) or from the target and prediction (sym)

best_possible  = sum of diagonal MBR penalties over TARGET nodes
                (score of a perfect prediction reproducing the target exactly)
```

for RTBS ​​target:

```
worst_possible = n_target * UNMATCHED_PENALTY
(counts each missing target from prediction)

```
Unmatched prediction nodes are NOT penalised explicitly — their cost appears implicitly because target nodes they "replaced" are counted as missing.

for RTBS symmetrical:

```
worst_possible = (n_target + n_prediction) * UNMATCHED_PENALTY/2
(counts each missing target node and each additional node from prediction)
```
The mismatch cost is distributed equally across the nodes from the target and prediction.
This measure better represents the structural similarity of two structures, but it does not allow for direct comparison of different predictions of a single target based on it, because it depends on the size of the prediction tree's structural elements.

---
## RTBS interpretation
RTBS after normalization takes values ​​from 0 to 1. Interpretation:
- Values ​​close to 1 indicate very close structural similarity.
- Values ​​close to 0 indicate very little or no similarity.

---

## Command-line options

| Argument           | Description                                                                         |
| ------------------ | ----------------------------------------------------------------------------------- |
| `target`           | Target RNA structure (`.pdb`)                                                       |
| `prediction`       | Single prediction structure (`.pdb`). Do not use together with `-p` / `--pred_dir`. |
| `-p`, `--pred_dir` | Directory containing prediction PDB files.                                          |
| `-o`, `--out_dir`  | Working output directory (default: `StruS_out`)                                     |
| `-m`, `--mbr`      | Path to the MBR JSON matrix (default: `mbr_matrix.json`)                            |

---

## Input files

The following arguments are always required:

- reference RNA structure (`--target`)
- one prediction or a directory containing multiple predictions 

The target structure must contain 3D atomic coordinates because they are required for RMSD calculation.

Dot-Bracket and BPSEQ files contain only secondary structure information. They are used exclusively to identify structural motifs and do **not** replace the target PDB structure.

---
## Output files

The output directory is `StruS_out` by default and can be changed using the `--out_dir` argument.

### Directory structure

```text
StruS_out/
├── annotated/
│   ├── target/
│   │   └── <target>.json
│   └── <prediction_1>.json
│   └── <prediction_2>.json
│   └── ...
│
├── converted/
│   ├── target/
│   │   └── <target>.json
│   └── <prediction_1>.json
│   └── <prediction_2>.json
│   └── ...
│
└── RTBS_results/
    ├── ranking.txt
    └── visualizations/
        ├── <prediction_1>_vs_<target>.png
        ├── <prediction_1>_vs_<target>.json
        ├── <prediction_2>_vs_<target>.png
        ├── <prediction_2>_vs_<target>.json
        └── ...
```

### `annotated/`

Contains the JSON files generated by the **RNApolis annotator** for the target structure and all prediction structures.

The target annotation is stored separately in the `target/` subdirectory.

### `converted/`

Contains the structural graphs.

Each JSON file represents the structural elements of the corresponding RNA structure, including their relationships in the tree.

The target graph is stored separately in the `target/` subdirectory.

### `RTBS_results/`

Contains the RTBS results.

#### `ranking.txt`

In directory mode, this file contains the RTBS ranking of all predictions:

```text
Prediction                                 RTBS_tgt   RTBS_sym
----------------------------------------------------------------
<prediction_1>                                0.9364    0.4164
<prediction_2>                                0.9243    0.3565
```

The columns are:

* `Prediction` — prediction structure name.
* `RTBS_tgt` — target-oriented RTBS measure.
* `RTBS_sym` — symmetric RTBS measure.

The measures are described in detail in the following section.

#### `visualizations/`

Contains a visualization for each prediction showing the structural-tree matching between the prediction and the target, including the assigned penalties.

For each prediction, two files are generated:

* `.png` — graphical visualization of the matching and penalties.
* `.json` — detailed information used to generate the visualization, including matched and unmatched nodes, BEAR values, and assigned penalties.

### Single-prediction mode

When a single prediction is provided instead of a prediction directory, the output structure is slightly different:

```text
StruS_out/
├── annotated/
│   └── target/
│       └── <target>.json
│   └── <prediction>.json
│
├── converted/
│   └── target/
│       └── <target>.json
│   └── <prediction>.json
│
└── RTBS_results/
    ├── <prediction>_vs_<target>.png
    └── <prediction>_vs_<target>.json
```

In this mode:

* no `visualizations/` subdirectory is created;
* visualization files are stored directly in `RTBS_results/`;
* no `ranking.txt` file is generated.
