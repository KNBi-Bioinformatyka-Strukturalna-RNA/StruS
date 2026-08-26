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
git clone https://github.com/KNBi-Bioinformatyka-Strukturalna-RNA/StruS StruS
cd StruS
git submodule init
git submodule update --recursive
chmod +x install_strus.sh StruS
./install_strus.sh
```

#### [Windows]

Clone the repo and create virtual environment.

```pwsh
git clone https://github.com/KNBi-Bioinformatyka-Strukturalna-RNA/StruS StruS
cd StruS
git submodule init
git submodule update --recursive
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
- automatic motif detection using **FR3D** (fr3d-python, `latest` branch)
- motif extraction directly from **Dot-Bracket** or **BPSEQ**
- reuse of previously generated motif lists
- automatic filtering of predictions using global TM-score
- all-atom local RMSD calculation for every motif, with an adjustable coverage threshold
- local TM-score calculation for every motif (where the motif is long enough)
- optional removal of isolated base pairs before motif detection
- support for multi-chain targets and predictions, including fully manual chain definitions
- CSV export of detailed and summary results
- normalized RMSD (nRMSD = RMSD / √(number of atoms used)) alongside raw RMSD, to account for motifs of different sizes
- full descriptive statistics (mean, std, min, Q1, median, Q3, max) for RMSD, nRMSD and TM-score, both per individual motif and aggregated by motif type
- optional pseudoknot-aware decomposition, separating pseudoknot-forming stems into their own category instead of forcing them into the standard stem/hairpin/loop scheme

## Assumptions

structRMSD assumes that:

- each execution analyses **one reference RNA structure (target)**;
- residue numbering is identical between the target and every prediction;
- corresponding residues represent the same nucleotides.

The program matches atoms using

```
(residue_id, atom_name)
```

and therefore does **not** perform sequence alignment or residue renumbering automatically.

If residue numbering differs between the target and predictions, RMSD values may be incorrect or motifs may be skipped because matching atoms cannot be found.

**Multi-chain targets and predictions are supported, and are never reduced automatically.** By default, the full target and the full prediction are always compared as-is - motifs (and, for `--dbn`, whole motif ranges) that don't fit within a given prediction's residue range are simply skipped individually rather than causing an error. If either the target or a prediction has more than one chain, the program prints which chains it found and reminds you that you can select just one.

To restrict the target to a single chain, use `--target-chain`. To restrict predictions to a single chain (applied uniformly to every prediction in the run), use `--prediction-chain`. Both accept the original chain letter from the file - or, if `--chain-mapping` was also used, the newly assigned chain letter from that remapping (chain selection is always applied after chain mapping, on its result).

To analyse multiple targets, run the program separately for each target.

---

## Input files

The following arguments are always required:

- reference RNA structure (`--target`)
- one prediction or a directory containing multiple predictions (`--prediction`)

The target structure must contain 3D atomic coordinates because they are required for RMSD calculation.

Dot-Bracket and BPSEQ files contain only secondary structure information. They are used exclusively to identify structural motifs and do **not** replace the target PDB structure.

---

## Motif sources

Structural motifs can be obtained from five different sources. Exactly one source must be selected; if none is given, the program defaults to **FR3D**.

## 1. FR3D (default)

If no other motif source is given, motifs are generated automatically using **FR3D** (`fr3d-python`, `latest` branch).

```
target.pdb
      │
      ▼
FR3D (NA_pairwise_interactions)
      │
      ▼
basepairs list
      │
      ▼
keep only canonical cWW A-U / G-C / G-U pairs
      │
      ▼
BPSEQ
      │
      ▼
motif list
```

## 2. Annotator 

```
--annotator
```

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
    BPSEQ
      │
      ▼
  motif list
```

---

## 3. Previously generated motif list

Previously generated motif lists can be reused.

```
--motif-tree target.structure_tree.json
```

No motif detection is performed. `--remove-isolated` has no effect here

---

## 4. Dot-Bracket file

Motifs can be extracted directly from a Dot-Bracket file.

```
--dbn target.dbn
```

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

Multi-chain Dot-Bracket files are supported: chains are detected and, if `--target-chain` is given, motifs spanning more than one chain are dropped and the rest are renumbered to that chain. Without `--target-chain`, the full, multi-chain target is used and cross-chain motifs are kept as-is.

---

## 5. BPSEQ file

Motifs can also be extracted from a BPSEQ file.

```
--bpseq target.bpseq
```

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

An externally supplied BPSEQ file describes one fixed numbering of the target as a whole. `--target-chain` has no effect with this source.

---

## Removing isolated base pairs

```
--remove-isolated
```

A single, unstacked base pair (a "stem" of length 1) is often geometric noise rather than a real structural element. With `--remove-isolated`, such pairs are dropped before the target is decomposed into motifs, and the decomposition is recomputed on the cleaned-up pairing table - neighbouring motifs are merged or resized accordingly.

This works for `--annotator`, `--dbn`, `--bpseq` and `--fr3d` (all four are internally converted to a `BpSeq` object before motif extraction). It has no effect with `--motif-tree`, since an already-computed motif list cannot be recomputed - a warning is printed instead.

---

## Manual chain mapping

```
--chain-mapping "t:A21-40;A1-20;;prediction_1.pdb:A1-20;B1-20"
```

lets you define the chain layout explicitly, per file.

Syntax:

- `;;` separates the definition blocks for different files.
- Each block starts with a file identifier followed by `:` - `t` for the target, or the exact filename for a prediction.
- Within a block, `;` separates logical chains, numbered in the order given (first one becomes the new chain A, second one chain B, and so on).
- Within a single logical chain, `+` joins segments (even from different original chains) into one, in the given order.
- Each segment is written as `<original_chain><first>-<last>` (e.g. `A21-40`).

Rules:

- A file not mentioned in `--chain-mapping` at all is left completely untouched.
- A file that **is** mentioned has its chain layout **fully replaced** by the given definition - anything not explicitly listed for that file is dropped, even if it belongs to a chain that was partially mentioned.
- Predictions not covered by `--chain-mapping` still go through the normal, automatic position-based matching against the (possibly remapped) target - there is no automatic verification that their chain order actually corresponds to the target's manual definition.
- `--target-chain` and `--prediction-chain` are applied *after* `--chain-mapping` - they refer to the newly created chain letters (A, B, ...), not the original ones.

## Workflow

```
                                 target.pdb
                                     │
           ┌───────────┬─────────────┼──────────────┬──────────────┐
           │           │             │              │              │
         FR3D      annotator      DBN/BPSEQ   structure_tree.json  │
           │           │             │              │              │
           └───────────┴─────────────┴──────────────┴──────────────┘
                                     │
                     (optional: remove isolated pairs)
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
                 motif RMSD + motif TM-score
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

If the fraction of matched atoms is below the coverage threshold, RMSD is not calculated for that motif. The threshold defaults to 90% and can be changed with `--min-coverage` (e.g. `--min-coverage 0.5`).

---

## Motif TM-score calculation

In addition to RMSD, the program also attempts to compute a local TM-score for every motif, using USalign on just that motif's residues (extracted to a temporary fragment PDB for both target and prediction).

This requires at least 3 residues in the motif - USalign cannot produce a usable result on shorter fragments. Motifs with fewer than 3 residues have an `n/a` `motif_tm_score` value, while `motif_rmsd` is still calculated normally for them.

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

## USalign

The repository includes a precompiled `USalign` binary, but it is compiled for **Linux**. On **macOS** or **Windows** this binary will not run - `USalign.cpp` must be compiled locally instead (or let the program do it automatically, see below).

If `USalign` is not found (as an executable in `PATH`, or as `./USalign` in the current directory), the program downloads `USalign.cpp` from the official source and compiles it automatically with `g++` on first run. This requires an internet connection and a working C++ compiler. The resulting `./USalign` binary is then reused on later runs.

A specific binary can also be provided manually with `--usalign-bin`.

---

## Command-line options

| Argument | Description |
|----------|-------------|
| `--target` | reference RNA structure (required) |
| `--prediction` | single prediction file or a directory of predictions (required) |
| `--target-chain` | chain to reduce the target to, for multi-chain targets |
| `--prediction-chain` | chain to reduce every prediction to, for multi-chain predictions |
| `--chain-mapping` | manual, per-file chain layout definition |
| `--motif-tree` | previously generated motif list |
| `--dbn` | Dot-Bracket secondary structure |
| `--bpseq` | BPSEQ secondary structure |
| `--annotator` | detect motifs using the rnapolis annotator |
| `--fr3d` | detect motifs using FR3D (default if no other source is given) |
| `--remove-isolated` | drop isolated base pairs before motif detection |
| `--decompose-pseudoknot-free` | separate pseudoknot-forming stems into their own category before the rest of the target is decomposed |
| `--tm-threshold` | minimum accepted TM-score (default: 0.45) |
| `--min-coverage` | minimum fraction of matched atoms required to compute a motif's RMSD (default: 0.9) |
| `--usalign-bin` | path to USalign executable (auto-detected/downloaded if not given) |
| `--out-per-motif` | output CSV containing per-motif results |
| `--out-summary` | output CSV containing summary statistics |
| `--out-by-type` | output CSV containing summary statistics aggregated by motif type |

Exactly one of `--motif-tree`, `--dbn`, `--bpseq`, `--annotator`, `--fr3d` may be given at a time.
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
- `motif_tm_score` - TM-score of this motif alone 
- `motif_rmsd` - RMSD of this motif alone
- `motif_nrmsd` - RMSD normalized by the square root of the number of atoms used

Only predictions that passed the global TM-score filter appear in this file.

---

## `motif_summary.csv`

Each row corresponds to one motif.

Columns:

- `motif_id`
- `motif_type`
- `residue_range`
- `n_predictions` - number of predictions for which RMSD was actually calculated for this motif
- for each of `rmsd`, `nrmsd` and `tm_score`: `mean_*`, `std_*`, `min_*`, `q1_*`, `median_*`, `q3_*`, `max_*` (`N/A` if no prediction produced a result for this motif)

---

## `motif_type_summary.csv`

Each row corresponds to one motif **type** (e.g. all stems across the whole target pooled together, regardless of which specific stem they came from).

Columns:

- `motif_type`
- `n_values` - number of (motif, prediction) pairs pooled into this row
- the same set of statistics as `motif_summary.csv` (mean/std/min/Q1/median/Q3/max) for `rmsd`, `nrmsd` and `tm_score`

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
| `pseudoknot_stem` | A stem-forming base pair that crosses another pair (a pseudoknot); only reported separately when `--decompose-pseudoknot-free` is used |

Loop motifs are classified automatically according to the number and length of their strands.

---

## Notes

- Structural motifs can be obtained from FR3D, annotator, Dot-Bracket, BPSEQ, or a previously generated motif list.
- If no motif source is given, the program defaults to FR3D.
- Dot-Bracket and BPSEQ files are used **only to identify structural motifs**.
- RMSD is always calculated from 3D atomic coordinates stored in the target and prediction structures.
- The program assumes consistent residue numbering between the target and all predictions.
- The program performs local superposition independently for every motif using the Kabsch algorithm implemented in BioPython.
- All atom **types** are included in the RMSD calculation (all-atom RMSD, not restricted to backbone or any other subset) - but for any given motif, only the atoms present in **both** the target and the prediction are actually superposed.
- Predictions whose global TM-score is below the selected threshold are excluded before motif-level analysis.
- If a prediction's total residue count doesn't match the target's, motifs whose residues extend beyond the prediction's range are skipped individually (with a message) rather than skipping the whole prediction.


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
