#!/usr/bin/env python3

import argparse
import subprocess
import sys
from pathlib import Path
from structRMSD import struct_rmsd_main
from config import (PYTHON_BIN, ANNOTATOR, CONVERTER, RTBS, MBR, FR3D, FR3D_TO_DBN, CONVERTER_FROM_DBN, ANNOTATOR_TO_DBN)

def parse_args():
    parser = argparse.ArgumentParser(
        prog="StruS.py",
        description="""
StruS - Structural RNA evaluation tool
Runs RTBS and/or structRMSD.

RNA Tree BEAR Similarity (RTBS) with RNA Tree Penalty Visualizer:
  Calculates two versions of the RTBS measure,
  which determines the structural similarity of structures based on a tree structure,
  and creates a visualization of penalties for structural differences between structures
  for each node of the tree of structural elements.
  Compares a target RNA secondary structure tree with a predicted tree,
  annotates each matched node of the prediction with a penalty score
  derived from the real MBR (Matrix of BEAR-encoded RNA secondary structures)
  from Mattei et al. 2014 (doi:10.1093/nar/gku283).

structRMSD:
  Evaluates the local structural accuracy of predicted RNA 3D models by comparing them to a target structure motif by motif, 
  rather than as a whole. Predictions are first filtered by their global TM-score, 
  computed via USalign (https://github.com/pylelab/USalign), rejecting those below a chosen similarity threshold (default 0.45). 
  For each retained prediction, every structural motif of the target — stems, hairpins, internal loops/bulges, 
  multi-branch junctions, and single strands, 
  identified via the rnapolis annotator (or supplied directly as a Dot-Bracket/BPSEQ secondary structure) — 
  is extracted by residue range and locally superposed onto its counterpart in the prediction using the Kabsch algorithm 
  (all-atom, Bio.PDB.Superimposer), yielding a per-motif RMSD. 
  Results are aggregated across all retained predictions into a mean and standard deviation RMSD for each motif, 
  highlighting which structural elements are reliably modeled and which are not.

Usage:
  Single prediction RTBS:
        StruS RTBS target.pdb prediction.pdb

  Multiple predictions structRMSD:
        StruS structRMSD target.pdb -p predictions/

  Run both tools with different output folder:
        StruS BOTH target.pdb prediction.pdb -o results/

"""
    )

    parser.add_argument("tool", choices=["RTBS", "structRMSD", "BOTH"], help="Tool to run.")
    parser.add_argument("target", help="Target RNA structure (.pdb)")
    parser.add_argument("prediction", nargs="?", default=None, help="Single prediction (.pdb). Do not use together with -p.")
    parser.add_argument("-p", "--pred_dir", default=None, help="Directory containing prediction PDB files.")
    parser.add_argument("-o", "--out_dir", default="StruS_out", help="Working output directory (default: StruS_out)")
    parser.add_argument("-m", "--mbr", default=MBR, help="Path to MBR JSON matrix (default: mbr_matrix.json next to this script) for RTBS.")
    parser.add_argument("-t", "--threshold", default=2, help="Number of matching nodes in the subtree to match the prediction and target in RTBS.")
    parser.add_argument("--annotator", action="store_true", help="Use the rnapolis annotator to find structural elements.")
    parser.add_argument("--fr3d", action="store_true", help="Use FR3D to find structural elements.")
    parser.add_argument("--dbn_rtbs", action="store_true", help="Use existing .dbn files to calculate RTBS.")
    parser.add_argument("--check_sequence_always", action="store_true", help="Require sequence compatibility when matching subtrees in RTBS, not only when matching remaining single nodes. Prevents accidental matching of similar nodes. Only possible for target comparison and prediction for the same sequence.")
    parser.add_argument("--remove_pseudoknots", action="store_true", help="Remove pseudoknots before processing. All characters other than '(' and ')' are permanently replaced with '.'.")
    parser.add_argument("--motif-tree", type=Path, default=None)
    parser.add_argument("--dbn", type=Path, default=None)
    parser.add_argument("--bpseq", type=Path, default=None)
    parser.add_argument("--tm-threshold", type=float, default=0.45)
    parser.add_argument("--target-chain", type=str, default=None)
    parser.add_argument("--usalign-bin", type=str, default=None)
    parser.add_argument("--out-per-motif", type=Path, default=Path("per_motif_rmsd.csv"))
    parser.add_argument("--out-summary", type=Path, default=Path("motif_summary.csv"))
    args = parser.parse_args()

    if sum([args.annotator, args.fr3d, args.dbn_rtbs]) > 1:
        parser.error("Use only one of --annotator, --fr3d, --dbn_rtbs for RTBS.")

    return args


def run_command(cmd, quiet=True):
    print("Running:")
    print(" ".join(map(str, cmd)))
    if quiet:
        result = subprocess.run(cmd,stdout=subprocess.DEVNULL)
    else:
        result = subprocess.run(cmd)

    if result.returncode != 0:
        raise RuntimeError(f"Command failed with exit code {result.returncode}")


def check_pdb(path):
    path = Path(path)
    if not path.exists():
        raise FileNotFoundError(path)
    if path.suffix.lower() != ".pdb":
        raise ValueError(f"{path} is not a .pdb file")
    return path


def check_dbn(path):
    path = Path(path)
    if not path.exists():
        raise FileNotFoundError(path)
    if path.suffix.lower() != ".dbn":
        raise ValueError(f"{path} is not a .dbn file")
    return path


def annotate_pdb(pdb_file, output_dir):
    output_dir.mkdir(parents=True, exist_ok=True)
    out_json = output_dir / f"{pdb_file.stem}.json"
    cmd = [PYTHON_BIN, ANNOTATOR, str(pdb_file), "--json", str(out_json)]
    run_command(cmd)
    print("Annotated")
    return out_json


def convert_annotation(json_file, output_dir):
    output_dir.mkdir(parents=True, exist_ok=True)
    out_name = output_dir / f"{json_file.stem}.json"
    cmd = [PYTHON_BIN, CONVERTER, str(json_file), "-o", str(out_name)]
    run_command(cmd)
    print("Converted")
    return out_name


def run_fr3d(pdb_file, output_dir):
    output_dir.mkdir(parents=True, exist_ok=True)
    cmd = [PYTHON_BIN, FR3D, str("./"+str(pdb_file)), "-o", str(output_dir)]
    run_command(cmd, quiet=True)
    print("FR3D annotated")
    basepair_file = output_dir / f"{pdb_file.stem}_basepair.txt"
    if not basepair_file.exists():
        raise FileNotFoundError(f"Expected FR3D output not found: {basepair_file}")
    return basepair_file


def fr3d_to_dbn(basepair_file, pdb_file, output_dir):
    output_dir.mkdir(parents=True, exist_ok=True)
    out_file = output_dir / f"{pdb_file.stem}.dbn"
    cmd = [PYTHON_BIN, FR3D_TO_DBN, str(basepair_file), str(pdb_file), "-o", str(out_file)]
    run_command(cmd)
    print("Converted FR3D basepairs to DBN")
    return out_file


def annotator_to_dbn(json_file, output_dir):
    output_dir.mkdir(parents=True, exist_ok=True)
    out_file = output_dir / f"{json_file.stem}.dbn"
    cmd = [PYTHON_BIN, ANNOTATOR_TO_DBN, str(json_file), "-o", str(out_file)]
    run_command(cmd)
    print("Converted annotator output to DBN")
    return out_file


def convert_from_dbn(dbn_file, output_dir, remove_pseudoknots=False):
    output_dir.mkdir(parents=True, exist_ok=True)
    out_name = output_dir / f"{Path(dbn_file).stem}.json"
    if remove_pseudoknots:
        cmd = [PYTHON_BIN, CONVERTER_FROM_DBN, str(dbn_file), "--output", str(out_name), "--remove_pseudoknots"]
    else:
        cmd = [PYTHON_BIN, CONVERTER_FROM_DBN, str(dbn_file), "--output", str(out_name)]
    run_command(cmd, quiet=False)
    print("Converted DBN")
    return out_name


def run_rtbs(target_json, prediction_jsons, pred_dir, output_dir, mbr, check_sequence_always=False, treshold=2):
    output_dir.mkdir(parents=True, exist_ok=True)
    if pred_dir is None:
        pred_json = prediction_jsons[0]
        out_file = output_dir / pred_json.stem
        cmd = [PYTHON_BIN, RTBS, str(target_json), str(pred_json), "-o", str(out_file), "-m", str(mbr), "-t", str(treshold)]
    else:
        cmd = [PYTHON_BIN, RTBS, str(target_json), "-p", str(prediction_jsons[0].parent), "-o", str(output_dir), "-m", str(mbr), "-t", str(treshold)]
    if check_sequence_always:
        cmd.append("--check_sequence_always")
    run_command(cmd, quiet=False)

def run_structRMSD(args, output_dir):
    output_dir.mkdir(parents=True, exist_ok=True)
    struct_rmsd_main(output_dir, args)


def execute_rtbs(target, predictions, pred_dir, workdir, mbr, source="fr3d", check_sequence_always=False, remove_pseudoknots = False, threshold = 2):
    print("\n=====RTBS=====")
    rtbs_dir = workdir / "RTBS_results"

    if source == "annotator":
        annotator_dir = workdir / "annotated_rnapolis"
        dbn_dir = workdir / "dbn_annotator"
        converted_dir = workdir / "converted_annotator"
        annotator_target_dir = annotator_dir / "target"
        dbn_target_dir = dbn_dir / "target"
        converted_target_dir = converted_dir / "target"

        print("\nAnnotating target:")
        target_json = annotate_pdb(target, annotator_target_dir)
        target_dbn = annotator_to_dbn(target_json, dbn_target_dir)
        converted_target = convert_from_dbn(target_dbn, converted_target_dir, remove_pseudoknots)

        prediction_jsons = []
        for pred in predictions:
            print(f"\nAnnotating prediction {pred.name}:")
            pred_json = annotate_pdb(pred, annotator_dir)
            pred_dbn = annotator_to_dbn(pred_json, dbn_dir)
            converted = convert_from_dbn(pred_dbn, converted_dir, remove_pseudoknots)
            prediction_jsons.append(converted)

    elif source == "fr3d":
        annotated_dir = workdir / "annotated_fr3d"
        dbn_dir = workdir / "dbn_fr3d"
        converted_dir = workdir / "converted_fr3d"
        annotated_target_dir = annotated_dir / "target"
        dbn_target_dir = dbn_dir / "target"
        converted_target_dir = converted_dir / "target"

        print("\nRunning FR3D on target:")
        target_bp = run_fr3d(target, annotated_target_dir)
        target_dbn = fr3d_to_dbn(target_bp, target, dbn_target_dir)
        converted_target = convert_from_dbn(target_dbn, converted_target_dir, remove_pseudoknots)

        prediction_jsons = []
        for pred in predictions:
            print(f"\nRunning FR3D on prediction {pred.name}:")
            pred_bp = run_fr3d(pred, annotated_dir)
            pred_dbn = fr3d_to_dbn(pred_bp, pred, dbn_dir)
            converted = convert_from_dbn(pred_dbn, converted_dir, remove_pseudoknots)
            prediction_jsons.append(converted)

    elif source == "dbn":
        converted_dir = workdir / "converted_dbn"
        converted_target_dir = converted_dir / "target"

        print("\nConverting target DBN:")
        converted_target = convert_from_dbn(target, converted_target_dir, remove_pseudoknots)

        prediction_jsons = []
        for pred in predictions:
            print(f"\nConverting prediction DBN {pred.name}:")
            converted = convert_from_dbn(pred, converted_dir, remove_pseudoknots)
            prediction_jsons.append(converted)

    else:
        raise ValueError(f"Unknown RTBS source: {source}")

    print("\nCalculating RTBS:")
    run_rtbs(converted_target, prediction_jsons, pred_dir, rtbs_dir, mbr, check_sequence_always=check_sequence_always, treshold=threshold)
    print("RTBS calculated")


def execute_struct_rmsd(args, workdir):
    print("\n=====structRMSD=====")
    rtbs_dir = workdir / "structRMSD_results"
    run_structRMSD(args, rtbs_dir)


def main():
    args = parse_args()
    workdir = Path(args.out_dir)
    workdir.mkdir(parents=True, exist_ok=True)

    if args.annotator:
        rtbs_source = "annotator"
    elif args.dbn_rtbs:
        rtbs_source = "dbn"
    else:
        rtbs_source = "fr3d"

    if args.tool == "BOTH" and rtbs_source == "dbn":
        raise RuntimeError("--dbn_rtbs cannot be combined with BOTH: structRMSD needs .pdb structures while RTBS --dbn_rtbs needs .dbn files.")

    dbn_mode = args.tool == "RTBS" and rtbs_source == "dbn"
    check_file = check_dbn if dbn_mode else check_pdb
    glob_pattern = "*.dbn" if dbn_mode else "*.pdb"

    target = check_file(args.target)
    mbr_path = Path(args.mbr)
    if not mbr_path.is_file():
        raise FileNotFoundError(f"MBR matrix file not found: {mbr_path}")
    pred_dir = None 

    if args.pred_dir:
        pred_dir = Path(args.pred_dir)
        if not pred_dir.exists():
            raise FileNotFoundError(pred_dir)

        predictions = sorted(pred_dir.glob(glob_pattern))
        if not predictions:
            raise RuntimeError(f"No {glob_pattern} files found in prediction directory")

    elif args.prediction:
        predictions = [check_file(args.prediction)]

    else:
        raise RuntimeError("Provide a prediction file or -p prediction_folder")

    if args.tool == "BOTH":
        execute_rtbs(target, predictions, pred_dir, workdir, args.mbr, source=rtbs_source, check_sequence_always=args.check_sequence_always, remove_pseudoknots = args.remove_pseudoknots, threshold=args.threshold)
        execute_struct_rmsd(args, workdir)

    elif args.tool == "RTBS":
        execute_rtbs(target, predictions, pred_dir, workdir, args.mbr, source=rtbs_source, check_sequence_always=args.check_sequence_always, remove_pseudoknots = args.remove_pseudoknots, threshold=args.threshold)

    elif args.tool == "structRMSD":
        execute_struct_rmsd(args, workdir)


if __name__ == "__main__":
    main()