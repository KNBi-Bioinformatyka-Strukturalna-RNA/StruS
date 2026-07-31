#!/usr/bin/env python3

import argparse
import subprocess
from pathlib import Path
from structRMSD import struct_rmsd_main

ANNOTATOR = "../../rnapolis-py/src/rnapolis/annotator.py"
CONVERTER = "annotation_converter.py"
RTBS = "RTBS.py"

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
  blebleble

Usage:
  Single prediction RTBS:
    python3 StruS.py RTBS target.pdb prediction.pdb

  Multiple predictions structRMSD:
    python3 StruS.py structRMSD target.pdb -p predictions/

  Run both tools with different output folder:
    python3 StruS.py target.pdb prediction.pdb -o results/

"""
    )

    parser.add_argument("tool", nargs="?", choices=["RTBS", "structRMSD"], help="Tool to run. If omitted, both tools are executed.")
    parser.add_argument("target", help="Target RNA structure (.pdb)")
    parser.add_argument("prediction", nargs="?", default=None, help="Single prediction (.pdb). Do not use together with -p.")
    parser.add_argument("-p", "--pred_dir", default=None, help="Directory containing prediction PDB files.")
    parser.add_argument("-o", "--out_dir", default="StruS_out", help="Working output directory (default: StruS_out)")
    parser.add_argument("-m", "--mbr", default="mbr_matrix.json", help="Path to MBR JSON matrix (default: mbr_matrix.json next to this script)")
    parser.add_argument("--motif-tree", type=Path, default=None)
    parser.add_argument("--dbn", type=Path, default=None)
    parser.add_argument("--bpseq", type=Path, default=None)
    parser.add_argument("--tm-threshold", type=float, default=0.45)
    parser.add_argument("--usalign-bin", type=str, default=None)
    parser.add_argument("--out-per-motif", type=Path, default=Path("per_motif_rmsd.csv"))
    parser.add_argument("--out-summary", type=Path, default=Path("motif_summary.csv"))
    return parser.parse_args()


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


def annotate_pdb(pdb_file, output_dir):
    output_dir.mkdir(parents=True, exist_ok=True)
    out_json = output_dir / f"{pdb_file.stem}.json"
    cmd = ["python3", ANNOTATOR, str(pdb_file), "--json", str(out_json)]
    run_command(cmd)
    print("Annotated")
    return out_json


def convert_annotation(json_file, output_dir):
    output_dir.mkdir(parents=True, exist_ok=True)
    out_name = output_dir / f"{json_file.stem}.json"
    cmd = ["python3", CONVERTER, str(json_file), "-o", str(out_name)]
    run_command(cmd)
    print("Converted")
    return out_name


def run_rtbs(target_json, prediction_jsons, pred_dir, output_dir, mbr):
    output_dir.mkdir(parents=True, exist_ok=True)
    if pred_dir is None:
        pred_json = prediction_jsons[0]
        out_file = output_dir / pred_json.stem
        cmd = ["python3", RTBS, str(target_json), str(pred_json), "-o", str(out_file), "-m", str(mbr)]
        run_command(cmd, quiet=False)
    else:
        cmd = ["python3", RTBS, str(target_json), "-p", str(prediction_jsons[0].parent), "-o", str(output_dir), "-m", str(mbr)]
        run_command(cmd, quiet=False)

def run_structRMSD(args, output_dir):
    output_dir.mkdir(parents=True, exist_ok=True)
    struct_rmsd_main(output_dir, args)


def execute_rtbs(target, predictions, pred_dir, workdir, mbr):
    print("\n=====RTBS=====")
    annotator_dir = workdir / "annotated"
    converted_dir = workdir / "converted"
    annotator_target_dir = annotator_dir / "target"
    converted_target_dir = converted_dir / "target"
    rtbs_dir = workdir / "RTBS_results"
    print("\nAnnotating target:")
    target_json = annotate_pdb(target, annotator_target_dir)
    converted_target = convert_annotation(target_json, converted_target_dir)
    prediction_jsons = []

    for pred in predictions:
        print(f"\nAnnotating prediction {pred.name}:")
        pred_json = annotate_pdb(pred, annotator_dir)
        converted = convert_annotation(pred_json, converted_dir)
        prediction_jsons.append(converted)
    
    print("\nCalculating RTBS:")
    run_rtbs(converted_target, prediction_jsons, pred_dir, rtbs_dir, mbr)
    print("RTBS calculated")


def execute_struct_rmsd(args, workdir):
    print("\n=====structRMSD=====")
    rtbs_dir = workdir / "structRMSD_results"
    run_structRMSD(args, rtbs_dir)



def main():
    args = parse_args()
    workdir = Path(args.out_dir)
    workdir.mkdir(parents=True, exist_ok=True)
    target = check_pdb(args.target)
    mbr_path = Path(args.mbr)
    if not mbr_path.is_file():
        raise FileNotFoundError(f"MBR matrix file not found: {mbr_path}")
    pred_dir = None 

    if args.pred_dir:
        pred_dir = Path(args.pred_dir)
        if not pred_dir.exists():
            raise FileNotFoundError(pred_dir)

        predictions = sorted(pred_dir.glob("*.pdb"))
        if not predictions:
            raise RuntimeError("No pdb files found in prediction directory")

    elif args.prediction:
        predictions = [check_pdb(args.prediction)]

    else:
        raise RuntimeError("Provide prediction.pdb or -p prediction_folder")

    if args.tool is None:
        execute_rtbs(target, predictions, pred_dir, workdir, args.mbr)
        execute_struct_rmsd(args, workdir)

    elif args.tool == "RTBS":
        execute_rtbs(target, predictions, pred_dir, workdir, args.mbr)

    elif args.tool == "structRMSD":
        execute_struct_rmsd(args, workdir)


if __name__ == "__main__":
    main()
