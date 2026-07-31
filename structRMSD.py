import argparse
import csv
import json
import os
import os
import shutil
import statistics
import subprocess
import urllib.request
from pathlib import Path

from Bio.PDB import PDBParser, Superimposer
from rnapolis.common import BpSeq, DotBracket


def parse_args():
    parser = argparse.ArgumentParser(
        description="RMSD motywow strukturalnych RNA (target vs predykcje), "
                     "z filtrowaniem po globalnym TM-score."
    )
    parser.add_argument("--target", required=True, type=Path)
    parser.add_argument("--prediction", required=True, type=Path)
    parser.add_argument("--motif-tree", type=Path, default=None)
    parser.add_argument("--dbn", type=Path, default=None)
    parser.add_argument("--bpseq", type=Path, default=None)
    parser.add_argument("--tm-threshold", type=float, default=0.45)
    parser.add_argument("--usalign-bin", type=str, default=None)
    parser.add_argument("--out-per-motif", type=Path, default=Path("per_motif_rmsd.csv"))
    parser.add_argument("--out-summary", type=Path, default=Path("motif_summary.csv"))
    return parser.parse_args()


def load_elements(data: dict) -> list[dict]:
    elements = []
    eid = 1

    def add_element(name, strands):
        nonlocal eid
        elements.append({
            "id": eid,
            "name": name,
            "strands": strands,
        })
        eid += 1

    for i, s in enumerate(data.get("stems", []), 1):
        add_element(
            f"Stem {i}",
            [
                {
                    "first": s["strand5p"]["first"],
                    "last": s["strand5p"]["last"],
                    "sequence": s["strand5p"]["sequence"],
                    "structure": s["strand5p"]["structure"],
                },
                {
                    "first": s["strand3p"]["first"],
                    "last": s["strand3p"]["last"],
                    "sequence": s["strand3p"]["sequence"],
                    "structure": s["strand3p"]["structure"],
                },
            ],
        )

    for i, s in enumerate(data.get("single_strands", []), 1):
        add_element(
            f"SingleStrand {i}",
            [{
                "first": s["strand"]["first"],
                "last": s["strand"]["last"],
                "sequence": s["strand"]["sequence"],
                "structure": s["strand"]["structure"],
            }],
        )

    for i, h in enumerate(data.get("hairpins", []), 1):
        add_element(
            f"Hairpin {i}",
            [{
                "first": h["strand"]["first"],
                "last": h["strand"]["last"],
                "sequence": h["strand"]["sequence"],
                "structure": h["strand"]["structure"],
            }],
        )

    for i, h in enumerate(data.get("loops", []), 1):
        add_element(
            f"Loop {i}",
            [
                {
                    "first": s["first"],
                    "last": s["last"],
                    "sequence": s["sequence"],
                    "structure": s["structure"],
                }
                for s in h["strands"]
            ],
        )

    return elements


def load_elements_from_bpseq(bpseq: "BpSeq") -> list[dict]:
    elements = []
    eid = 1

    def add_element(name, strands):
        nonlocal eid
        elements.append({"id": eid, "name": name, "strands": strands})
        eid += 1

    def strand_dict(strand):
        return {
            "first": strand.first,
            "last": strand.last,
            "sequence": strand.sequence,
            "structure": strand.structure,
        }

    stems, single_strands, hairpins, loops = bpseq.elements

    for i, s in enumerate(stems, 1):
        add_element(f"Stem {i}", [strand_dict(s.strand5p), strand_dict(s.strand3p)])

    for i, s in enumerate(single_strands, 1):
        add_element(f"SingleStrand {i}", [strand_dict(s.strand)])

    for i, h in enumerate(hairpins, 1):
        add_element(f"Hairpin {i}", [strand_dict(h.strand)])

    for i, l in enumerate(loops, 1):
        add_element(f"Loop {i}", [strand_dict(s) for s in l.strands])

    return elements


def run_annotator(target_pdb: Path, output_json: Path) -> Path:
    try:
        subprocess.run(
            ["annotator", "--json", str(output_json), "--extended", str(target_pdb)],
            check=True,
            capture_output=True,
            text=True,
        )
    except subprocess.CalledProcessError as e:
        print(f"annotator zakonczyl sie bledem:\n{e.stderr}")
        raise
    return output_json


def get_motif_tree(
    target_pdb: Path,
    motif_tree_path: Path | None,
    dbn_path: Path | None = None,
    bpseq_path: Path | None = None,
) -> list[dict]:
    if motif_tree_path is not None:
        with open(motif_tree_path) as f:
            return json.load(f)

    if dbn_path is not None:
        bpseq = BpSeq.from_dotbracket(DotBracket.from_file(str(dbn_path)))
    elif bpseq_path is not None:
        bpseq = BpSeq.from_file(str(bpseq_path))
    else:
        bpseq = None

    if bpseq is not None:
        elements = load_elements_from_bpseq(bpseq)
    else:
        raw_json_path = Path(target_pdb).with_suffix(".annotator.json")
        run_annotator(target_pdb, raw_json_path)

        with open(raw_json_path) as f:
            raw_data = json.load(f)

        elements = load_elements(raw_data)

    tree_path = Path(target_pdb).with_suffix(".structure_tree.json")
    with open(tree_path, "w") as f:
        json.dump(elements, f, indent=2)
    print(f"Zapisano liste motywow: {tree_path}")

    return elements


def motif_type_from_name(motif: dict) -> str:
    prefix = motif["name"].split(" ")[0]

    if prefix == "Stem":
        return "stem"
    if prefix == "Hairpin":
        return "hairpin"
    if prefix == "SingleStrand":
        return "single_strand"

    if prefix == "Loop":
        n_strands = len(motif["strands"])
        if n_strands >= 3:
            return f"junction_{n_strands}way"
        if n_strands == 2:
            lengths = [s["last"] - s["first"] + 1 for s in motif["strands"]]
            if min(lengths) <= 0:
                return "bulge"
            return "internal_loop"
        return "loop_unexpected_strand_count"

    return prefix.lower()


def load_structure(pdb_path: Path, chain_id: str = "A"):
    parser = PDBParser(QUIET=True)
    structure = parser.get_structure(pdb_path.stem, str(pdb_path))
    model = next(iter(structure))
    return model[chain_id]


def list_prediction_files(predictions_arg: Path) -> list[Path]:
    if predictions_arg.is_dir():
        return sorted(predictions_arg.glob("*.pdb"))
    return [predictions_arg]


def find_usalign_binary(usalign_bin: str | None) -> str:
    if usalign_bin:
        return usalign_bin

    found = shutil.which("USalign")
    if found:
        return found

    if Path("USalign").exists():
        return "./USalign"

    print("USalign nie znaleziony - pobieram zrodlo i kompiluje (g++)...")
    if not Path("USalign.cpp").exists():
        urllib.request.urlretrieve(
            "https://zhanggroup.org/US-align/bin/module/USalign.cpp", "USalign.cpp"
        )
    subprocess.run(
        ["g++", "-static", "-O3", "-ffast-math", "-o", "USalign", "USalign.cpp"],
        check=True,
    )
    return "./USalign"


def compute_tm_score(target_pdb: Path, prediction_pdb: Path, usalign_bin: str) -> float:
    result = subprocess.run(
        [usalign_bin, str(target_pdb), str(prediction_pdb), "-outfmt", "2"],
        capture_output=True,
        text=True,
        check=True,
    )

    lines = result.stdout.strip().split("\n")
    if len(lines) < 2:
        raise RuntimeError(
            f"Nieoczekiwany output USalign dla {target_pdb} vs {prediction_pdb}: "
            f"{result.stdout!r}"
        )

    tm_score = float(lines[1].split("\t")[2])
    return tm_score


def filter_predictions_by_tm_score(
    target_pdb: Path,
    prediction_paths: list[Path],
    threshold: float,
    usalign_bin: str,
) -> list[tuple[Path, float]]:
    passing = []
    for prediction_pdb in prediction_paths:
        tm_score = compute_tm_score(target_pdb, prediction_pdb, usalign_bin)
        if tm_score >= threshold:
            passing.append((prediction_pdb, tm_score))
        else:
            print(
                f"[pominieto] {prediction_pdb.name}: TM-score={tm_score:.3f} "
                f"< {threshold}"
            )
    return passing


def collect_motif_residue_ids(motif: dict) -> list[int]:
    residue_ids = set()
    for strand in motif["strands"]:
        residue_ids.update(range(strand["first"], strand["last"] + 1))
    return sorted(residue_ids)


def get_atoms_by_residue(chain, residue_ids: list[int]) -> dict[tuple[int, str], object]:
    atoms = {}
    for res_id in residue_ids:
        if res_id not in chain:
            continue
        residue = chain[res_id]
        for atom in residue:
            atoms[(res_id, atom.get_name())] = atom
    return atoms


def compute_motif_rmsd(
    target_chain, prediction_chain, motif: dict, min_coverage: float = 0.9
) -> float | None:
    residue_ids = collect_motif_residue_ids(motif)
    target_atoms = get_atoms_by_residue(target_chain, residue_ids)
    prediction_atoms = get_atoms_by_residue(prediction_chain, residue_ids)

    if not target_atoms:
        return None

    common_keys = sorted(set(target_atoms) & set(prediction_atoms))
    coverage = len(common_keys) / len(target_atoms)
    if coverage < min_coverage:
        print(
            f"[pominieto motyw] {motif['name']}: pokrycie atomow "
            f"{coverage:.0%} < {min_coverage:.0%}"
        )
        return None

    fixed = [target_atoms[key] for key in common_keys]
    moving = [prediction_atoms[key] for key in common_keys]

    superimposer = Superimposer()
    superimposer.set_atoms(fixed, moving)
    return superimposer.rms


def format_residue_range(motif: dict) -> str:
    return ",".join(f"{s['first']}-{s['last']}" for s in motif["strands"])


def process_target(
    target_pdb: Path,
    motif_tree: list[dict],
    passing_predictions: list[tuple[Path, float]],
) -> list[dict]:

    target_chain = load_structure(target_pdb)
    records = []

    for prediction_path, tm_score in passing_predictions:
        prediction_chain = load_structure(prediction_path)

        for motif in motif_tree:
            rmsd = compute_motif_rmsd(target_chain, prediction_chain, motif)
            records.append({
                "motif_id": motif["id"],
                "motif_type": motif_type_from_name(motif),
                "residue_range": format_residue_range(motif),
                "prediction_file": prediction_path.name,
                "tm_score": tm_score,
                "rmsd": rmsd,
            })

    return records


def aggregate_stats(records: list[dict]) -> list[dict]:
    """
    Grupuje po motif_id, liczy: n_predictions (tylko te z rmsd != None),
    mean_rmsd, std_rmsd. Zwraca liste rekordow zbiorczych, po jednym per
    motyw (posortowane po motif_id).
    """
    by_motif: dict[int, list[dict]] = {}
    for record in records:
        by_motif.setdefault(record["motif_id"], []).append(record)

    summary = []
    for motif_id in sorted(by_motif):
        group = by_motif[motif_id]
        rmsds = [r["rmsd"] for r in group if r["rmsd"] is not None]

        if rmsds:
            mean_rmsd = statistics.mean(rmsds)
            std_rmsd = statistics.pstdev(rmsds) if len(rmsds) > 1 else 0.0
        else:
            mean_rmsd = None
            std_rmsd = None

        summary.append({
            "motif_id": motif_id,
            "motif_type": group[0]["motif_type"],
            "residue_range": group[0]["residue_range"],
            "n_predictions": len(rmsds),
            "mean_rmsd": mean_rmsd,
            "std_rmsd": std_rmsd,
        })

    return summary


def write_per_motif_csv(records: list[dict], out_path: Path):
    fieldnames = [
        "motif_id", "motif_type", "residue_range",
        "prediction_file", "tm_score", "rmsd",
    ]
    with open(out_path, "w", newline="") as f:
        writer = csv.DictWriter(f, fieldnames=fieldnames)
        writer.writeheader()
        writer.writerows(records)


def write_summary_csv(summary: list[dict], out_path: Path):
    fieldnames = [
        "motif_id", "motif_type", "residue_range",
        "n_predictions", "mean_rmsd", "std_rmsd",
    ]
    with open(out_path, "w", newline="") as f:
        writer = csv.DictWriter(f, fieldnames=fieldnames)
        writer.writeheader()
        writer.writerows(summary)


def main():
    # args = kwargs
    args = parse_args()
    prediction = args.prediction if args.prediction else args.pred_dir
    motif_tree = get_motif_tree(args.target, args.motif_tree, args.dbn, args.bpseq)
    usalign_bin = find_usalign_binary(args.usalign_bin)

    prediction_paths = list_prediction_files(prediction)
    passing_predictions = filter_predictions_by_tm_score(
        args.target, prediction_paths, args.tm_threshold, usalign_bin
    )
    print(
        f"{len(passing_predictions)}/{len(prediction_paths)} predykcji "
        f"przeszlo prog TM-score >= {args.tm_threshold}"
    )

    records = process_target(args.target, motif_tree, passing_predictions)
    summary = aggregate_stats(records)

    write_per_motif_csv(records, args.out_per_motif)
    write_summary_csv(summary, args.out_summary)
    print(f"Zapisano {args.out_per_motif} i {args.out_summary}")


def struct_rmsd_main(workdir, kwargs):
    args = kwargs
    prediction = Path(args.prediction) if args.prediction else Path(args.pred_dir)
    motif_tree = get_motif_tree(args.target, args.motif_tree, args.dbn, args.bpseq)
    usalign_bin = find_usalign_binary(args.usalign_bin)

    prediction_paths = list_prediction_files(prediction)
    passing_predictions = filter_predictions_by_tm_score(
        args.target, prediction_paths, args.tm_threshold, usalign_bin
    )
    print(
        f"{len(passing_predictions)}/{len(prediction_paths)} predykcji "
        f"przeszlo prog TM-score >= {args.tm_threshold}"
    )

    records = process_target(Path(args.target), motif_tree, passing_predictions)
    summary = aggregate_stats(records)

    out_per_motif_path = os.path.join(workdir, args.out_per_motif)
    out_summary_path = os.path.join(workdir, args.out_summary)
    write_per_motif_csv(records, out_per_motif_path)
    write_summary_csv(summary, out_summary_path)
    print(f"Zapisano {out_per_motif_path} i {out_summary_path}")


if __name__ == "__main__":
    main()
