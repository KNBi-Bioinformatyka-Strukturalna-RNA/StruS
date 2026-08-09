import argparse
import csv
import json
import os
import shutil
import statistics
import subprocess
import tempfile
import urllib.request
from copy import deepcopy
from pathlib import Path
from tempfile import NamedTemporaryFile

from Bio.PDB import PDBIO, PDBParser, Select, Superimposer
from rnapolis.common import BpSeq, DotBracket
from config import ANNOTATOR, PYTHON_BIN

MIN_RESIDUES_FOR_TM_SCORE = 3

CANONICAL_BASE_PAIRS = {
    frozenset({"A", "U"}),
    frozenset({"G", "C"}),
    frozenset({"G", "U"}),
}


def validate_motif_source(args) -> None:
    motif_sources_given = [
        getattr(args, "motif_tree", None) is not None,
        getattr(args, "dbn", None) is not None,
        getattr(args, "bpseq", None) is not None,
        getattr(args, "annotator", False),
        getattr(args, "fr3d", False),
    ]
    n = sum(motif_sources_given)
    if n == 0:
        args.fr3d = True
        print("No motif source given - defaulting to --fr3d.")
        return
    if n > 1:
        raise SystemExit(
            "error: more than one motif source given - pick exactly one "
            "of: --motif-tree, --dbn, --bpseq, --annotator, --fr3d"
        )


def parse_args():
    parser = argparse.ArgumentParser(description="RMSD of RNA structural motifs (target vs predictions), with global TM-score filtering.")
    parser.add_argument("--target", required=True, type=Path)
    parser.add_argument("--prediction", required=True, type=Path)
    parser.add_argument("--target-chain", type=str, default=None)
    parser.add_argument("--motif-tree", type=Path, default=None)
    parser.add_argument("--dbn", type=Path, default=None)
    parser.add_argument("--bpseq", type=Path, default=None)
    parser.add_argument("--annotator", action="store_true")
    parser.add_argument("--fr3d", action="store_true")
    parser.add_argument("--tm-threshold", type=float, default=0.45)
    parser.add_argument("--usalign-bin", type=str, default=None)
    parser.add_argument("--out-per-motif", type=Path, default=Path("per_motif_rmsd.csv"))
    parser.add_argument("--out-summary", type=Path, default=Path("motif_summary.csv"))
    args = parser.parse_args()

    try:
        validate_motif_source(args)
    except SystemExit as e:
        parser.error(str(e))

    return args


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


def run_command(cmd, quiet=True):
    print("Running:")
    print(" ".join(map(str, cmd)))
    if quiet:
        result = subprocess.run(cmd, stdout=subprocess.DEVNULL)
    else:
        result = subprocess.run(cmd)

    if result.returncode != 0:
        raise RuntimeError(f"Command failed with exit code {result.returncode}")


def run_annotator(target_pdb: Path, output_json: Path) -> Path:
    cmd = [
        PYTHON_BIN,
        ANNOTATOR,
        str(target_pdb),
        "--json",
        str(output_json),
        "--extended",
    ]
    run_command(cmd, quiet=False)
    return output_json


def merge_multichain_dbn(dbn_path: Path) -> tuple[Path, dict[str, tuple[int, int]]]:
    chain_ranges = {}
    sequences = []
    structures = []
    current_start = 1

    with open(dbn_path) as f:
        lines = [line.strip() for line in f if line.strip()]

    i = 0
    while i < len(lines):
        if not lines[i].startswith(">"):
            raise ValueError(f"Expected chain header, got: {lines[i]}")

        chain_name = lines[i].replace(">strand_", "").replace(">", "")

        if i + 2 >= len(lines):
            raise ValueError(f"Incomplete DBN record for {chain_name}")

        sequence = lines[i + 1]
        structure = lines[i + 2]

        sequences.append(sequence)
        structures.append(structure)

        end = current_start + len(sequence) - 1
        chain_ranges[chain_name] = (current_start, end)
        current_start = end + 1

        i += 3

    if len(sequences) == 1:
        return dbn_path, chain_ranges

    merged = NamedTemporaryFile(suffix=".dbn", delete=False, mode="w")
    merged.write(">RNA\n")
    merged.write("".join(sequences) + "\n")
    merged.write("".join(structures) + "\n")
    merged.close()

    return Path(merged.name), chain_ranges


def restrict_to_target_chain(
    elements: list[dict],
    chain_ranges: dict[str, tuple[int, int]],
    chain_id: str,
) -> list[dict]:
    if chain_id not in chain_ranges:
        raise ValueError(
            f"Chain {chain_id!r} not found. Available: {list(chain_ranges)}"
        )

    first, last = chain_ranges[chain_id]
    filtered = []

    for motif in elements:
        strand_ranges = [(s["first"], s["last"]) for s in motif["strands"]]
        inside = [
            first <= s_first <= last and first <= s_last <= last
            for s_first, s_last in strand_ranges
        ]

        if all(inside):
            motif = deepcopy(motif)
            for strand in motif["strands"]:
                strand["first"] -= first - 1
                strand["last"] -= first - 1
            filtered.append(motif)
        elif any(inside):
            print(
                f"[skipped] {motif['name']}: spans multiple chains "
                f"(cannot compare only chain {chain_id})"
            )

    return filtered


def run_fr3d(target_pdb: Path) -> Path:
    target_dir = target_pdb.parent if target_pdb.parent != Path("") else Path(".")
    stem = target_pdb.stem

    cmd = [
        PYTHON_BIN, "-m", "fr3d.classifiers.NA_pairwise_interactions",
        "-i", str(target_dir), "-o", str(target_dir), "-c", "basepair", stem,
    ]
    run_command(cmd, quiet=False)

    basepairs_path = target_dir / f"{stem}_basepair.txt"
    if not basepairs_path.exists():
        raise RuntimeError(
            f"FR3D did not produce the expected output file: {basepairs_path}"
        )
    return basepairs_path


def parse_fr3d_basepairs(basepairs_path: Path) -> list[tuple[dict, dict, str]]:
    seen = set()
    pairs = []

    with open(basepairs_path) as f:
        for line in f:
            line = line.strip()
            if not line:
                continue

            unit_id1, interaction_type, unit_id2, _crossing = line.split("\t")
            parts1 = unit_id1.split("|")
            parts2 = unit_id2.split("|")

            nt1 = {"chain": parts1[2], "base": parts1[3], "number": int(parts1[4])}
            nt2 = {"chain": parts2[2], "base": parts2[3], "number": int(parts2[4])}

            key = tuple(sorted([
                (nt1["chain"], nt1["number"]),
                (nt2["chain"], nt2["number"]),
            ]))
            if key in seen:
                continue
            seen.add(key)

            if (nt1["chain"], nt1["number"]) > (nt2["chain"], nt2["number"]):
                nt1, nt2 = nt2, nt1

            pairs.append((nt1, nt2, interaction_type))

    return pairs


def filter_canonical_pairs(pairs: list[tuple[dict, dict, str]]) -> list[tuple[dict, dict]]:
    canonical = []
    for nt1, nt2, interaction_type in pairs:
        if interaction_type != "cWW":
            continue
        if frozenset({nt1["base"], nt2["base"]}) not in CANONICAL_BASE_PAIRS:
            continue
        canonical.append((nt1, nt2))
    return canonical


def build_bpseq_from_fr3d(
    canonical_pairs: list[tuple[dict, dict]],
    index_to_residue: dict[int, object],
    out_path: Path,
) -> Path:
    chain_local_to_global = {
        (residue.parent.id, residue.get_id()[1]): idx
        for idx, residue in index_to_residue.items()
    }

    n = len(index_to_residue)
    partner = [0] * (n + 1)

    for nt1, nt2 in canonical_pairs:
        key1 = (nt1["chain"], nt1["number"])
        key2 = (nt2["chain"], nt2["number"])

        if key1 not in chain_local_to_global or key2 not in chain_local_to_global:
            print(f"[fr3d] skipped pair {key1}-{key2}: not found in target structure")
            continue

        g1 = chain_local_to_global[key1]
        g2 = chain_local_to_global[key2]

        if partner[g1] != 0 or partner[g2] != 0:
            print(f"[fr3d] skipped pair {key1}-{key2}: residue already paired (conflict)")
            continue

        partner[g1] = g2
        partner[g2] = g1

    with open(out_path, "w") as f:
        for idx in range(1, n + 1):
            base = index_to_residue[idx].get_resname().strip()
            f.write(f"{idx} {base} {partner[idx]}\n")

    return out_path


def get_motif_tree(
    target_pdb: Path,
    motif_tree_path: Path | None,
    dbn_path: Path | None = None,
    bpseq_path: Path | None = None,
    use_annotator: bool = False,
    use_fr3d: bool = False,
    target_chain: str | None = None,
) -> list[dict]:
    if motif_tree_path is not None:
        with open(motif_tree_path) as f:
            return json.load(f)

    if dbn_path is not None:
        dbn_path, chain_ranges = merge_multichain_dbn(dbn_path)
        print("Detected DBN chain ranges:")
        print(chain_ranges)

        bpseq = BpSeq.from_dotbracket(DotBracket.from_file(str(dbn_path)))
        elements = load_elements_from_bpseq(bpseq)

        if target_chain is not None:
            elements = restrict_to_target_chain(elements, chain_ranges, target_chain)
    elif bpseq_path is not None:
        if target_chain is not None:
            print(f"--target-chain {target_chain!r} has no effect with --bpseq")
        bpseq = BpSeq.from_file(str(bpseq_path))
        elements = load_elements_from_bpseq(bpseq)
    elif use_annotator:
        raw_json_path = Path(target_pdb).with_suffix(".annotator.json")
        run_annotator(target_pdb, raw_json_path)

        with open(raw_json_path) as f:
            raw_data = json.load(f)

        elements = load_elements(raw_data)
    elif use_fr3d:
        basepairs_path = run_fr3d(target_pdb)
        raw_pairs = parse_fr3d_basepairs(basepairs_path)
        canonical_pairs = filter_canonical_pairs(raw_pairs)

        _model, index_to_residue = load_structure(target_pdb)
        fr3d_bpseq_path = Path(target_pdb).with_suffix(".fr3d.bpseq")
        build_bpseq_from_fr3d(canonical_pairs, index_to_residue, fr3d_bpseq_path)
        print(f"Built BPSEQ from FR3D output: {fr3d_bpseq_path}")

        bpseq = BpSeq.from_file(str(fr3d_bpseq_path))
        elements = load_elements_from_bpseq(bpseq)
    else:
        raise ValueError("get_motif_tree: no motif source given")

    tree_path = Path(target_pdb).with_suffix(".structure_tree.json")
    with open(tree_path, "w") as f:
        json.dump(elements, f, indent=2)
    print(f"Saved motif list: {tree_path}")

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


def get_chain_ids(pdb_path: Path) -> list[str]:
    parser = PDBParser(QUIET=True)
    structure = parser.get_structure(pdb_path.stem, str(pdb_path))
    model = next(iter(structure))
    return [chain.id for chain in model]


def resolve_target_pdb(
    target_pdb: Path,
    target_chain: str | None,
    prediction_paths: list[Path],
    use_bpseq: bool = False,
    use_motif_tree: bool = False,
) -> tuple[Path, str | None]:
    if use_bpseq or use_motif_tree:
        return target_pdb, target_chain

    target_chain_ids = get_chain_ids(target_pdb)

    if target_chain is not None:
        if target_chain not in target_chain_ids:
            raise ValueError(
                f"--target-chain {target_chain!r} not found in {target_pdb} "
                f"- available chains: {target_chain_ids}"
            )
        if len(target_chain_ids) == 1:
            return target_pdb, target_chain
    else:
        if len(target_chain_ids) == 1:
            return target_pdb, target_chain

        prediction_chain_counts = {
            len(get_chain_ids(prediction_path))
            for prediction_path in prediction_paths
        }

        if prediction_chain_counts == {len(target_chain_ids)}:
            return target_pdb, target_chain

        if prediction_chain_counts == {1}:
            target_chain = target_chain_ids[0]
            print(
                f"No --target-chain specified. Target has "
                f"{len(target_chain_ids)} chains and predictions are "
                f"single-chain. Using first target chain "
                f"{target_chain!r} by default."
            )
        else:
            print(
                f"Target has {len(target_chain_ids)} chains, but "
                f"prediction chain counts are not uniform "
                f"({sorted(prediction_chain_counts)}) - falling back to "
                f"the full target."
            )
            return target_pdb, target_chain

    class _ChainSelector(Select):
        def accept_chain(self, chain):
            return chain.id == target_chain

    parser = PDBParser(QUIET=True)
    structure = parser.get_structure(target_pdb.stem, str(target_pdb))
    model = next(iter(structure))

    out_path = target_pdb.with_suffix(f".chain{target_chain}.pdb")
    io = PDBIO()
    io.set_structure(model)
    io.save(str(out_path), _ChainSelector())

    print(f"Using target chain {target_chain!r}: {out_path}")
    return out_path, target_chain


def load_structure(pdb_path: Path):
    parser = PDBParser(QUIET=True)
    structure = parser.get_structure(pdb_path.stem, str(pdb_path))
    model = next(iter(structure))

    index_to_residue = {}
    global_index = 1
    for chain in model:
        for residue in chain:
            index_to_residue[global_index] = residue
            global_index += 1

    return model, index_to_residue


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

    print("USalign not found - downloading source and compiling (g++)...")
    if not Path("USalign.cpp").exists():
        urllib.request.urlretrieve(
            "https://zhanggroup.org/US-align/bin/module/USalign.cpp", "USalign.cpp"
        )
    subprocess.run(
        ["g++", "-O3", "-ffast-math", "-o", "USalign", "USalign.cpp"],
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
            f"Unexpected USalign output for {target_pdb} vs {prediction_pdb}: "
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
                f"[skipped] {prediction_pdb.name}: TM-score={tm_score:.3f} "
                f"< {threshold}"
            )
    return passing


def collect_motif_residue_ids(motif: dict) -> list[int]:
    residue_ids = set()
    for strand in motif["strands"]:
        residue_ids.update(range(strand["first"], strand["last"] + 1))
    return sorted(residue_ids)


def get_atoms_by_residue(
    index_to_residue: dict[int, object], residue_ids: list[int]
) -> dict[tuple[int, str], object]:
    atoms = {}
    for res_id in residue_ids:
        residue = index_to_residue.get(res_id)
        if residue is None:
            continue
        for atom in residue:
            atoms[(res_id, atom.get_name())] = atom
    return atoms


class _MotifResidueSelector(Select):
    def __init__(self, residue_keys: set[tuple]):
        self.residue_keys = residue_keys

    def accept_residue(self, residue):
        return (residue.parent.id, residue.get_id()) in self.residue_keys


def write_motif_fragment(
    model, index_to_residue: dict[int, object], residue_ids: list[int], out_path: Path
) -> Path:
    residue_keys = {
        (index_to_residue[res_id].parent.id, index_to_residue[res_id].get_id())
        for res_id in residue_ids
        if res_id in index_to_residue
    }
    io = PDBIO()
    io.set_structure(model)
    io.save(str(out_path), _MotifResidueSelector(residue_keys))
    return out_path


def compute_motif_tm_score(
    target_model, target_index, prediction_model, prediction_index,
    motif: dict, usalign_bin: str
) -> float | None:
    residue_ids = collect_motif_residue_ids(motif)

    if len(residue_ids) < MIN_RESIDUES_FOR_TM_SCORE:
        print(
            f"[motif skipped - TM-score] {motif['name']}: too short for "
            f"TM-score ({len(residue_ids)} residues, USalign needs at "
            f"least {MIN_RESIDUES_FOR_TM_SCORE})"
        )
        return None

    with tempfile.TemporaryDirectory() as tmp_dir:
        tmp_dir = Path(tmp_dir)
        target_fragment = tmp_dir / "target_fragment.pdb"
        prediction_fragment = tmp_dir / "prediction_fragment.pdb"

        write_motif_fragment(target_model, target_index, residue_ids, target_fragment)
        write_motif_fragment(prediction_model, prediction_index, residue_ids, prediction_fragment)

        try:
            return compute_tm_score(target_fragment, prediction_fragment, usalign_bin)
        except (subprocess.CalledProcessError, RuntimeError) as e:
            print(f"[motif skipped - TM-score] {motif['name']}: USalign failed ({e})")
            return None


def compute_motif_rmsd(
    target_index, prediction_index, motif: dict, min_coverage: float = 0.9
) -> float | None:
    residue_ids = collect_motif_residue_ids(motif)
    target_atoms = get_atoms_by_residue(target_index, residue_ids)
    prediction_atoms = get_atoms_by_residue(prediction_index, residue_ids)

    if not target_atoms:
        print(
            f"[motif skipped] {motif['name']}: none of its residues exist "
            f"in the target representation being used for this comparison"
        )
        return None

    common_keys = sorted(set(target_atoms) & set(prediction_atoms))
    coverage = len(common_keys) / len(target_atoms)
    if coverage < min_coverage:
        print(
            f"[motif skipped] {motif['name']}: atom coverage "
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
    usalign_bin: str,
) -> list[dict]:

    target_model, target_index = load_structure(target_pdb)
    records = []

    for prediction_path, _global_tm_score in passing_predictions:
        prediction_model, prediction_index = load_structure(prediction_path)

        if len(prediction_index) != len(target_index):
            print(
                f"[warning] {prediction_path.name}: target contains "
                f"{len(target_index)} residues while prediction contains "
                f"{len(prediction_index)}. Motifs extending beyond the "
                f"prediction residue range will be skipped."
            )
            prediction_last = len(prediction_index)
        else:
            prediction_last = None

        for motif in motif_tree:
            if prediction_last is not None:
                motif_last = max(strand["last"] for strand in motif["strands"])
                if motif_last > prediction_last:
                    print(
                        f"[skipped] {motif['name']}: residues extend beyond "
                        f"prediction ({motif_last} > {prediction_last})"
                    )
                    continue

            motif_rmsd = compute_motif_rmsd(target_index, prediction_index, motif)
            motif_tm_score = compute_motif_tm_score(
                target_model, target_index, prediction_model, prediction_index,
                motif, usalign_bin
            )

            records.append({
                "motif_id": motif["id"],
                "motif_type": motif_type_from_name(motif),
                "residue_range": format_residue_range(motif),
                "prediction_file": prediction_path.name,
                "motif_tm_score": motif_tm_score,
                "motif_rmsd": motif_rmsd,
            })

    return records


def aggregate_stats(records: list[dict]) -> list[dict]:
    by_motif: dict[int, list[dict]] = {}
    for record in records:
        by_motif.setdefault(record["motif_id"], []).append(record)

    summary = []
    for motif_id in sorted(by_motif):
        group = by_motif[motif_id]
        rmsds = [r["motif_rmsd"] for r in group if r["motif_rmsd"] is not None]

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
        "prediction_file", "motif_tm_score", "motif_rmsd",
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
    args = parse_args()

    prediction = args.prediction if args.prediction else args.pred_dir
    prediction_paths = list_prediction_files(prediction)

    target_pdb, target_chain = resolve_target_pdb(
        args.target,
        getattr(args, "target_chain", None),
        prediction_paths,
        use_bpseq=args.bpseq is not None,
        use_motif_tree=args.motif_tree is not None,
    )

    motif_tree = get_motif_tree(
        target_pdb, args.motif_tree, args.dbn, args.bpseq,
        args.annotator, args.fr3d, target_chain,
    )

    usalign_bin = find_usalign_binary(args.usalign_bin)

    passing_predictions = filter_predictions_by_tm_score(
        target_pdb, prediction_paths, args.tm_threshold, usalign_bin
    )
    print(
        f"{len(passing_predictions)}/{len(prediction_paths)} predictions "
        f"passed the TM-score threshold >= {args.tm_threshold}"
    )

    records = process_target(target_pdb, motif_tree, passing_predictions, usalign_bin)
    summary = aggregate_stats(records)

    write_per_motif_csv(records, args.out_per_motif)
    write_summary_csv(summary, args.out_summary)
    print(f"Saved {args.out_per_motif} and {args.out_summary}")


def struct_rmsd_main(workdir, kwargs):
    args = kwargs
    validate_motif_source(args)

    prediction = Path(args.prediction) if args.prediction else Path(args.pred_dir)
    prediction_paths = list_prediction_files(prediction)

    target_pdb, target_chain = resolve_target_pdb(
        Path(args.target),
        getattr(args, "target_chain", None),
        prediction_paths,
        use_bpseq=args.bpseq is not None,
        use_motif_tree=args.motif_tree is not None,
    )

    motif_tree = get_motif_tree(
        target_pdb, args.motif_tree, args.dbn, args.bpseq,
        args.annotator, args.fr3d, target_chain,
    )

    usalign_bin = find_usalign_binary(args.usalign_bin)

    passing_predictions = filter_predictions_by_tm_score(
        target_pdb, prediction_paths, args.tm_threshold, usalign_bin
    )
    print(
        f"{len(passing_predictions)}/{len(prediction_paths)} predictions "
        f"passed the TM-score threshold >= {args.tm_threshold}"
    )

    records = process_target(target_pdb, motif_tree, passing_predictions, usalign_bin)
    summary = aggregate_stats(records)

    out_per_motif_path = os.path.join(workdir, args.out_per_motif)
    out_summary_path = os.path.join(workdir, args.out_summary)

    write_per_motif_csv(records, out_per_motif_path)
    write_summary_csv(summary, out_summary_path)
    print(f"Saved {out_per_motif_path} and {out_summary_path}")


if __name__ == "__main__":
    main()