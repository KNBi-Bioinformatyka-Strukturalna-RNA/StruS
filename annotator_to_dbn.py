#!/usr/bin/env python3

import argparse
import json
from pathlib import Path
from dataclasses import dataclass


@dataclass
class Chain:
    chain_id: str
    sequence: str
    structure: str


def parse_dot_bracket(db: str) -> list[Chain]:
    lines = [line for line in db.splitlines() if line.strip() != ""]

    if len(lines) < 3:
        raise ValueError(
            "Unexpected dot_bracket format. "
            "Expected at least 3 lines: >strand, sequence, structure."
        )

    chains: list[Chain] = []
    i = 0
    while i < len(lines):
        header = lines[i].strip()
        if not header.startswith(">"):
            raise ValueError(
                f"Expected a header line starting with '>' at line {i + 1}, "
                f"got: {header!r}"
            )
        if i + 2 >= len(lines):
            raise ValueError(
                f"Truncated block for header {header!r}: "
                "missing sequence and/or structure line."
            )

        chain_id = header[1:].strip()
        sequence = lines[i + 1].strip()
        structure = lines[i + 2].strip()

        if len(sequence) != len(structure):
            raise ValueError(
                f"Sequence and structure have different lengths for chain "
                f"{chain_id!r}: {len(sequence)} vs {len(structure)}"
            )

        chains.append(Chain(chain_id, sequence, structure))
        i += 3

    return chains


def chains_from_bpseq_index(sequence: str, structure: str, bpseq_index: dict) -> list[Chain] | None:
    n = len(sequence)
    if len(bpseq_index) != n:
        return None

    chains: list[Chain] = []
    current_chain_id = None
    current_seq: list[str] = []
    current_struct: list[str] = []

    for position in range(1, n + 1):
        info = bpseq_index.get(str(position))
        if info is None or "auth" not in info or "chain" not in info["auth"]:
            return None
        chain_id = info["auth"]["chain"]

        if chain_id != current_chain_id:
            if current_chain_id is not None:
                chains.append(Chain(current_chain_id, "".join(current_seq), "".join(current_struct)))
            current_chain_id = chain_id
            current_seq = []
            current_struct = []

        current_seq.append(sequence[position - 1])
        current_struct.append(structure[position - 1])

    if current_chain_id is not None:
        chains.append(Chain(current_chain_id, "".join(current_seq), "".join(current_struct)))

    seen_counts: dict[str, int] = {}
    for chain in chains:
        seen_counts[chain.chain_id] = seen_counts.get(chain.chain_id, 0) + 1
    running_counts: dict[str, int] = {}
    for chain in chains:
        if seen_counts[chain.chain_id] > 1:
            running_counts[chain.chain_id] = running_counts.get(chain.chain_id, 0) + 1
            chain.chain_id = f"{chain.chain_id}_seg{running_counts[chain.chain_id]}"

    return chains


def diagnose_empty_dot_bracket(data: dict, input_path: Path) -> str:
    num_base_pairs = len(data.get("base_pairs", []) or [])
    num_stackings = len(data.get("stackings", []) or [])
    bpseq_sequence = (data.get("bpseq", {}) or {}).get("sequence", "")

    lines = [
        f"RNAPolis annotator returned an empty 'dot_bracket' for '{input_path}'.",
        f"  base_pairs found:    {num_base_pairs}",
        f"  stackings found:     {num_stackings}",
        f"  bpseq.sequence:      {bpseq_sequence!r}",
    ]

    if num_stackings > 0 and not bpseq_sequence:
        lines.append(
            "This combination (pairwise interactions detected, but no "
            "sequence/bpseq at all) usually means RNAPolis could not "
            "reconstruct an ordered nucleotide chain from the backbone "
            "geometry (e.g. missing or degenerate P/O3'/O5' atoms, broken "
            "or overlapping residues) for this predicted structure. "
            "Try to use FR3D instead of annotator. "
            "This can also be a property of the predicted 3D structure itself "
            "- check the geometry/backbone connectivity of the prediction."
        )
    elif num_base_pairs == 0:
        lines.append(
            "No canonical base pairs were detected in this structure, so "
            "there is no secondary structure to convert."
        )

    return "\n".join(lines)


def convert_annotator_to_dotbracket(input_path: Path, output_path: Path):
    with open(input_path, "r") as f:
        data = json.load(f)

    db = data.get("dot_bracket", "") or ""
    if not db.strip():
        raise ValueError(diagnose_empty_dot_bracket(data, input_path))

    rnapolis_chains = parse_dot_bracket(db)

    sequence = "".join(chain.sequence for chain in rnapolis_chains)
    structure = "".join(chain.structure for chain in rnapolis_chains)

    bpseq_index = data.get("bpseq_index") or {}
    chains = chains_from_bpseq_index(sequence, structure, bpseq_index)
    used_auth_chains = chains is not None
    if chains is None:
        chains = rnapolis_chains

    structure_id = input_path.stem

    with open(output_path, "w") as f:
        for chain in chains:
            f.write(f">{chain.chain_id}\n")
            f.write(f"{chain.sequence}\n")
            f.write(f"{chain.structure}\n")

    total_length = sum(len(c.sequence) for c in chains)

    print(f"ID:          {structure_id}")
    print(f"Chain ids:   {'from bpseq_index (original PDB chain letters)' if used_auth_chains else 'from RNAPolis dot_bracket blocks (fallback - bpseq_index unavailable/incomplete)'}")
    print(f"Chains:      {len(chains)}")
    for chain in chains:
        print(f"  - {chain.chain_id}: length={len(chain.sequence)}")
        print(f"      Sequence:  {chain.sequence}")
        print(f"      Structure: {chain.structure}")
    print(f"Total length: {total_length}")
    print(f"Saved:        {output_path}")


def main():
    parser = argparse.ArgumentParser(
        description=(
            "Convert RNAPolis annotator JSON to full dot-bracket format "
            "(supports single- and multi-chain RNA structures)."
        )
    )

    parser.add_argument(
        "input",
        type=Path,
        help="RNAPolis annotator JSON file"
    )

    parser.add_argument(
        "-o",
        "--output",
        type=Path,
        required=True,
        help="Output dot-bracket file"
    )

    args = parser.parse_args()

    try:
        convert_annotator_to_dotbracket(args.input, args.output)
    except (ValueError, KeyError, OSError, json.JSONDecodeError) as e:
        parser.error(str(e))


if __name__ == "__main__":
    main()