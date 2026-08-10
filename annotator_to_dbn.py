#!/usr/bin/env python3

import argparse
import json
from pathlib import Path


def convert_annotator_to_dotbracket(input_path: Path, output_path: Path):
    with open(input_path, "r") as f:
        data = json.load(f)
    db = data["dot_bracket"]
    lines = db.splitlines()

    if len(lines) < 3:
        raise ValueError(
            "Unexpected dot_bracket format. "
            "Expected at least 3 lines: >strand, sequence, structure."
        )

    sequence = lines[1].strip()
    structure = lines[2].strip()

    if len(sequence) != len(structure):
        raise ValueError(
            f"Sequence and structure have different lengths: "
            f"{len(sequence)} vs {len(structure)}"
        )
    structure_id = input_path.stem

    with open(output_path, "w") as f:
        f.write(f">{structure_id}\n")
        f.write(f"{sequence}\n")
        f.write(f"{structure}\n")

    print(f"ID:        {structure_id}")
    print(f"Length:    {len(sequence)}")
    print(f"Sequence:  {sequence}")
    print(f"Structure: {structure}")
    print(f"Saved:     {output_path}")


def main():
    parser = argparse.ArgumentParser(
        description=(
            "Convert RNAPolis annotator JSON to full dot-bracket format."
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

    convert_annotator_to_dotbracket(args.input, args.output)


if __name__ == "__main__":
    main()
