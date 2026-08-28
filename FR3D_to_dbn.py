import argparse
from dataclasses import dataclass
from pathlib import Path
from Bio.PDB import PDBParser


# This script converts FR3D cWW base-pair annotations into an extended
# dot-bracket structure. Different bracket types are used for pseudoknots.
# Multi-chain structures are written as one '>chain_id' block per PDB
# chain, so downstream tools (e.g. dbn_to_structure_graph.py) can pick up
# the original chain letters and support inter-chain base pairs.

parser = argparse.ArgumentParser(
    description="Convert FR3D cWW base-pair annotations to extended, multi-chain DBN format"
)
parser.add_argument("fr3d", help="Input FR3D annotation file")
parser.add_argument("pdb", help="Input PDB file")
parser.add_argument("-o", "--output", default="structure.dbn", help="Output DBN file")
args = parser.parse_args()

FR3D = Path(args.fr3d)
PDB = Path(args.pdb)
OUTPUT = Path(args.output)

BRACKET_TYPES = [("(", ")"), ("[", "]"), ("{", "}"), ("<", ">")]
BRACKET_TYPES += [(chr(ord("A") + i), chr(ord("a") + i)) for i in range(26)]

MODIFIED_NUCLEOTIDE_PARENTS = {
    # Adenosine derivatives
    "1MA": "A", "6IA": "A", "6MZ": "A", "A2M": "A", "MA6": "A", "MIA": "A", "I": "A", "MAD": "A", "M2A": "A",
    # Cytidine derivatives
    "OMC": "C", "5MC": "C", "M5C": "C", "4OC": "C", "AGM": "C", "CCC": "C",
    # Guanosine derivatives
    "2MG": "G", "7MG": "G", "M2G": "G", "OMG": "G", "YG": "G", "YYG": "G", "1MG": "G", "QUO": "G", "G7M": "G", "G46": "G",
    # Uridine derivatives
    "H2U": "U", "PSU": "U", "OMU": "U", "5MU": "U", "4SU": "U", "T": "U", "UR3": "U", "5BU": "U", "70U": "U", "DHU": "U", "PGP": "U",
}


@dataclass
class Chain:
    chain_id: str
    sequence: str
    structure: str


def resname_to_base(resname):
    if resname in {"A", "C", "G", "U"}:
        return resname
    return MODIFIED_NUCLEOTIDE_PARENTS.get(resname)


def read_pdb_sequence(path):
    if not path.is_file():
        raise ValueError(f"PDB file does not exist: {path}")
    try:
        parser = PDBParser(QUIET=True)
        structure = parser.get_structure("structure", path)
    except Exception as e:
        raise ValueError(f"Cannot parse PDB file '{path}': {e}")

    residues = []
    unmapped_resnames = set()
    for model in structure:
        for chain in model:
            for residue in chain:
                hetflag, resseq, icode = residue.id
                resname = residue.resname.strip().upper()
                base = resname_to_base(resname)
                if base is not None:
                    residues.append((chain.id, resseq, icode, base))
                    continue
                if hetflag.strip():
                    unmapped_resnames.add(resname)
        break

    if unmapped_resnames:
        print(
            "Warning: encountered HETATM residues not in the standard/modified "
            f"ribonucleotide tables, skipped: {sorted(unmapped_resnames)}. "
            "If FR3D reports a base pair involving one of these, add it to "
            "MODIFIED_NUCLEOTIDE_PARENTS in this script."
        )
    if not residues:
        raise ValueError(f"No RNA residues (A, C, G, U, or recognized modified ribonucleotides) found in PDB file '{path}'")
    return residues


def read_fr3d(path, residues):
    if not path.is_file():
        raise ValueError(f"FR3D file does not exist: {path}")

    residue_positions = {}
    for index, (chain, resseq, icode, resname) in enumerate(residues, start=1):
        key = (chain, resseq)
        if key in residue_positions:
            raise ValueError(f"Ambiguous PDB residue numbering: chain {chain}, residue {resseq} occurs more than once")
        residue_positions[key] = index

    pairs = set()
    try:
        with path.open() as f:
            for line_number, line in enumerate(f, start=1):
                line = line.strip()
                if not line:
                    continue
                fields = line.split()
                if len(fields) < 3:
                    raise ValueError(f"Invalid FR3D line {line_number}: expected at least 3 fields")
                residue1 = fields[0]
                interaction = fields[1]
                residue2 = fields[2]
                if interaction != "cWW":
                    continue
                try:
                    parts1 = residue1.split("|")
                    parts2 = residue2.split("|")
                    if len(parts1) < 5 or len(parts2) < 5:
                        raise ValueError
                    chain1 = parts1[2]
                    resseq1 = int(parts1[4])
                    chain2 = parts2[2]
                    resseq2 = int(parts2[4])
                except (IndexError, ValueError):
                    raise ValueError(f"Invalid residue identifier in FR3D line {line_number}: {line}")

                key1 = (chain1, resseq1)
                key2 = (chain2, resseq2)
                if key1 not in residue_positions:
                    raise ValueError(f"FR3D residue {residue1} was not found in PDB")
                if key2 not in residue_positions:
                    raise ValueError(f"FR3D residue {residue2} was not found in PDB")

                pos1 = residue_positions[key1]
                pos2 = residue_positions[key2]
                if pos1 == pos2:
                    raise ValueError(f"Self-pair detected at position {pos1}")
                pairs.add(tuple(sorted((pos1, pos2))))
    except OSError as e:
        raise ValueError(f"Cannot read FR3D file '{path}': {e}")

    return sorted(pairs)


def pairs_cross(pair1, pair2):
    left1, right1 = pair1
    left2, right2 = pair2
    return left1 < left2 < right1 < right2 or left2 < left1 < right2 < right1


def assign_bracket_types(pairs):
    assignments = {}
    for pair in sorted(pairs):
        used_types = set()
        for previous_pair, bracket_type in assignments.items():
            if pairs_cross(pair, previous_pair):
                used_types.add(bracket_type)
        for bracket_type in range(len(BRACKET_TYPES)):
            if bracket_type not in used_types:
                assignments[pair] = bracket_type
                break
        else:
            raise ValueError(f"Too many mutually crossing base pairs. Maximum supported bracket types: {len(BRACKET_TYPES)}")
    return assignments


def build_dbn(residues, pairs):
    structure = ["."] * len(residues)
    assignments = assign_bracket_types(pairs)
    for left, right in pairs:
        bracket_type = assignments[(left, right)]
        opening, closing = BRACKET_TYPES[bracket_type]
        if structure[left - 1] != ".":
            raise ValueError(f"Position {left} is involved in multiple base pairs")
        if structure[right - 1] != ".":
            raise ValueError(f"Position {right} is involved in multiple base pairs")
        structure[left - 1] = opening
        structure[right - 1] = closing

    sequence = "".join(residue[3] for residue in residues)
    return sequence, "".join(structure), assignments


def split_by_chain(residues, sequence, structure):
    if len(residues) != len(sequence) or len(residues) != len(structure):
        raise ValueError("Internal error: residues/sequence/structure length mismatch")

    chains = []
    current_chain_id = None
    current_seq = []
    current_struct = []
    for index, residue in enumerate(residues):
        chain_id = residue[0]
        if chain_id != current_chain_id:
            if current_chain_id is not None:
                chains.append(Chain(current_chain_id, "".join(current_seq), "".join(current_struct)))
            current_chain_id = chain_id
            current_seq = []
            current_struct = []
        current_seq.append(sequence[index])
        current_struct.append(structure[index])
    if current_chain_id is not None:
        chains.append(Chain(current_chain_id, "".join(current_seq), "".join(current_struct)))

    seen_counts = {}
    for chain in chains:
        seen_counts[chain.chain_id] = seen_counts.get(chain.chain_id, 0) + 1
    running_counts = {}
    for chain in chains:
        if seen_counts[chain.chain_id] > 1:
            running_counts[chain.chain_id] = running_counts.get(chain.chain_id, 0) + 1
            chain.chain_id = f"{chain.chain_id}_seg{running_counts[chain.chain_id]}"

    return chains


def main():
    try:
        residues = read_pdb_sequence(PDB)
        pairs = read_fr3d(FR3D, residues)
        sequence, structure, assignments = build_dbn(residues, pairs)
        chains = split_by_chain(residues, sequence, structure)

        if OUTPUT.exists() and OUTPUT.is_dir():
            raise ValueError(f"Output path is a directory: {OUTPUT}")

        OUTPUT.parent.mkdir(parents=True, exist_ok=True)
        with OUTPUT.open("w") as f:
            for chain in chains:
                f.write(f">{chain.chain_id}\n")
                f.write(f"{chain.sequence}\n")
                f.write(f"{chain.structure}\n")

    except (ValueError, OSError) as e:
        parser.error(str(e))

    used_brackets = sorted(set(assignments.values()))
    print(f"Saved {OUTPUT}")
    print(f"RNA residues: {len(residues)}")
    print(f"Chains: {len(chains)}")
    for chain in chains:
        print(f"  - {chain.chain_id}: length={len(chain.sequence)}")
    print(f"Retained cWW pairs: {len(pairs)}")
    print(f"Bracket types used: {len(used_brackets)}")


if __name__ == "__main__":
    main()