import argparse
import json
from pathlib import Path

parser = argparse.ArgumentParser(
    description="Convert an extended DBN structure into an RNA structural graph annotated with BEAR secondary-structure-element types."
)
parser.add_argument("input", help="Input DBN file")
parser.add_argument(
    "-o",
    "--output",
    default="structure_tree.json",
    help="Output JSON file"
)

args = parser.parse_args()

INPUT = Path(args.input)
OUTPUT = Path(args.output)

BRACKET_PAIRS = [
    ("(", ")"),
    ("[", "]"),
    ("{", "}"),
    ("<", ">"),
]

BRACKET_PAIRS += [
    (chr(ord("A") + i), chr(ord("a") + i))
    for i in range(26)
]

OPEN_TO_CLOSE = dict(BRACKET_PAIRS)
CLOSE_TO_OPEN = {
    close: opening
    for opening, close in BRACKET_PAIRS
}


def read_dbn(path):
    if not path.is_file():
        raise ValueError(f"Input file does not exist: {path}")

    if path.suffix.lower() != ".dbn":
        raise ValueError(f"Input file must have a .dbn extension: {path}")

    try:
        with path.open() as f:
            lines = [line.strip() for line in f if line.strip()]
    except OSError as e:
        raise ValueError(f"Cannot read input file '{path}': {e}")

    if len(lines) != 3:
        raise ValueError(
            f"Invalid DBN file '{path}': expected exactly 3 non-empty lines "
            f"(header, sequence, structure), found {len(lines)}"
        )

    if not lines[0].startswith(">"):
        raise ValueError(
            f"Invalid DBN header in '{path}': first line must start with '>'"
        )

    sequence = lines[1].upper()
    structure = lines[2]

    invalid_bases = set(sequence) - set("ACGUT")
    if invalid_bases:
        raise ValueError(
            f"Invalid nucleotide(s) in '{path}': "
            f"{', '.join(sorted(invalid_bases))}"
        )

    if len(sequence) != len(structure):
        raise ValueError(
            f"Sequence and structure have different lengths in '{path}': "
            f"{len(sequence)} != {len(structure)}"
        )

    return sequence, structure


def get_pairs(structure):
    stacks = {
        opening: []
        for opening, _ in BRACKET_PAIRS
    }

    pairs = []

    for position, char in enumerate(structure, start=1):
        if char == ".":
            continue

        if char in OPEN_TO_CLOSE:
            stacks[char].append(position)
            continue

        if char in CLOSE_TO_OPEN:
            opening = CLOSE_TO_OPEN[char]

            if not stacks[opening]:
                raise ValueError(
                    f"Invalid dot-bracket structure: unmatched "
                    f"'{char}' at position {position}"
                )

            left = stacks[opening].pop()

            pairs.append(
                (left, position, opening)
            )
            continue

        raise ValueError(
            f"Invalid dot-bracket structure: unsupported character "
            f"'{char}' at position {position}"
        )

    unmatched = []

    for opening, stack in stacks.items():
        unmatched.extend(
            (opening, position)
            for position in stack
        )

    if unmatched:
        details = ", ".join(
            f"'{opening}' at {position}"
            for opening, position in sorted(
                unmatched,
                key=lambda x: x[1]
            )
        )

        raise ValueError(
            f"Invalid dot-bracket structure: unmatched opening brackets: "
            f"{details}"
        )

    return sorted(
        pairs,
        key=lambda pair: pair[0]
    )


def find_stems(pairs):
    pair_set = set(pairs)
    used = set()
    stems = []

    for left, right, bracket in pairs:
        pair = (left, right, bracket)

        if pair in used:
            continue

        current = [pair]
        used.add(pair)

        left_pos = left
        right_pos = right

        while (
            left_pos + 1,
            right_pos - 1,
            bracket
        ) in pair_set:
            left_pos += 1
            right_pos -= 1

            next_pair = (
                left_pos,
                right_pos,
                bracket
            )

            current.append(next_pair)
            used.add(next_pair)

        stems.append({
            "left_first": current[0][0],
            "left_last": current[-1][0],
            "right_first": current[-1][1],
            "right_last": current[0][1],
            "bracket": bracket,
        })

    stems.sort(
        key=lambda stem: (
            stem["left_first"],
            stem["right_last"]
        )
    )

    return stems


def stems_cross(stem1, stem2):
    a = stem1["left_first"]
    b = stem1["right_last"]
    c = stem2["left_first"]
    d = stem2["right_last"]

    return (
        a < c < b < d
        or
        c < a < d < b
    )


def stem_contains(parent, child):
    return (
        parent["left_first"] < child["left_first"]
        and
        child["right_last"] < parent["right_last"]
    )


def assign_loop_children(stems):
    for i, child in enumerate(stems):
        candidates = []

        child_start = child["left_first"]
        child_end = child["right_last"]

        for j, parent in enumerate(stems):
            if i == j:
                continue

            parent_start = parent["left_first"]
            parent_end = parent["right_last"]

            nested = (
                parent_start < child_start
                and
                child_end < parent_end
            )

            crossing = (
                parent_start < child_start < parent_end < child_end
                or
                child_start < parent_start < child_end < parent_end
            )

            if nested or crossing:
                candidates.append(
                    (
                        abs(parent_start - child_start),
                        parent_end - parent_start,
                        j
                    )
                )

        if candidates:
            _, _, parent_index = min(candidates)
            child["loop_parent"] = parent_index
        else:
            child["loop_parent"] = None

    return stems


def make_strand(first, last, sequence, structure):
    if first > last:
        return None

    if first < 1 or last > len(sequence):
        raise ValueError(
            f"Invalid strand range: {first}-{last}"
        )

    return {
        "first": first,
        "last": last,
        "sequence": sequence[first - 1:last],
        "structure": structure[first - 1:last],
    }


def add_element(
    elements,
    name,
    bear_type,
    strands,
    branching=False,
    extra=None
):
    element_id = len(elements) + 1

    element = {
        "id": element_id,
        "name": name,
        "bear_type": bear_type,
        "branching": branching,
        "parent": None,
        "children": [],
        "strands": strands,
    }

    if extra:
        element.update(extra)

    elements.append(element)

    return element_id


def link(elements, parent_id, child_id):
    if child_id not in elements[parent_id - 1]["children"]:
        elements[parent_id - 1]["children"].append(child_id)

    if elements[child_id - 1]["parent"] is None:
        elements[child_id - 1]["parent"] = parent_id


def classify_loop_for_stem(
    stem_index,
    stems,
    sequence,
    structure
):
    stem = stems[stem_index]

    loop_start = stem["left_last"] + 1
    loop_end = stem["right_first"] - 1

    children = [
        i
        for i, child in enumerate(stems)
        if child["loop_parent"] == stem_index
    ]

    children.sort(
        key=lambda i: stems[i]["left_first"]
    )

    if not children:
        strand = make_strand(
            loop_start,
            loop_end,
            sequence,
            structure
        )

        return {
            "type": "Hairpin",
            "branching": False,
            "strands": [strand] if strand else [],
            "children": [],
        }

    crossing_children = [
        i
        for i in children
        if stems_cross(
            stem,
            stems[i]
        )
    ]

    if len(children) == 1 and not crossing_children:
        child = stems[children[0]]

        left_start = loop_start
        left_end = child["left_first"] - 1

        right_start = child["right_last"] + 1
        right_end = loop_end

        left_length = max(
            0,
            left_end - left_start + 1
        )

        right_length = max(
            0,
            right_end - right_start + 1
        )

        if left_length == 0 and right_length == 0:
            return {
                "type": "Stack",
                "branching": False,
                "children": children,
                "segments": [],
            }

        if left_length > 0 and right_length > 0:
            return {
                "type": "InternalLoop",
                "branching": False,
                "children": children,
                "segments": [
                    (
                        "LEFT",
                        left_start,
                        left_end
                    ),
                    (
                        "RIGHT",
                        right_start,
                        right_end
                    ),
                ],
            }

        if left_length > 0:
            return {
                "type": "BulgeLeft",
                "branching": False,
                "children": children,
                "segments": [
                    (
                        "LEFT",
                        left_start,
                        left_end
                    )
                ],
            }

        return {
            "type": "BulgeRight",
            "branching": False,
            "children": children,
            "segments": [
                (
                    "RIGHT",
                    right_start,
                    right_end
                )
            ],
        }

    segments = []
    current = loop_start

    for child_index in children:
        child = stems[child_index]

        child_left = child["left_first"]
        child_right = child["right_last"]

        if child_left > loop_end:
            continue

        segment_end = min(
            child_left - 1,
            loop_end
        )

        if current <= segment_end:
            segments.append(
                (current, segment_end)
            )

        current = max(
            current,
            child_right + 1
        )

        if current > loop_end:
            break

    if current <= loop_end:
        segments.append(
            (current, loop_end)
        )

    return {
        "type": "Junction",
        "branching": True,
        "children": children,
        "segments": segments,
    }


def chain_top_level(elements):
    top_level = [
        element
        for element in elements
        if element["parent"] is None
    ]

    top_level.sort(
        key=lambda element: (
            min(
                strand["first"]
                for strand in element["strands"]
            )
            if element["strands"]
            else 0
        )
    )

    for previous, current in zip(
        top_level,
        top_level[1:]
    ):
        link(
            elements,
            previous["id"],
            current["id"]
        )

    return elements


def build_structure_graph(sequence, structure):
    pairs = get_pairs(structure)

    if not pairs:
        if not sequence:
            return []

        return [{
            "id": 1,
            "name": "SingleStrand 1",
            "bear_type": "Unclassified",
            "branching": False,
            "parent": None,
            "children": [],
            "strands": [
                make_strand(
                    1,
                    len(sequence),
                    sequence,
                    structure
                )
            ],
        }]

    stems = find_stems(pairs)
    stems = assign_loop_children(stems)

    elements = []
    stem_element_id = {}
    stem_branching = [False] * len(stems)

    for index, stem in enumerate(stems):
        stem_element_id[index] = add_element(
            elements,
            f"Stem {index + 1}",
            "Stem",
            [
                make_strand(
                    stem["left_first"],
                    stem["left_last"],
                    sequence,
                    structure
                ),
                make_strand(
                    stem["right_first"],
                    stem["right_last"],
                    sequence,
                    structure
                ),
            ]
        )

    counters = {
        "Hairpin": 0,
        "BulgeLeft": 0,
        "BulgeRight": 0,
        "LeftInternalLoop": 0,
        "RightInternalLoop": 0,
        "Junction": 0,
    }

    for index in range(len(stems)):
        motif = classify_loop_for_stem(
            index,
            stems,
            sequence,
            structure
        )

        motif_type = motif["type"]

        if motif_type == "Stack":
            child_index = motif["children"][0]

            link(
                elements,
                stem_element_id[index],
                stem_element_id[child_index]
            )

            continue

        if motif_type == "Junction":
            stem_branching[index] = True

            for child_index in motif["children"]:
                stem_branching[child_index] = True

            counters["Junction"] += 1

            junction_id = add_element(
                elements,
                f"Junction {counters['Junction']}",
                "Junction",
                [
                    make_strand(
                        segment[0],
                        segment[1],
                        sequence,
                        structure
                    )
                    for segment in motif["segments"]
                    if segment[0] <= segment[1]
                ],
                branching=True
            )

            link(
                elements,
                stem_element_id[index],
                junction_id
            )

            for child_index in motif["children"]:
                link(
                    elements,
                    junction_id,
                    stem_element_id[child_index]
                )

            continue

        if motif_type == "Hairpin":
            counters["Hairpin"] += 1

            loop_id = add_element(
                elements,
                f"Hairpin {counters['Hairpin']}",
                "LOOP",
                motif["strands"]
            )

            link(
                elements,
                stem_element_id[index],
                loop_id
            )

            continue

        if motif_type == "BulgeLeft":
            counters["BulgeLeft"] += 1

            segment = motif["segments"][0]

            strand = make_strand(
                segment[1],
                segment[2],
                sequence,
                structure
            )

            bulge_id = add_element(
                elements,
                f"BulgeLeft {counters['BulgeLeft']}",
                "BULGELEFT",
                [strand]
            )

            link(
                elements,
                stem_element_id[index],
                bulge_id
            )

            child_index = motif["children"][0]

            link(
                elements,
                bulge_id,
                stem_element_id[child_index]
            )

            continue

        if motif_type == "BulgeRight":
            counters["BulgeRight"] += 1

            segment = motif["segments"][0]

            strand = make_strand(
                segment[1],
                segment[2],
                sequence,
                structure
            )

            bulge_id = add_element(
                elements,
                f"BulgeRight {counters['BulgeRight']}",
                "BULGERIGHT",
                [strand]
            )

            link(
                elements,
                stem_element_id[index],
                bulge_id
            )

            child_index = motif["children"][0]

            link(
                elements,
                bulge_id,
                stem_element_id[child_index]
            )

            continue

        if motif_type == "InternalLoop":
            child_index = motif["children"][0]

            left_segment = motif["segments"][0]
            right_segment = motif["segments"][1]

            left_strand = make_strand(
                left_segment[1],
                left_segment[2],
                sequence,
                structure
            )

            right_strand = make_strand(
                right_segment[1],
                right_segment[2],
                sequence,
                structure
            )

            counters["LeftInternalLoop"] += 1

            left_id = add_element(
                elements,
                f"LeftInternalLoop {counters['LeftInternalLoop']}",
                "LEFTINTERNALLOOP",
                [left_strand]
            )

            counters["RightInternalLoop"] += 1

            right_id = add_element(
                elements,
                f"RightInternalLoop {counters['RightInternalLoop']}",
                "RIGHTINTERNALLOOP",
                [right_strand]
            )

            link(
                elements,
                stem_element_id[index],
                left_id
            )

            link(
                elements,
                left_id,
                right_id
            )

            link(
                elements,
                right_id,
                stem_element_id[child_index]
            )

            continue

    for index, is_branching in enumerate(stem_branching):
        if is_branching:
            stem_id = stem_element_id[index]

            elements[stem_id - 1]["bear_type"] = "STEM_branch"
            elements[stem_id - 1]["branching"] = True

    paired_positions = {
        position
        for left, right, _ in pairs
        for position in (left, right)
    }

    regions = []
    start = None

    for position in range(
        1,
        len(structure) + 1
    ):
        if position not in paired_positions:
            if start is None:
                start = position

        elif start is not None:
            regions.append(
                (start, position - 1)
            )
            start = None

    if start is not None:
        regions.append(
            (start, len(structure))
        )

    single_id = 1

    for first, last in regions:
        inside_structure = any(
            stem["left_first"] < first
            and last < stem["right_last"]
            for stem in stems
        )

        if inside_structure:
            continue

        add_element(
            elements,
            f"SingleStrand {single_id}",
            "Unclassified",
            [
                make_strand(
                    first,
                    last,
                    sequence,
                    structure
                )
            ]
        )

        single_id += 1

    return chain_top_level(elements)


def main():
    try:
        sequence, structure = read_dbn(INPUT)

        graph = build_structure_graph(
            sequence,
            structure
        )

        if OUTPUT.exists() and OUTPUT.is_dir():
            raise ValueError(
                f"Output path is a directory: {OUTPUT}"
            )

        OUTPUT.parent.mkdir(
            parents=True,
            exist_ok=True
        )

        with OUTPUT.open("w") as f:
            json.dump(
                graph,
                f,
                indent=2
            )

    except (
        ValueError,
        OSError,
        json.JSONDecodeError
    ) as e:
        parser.error(str(e))

    print(f"Saved {OUTPUT}")


if __name__ == "__main__":
    main()