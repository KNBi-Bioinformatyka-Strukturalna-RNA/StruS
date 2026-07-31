import json
from pathlib import Path
import argparse

#This is a simple script that converts the output of the annotator tool into a tree structure.
parser = argparse.ArgumentParser(description="Json annotation to tree structure converter")
parser.add_argument("input", help="Input JSON file from annotator with --json")
parser.add_argument(
    "-o", "--output",
    default="structure_tree.json",
    help="Output JSON file"
)

args = parser.parse_args()

INPUT = args.input
OUTPUT = args.output

def load_elements(data):
    elements = []
    eid = 1

    def add_element(name, strands):
        nonlocal eid
        firsts = [s["first"] for s in strands]
        lasts  = [s["last"] for s in strands]

        elements.append({
            "id": eid,
            "name": name,
            "_first": firsts,
            "_last": lasts,
            "parent": None,
            "children": [],
            "strands": strands
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
                    "structure": s["strand5p"]["structure"]
                },
                {
                    "first": s["strand3p"]["first"],
                    "last": s["strand3p"]["last"],
                    "sequence": s["strand3p"]["sequence"],
                    "structure": s["strand3p"]["structure"]
                }
            ]
        )

    for i, s in enumerate(data.get("single_strands", []), 1):
        add_element(
            f"SingleStrand {i}",
            [{
                "first": s["strand"]["first"],
                "last": s["strand"]["last"],
                "sequence": s["strand"]["sequence"],
                "structure": s["strand"]["structure"]
            }]
        )

    for i, h in enumerate(data.get("hairpins", []), 1):
        add_element(
            f"Hairpin {i}",
            [{
                "first": h["strand"]["first"],
                "last": h["strand"]["last"],
                "sequence": h["strand"]["sequence"],
                "structure": h["strand"]["structure"]
            }]
        )
    
    for i, h in enumerate(data.get("loops", []), 1):
        add_element(
            f"Loop {i}",
            [{
                "first": s["first"],
                "last": s["last"],
                "sequence": s["sequence"],
                "structure": s["structure"]
                }
                for s in h["strands"]
            ]
    )

    return elements

def build_tree(elements):
    for child in elements:
        candidates = []
        for parent in elements:
            if parent["id"] == child["id"]:
                continue
            if any(
                c_first == p_last
                for c_first in child["_first"]
                for p_last in parent["_last"]
            ):
                candidates.append(parent)

        if candidates:
            parent = min(candidates, key=lambda p: min(p["_first"]))
            child["parent"] = parent["id"]
            parent["children"].append(child["id"])
    
    # Usuwamy pomocnicze pola
    for e in elements:
        del e["_first"]
        del e["_last"]

    return elements

if __name__ == "__main__":
    with open(INPUT) as f:
        data = json.load(f)

    elements = load_elements(data)
    tree = build_tree(elements)

    with open(OUTPUT, "w") as f:
        json.dump(tree, f, indent=2)

    print(f"Zapisano {OUTPUT}")
