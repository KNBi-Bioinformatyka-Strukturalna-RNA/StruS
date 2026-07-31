"""
RNA Tree BEAR Similarity with RNA Tree Penalty Visualizer
============================
Compares a target RNA secondary structure tree with a predicted tree,
annotates each matched node of the prediction with a penalty score
derived from the real MBR (Matrix of BEAR-encoded RNA secondary structures)
from Mattei et al. 2014 (doi:10.1093/nar/gku283).
Calculates two versions of the RTBS measure,
which determines the structural similarity of structures based on a tree structure,
and creates a visualization of penalties for structural differences between structures
for each node of the tree of structural elements.

Penalty = -MBR_score(BEAR_char_target, BEAR_char_pred)
  → negative MBR scores become positive penalties (bad)
  → positive MBR scores become negative penalties (reward, node turns green)

Normalisation:
  best_possible  = sum of diagonal MBR penalties over TARGET nodes
                   (score of a perfect prediction reproducing the target exactly)
  
  for RTBS ​​target:        
    worst_possible = n_target * UNMATCHED_PENALTY
                   (every target node missing from prediction)
    Unmatched prediction nodes are NOT penalised explicitly — their cost appears
    implicitly because target nodes they "replaced" are counted as missing.

  for RTBS symmetrical:
    worst_possible = (n_target + n_prediction) * UNMATCHED_PENALTY/2
    The mismatch cost is distributed equally across the nodes from the target and prediction.
    This measure better represents the structural similarity of two structures,
    but it does not allow for direct comparison of different predictions of a single target based on it,
    because it depends on the size of the prediction tree's structural elements.
  
  sum_penalty = the sum of penalties from all matched nodes and the penalty for mismatches
  from the target (tgt) or from the target and prediction (sym)

  RTBS = 1 - (sum_penalty − best_possible) / (worst_possible − best_possible)

Interpretation:
  Values ​​close to 1 indicate very close structural similarity,
  values ​​close to 0 indicate very little or no similarity.

Usage:
  python rna_tree_penalty_viz.py target.json prediction.json [-o out.png] [--show]
  python rna_tree_penalty_viz.py target.json -p predictions/ [-o out_dir/] [--show]

This script requires json files from annotation_converter.py.
"""

import json
import argparse
import re
import csv
import pathlib
import sys
import matplotlib.pyplot as plt
import matplotlib.cm as cm
import matplotlib.colors as mcolors
import numpy as np


#arguments
parser = argparse.ArgumentParser(
    description="Visualise RNA prediction tree(s) with real MBR penalty scores"
)
parser.add_argument("target",
    help="Target structure JSON")
parser.add_argument("prediction", nargs="?", default=None,
    help="Single predicted structure JSON (omit when using -p)")
parser.add_argument("-p", "--pred_dir", default=None,
    help="Folder of predicted structure JSONs (alternative to single prediction)")
parser.add_argument("-m", "--mbr", default=None,
    help="Path to MBR CSV matrix (default: mbr_matrix.json next to this script)")
parser.add_argument("-o", "--output", default=None,
    help="Output path: image file for single prediction, folder for -p mode")
parser.add_argument("--show", action="store_true",
    help="Display each plot interactively")
args = parser.parse_args()

#Validate argument combinations
if args.pred_dir is None and args.prediction is None:
    parser.error("Provide either a single prediction file or -p/--pred_dir folder.")
if args.pred_dir is not None and args.prediction is not None:
    parser.error("Cannot use both a single prediction file and -p/--pred_dir.")

#BEAR alphabet
BEAR_ALPHABET = {
    "stem":               list("abcdefghi="),
    "stem_branch":        list("ABCDEFGHIJ"),
    "hairpin":            list("jklmnopqrstuvwxyz^"),
    "loop":               list("jklmnopqrstuvwxyz^"),
    "internalloop_left":  ['!', '"', '#', '$', '%', '&', "'", '(', ')', '+'],
    "internalloop_right": ['2','3','4','5','6','7','8','9','0','>'],
    "bulge":              ['['],
    "singlestrand":       [':'],
}

def name_to_bear_chars(name: str) -> list:
    n = re.sub(r'\s*\d+$', '', name).strip().lower()
    if n.startswith("stem"):    return BEAR_ALPHABET["stem"]
    if n.startswith("hairpin"): return BEAR_ALPHABET["hairpin"]
    if n.startswith("loop"):    return BEAR_ALPHABET["loop"]
    if n.startswith("single"):  return BEAR_ALPHABET["singlestrand"]
    if n.startswith("bulge"):   return ['[', ']']
    if n.startswith("internal"):return BEAR_ALPHABET["internalloop_left"]
    return [':']

def get_structural_length(node: dict, name: str) -> int:
    strands = node.get("strands", [])
    n = re.sub(r'\s*\d+$', '', name).strip().lower()
    if not strands:
        return 1
    if n.startswith("stem"):
        return max(strands[0].get("last", 0) - strands[0].get("first", 0) + 1, 1)
    total = sum(s.get("last", 0) - s.get("first", 0) + 1 for s in strands)
    return max(total, 1)

def node_to_bear(node: dict) -> str:
    name  = node.get("name", "")
    chars = name_to_bear_chars(name)
    length = get_structural_length(node, name)
    idx = max(0, min(length - 1, len(chars) - 1))
    return chars[idx]

#Load MBR matrix
def load_mbr_from_csv(csv_path: str) -> dict:
    rows = []
    with open(csv_path, encoding='utf-8') as f:
        for row in csv.reader(f):
            rows.append(row)
    header = rows[0]
    chars = [c for c in header if c.strip() != '']
    mbr = {}
    for row in rows[1:]:
        if not row or not row[0].strip():
            break
        row_char = row[0]
        values = []
        for v in row[1:len(chars)+1]:
            v = v.strip().strip('"').replace(',', '.')
            try:    values.append(float(v))
            except: values.append(0.0)
        mbr[row_char] = dict(zip(chars, values))
    return mbr

def load_mbr_from_json(json_path: str) -> dict:
    with open(json_path) as f:
        return json.load(f)

script_dir = pathlib.Path(__file__).parent
mbr_json   = script_dir / "mbr_matrix.json"
mbr_csv    = args.mbr

if mbr_csv and pathlib.Path(mbr_csv).exists():
    MBR = load_mbr_from_csv(mbr_csv)
    print(f"Loaded MBR from CSV: {mbr_csv}")
elif mbr_json.exists():
    MBR = load_mbr_from_json(str(mbr_json))
    print(f"Loaded MBR from JSON cache: {mbr_json}")
else:
    raise FileNotFoundError(
        "MBR matrix not found. Provide --mbr path/to/matrix.csv "
        "or place mbr_matrix.json next to this script."
    )

MAX_MBR_SCORE     = max(MBR[r][c] for r in MBR for c in MBR[r] if MBR[r][c] is not None)
MIN_MBR_SCORE     = min(MBR[r][c] for r in MBR for c in MBR[r] if MBR[r][c] is not None)
UNMATCHED_PENALTY = -MIN_MBR_SCORE   # ≈ 3.68

def mbr_penalty(node_target: dict, node_pred: dict) -> tuple:
    bt = node_to_bear(node_target)
    bp = node_to_bear(node_pred)
    if bt in MBR and bp in MBR[bt]:
        score = MBR[bt][bp]
    elif bp in MBR and bt in MBR[bp]:
        score = MBR[bp][bt]
    else:
        score = MIN_MBR_SCORE
    return -score, bt, bp

# Load / cycle-break helpers
def load_tree(path) -> dict:
    with open(path, encoding="utf-8") as f:
        data = json.load(f)
    return {n["id"]: n for n in data}

def find_root(nodes):
    no_parent = [nid for nid, n in nodes.items() if n.get("parent") is None]
    if len(no_parent) == 1:
        return no_parent[0]
    kids = set()
    for n in nodes.values():
        kids.update(n.get("children", []))
    non_kids = [nid for nid in nodes if nid not in kids]
    if len(non_kids) == 1:
        return non_kids[0]
    def min_first(nid):
        return min((s.get("first", 9999) for s in nodes[nid].get("strands", [])), default=9999)
    return min(nodes.keys(), key=min_first)

def break_cycles(root, nodes):
    clean = {nid: dict(n) for nid, n in nodes.items()}
    for nid in clean:
        clean[nid]["children"] = list(clean[nid].get("children", []))
    visited, queue = set(), [root]
    while queue:
        nid = queue.pop(0)
        if nid not in clean or nid in visited:
            continue
        visited.add(nid)
        valid = []
        for ch in clean[nid]["children"]:
            if ch not in visited and ch in clean:
                valid.append(ch)
                queue.append(ch)
        clean[nid]["children"] = valid
    return {nid: clean[nid] for nid in visited}

#Local subtree matching with additional checking of the remaining nodes
def normalize_name(name: str) -> str:
    return re.sub(r'\s*\d+$', '', name).strip().lower()


#Nodes are considered compatible if the nucleotide range of one is within the range of the other.
def nucleotide_range_match(node_a, node_b):
    strands_a = node_a.get("strands", [])
    strands_b = node_b.get("strands", [])
    if not strands_a or not strands_b:
        return False
    ranges_a = [
        (s.get("first"), s.get("last"))
        for s in strands_a
        if s.get("first") is not None and s.get("last") is not None
    ]
    ranges_b = [
        (s.get("first"), s.get("last"))
        for s in strands_b
        if s.get("first") is not None and s.get("last") is not None
    ]
    for a1, a2 in ranges_a:
        for b1, b2 in ranges_b:
            if ((a1 <= b1 and a2 >= b2) or (b1 <= a1 and b2 >= a2)):
                return True
    return False



def match_subtree(nodes_a, root_a, nodes_b, root_b, memo=None, visiting=None):
    if memo is None:
        memo = {}
    if visiting is None:
        visiting = set()
    state = (root_a, root_b)
    if state in memo:
        return memo[state]
    if state in visiting:
        return 0, []
    na = nodes_a[root_a]
    nb = nodes_b[root_b]
    if normalize_name(na["name"]) != normalize_name(nb["name"]):
        return 0, []
    visiting.add(state)
    mapping = [(root_a, root_b)]
    score = 1
    children_a = na.get("children", [])
    children_b = nb.get("children", [])

    if children_a and children_b:
        candidates = []
        for ca in children_a:
            for cb in children_b:
                sc, mp = match_subtree(nodes_a, ca, nodes_b, cb, memo, visiting)
                if sc > 0:
                    candidates.append((sc, mp, ca, cb))
        candidates.sort(key=lambda x: x[0], reverse=True)
        used_a = set()
        used_b = set()
        for sc, mp, ca, cb in candidates:
            if ca in used_a or cb in used_b:
                continue
            score += sc
            mapping.extend(mp)
            used_a.add(ca)
            used_b.add(cb)
    visiting.remove(state)

    seen_a = set()
    seen_b = set()
    clean = []
    for a, b in mapping:
        if a not in seen_a and b not in seen_b:
            clean.append((a, b))
            seen_a.add(a)
            seen_b.add(b)
    result = (len(clean), clean)
    memo[state] = result

    return result



def find_best_mapping(nodes_t, nodes_p):
    memo = {}
    all_matches = []
    for nid_t in nodes_t:
        for nid_p in nodes_p:
            score, mapping = match_subtree(nodes_t, nid_t, nodes_p, nid_p, memo)
            if score >= 2:
                all_matches.append((score, mapping))

    all_matches.sort(key=lambda x: x[0], reverse=True)
    best_mapping = []
    used_t = set()
    used_p = set()

    for score, mapping in all_matches:
        if any(t in used_t or p in used_p for t, p in mapping):
            continue
        best_mapping.extend(mapping)
        for t, p in mapping:
            used_t.add(t)
            used_p.add(p)

    #target_leaves = [nid for nid, node in nodes_t.items() if not node.get("children")]
    #prediction_leaves = [nid for nid, node in nodes_p.items() if not node.get("children")]
    #for lt in target_leaves:
    #    if lt in used_t:
    #        continue
    #    for lp in prediction_leaves:
    #        if lp in used_p:
    #            continue
    #        if normalize_name(nodes_t[lt]["name"]) != normalize_name(nodes_p[lp]["name"]):
    #            continue
    #        if leaf_nucleotide_match(nodes_t[lt], nodes_p[lp]):
    #            best_mapping.append((lt, lp))
    #            used_t.add(lt)
    #            used_p.add(lp)
    #            break

    remaining_t = [nid for nid in nodes_t if nid not in used_t]
    remaining_p = [nid for nid in nodes_p if nid not in used_p]
    for nt in remaining_t:
        for np in remaining_p:
            if np in used_p:
                continue
            if normalize_name(nodes_t[nt]["name"]) != normalize_name(nodes_p[np]["name"]):
                continue
            if nucleotide_range_match(nodes_t[nt], nodes_p[np]):
                best_mapping.append((nt, np))
                used_t.add(nt)
                used_p.add(np)
                break

    return best_mapping

# Load target (once — shared across all predictions)
nodes_t_raw = load_tree(args.target)
root_t_raw  = find_root(nodes_t_raw)
nodes_t     = break_cycles(root_t_raw, nodes_t_raw)
root_t      = find_root(nodes_t)
n_dropped_t = len(nodes_t_raw) - len(nodes_t)
if n_dropped_t:
    print(f"WARNING: target had cyclic edges; {n_dropped_t} back-edge(s) removed.")

target_name = pathlib.Path(args.target).stem
n_target    = len(nodes_t)

# best_possible: computed once over target nodes, shared by all predictions
best_possible = 0.0
for nid_t in nodes_t:
    bt   = node_to_bear(nodes_t[nid_t])
    diag = MBR.get(bt, {}).get(bt, None)
    best_possible += (-diag if diag is not None else -MAX_MBR_SCORE)

# worst_possible: every target node unmatched (n_target × UNMATCHED_PENALTY)
worst_possible = n_target * UNMATCHED_PENALTY
penalty_range  = worst_possible - best_possible

#Core processing function
H_SPACING = 3.2
V_SPACING = 2.5

def process_prediction(pred_path: pathlib.Path, out_path: pathlib.Path):
    #load & clean prediction tree
    nodes_p_raw = load_tree(pred_path)
    root_p_raw  = find_root(nodes_p_raw)
    nodes_p     = break_cycles(root_p_raw, nodes_p_raw)
    root_p      = find_root(nodes_p)
    n_dropped_p = len(nodes_p_raw) - len(nodes_p)
    if n_dropped_p:
        print(f"  WARNING: prediction had cyclic edges; {n_dropped_p} back-edge(s) removed.")

    pred_name = pred_path.stem
    n_pred    = len(nodes_p)

    #mapping
    best_mapping   = find_best_mapping(nodes_t, nodes_p)
    pred_to_target = {b: a for a, b in best_mapping}
    target_to_pred = {a: b for a, b in best_mapping}
    n_matched           = len(best_mapping)
    n_unmatched_pred    = n_pred   - n_matched
    n_unmatched_target  = n_target - n_matched

    #penalties per prediction node
    node_info = {}
    for nid_p in nodes_p:
        if nid_p in pred_to_target:
            nid_t = pred_to_target[nid_p]
            pen, bt, bp = mbr_penalty(nodes_t[nid_t], nodes_p[nid_p])
            node_info[nid_p] = {"penalty": pen, "bear_t": bt, "bear_p": bp,
                                 "matched": True, "target_id": nid_t}
        else:
            bp = node_to_bear(nodes_p[nid_p])
            node_info[nid_p] = {"penalty": UNMATCHED_PENALTY, "bear_t": "—",
                                 "bear_p": bp, "matched": False}

    unmatched_target_ids = [nid for nid in nodes_t if nid not in target_to_pred]

    # total penalty: MBR costs for matched pairs
    #              + UNMATCHED_PENALTY for every unmatched prediction node
    #              + UNMATCHED_PENALTY for every unmatched target node
    pred_matched_penalty = sum(
        info["penalty"] for info in node_info.values() if info["matched"]
    )
    total_pen = (pred_matched_penalty
                 + n_unmatched_pred   * UNMATCHED_PENALTY
                 + n_unmatched_target * UNMATCHED_PENALTY)
    avg_pen = total_pen / (n_pred + n_target - n_matched)

    # Separate total for target-anchored metric - excludes unmatched prediction
    # nodes so that predictions of different sizes are directly comparable.
    total_pen_target = pred_matched_penalty + n_unmatched_target * UNMATCHED_PENALTY
    avg_pen_target   = total_pen_target / n_target

    # Metric 1: norm_penalty  (symmetric, sensitive to prediction size)
    #   worst = (n_target + n_pred) x UNMATCHED_PENALTY
    #   Extra nodes in prediction raise the ceiling - larger predictions are
    #   penalised more for any unmatched node.
    worst_sym    = (n_target + n_pred) * UNMATCHED_PENALTY
    range_sym    = worst_sym - best_possible
    norm_penalty = (total_pen - best_possible) / range_sym if range_sym > 0 else 0.0
    norm_penalty = 1.0-float(np.clip(norm_penalty, 0.0, 1.0))

    # Metric 2: norm_penalty_target  (target-anchored, independent of prediction size)
    #   Uses total_pen_target (no unmatched pred nodes) and worst = n_target x UNMATCHED.
    #   Both numerator and denominator are fixed to target size, so predictions
    #   with extra or fewer nodes than the target are fully comparable.
    worst_target        = n_target * UNMATCHED_PENALTY
    range_target        = worst_target - best_possible
    norm_penalty_target = (total_pen_target - best_possible) / range_target if range_target > 0 else 0.0
    norm_penalty_target = 1.0-float(np.clip(norm_penalty_target, 0.0, 1.0))

    #layout
    def layout_tree(node_id, nodes, x=0, y=0, pos=None, visited=None):
        if pos     is None: pos     = {}
        if visited is None: visited = set()
        if node_id in visited: return x, x, pos
        visited.add(node_id)
        node     = nodes[node_id]
        children = [c for c in node.get("children", []) if c in nodes]
        if not children:
            pos[node_id] = (x, y); return x, x, pos
        left = right = None
        cur_x = x
        for child in children:
            l, r, _ = layout_tree(child, nodes, cur_x, y - V_SPACING, pos, set(visited))
            cur_x = r + H_SPACING
            left  = l if left  is None else min(left,  l)
            right = r if right is None else max(right, r)
        cx = (left + right) / 2
        pos[node_id] = (cx, y)
        return left, right, pos

    _, _, pos_p = layout_tree(root_p, nodes_p)

    #colour mapping
    all_penalties = [info["penalty"] for info in node_info.values()]
    p_min = min(all_penalties)
    p_max = max(all_penalties)

    GLOBAL_PENALTY_MIN = -MAX_MBR_SCORE
    GLOBAL_PENALTY_MAX = -MIN_MBR_SCORE

    safe_vmin = GLOBAL_PENALTY_MIN
    safe_vmax = GLOBAL_PENALTY_MAX
    vcenter   = 0.0
    #if vcenter <= safe_vmin: vcenter = safe_vmin + abs(safe_vmin) * 0.01 + 1e-6
    #if vcenter >= safe_vmax: vcenter = safe_vmax - abs(safe_vmax) * 0.01 - 1e-6

    cmap = mcolors.LinearSegmentedColormap.from_list(
        "penalty", ["#2ecc71", "#f5f5f5", "#e74c3c"])
    norm = mcolors.TwoSlopeNorm(vmin=safe_vmin, vcenter=vcenter, vmax=safe_vmax)

    def penalty_color(p):
        return cmap(norm(np.clip(p, safe_vmin, safe_vmax)))

    #draw
    xs = [x for x, y in pos_p.values()]
    ys = [y for x, y in pos_p.values()]
    x_min, x_max = min(xs), max(xs)
    y_min, y_max = min(ys), max(ys)
    x_pad        = max(H_SPACING * 1.5, 3.0)
    y_pad_top    = V_SPACING * 0.6
    y_pad_bottom = V_SPACING * 1.2

    width  = max(12, (x_max - x_min + 2 * x_pad) * 0.9)
    height = max(7,  (y_max - y_min + y_pad_top + y_pad_bottom) * 0.85)
    fig, ax = plt.subplots(figsize=(width, height))
    ax.axis("off")
    ax.set_xlim(x_min - x_pad,       x_max + x_pad)
    ax.set_ylim(y_min - y_pad_bottom, y_max + y_pad_top)

    # edges
    for nid, (x, y) in pos_p.items():
        for child_id in nodes_p[nid].get("children", []):
            if child_id in pos_p:
                cx, cy = pos_p[child_id]
                ax.plot([x, cx], [y, cy], color="#95a5a6", lw=1.2, zorder=1, alpha=0.7)

    # nodes
    label_offset = V_SPACING * 0.1
    for nid, (x, y) in pos_p.items():
        info = node_info[nid]
        pen  = info["penalty"]
        col  = penalty_color(pen)
        name = nodes_p[nid]["name"]
        bt   = info["bear_t"]
        bp   = info["bear_p"]

        ax.text(x, y, name, ha="center", va="center", fontsize=8,
                fontweight="bold", zorder=5,
                bbox=dict(facecolor=col, edgecolor="#2c3e50",
                          boxstyle="round,pad=0.4", alpha=0.93, linewidth=1.1))

        if info["matched"]:
            tname = nodes_t[info["target_id"]]["name"]
            lbl = f"pen={pen:+.2f}  BEAR: {bt}→{bp}  ({tname})"
        else:
            lbl = f"pen={UNMATCHED_PENALTY:+.2f}  UNMATCHED  BEAR: {bp}"

        ax.text(x, y - label_offset, lbl,
                ha="center", va="top", fontsize=8, color="#34495e",
                zorder=6, style="italic")

    # colorbar
    fig.subplots_adjust(left=0.01, right=0.87, top=0.93, bottom=0.04)
    cbar_ax = fig.add_axes([0.89, 0.20, 0.015, 0.55])
    sm = cm.ScalarMappable(cmap=cmap, norm=norm)
    sm.set_array([])
    cbar = fig.colorbar(sm, cax=cbar_ax)
    cbar.set_label("Penalty  (−MBR score)", fontsize=8)
    tick_vals = sorted({safe_vmin, safe_vmax, 0.0})
    cbar.set_ticks(tick_vals)
    cbar.set_ticklabels([f"{v:+.2f}" for v in tick_vals])
    #tick_vals = sorted({
    #    round(p_min, 2),
    #    round(p_max, 2),
    #    0.0,
    #    round(UNMATCHED_PENALTY, 2),
    #})
    #tick_vals = [t for t in tick_vals if safe_vmin <= t <= safe_vmax]
    #cbar.set_ticks(tick_vals)
    #cbar.set_ticklabels([f"{v:+.2f}" for v in tick_vals], fontsize=7)
    cbar_ax.text(0.5,  1.01, "worst\nmatch", ha="center", va="bottom",
                 fontsize=9, color="#e74c3c", fontweight="bold",
                 transform=cbar_ax.transAxes)
    cbar_ax.text(0.5, -0.01, "best\nmatch", ha="center", va="top",
                 fontsize=9, color="#2ecc71", fontweight="bold",
                 transform=cbar_ax.transAxes)

    # info box
    info_txt = (
        f"Target:     {target_name}  ({n_target} nodes)\n"
        f"Prediction: {pred_name}  ({n_pred} nodes)\n"
        f"  Unmatched pred:   {n_unmatched_pred}\n"
        f"  Unmatched target: {n_unmatched_target}\n"
        f"Sum penalty (sym):  {total_pen:+.2f}  (avg/node: {avg_pen:+.2f})\n"
        f"Sum penalty (tgt):  {total_pen_target:+.2f}  (avg/node: {avg_pen_target:+.2f})\n"
        f"RTBS (sym): {norm_penalty:.4f}  (1=best, 0=worst)\n"
        f"RTBS (tgt): {norm_penalty_target:.4f}  (1=best, 0=worst)\n"
        f"  best={best_possible:.2f}  worst_sym={worst_sym:.2f}  worst_tgt={worst_target:.2f}"
    )
    ax.text(0.01, 0.99, info_txt,
            transform=ax.transAxes, va="top", ha="left", fontsize=8,
            family="monospace",
            bbox=dict(facecolor="white", edgecolor="#bdc3c7",
                      boxstyle="round,pad=0.5", alpha=0.88))

    fig.suptitle(
        f"RNA Tree BEAR Similarity (RTBS) — Penalty Visualisation\n"
        f"Target: {target_name}   |   Prediction: {pred_name}",
        fontsize=11, y=0.99)

    plt.savefig(out_path, dpi=150, bbox_inches="tight")
    if args.show:
        plt.show()
    plt.close(fig)

    # JSON report
    report_pred = []
    for nid, info in node_info.items():
        entry = {
            "pred_node_id":   nid,
            "pred_node_name": nodes_p[nid]["name"],
            "bear_pred":      info["bear_p"],
            "penalty":        round(info["penalty"], 4),
            "matched":        info["matched"],
        }
        if info["matched"]:
            tid = info["target_id"]
            entry["target_node_id"]   = tid
            entry["target_node_name"] = nodes_t[tid]["name"]
            entry["bear_target"]      = info["bear_t"]
        report_pred.append(entry)
    report_pred.sort(key=lambda x: -x["penalty"])

    report_target_unmatched = [
        {"target_node_id":   nid,
         "target_node_name": nodes_t[nid]["name"],
         "bear_target":      node_to_bear(nodes_t[nid]),
         "penalty":          round(UNMATCHED_PENALTY, 4),
         "matched":          False}
        for nid in unmatched_target_ids
    ]

    summary = {
        "target":                    target_name,
        "prediction":                pred_name,
        "n_target_nodes":            n_target,
        "n_pred_nodes":              n_pred,
        "n_unmatched_pred":          n_unmatched_pred,
        "n_unmatched_target":        n_unmatched_target,
        "sum_penalty_sym":            round(total_pen, 4),
        "sum_penalty_target":         round(total_pen_target, 4),
        "avg_penalty_sym":            round(avg_pen, 4),
        "avg_penalty_target":         round(avg_pen_target, 4),
        "RTBS_sym":          round(norm_penalty, 6),
        "RTBS_target":       round(norm_penalty_target, 6),
        "best_possible":             round(best_possible, 4),
        "worst_possible_sym":        round(worst_sym, 4),
        "worst_possible_target":     round(worst_target, 4),
        "unmatched_penalty":         round(UNMATCHED_PENALTY, 4),
        "prediction_nodes":          report_pred,
        "unmatched_target_nodes":    report_target_unmatched,
    }

    rpath = out_path.with_suffix(".json")
    with open(rpath, "w") as f:
        json.dump(summary, f, indent=2)

    print(f"  Saved: {out_path}  |  RTBS_sym={norm_penalty:.4f}  RTBS_tgt={norm_penalty_target:.4f}")
    return norm_penalty, norm_penalty_target, summary

def validate_prediction_file(path):
    try:
        with open(path, encoding="utf-8") as f:
            data = json.load(f)
        if not isinstance(data, list) or len(data) == 0:
            return False, "empty or invalid JSON structure"
        return True, None
    except Exception as e:
        return False, str(e)
    
#Run: single file or folder
if args.pred_dir:
    #folder mode
    pred_dir = pathlib.Path(args.pred_dir)
    if not pred_dir.is_dir():
        sys.exit(f"ERROR: {pred_dir} is not a directory.")

    pred_files = sorted(pred_dir.glob("*.json"))
    if not pred_files:
        sys.exit(f"ERROR: no JSON files found in {pred_dir}.")

    out_dir = pathlib.Path(args.output) if args.output else pathlib.Path(f"{target_name}_results")
    out_dir.mkdir(parents=True, exist_ok=True)
    print(f"Output folder: {out_dir}")
    print(f"Processing {len(pred_files)} prediction(s) against target: {target_name}\n")

    all_scores = []
    failed_files = []
    for pred_file in pred_files:
        print(f"[{pred_file.name}]")
        valid, reason = validate_prediction_file(pred_file)
        if not valid:
            print(f"  WARNING: skipped {pred_file.name} ({reason})")
            failed_files.append(
                {
                    "prediction": pred_file.name,
                    "error": reason,
                    "type": "ValidationError"
                }
            )
            continue
        try:
            vis_dir = out_dir / "visualizations"
            vis_dir.mkdir(parents=True, exist_ok=True)
            out_img = vis_dir / f"{pred_file.stem}_vs_{target_name}.png"
            norm_pen, norm_pen_tgt, _ = process_prediction(pred_file, out_img)
            all_scores.append((pred_file.stem, norm_pen, norm_pen_tgt))
        except Exception as e:

            print(f"  WARNING: skipped {pred_file.name} "
                f"({type(e).__name__}: {e})"
            )

            failed_files.append(
                {
                    "prediction": pred_file.name,
                    "error": str(e),
                    "type": type(e).__name__
                }
            )
            continue

    #summary ranking (sorted by target-anchored metric)
    all_scores.sort(key=lambda x: x[2], reverse=True)
    summary_path = out_dir / "ranking.txt"
    with open(summary_path, "w", encoding="utf-8") as f:
        f.write(f"{'Prediction':<40} {'RTBS_tgt':>10} {'RTBS_sym':>10}\n")
        f.write("-" * 64 + "\n")
        for name, sc_sym, sc_tgt in all_scores:
            f.write(
                f"{name:<40}"
                f"{sc_tgt:>10.4f}"
                f"{sc_sym:>10.4f}\n"
            )
    print(f"\nRanking saved: {summary_path}")

    if failed_files:
        failed_path = out_dir / "failed_predictions.txt"
        with open(failed_path, "w", encoding="utf-8") as f:
            f.write(f"{'Prediction':<40} {'Type':<25} Error\n")
            f.write("-" * 100 + "\n")
            for item in failed_files:
                f.write(
                f"{item['prediction']:<40} "
                f"{item['type']:<25} "
                f"{item['error']}\n"
            )
        print(
            f"\nSkipped {len(failed_files)} invalid prediction(s). "
            f"Details: {failed_path}"
        )

else:
    #single file mode
    pred_path = pathlib.Path(args.prediction)
    out_path  = pathlib.Path(args.output) if args.output else pathlib.Path("penalty_tree.png")
    print(f"Processing: {pred_path.name}")
    process_prediction(pred_path, out_path)