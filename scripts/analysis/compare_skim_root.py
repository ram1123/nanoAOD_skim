#!/usr/bin/env python3
import argparse
import csv
import math
import os
import sys

import numpy as np

try:
    import uproot
except ImportError as exc:
    raise SystemExit(
        "ERROR: This script requires `uproot`. "
        "Install it in your environment or run inside a setup that already provides it."
    ) from exc


DEFAULT_SUMMARY_BRANCHES = [
    "massZ1",
    "pTZ1",
    "pTL1",
    "pTL2",
    "pT_MET",
    "HZZ2l2qNu_nJets",
    "HZZ2l2qNu_nMediumBtagJets",
    "HZZ2l2nu_ifVBF",
    "HZZ2l2qNu_isELE",
    "passZZ2l2nuSelection",
    "Triggers_HZZ2l2nu_SingleLep",
    "HZZ2l2qNu_cutOppositeChargeFlag",
]

DEFAULT_VALUE_DIFF_BRANCHES = [
    "HZZ2l2qNu_cutOppositeChargeFlag",
    "HZZ2l2qNu_isELE",
    "passZZ2l2nuSelection",
    "massZ1",
    "pTZ1",
    "pTL1",
    "pTL2",
    "pT_MET",
    "HZZ2l2qNu_nJets",
    "HZZ2l2qNu_nMediumBtagJets",
    "HZZ2l2nu_ifVBF",
    "Triggers_HZZ2l2nu_SingleLep",
]


def parse_args():
    parser = argparse.ArgumentParser(
        description="Compare two skim ROOT files at the tree, branch, and event level."
    )
    parser.add_argument("file1", help="First ROOT file")
    parser.add_argument("file2", help="Second ROOT file")
    parser.add_argument(
        "--tree",
        default="Events",
        help="Primary tree to compare in detail (default: Events)",
    )
    parser.add_argument(
        "--max-branches",
        type=int,
        default=20,
        help="Maximum number of branch names to print in added/removed lists",
    )
    parser.add_argument(
        "--event-id-branches",
        nargs=3,
        default=["run", "luminosityBlock", "event"],
        metavar=("RUN", "LUMI", "EVENT"),
        help="Branches used to build event identifiers",
    )
    parser.add_argument(
        "--summary-branches",
        nargs="*",
        default=DEFAULT_SUMMARY_BRANCHES,
        help="Scalar branches for quick numeric summaries",
    )
    parser.add_argument(
        "--max-event-diff",
        type=int,
        default=10,
        help="Maximum number of unique event ids to print per file",
    )
    parser.add_argument(
        "--value-diff-branches",
        nargs="*",
        default=DEFAULT_VALUE_DIFF_BRANCHES,
        help="Scalar branches to compare event-by-event for common events",
    )
    parser.add_argument(
        "--max-value-diff-events",
        type=int,
        default=20,
        help="Maximum number of differing common events to print",
    )
    parser.add_argument(
        "--csv-output",
        default=None,
        help="Optional CSV file to write event-by-event branch differences and non-common events",
    )
    return parser.parse_args()


def fail(message):
    print(f"ERROR: {message}", file=sys.stderr)
    sys.exit(1)


def list_trees(root_file):
    trees = {}
    for key, obj in root_file.items():
        if getattr(obj, "classname", "").startswith("TTree"):
            trees[key.split(";")[0]] = obj
    return trees


def print_header(title):
    print(f"\n{title}")
    print("-" * len(title))


def format_list(values, max_items):
    if not values:
        return "None"
    shown = values[:max_items]
    suffix = "" if len(values) <= max_items else f" ... (+{len(values) - max_items} more)"
    return ", ".join(shown) + suffix


def compare_tree_inventory(file1_trees, file2_trees, max_branches):
    names1 = set(file1_trees)
    names2 = set(file2_trees)
    common = sorted(names1 & names2)
    only1 = sorted(names1 - names2)
    only2 = sorted(names2 - names1)

    print_header("Tree Inventory")
    print(f"Common trees ({len(common)}): {format_list(common, max_branches)}")
    print(f"Only in file1 ({len(only1)}): {format_list(only1, max_branches)}")
    print(f"Only in file2 ({len(only2)}): {format_list(only2, max_branches)}")


def compare_branch_inventory(tree1, tree2, max_branches):
    branches1 = set(tree1.keys())
    branches2 = set(tree2.keys())
    common = sorted(branches1 & branches2)
    only1 = sorted(branches1 - branches2)
    only2 = sorted(branches2 - branches1)

    print_header(f"Branch Inventory: {tree1.name}")
    print(f"Entries: file1={tree1.num_entries}, file2={tree2.num_entries}")
    print(f"Common branches ({len(common)}): {len(common)}")
    print(f"Only in file1 ({len(only1)}): {format_list(only1, max_branches)}")
    print(f"Only in file2 ({len(only2)}): {format_list(only2, max_branches)}")
    return common


def branch_is_scalar(branch):
    interpretation = getattr(branch, "interpretation", None)
    if interpretation is None:
        return False
    return "AsJagged" not in interpretation.__class__.__name__


def summarize_numeric_branch(tree, branch_name):
    branch = tree[branch_name]
    if not branch_is_scalar(branch):
        return None

    array = branch.array(library="np")
    if array.dtype.kind not in "biuf":
        return None
    if array.size == 0:
        return {"count": 0, "mean": math.nan, "min": math.nan, "max": math.nan}

    finite = array[np.isfinite(array)] if array.dtype.kind == "f" else array
    if finite.size == 0:
        return {"count": int(array.size), "mean": math.nan, "min": math.nan, "max": math.nan}

    return {
        "count": int(array.size),
        "mean": float(np.mean(finite)),
        "min": float(np.min(finite)),
        "max": float(np.max(finite)),
    }


def print_summary_comparison(tree1, tree2, branch_names):
    print_header(f"Scalar Summaries: {tree1.name}")
    printed = False
    for branch_name in branch_names:
        if branch_name not in tree1 or branch_name not in tree2:
            continue
        summary1 = summarize_numeric_branch(tree1, branch_name)
        summary2 = summarize_numeric_branch(tree2, branch_name)
        if summary1 is None or summary2 is None:
            continue
        printed = True
        print(
            f"{branch_name}: "
            f"file1(mean={summary1['mean']:.6g}, min={summary1['min']:.6g}, max={summary1['max']:.6g}, n={summary1['count']}) | "
            f"file2(mean={summary2['mean']:.6g}, min={summary2['min']:.6g}, max={summary2['max']:.6g}, n={summary2['count']})"
        )
    if not printed:
        print("No comparable scalar numeric summary branches were found.")


def load_event_ids(tree, event_id_branches):
    missing = [branch for branch in event_id_branches if branch not in tree]
    if missing:
        return None, missing

    arrays = tree.arrays(event_id_branches, library="np")
    event_ids = set(zip(arrays[event_id_branches[0]], arrays[event_id_branches[1]], arrays[event_id_branches[2]]))
    return event_ids, None


def print_event_overlap(tree1, tree2, event_id_branches, max_event_diff):
    print_header(f"Event Overlap: {tree1.name}")
    ids1, missing1 = load_event_ids(tree1, event_id_branches)
    ids2, missing2 = load_event_ids(tree2, event_id_branches)

    if missing1 or missing2:
        missing = missing1 if missing1 else missing2
        print(f"Skipping event-id comparison because branches are missing: {', '.join(missing)}")
        return

    common = ids1 & ids2
    only1 = sorted(ids1 - ids2)[:max_event_diff]
    only2 = sorted(ids2 - ids1)[:max_event_diff]

    print(f"Events in file1: {len(ids1)}")
    print(f"Events in file2: {len(ids2)}")
    print(f"Common event ids: {len(common)}")
    print(f"Only in file1: {len(ids1 - ids2)}")
    for event_id in only1:
        print(f"  file1-only event: {event_id}")
    print(f"Only in file2: {len(ids2 - ids1)}")
    for event_id in only2:
        print(f"  file2-only event: {event_id}")


def branch_arrays_for_events(tree, branch_names, event_id_branches):
    requested = list(dict.fromkeys(event_id_branches + branch_names))
    available = [branch for branch in requested if branch in tree]
    arrays = tree.arrays(available, library="np")

    event_map = {}
    for idx in range(tree.num_entries):
        event_id = tuple(arrays[branch][idx] for branch in event_id_branches)
        values = {}
        for branch in branch_names:
            if branch in arrays:
                values[branch] = arrays[branch][idx]
        event_map[event_id] = values

    missing = [branch for branch in branch_names if branch not in arrays]
    return event_map, missing


def values_differ(value1, value2):
    if isinstance(value1, np.generic):
        value1 = value1.item()
    if isinstance(value2, np.generic):
        value2 = value2.item()

    if isinstance(value1, float) or isinstance(value2, float):
        if math.isnan(value1) and math.isnan(value2):
            return False
        return not math.isclose(value1, value2, rel_tol=1e-9, abs_tol=1e-9)
    return value1 != value2


def normalize_value(value):
    if isinstance(value, np.generic):
        value = value.item()
    return value


def compute_delta(value1, value2):
    value1 = normalize_value(value1)
    value2 = normalize_value(value2)

    if isinstance(value1, bool) or isinstance(value2, bool):
        return int(value2) - int(value1)
    if isinstance(value1, (int, float)) and isinstance(value2, (int, float)):
        if isinstance(value1, float) and math.isnan(value1):
            return ""
        if isinstance(value2, float) and math.isnan(value2):
            return ""
        return value2 - value1
    return ""


def print_value_differences(
    tree1,
    tree2,
    event_id_branches,
    branch_names,
    max_value_diff_events,
    csv_output=None,
):
    print_header(f"Event-by-Event Value Differences: {tree1.name}")

    ids1, missing1 = load_event_ids(tree1, event_id_branches)
    ids2, missing2 = load_event_ids(tree2, event_id_branches)
    if missing1 or missing2:
        missing = missing1 if missing1 else missing2
        print(f"Skipping value-diff comparison because event-id branches are missing: {', '.join(missing)}")
        return

    event_map1, branch_missing1 = branch_arrays_for_events(tree1, branch_names, event_id_branches)
    event_map2, branch_missing2 = branch_arrays_for_events(tree2, branch_names, event_id_branches)

    branch_missing = sorted(set(branch_missing1 + branch_missing2))
    if branch_missing:
        print(f"Skipped missing branches: {', '.join(branch_missing)}")

    csv_rows = []
    diff_count = 0
    for sync_instance, event_id in enumerate(sorted(ids1 & ids2), start=1):
        branch_diffs = []
        row = {
            "run": event_id[0],
            "luminosityBlock": event_id[1],
            "event": event_id[2],
            "_sync_instance": sync_instance,
        }
        for branch in branch_names:
            if branch not in event_map1[event_id] or branch not in event_map2[event_id]:
                continue
            value1 = event_map1[event_id][branch]
            value2 = event_map2[event_id][branch]
            norm1 = normalize_value(value1)
            norm2 = normalize_value(value2)
            row[f"{branch}_1"] = norm1
            row[f"{branch}_2"] = norm2
            row[f"delta_{branch}"] = compute_delta(norm1, norm2)
            if values_differ(value1, value2):
                branch_diffs.append((branch, norm1, norm2))

        if branch_diffs:
            diff_count += 1
            row["event_status"] = "different_common_event"
            csv_rows.append(row)
            if diff_count <= max_value_diff_events:
                print(f"Event {event_id}")
                for branch, value1, value2 in branch_diffs:
                    print(f"  {branch}: file1={value1} file2={value2}")

    print(f"Events with at least one requested branch difference: {diff_count}")
    if diff_count > max_value_diff_events:
        print(f"Printed first {max_value_diff_events} differing events.")

    file1_only_ids = sorted(ids1 - ids2)
    file2_only_ids = sorted(ids2 - ids1)

    if file1_only_ids:
        print(f"Adding {len(file1_only_ids)} file1-only events to the CSV output.")
    for sync_instance, event_id in enumerate(file1_only_ids, start=1):
        row = {
            "run": event_id[0],
            "luminosityBlock": event_id[1],
            "event": event_id[2],
            "_sync_instance": f"file1_only_{sync_instance}",
            "event_status": "only_in_file1",
        }
        for branch in branch_names:
            if branch in branch_missing:
                continue
            value1 = event_map1.get(event_id, {}).get(branch, "")
            row[f"{branch}_1"] = normalize_value(value1) if value1 != "" else ""
            row[f"{branch}_2"] = ""
            row[f"delta_{branch}"] = ""
        csv_rows.append(row)

    if file2_only_ids:
        print(f"Adding {len(file2_only_ids)} file2-only events to the CSV output.")
    for sync_instance, event_id in enumerate(file2_only_ids, start=1):
        row = {
            "run": event_id[0],
            "luminosityBlock": event_id[1],
            "event": event_id[2],
            "_sync_instance": f"file2_only_{sync_instance}",
            "event_status": "only_in_file2",
        }
        for branch in branch_names:
            if branch in branch_missing:
                continue
            value2 = event_map2.get(event_id, {}).get(branch, "")
            row[f"{branch}_1"] = ""
            row[f"{branch}_2"] = normalize_value(value2) if value2 != "" else ""
            row[f"delta_{branch}"] = ""
        csv_rows.append(row)

    if csv_output:
        fieldnames = ["run", "luminosityBlock", "event", "_sync_instance", "event_status"]
        for branch in branch_names:
            if branch in branch_missing:
                continue
            fieldnames.extend([f"{branch}_1", f"{branch}_2", f"delta_{branch}"])
        with open(csv_output, "w", newline="") as handle:
            writer = csv.DictWriter(
                handle,
                fieldnames=fieldnames,
            )
            writer.writeheader()
            writer.writerows(csv_rows)
        print(f"CSV diff report written to: {csv_output}")


def main():
    args = parse_args()
    for path in [args.file1, args.file2]:
        if not os.path.isfile(path):
            fail(f"File not found: {path}")

    with uproot.open(args.file1) as file1, uproot.open(args.file2) as file2:
        trees1 = list_trees(file1)
        trees2 = list_trees(file2)
        compare_tree_inventory(trees1, trees2, args.max_branches)

        if args.tree not in trees1:
            fail(f"Tree `{args.tree}` not found in file1")
        if args.tree not in trees2:
            fail(f"Tree `{args.tree}` not found in file2")

        tree1 = trees1[args.tree]
        tree2 = trees2[args.tree]

        compare_branch_inventory(tree1, tree2, args.max_branches)
        print_summary_comparison(tree1, tree2, args.summary_branches)
        print_event_overlap(tree1, tree2, args.event_id_branches, args.max_event_diff)
        print_value_differences(
            tree1,
            tree2,
            args.event_id_branches,
            args.value_diff_branches,
            args.max_value_diff_events,
            args.csv_output,
        )


if __name__ == "__main__":
    main()
