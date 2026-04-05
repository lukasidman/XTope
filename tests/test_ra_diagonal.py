#!/usr/bin/env python3
"""
Benchmark the reduced-alphabet diagonal chain filter against the 500-sequence
benchmark dataset. Measures:
  1. Filter rate: what % of background pairs are eliminated
  2. False negative rate: which known family pairs are missed
  3. Comparison with current two-pass filter
  4. Tuning: sweep k and min_chain to find the sweet spot

Uses tag-stripped sequences from the benchmark CSV.
"""

import sys
import time
import csv
from pathlib import Path
from collections import defaultdict

# Add src to path
sys.path.insert(0, str(Path(__file__).parent / "src"))

from antigen_screener.ra_diagonal_filter import (
    RADiagonalFilter,
    reduce_sequence,
    N_GROUPS,
)
from antigen_screener.tag_stripper import strip_tag
from antigen_screener.kmer_filter import TwoPassFilter


# ---------------------------------------------------------------------------
# Load benchmark data
# ---------------------------------------------------------------------------
BENCHMARK_CSV = Path(__file__).parent / "benchmark_antigens.csv"

FAMILY_PREFIXES = ["CYTK_", "KINAS_", "ALLO_", "IGF_", "RECPT_", "GRAD_", "MOSAIC_", "CONV_"]

def load_benchmark():
    """Load and tag-strip benchmark sequences."""
    records = []
    with open(BENCHMARK_CSV, encoding="utf-8-sig", newline="") as f:
        reader = csv.reader(f, delimiter=";")
        header = next(reader)
        for row in reader:
            if len(row) < 3:
                continue
            name = row[1].strip()
            seq = row[2].strip()
            if not name or not seq:
                continue
            stripped = strip_tag(seq)["stripped"]
            records.append((name, stripped))
    return records


def get_family(name: str) -> str | None:
    """Extract family prefix from sequence name, or None for background."""
    for prefix in FAMILY_PREFIXES:
        if name.startswith(prefix):
            return prefix.rstrip("_")
    return None


def get_expected_pairs(records):
    """Build set of all expected within-family pairs."""
    families = defaultdict(list)
    for name, _ in records:
        fam = get_family(name)
        if fam:
            families[fam].append(name)

    expected = set()
    for fam, members in families.items():
        for i in range(len(members)):
            for j in range(i + 1, len(members)):
                a, b = sorted([members[i], members[j]])
                expected.add((a, b))
    return expected, families


# ---------------------------------------------------------------------------
# Run filter and collect results
# ---------------------------------------------------------------------------
def benchmark_filter(records, k, min_chain, label=""):
    """Run RADiagonalFilter with given params and return metrics."""
    filt = RADiagonalFilter(k=k, min_chain=min_chain)

    t0 = time.time()
    filt.add_batch(records)
    index_time = time.time() - t0

    expected_pairs, families = get_expected_pairs(records)
    n = len(records)
    total_pairs = n * (n - 1) // 2

    # Run all queries and collect passing pairs
    passing_pairs = set()
    t0 = time.time()
    for name, seq in records:
        passing = filt.query(seq, exclude_id=name)
        for target_id in passing:
            pair = tuple(sorted([name, target_id]))
            passing_pairs.add(pair)
    query_time = time.time() - t0

    # Classify results
    true_positives = passing_pairs & expected_pairs
    false_negatives = expected_pairs - passing_pairs
    background_passing = passing_pairs - expected_pairs
    background_total = total_pairs - len(expected_pairs)
    background_rejected = background_total - len(background_passing)

    filter_rate = background_rejected / background_total * 100 if background_total > 0 else 0
    sensitivity = len(true_positives) / len(expected_pairs) * 100 if expected_pairs else 0

    # Per-family sensitivity
    family_stats = {}
    for fam, members in families.items():
        fam_pairs = set()
        for i in range(len(members)):
            for j in range(i + 1, len(members)):
                fam_pairs.add(tuple(sorted([members[i], members[j]])))
        fam_tp = len(fam_pairs & passing_pairs)
        fam_total = len(fam_pairs)
        family_stats[fam] = (fam_tp, fam_total)

    return {
        "label": label,
        "k": k,
        "min_chain": min_chain,
        "index_time": index_time,
        "query_time": query_time,
        "total_pairs": total_pairs,
        "expected_pairs": len(expected_pairs),
        "passing_pairs": len(passing_pairs),
        "true_positives": len(true_positives),
        "false_negatives": len(false_negatives),
        "background_passing": len(background_passing),
        "background_rejected": background_rejected,
        "filter_rate": filter_rate,
        "sensitivity": sensitivity,
        "family_stats": family_stats,
        "fn_details": false_negatives,
    }


def benchmark_twopass(records):
    """Run existing TwoPassFilter for comparison."""
    filt = TwoPassFilter()

    t0 = time.time()
    filt.add_batch(records)
    index_time = time.time() - t0

    expected_pairs, families = get_expected_pairs(records)
    n = len(records)
    total_pairs = n * (n - 1) // 2

    passing_pairs = set()
    t0 = time.time()
    for name, seq in records:
        results = filt.query(seq, exclude_id=name)
        for target_id, _ in results:
            pair = tuple(sorted([name, target_id]))
            passing_pairs.add(pair)
    query_time = time.time() - t0

    true_positives = passing_pairs & expected_pairs
    false_negatives = expected_pairs - passing_pairs
    background_passing = passing_pairs - expected_pairs
    background_total = total_pairs - len(expected_pairs)
    background_rejected = background_total - len(background_passing)

    filter_rate = background_rejected / background_total * 100 if background_total > 0 else 0
    sensitivity = len(true_positives) / len(expected_pairs) * 100 if expected_pairs else 0

    family_stats = {}
    for fam, members in families.items():
        fam_pairs = set()
        for i in range(len(members)):
            for j in range(i + 1, len(members)):
                fam_pairs.add(tuple(sorted([members[i], members[j]])))
        fam_tp = len(fam_pairs & passing_pairs)
        fam_total = len(fam_pairs)
        family_stats[fam] = (fam_tp, fam_total)

    return {
        "label": "TwoPass (current)",
        "index_time": index_time,
        "query_time": query_time,
        "total_pairs": total_pairs,
        "expected_pairs": len(expected_pairs),
        "passing_pairs": len(passing_pairs),
        "true_positives": len(true_positives),
        "false_negatives": len(false_negatives),
        "background_passing": len(background_passing),
        "background_rejected": background_rejected,
        "filter_rate": filter_rate,
        "sensitivity": sensitivity,
        "family_stats": family_stats,
        "fn_details": false_negatives,
    }


# ---------------------------------------------------------------------------
# Detailed analysis of specific families
# ---------------------------------------------------------------------------
def analyze_igf_family(records):
    """Deep-dive into IGF family detection across parameter settings."""
    igf_records = [(n, s) for n, s in records if n.startswith("IGF_")]
    all_records_dict = {n: s for n, s in records}

    print("\n" + "=" * 70)
    print("  IGF FAMILY DEEP DIVE (8aa shared epitope: RPCHQFNV)")
    print("=" * 70)

    # Show reduced sequences for IGF members
    for name, seq in igf_records:
        reduced = reduce_sequence(seq)
        # Find the epitope region
        epitope = "RPCHQFNV"
        pos = seq.find(epitope)
        if pos >= 0:
            region = seq[max(0, pos-3):pos+len(epitope)+3]
            reduced_region = reduced[max(0, pos-3):pos+len(epitope)+3]
            print(f"\n  {name}:")
            print(f"    Epitope at pos {pos}: ...{region}...")
            print(f"    Reduced:             ...{reduced_region}...")

    # Test different k values
    print(f"\n  --- Chain detection across k values ---")
    for k in [3, 4, 5, 6]:
        filt = RADiagonalFilter(k=k, min_chain=2)
        filt.add_batch(records)
        print(f"\n  k={k}, min_chain=2:")
        for name, seq in igf_records:
            details = filt.query_with_details(seq, exclude_id=name)
            igf_hits = [(t, c, d) for t, c, d in details if t.startswith("IGF_")]
            bg_chains = [c for t, c, d in details if not any(t.startswith(p) for p in FAMILY_PREFIXES)]
            max_bg = max(bg_chains) if bg_chains else 0
            bg_above_1 = sum(1 for c in bg_chains if c >= 2)
            print(f"    {name}: IGF hits={[(t,c) for t,c,_ in igf_hits]}  "
                  f"max_bg_chain={max_bg}  bg_with_chain≥2={bg_above_1}")


# ---------------------------------------------------------------------------
# Main
# ---------------------------------------------------------------------------
def main():
    print("Loading benchmark dataset...")
    records = load_benchmark()
    print(f"Loaded {len(records)} sequences")

    expected_pairs, families = get_expected_pairs(records)
    print(f"Expected family pairs: {len(expected_pairs)}")
    for fam, members in sorted(families.items()):
        n = len(members)
        print(f"  {fam}: {n} members, {n*(n-1)//2} pairs")

    # --- Parameter sweep ---
    print("\n" + "=" * 70)
    print("  PARAMETER SWEEP: k × min_chain")
    print("=" * 70)

    configs = [
        (4, 2, "k=4 mc=2 (generous)"),
        (4, 3, "k=4 mc=3"),
        (5, 2, "k=5 mc=2 (proposed)"),
        (5, 3, "k=5 mc=3"),
        (6, 2, "k=6 mc=2"),
        (6, 3, "k=6 mc=3"),
        (3, 3, "k=3 mc=3"),
        (3, 4, "k=3 mc=4"),
    ]

    results = []
    for k, mc, label in configs:
        r = benchmark_filter(records, k=k, min_chain=mc, label=label)
        results.append(r)

    # Also run current two-pass for comparison
    tp_result = benchmark_twopass(records)

    # Print comparison table
    print(f"\n{'Label':<22} {'Filter%':>8} {'Sens%':>7} {'TP':>5} {'FN':>5} "
          f"{'BG pass':>8} {'Time':>7}")
    print("-" * 75)
    for r in results:
        print(f"{r['label']:<22} {r['filter_rate']:>7.1f}% {r['sensitivity']:>6.1f}% "
              f"{r['true_positives']:>5} {r['false_negatives']:>5} "
              f"{r['background_passing']:>8} {r['query_time']:>6.2f}s")
    print("-" * 75)
    print(f"{'TwoPass (current)':<22} {tp_result['filter_rate']:>7.1f}% "
          f"{tp_result['sensitivity']:>6.1f}% "
          f"{tp_result['true_positives']:>5} {tp_result['false_negatives']:>5} "
          f"{tp_result['background_passing']:>8} {tp_result['query_time']:>6.2f}s")

    # Per-family breakdown for key configs
    print("\n" + "=" * 70)
    print("  PER-FAMILY SENSITIVITY")
    print("=" * 70)

    key_results = [r for r in results if r["label"] in
                   ["k=5 mc=2 (proposed)", "k=4 mc=2 (generous)", "k=6 mc=2"]]
    key_results.append(tp_result)

    fam_names = sorted(families.keys())
    header = f"{'Family':<10}"
    for r in key_results:
        header += f" {r['label'][:16]:>16}"
    print(header)
    print("-" * (10 + 17 * len(key_results)))

    for fam in fam_names:
        row = f"{fam:<10}"
        for r in key_results:
            tp, total = r["family_stats"][fam]
            row += f" {tp:>5}/{total:<4} {'✓' if tp == total else '✗':>3}"
        print(row)

    # Show which pairs are missed by the proposed config
    proposed = [r for r in results if r["label"] == "k=5 mc=2 (proposed)"][0]
    if proposed["fn_details"]:
        print(f"\n  Pairs MISSED by k=5 mc=2:")
        for a, b in sorted(proposed["fn_details"]):
            fam_a = get_family(a) or "BG"
            fam_b = get_family(b) or "BG"
            print(f"    {a} <-> {b}  (families: {fam_a}, {fam_b})")

    # --- IGF deep dive ---
    analyze_igf_family(records)

    # --- Selectivity analysis ---
    print("\n" + "=" * 70)
    print("  SELECTIVITY ANALYSIS (k=5, mc=2)")
    print("=" * 70)

    filt = RADiagonalFilter(k=5, min_chain=2)
    filt.add_batch(records)

    # Distribution of passing targets per query
    pass_counts = []
    for name, seq in records:
        passing = filt.query(seq, exclude_id=name)
        pass_counts.append(len(passing))

    avg_pass = sum(pass_counts) / len(pass_counts)
    max_pass = max(pass_counts)
    min_pass = min(pass_counts)
    median_pass = sorted(pass_counts)[len(pass_counts) // 2]

    print(f"  Targets passing per query (out of {len(records)-1}):")
    print(f"    Mean:   {avg_pass:.1f} ({avg_pass/(len(records)-1)*100:.1f}%)")
    print(f"    Median: {median_pass} ({median_pass/(len(records)-1)*100:.1f}%)")
    print(f"    Min:    {min_pass}")
    print(f"    Max:    {max_pass}")
    print(f"    → Average rejection rate: {(1 - avg_pass/(len(records)-1))*100:.1f}%")

    # Distribution histogram
    bins = [0, 5, 10, 20, 50, 100, 200, 500]
    print(f"\n  Distribution of passing targets per query:")
    for i in range(len(bins) - 1):
        count = sum(1 for c in pass_counts if bins[i] <= c < bins[i + 1])
        bar = "█" * count
        print(f"    {bins[i]:>3}-{bins[i+1]:>3}: {count:>4} {bar}")
    count = sum(1 for c in pass_counts if c >= bins[-1])
    print(f"    {bins[-1]:>3}+  : {count:>4} {'█' * count}")


if __name__ == "__main__":
    main()
