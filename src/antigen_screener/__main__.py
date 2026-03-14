#!/usr/bin/env python3
"""
antigen-screener CLI
====================

Commands:
  run      - Run all-vs-all precomputation across the full database
  query    - Query a single antigen against precomputed results
  export   - Export results to CSV

Examples:
  python -m antigen_screener run --input antigens.csv --db results.db
  python -m antigen_screener run --input antigens.csv --db results.db --tag MHHHHHHGSSG
  python -m antigen_screener query --db results.db --id AG_001
  python -m antigen_screener query --db results.db --seq MKALVPVFAGLLLVAGLAAVHSQSLD
  python -m antigen_screener export --db results.db --output results.csv --min-score 1.5
"""

import argparse
import sys
import time
from pathlib import Path


def cmd_run(args):
    if args.backend == "vectorized":
        _cmd_run_vectorized(args)
    else:
        _cmd_run_kmer(args)


def _cmd_run_kmer(args):
    from antigen_screener.pipeline import run_pipeline

    run_pipeline(
        input_path=args.input,
        db_path=args.db,
        tag_sequence=args.tag,
        id_col=args.id_col,
        seq_col=args.seq_col,
        kmer_threshold=args.kmer_threshold,
        min_norm_score=args.min_score,
        min_aligned_len=args.min_aligned,
        matrix=args.matrix,
        resume=not args.no_resume,
    )


def _cmd_run_vectorized(args):
    import datetime
    from antigen_screener.db_loader import load_sequences
    from antigen_screener.tag_stripper import strip_tag, set_tag
    from antigen_screener.store import ResultsStore
    from antigen_screener.vectorized_sw import run_vectorized_pipeline

    print(f"\n{'='*60}")
    print(f"  Antigen Cross-Reactivity Screener (vectorized NumPy)")
    print(f"  Started: {datetime.datetime.now().strftime('%Y-%m-%d %H:%M:%S')}")
    print(f"  Backend: CPU (NumPy)")
    print(f"{'='*60}\n")

    if args.tag:
        set_tag(args.tag)
        print(f"  Tag set to: {args.tag}\n")

    print(f"[1/3] Loading sequences from: {args.input}")
    raw_records = load_sequences(args.input, id_col=args.id_col, seq_col=args.seq_col, verbose=True)
    total = len(raw_records)
    print(f"  Total antigens: {total:,}\n")

    print(f"[2/3] Stripping His6-ABP tags...")
    store = ResultsStore(args.db)
    stripped_seqs: dict[str, str] = {}
    tag_found_count = 0
    antigen_batch = []

    for antigen_id, sequence in raw_records:
        result = strip_tag(sequence)
        stripped = result["stripped"]
        found = result["tag_found"]
        if found:
            tag_found_count += 1
        antigen_batch.append((antigen_id, sequence, stripped, found))
        stripped_seqs[antigen_id] = stripped

    store.upsert_antigens_batch(antigen_batch)
    print(f"  Tags detected and removed in {tag_found_count:,} / {total:,} sequences\n")

    print(f"[3/3] Running vectorized all-vs-all alignment...")
    print(f"  Parameters: min_norm_score={args.min_score}")
    stats = run_vectorized_pipeline(
        sequences=stripped_seqs,
        store=store,
        min_norm_score=args.min_score,
        resume=not args.no_resume,
    )

    print(f"\n{'='*60}")
    print(f"  Run complete!")
    print(f"  Total time:          {stats['elapsed_seconds']:.1f}s")
    print(f"  Antigens screened:   {total:,}")
    print(f"  Similar pairs found: {stats['total_pairs']:,}")
    print(f"  Results saved to:    {args.db}")
    print(f"{'='*60}\n")

    store.set_meta("last_run", datetime.datetime.now().isoformat())
    store.set_meta("total_antigens", str(total))
    store.set_meta("total_pairs", str(stats["total_pairs"]))
    store.set_meta("backend", "vectorized_numpy")
    store.close()


def cmd_query(args):
    from antigen_screener.store       import ResultsStore
    from antigen_screener.tag_stripper import strip_tag, set_tag
    from antigen_screener.kmer_filter  import KmerIndex
    from antigen_screener.aligner      import batch_align

    db_path = Path(args.db)
    if not db_path.exists():
        print(f"Error: Database not found at {args.db}")
        print("Run `python -m antigen_screener run` first.")
        sys.exit(1)

    store = ResultsStore(args.db)

    # --- Query by ID (precomputed) ---
    if args.id:
        print(f"\nQuerying precomputed results for: {args.id}")
        results = store.query_similar(args.id, min_score=args.min_score, top_n=args.top_n)
        if not results:
            print("  No similar antigens found above threshold.")
        else:
            print(f"\n  {'Rank':<5} {'Partner ID':<25} {'Norm Score':<12} {'Aligned Len':<13} {'Q Len':<7} {'T Len'}")
            print(f"  {'-'*5} {'-'*25} {'-'*12} {'-'*13} {'-'*7} {'-'*6}")
            for i, r in enumerate(results, 1):
                print(f"  {i:<5} {r['partner_id']:<25} {r['normalized_score']:<12.3f} "
                      f"{r['aligned_region_len']:<13} {r['query_length']:<7} {r['target_length']}")
        store.close()
        return

    # --- Query by raw sequence (live alignment) ---
    if args.seq:
        if args.tag:
            set_tag(args.tag)
        stripped_info = strip_tag(args.seq)
        stripped = stripped_info["stripped"]
        print(f"\nQuery sequence (stripped): {stripped}")
        print(f"Tag found: {stripped_info['tag_found']}\n")

        # Load all stripped sequences from DB
        print("Loading database sequences...")
        all_seqs = store.get_all_stripped()
        print(f"Running live alignment against {len(all_seqs):,} sequences...\n")

        from antigen_screener.aligner import batch_align
        results = batch_align(
            query_seq=stripped,
            query_id="<query>",
            candidates=all_seqs,
            min_norm_score=args.min_score,
            min_aligned_len=args.min_aligned,
        )

        results = results[:args.top_n]
        if not results:
            print("  No similar antigens found above threshold.")
        else:
            print(f"  {'Rank':<5} {'Target ID':<25} {'Norm Score':<12} {'Aligned Len':<13} {'T Len'}")
            print(f"  {'-'*5} {'-'*25} {'-'*12} {'-'*13} {'-'*6}")
            for i, r in enumerate(results, 1):
                print(f"  {i:<5} {r.target_id:<25} {r.normalized_score:<12.3f} "
                      f"{r.aligned_region_len:<13} {r.target_length}")

        store.close()
        return

    print("Error: Provide either --id or --seq for query mode.")
    sys.exit(1)


def cmd_export(args):
    from antigen_screener.store import ResultsStore

    db_path = Path(args.db)
    if not db_path.exists():
        print(f"Error: Database not found at {args.db}")
        sys.exit(1)

    store = ResultsStore(args.db)
    output = args.output or str(db_path.with_suffix(".csv"))
    n = store.export_csv(output, min_score=args.min_score)
    print(f"Exported {n:,} similarity pairs to: {output}")
    store.close()


def cmd_stats(args):
    from antigen_screener.store import ResultsStore
    db_path = Path(args.db)
    if not db_path.exists():
        print(f"Error: Database not found at {args.db}")
        sys.exit(1)

    store = ResultsStore(args.db)
    print(f"\n  Database: {args.db}")
    print(f"  Antigens indexed:  {store.antigen_count():,}")
    print(f"  Similarity pairs:  {store.similarity_count():,}")
    print(f"  Last run:          {store.get_meta('last_run') or 'unknown'}")
    store.close()


def main():
    parser = argparse.ArgumentParser(
        prog="antigen-screener",
        description="Antibody cross-reactivity screener via antigen sequence similarity",
    )
    sub = parser.add_subparsers(dest="command", required=True)

    # ---- run ----
    p_run = sub.add_parser("run", help="Run all-vs-all precomputation")
    p_run.add_argument("--input",         required=True,  help="Path to CSV/TSV antigen database")
    p_run.add_argument("--db",            required=True,  help="Path to output SQLite database")
    p_run.add_argument("--tag",           default=None,   help="His6-ABP tag sequence to strip (overrides default)")
    p_run.add_argument("--id-col",        default=None,   help="Column name for antigen IDs")
    p_run.add_argument("--seq-col",       default=None,   help="Column name for sequences")
    p_run.add_argument("--kmer-threshold",type=float, default=0.04, help="Jaccard threshold for pre-filter (default: 0.04)")
    p_run.add_argument("--min-score",     type=float, default=1.0,  help="Min normalised SW score to record (default: 1.0)")
    p_run.add_argument("--min-aligned",   type=int,   default=8,    help="Min aligned region length in aa (default: 8)")
    p_run.add_argument("--matrix",        default="blosum62", choices=["blosum62","blosum45","blosum80"], help="Substitution matrix")
    p_run.add_argument("--no-resume",     action="store_true", help="Start fresh even if DB has partial results")
    p_run.add_argument("--backend",       default="kmer", choices=["kmer", "vectorized"],
                       help="Alignment backend: 'kmer' (k-mer filter + per-pair SW) or 'vectorized' (NumPy batched SW, no pre-filter)")
    p_run.set_defaults(func=cmd_run)

    # ---- query ----
    p_q = sub.add_parser("query", help="Query one antigen against precomputed results")
    p_q.add_argument("--db",          required=True, help="Path to SQLite database")
    p_q_excl = p_q.add_mutually_exclusive_group(required=True)
    p_q_excl.add_argument("--id",     help="Antigen ID to look up precomputed results")
    p_q_excl.add_argument("--seq",    help="Raw sequence to run live alignment")
    p_q.add_argument("--tag",         default=None, help="Tag to strip from --seq (if using live mode)")
    p_q.add_argument("--min-score",   type=float, default=1.0, help="Min normalised score (default: 1.0)")
    p_q.add_argument("--min-aligned", type=int,   default=8,   help="Min aligned region length (default: 8)")
    p_q.add_argument("--top-n",       type=int,   default=25,  help="Number of results to show (default: 25)")
    p_q.set_defaults(func=cmd_query)

    # ---- export ----
    p_e = sub.add_parser("export", help="Export results to CSV")
    p_e.add_argument("--db",        required=True, help="Path to SQLite database")
    p_e.add_argument("--output",    default=None,  help="Output CSV path (default: same name as db)")
    p_e.add_argument("--min-score", type=float, default=0.0, help="Min normalised score to include (default: 0.0)")
    p_e.set_defaults(func=cmd_export)

    # ---- stats ----
    p_s = sub.add_parser("stats", help="Show database statistics")
    p_s.add_argument("--db", required=True, help="Path to SQLite database")
    p_s.set_defaults(func=cmd_stats)

    args = parser.parse_args()
    t0 = time.perf_counter()
    args.func(args)
    elapsed = time.perf_counter() - t0
    if elapsed >= 1.0:
        print(f"\nCompleted in {elapsed:.3f}s")
    else:
        print(f"\nCompleted in {elapsed * 1000:.1f}ms")


if __name__ == "__main__":
    main()
