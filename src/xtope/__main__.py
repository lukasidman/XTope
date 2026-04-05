#!/usr/bin/env python3
"""
xtope CLI
====================

Commands:
  run      - Run all-vs-all precomputation across the full database
  query    - Query a single antigen against precomputed results
  export   - Export results to CSV

Examples:
  python -m xtope run --input antigens.csv --db results.db
  python -m xtope run --input antigens.csv --db results.db --tag MHHHHHHGSSG
  python -m xtope query --db results.db --id AG_001
  python -m xtope query --db results.db --seq MKALVPVFAGLLLVAGLAAVHSQSLD
  python -m xtope export --db results.db --output results.csv --max-evalue 0.001
"""

import argparse
import sys
import time
from pathlib import Path


def cmd_run(args):
    if args.backend == "kmer":
        _cmd_run_kmer(args)
    else:
        _cmd_run_vectorized(args)


def _cmd_run_kmer(args):
    import datetime
    from xtope.pipeline import run_pipeline

    print(f"\n{'='*60}")
    print(f"  XTope")
    print(f"  Started:    {datetime.datetime.now().strftime('%Y-%m-%d %H:%M:%S')}")
    print(f"  Backend:    kmer (quick screen)")
    print(f"{'='*60}\n")

    run_pipeline(
        input_path=args.input,
        db_path=args.db,
        tag_sequence=args.tag,
        id_col=args.id_col,
        seq_col=args.seq_col,
        kmer_threshold=args.kmer_threshold,
        max_evalue=args.max_evalue,
        min_aligned_len=args.min_aligned,
        matrix=args.matrix,
        resume=not args.no_resume,
    )


def _cmd_run_vectorized(args):
    import datetime
    from xtope.db_loader import load_sequences
    from xtope.tag_stripper import strip_tag, set_tag
    from xtope.store import ResultsStore
    from xtope.vectorized_sw import run_vectorized_pipeline
    from xtope.evalue import format_evalue

    prefilter_label = "off (exhaustive)" if args.no_prefilter else "RA diagonal + BLOSUM62 rescue"
    print(f"\n{'='*60}")
    print(f"  XTope")
    print(f"  Started:    {datetime.datetime.now().strftime('%Y-%m-%d %H:%M:%S')}")
    print(f"  Backend:    vectorized NumPy")
    print(f"  Pre-filter: {prefilter_label}")
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

    prefilter = not args.no_prefilter
    print(f"[3/3] Running vectorized all-vs-all alignment...")
    print(f"  Parameters: max_evalue={format_evalue(args.max_evalue)}, "
          f"prefilter={'RA diagonal + BLOSUM62 rescue' if prefilter else 'off (exhaustive)'}")
    stats = run_vectorized_pipeline(
        sequences=stripped_seqs,
        store=store,
        max_evalue=args.max_evalue,
        resume=not args.no_resume,
        prefilter=prefilter,
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
    from xtope.store       import ResultsStore
    from xtope.tag_stripper import strip_tag, set_tag
    from xtope.kmer_filter  import KmerIndex
    from xtope.aligner      import batch_align
    from xtope.evalue       import format_evalue, evalue_significance

    db_path = Path(args.db)
    if not db_path.exists():
        print(f"Error: Database not found at {args.db}")
        print("Run `python -m xtope run` first.")
        sys.exit(1)

    store = ResultsStore(args.db)

    # --- Query by ID (precomputed) ---
    if args.id:
        print(f"\nQuerying precomputed results for: {args.id}")
        results = store.query_similar(args.id, max_evalue=args.max_evalue, top_n=args.top_n)
        if not results:
            print("  No similar antigens found below E-value threshold.")
        else:
            print(f"\n  {'Rank':<5} {'Partner ID':<25} {'E-value':<12} "
                  f"{'Bit Score':<10} {'Significance':<16} {'Aligned':<9} "
                  f"{'PC Score':<10} {'Source'}")
            print(f"  {'-'*5} {'-'*25} {'-'*12} {'-'*10} {'-'*16} {'-'*9} "
                  f"{'-'*10} {'-'*14}")
            for i, r in enumerate(results, 1):
                ev = r['evalue']
                sig = evalue_significance(ev)
                pc  = r['pc_composite_score']
                src = r['detection_source'] or ""
                pc_str  = f"{pc:.3f}" if pc is not None else "  —   "
                print(f"  {i:<5} {r['partner_id']:<25} {format_evalue(ev):<12} "
                      f"{r['bit_score']:<10.1f} {sig:<16} "
                      f"{r['aligned_region_len']:<9} {pc_str:<10} {src}")

            # Show physicochemical detail for top hit if scores are present
            top = results[0]
            if top.get('pc_composite_score') is not None:
                print(f"\n  Physicochemical detail (top hit: {top['partner_id']}):")
                print(f"    Hydrophobicity corr (aligned):  {top['hydrophobicity_corr']:.3f}")
                print(f"    Charge corr (aligned):          {top['charge_corr']:.3f}")
                print(f"    pI difference:                  {top['pi_diff']:.2f}")
                print(f"    Hydrophobicity cross-corr:      {top['hydrophobicity_cross_corr']:.3f}")
                print(f"    Charge cross-corr:              {top['charge_cross_corr']:.3f}")
                print(f"    Binary H/P cross-corr:          {top['binary_hydro_cross_corr']:.3f}")
                print(f"    PC composite score:             {top['pc_composite_score']:.3f}  (≥0.55 = physicochemical hit)")
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
        db_size = sum(len(seq) for _, seq in all_seqs)
        print(f"Running live alignment against {len(all_seqs):,} sequences "
              f"({db_size:,} residues)...\n")

        from xtope.aligner import batch_align
        results = batch_align(
            query_seq=stripped,
            query_id="<query>",
            candidates=all_seqs,
            max_evalue=args.max_evalue,
            min_aligned_len=args.min_aligned,
            db_size=db_size,
        )

        results = results[:args.top_n]
        if not results:
            print("  No similar antigens found below E-value threshold.")
        else:
            print(f"  {'Rank':<5} {'Target ID':<25} {'E-value':<12} "
                  f"{'Bit Score':<11} {'Significance':<16} {'Aligned':<9} {'T Len'}")
            print(f"  {'-'*5} {'-'*25} {'-'*12} {'-'*11} {'-'*16} {'-'*9} {'-'*6}")
            for i, r in enumerate(results, 1):
                sig = evalue_significance(r.evalue)
                print(f"  {i:<5} {r.target_id:<25} {format_evalue(r.evalue):<12} "
                      f"{r.bit_score:<11.1f} {sig:<16} "
                      f"{r.aligned_region_len:<9} {r.target_length}")

        store.close()
        return

    print("Error: Provide either --id or --seq for query mode.")
    sys.exit(1)


def cmd_export(args):
    from xtope.store import ResultsStore
    from xtope.evalue import format_evalue

    db_path = Path(args.db)
    if not db_path.exists():
        print(f"Error: Database not found at {args.db}")
        sys.exit(1)

    store = ResultsStore(args.db)
    output = args.output or str(db_path.with_suffix(".csv"))
    n = store.export_csv(output, max_evalue=args.max_evalue)
    print(f"Exported {n:,} similarity pairs to: {output}")
    print(f"  E-value cutoff: {format_evalue(args.max_evalue)}")
    store.close()


def cmd_stats(args):
    from xtope.store import ResultsStore
    db_path = Path(args.db)
    if not db_path.exists():
        print(f"Error: Database not found at {args.db}")
        sys.exit(1)

    store = ResultsStore(args.db)
    print(f"\n  Database: {args.db}")
    print(f"  Antigens indexed:  {store.antigen_count():,}")
    print(f"  Total residues:    {store.total_residues():,}")
    print(f"  Similarity pairs:  {store.similarity_count():,}")
    print(f"  Last run:          {store.get_meta('last_run') or 'unknown'}")
    db_size = store.get_meta('db_size_residues')
    if db_size:
        print(f"  DB size (residues): {int(db_size):,}")
    store.close()


def cmd_help_scores(args):
    """Display E-value and bit-score interpretation guide."""
    print("""
  E-value & Bit-score Interpretation Guide
  =========================================

  E-value (Expected value)
  ------------------------
  The E-value is the number of alignments with this score (or better)
  that you would expect to see by chance in a database of this size.
  Lower E-values indicate more significant hits.

    E-value         Significance     Interpretation
    -------         ------------     --------------
    < 1e-10         Very high        Almost certainly share an epitope
    1e-10 - 1e-4    Strong           Likely cross-reactive
    1e-4  - 0.1     Moderate         Worth investigating
    0.1   - 1.0     Weak             Possibly coincidental
    > 1.0           Not significant  Probably random

  Bit-score
  ---------
  The bit-score normalises the raw SW alignment score into a common
  scale (bits of information) independent of sequence length, database
  size, and scoring system. Higher bit-scores indicate stronger matches.

    Bit-score       Interpretation
    ---------       --------------
    > 50            Very strong similarity
    30 - 50         Strong similarity
    20 - 30         Moderate similarity
    < 20            Weak or no similarity

  Unlike the old normalised score (raw_score / min_length), both E-value
  and bit-score account for sequence length and database composition,
  eliminating the length bias that inflated scores for short sequences.

  The E-value depends on database size: the same raw score produces a
  larger (worse) E-value in a bigger database, because there are more
  opportunities for chance matches. Bit-scores are database-independent.

  Karlin-Altschul parameters used:
    Matrix:     BLOSUM62
    Gap open:   10  (first gap position)
    Gap extend:  1  (each additional position)
""")


def main():
    parser = argparse.ArgumentParser(
        prog="xtope",
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
    p_run.add_argument("--max-evalue",    type=float, default=0.01, help="Max E-value to record a pair (default: 0.01)")
    p_run.add_argument("--min-aligned",   type=int,   default=8,    help="Min aligned region length in aa (default: 8)")
    p_run.add_argument("--matrix",        default="blosum62", choices=["blosum62","blosum45","blosum80"], help="Substitution matrix")
    p_run.add_argument("--no-resume",     action="store_true", help="Start fresh even if DB has partial results")
    p_run.add_argument("--backend",       default="vectorized", choices=["vectorized", "kmer"],
                       help="'vectorized' (default): NumPy batched SW, sensitive, catches contiguous epitopes. "
                            "'kmer': fast k-mer pre-filter + per-pair SW, good for quick screens on large datasets.")
    p_run.add_argument("--no-prefilter",  action="store_true",
                       help="(vectorized only) Disable the RA diagonal pre-filter and score every pair exhaustively. "
                            "More sensitive to gapped alignments; slower on large datasets.")
    p_run.set_defaults(func=cmd_run)

    # ---- query ----
    p_q = sub.add_parser("query", help="Query one antigen against precomputed results")
    p_q.add_argument("--db",          required=True, help="Path to SQLite database")
    p_q_excl = p_q.add_mutually_exclusive_group(required=True)
    p_q_excl.add_argument("--id",     help="Antigen ID to look up precomputed results")
    p_q_excl.add_argument("--seq",    help="Raw sequence to run live alignment")
    p_q.add_argument("--tag",         default=None, help="Tag to strip from --seq (if using live mode)")
    p_q.add_argument("--max-evalue",  type=float, default=0.01, help="Max E-value threshold (default: 0.01)")
    p_q.add_argument("--min-aligned", type=int,   default=8,    help="Min aligned region length (default: 8)")
    p_q.add_argument("--top-n",       type=int,   default=25,   help="Number of results to show (default: 25)")
    p_q.set_defaults(func=cmd_query)

    # ---- export ----
    p_e = sub.add_parser("export", help="Export results to CSV")
    p_e.add_argument("--db",         required=True, help="Path to SQLite database")
    p_e.add_argument("--output",     default=None,  help="Output CSV path (default: same name as db)")
    p_e.add_argument("--max-evalue", type=float, default=10.0, help="Max E-value to include (default: 10.0)")
    p_e.set_defaults(func=cmd_export)

    # ---- stats ----
    p_s = sub.add_parser("stats", help="Show database statistics")
    p_s.add_argument("--db", required=True, help="Path to SQLite database")
    p_s.set_defaults(func=cmd_stats)

    # ---- help-scores ----
    p_hs = sub.add_parser("help-scores", help="Show E-value and bit-score interpretation guide")
    p_hs.set_defaults(func=cmd_help_scores)

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
