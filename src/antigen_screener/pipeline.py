"""
Main pipeline orchestrator.
Runs the full all-vs-all antigen similarity screening pipeline:
  1. Load CSV/TSV database
  2. Strip His6-ABP tags from all sequences
  3. Build k-mer index for fast pre-filtering
  4. For each antigen, find candidates via k-mer filter, then score with SW
  5. Persist results to SQLite for future querying
"""

import time
import datetime
from pathlib import Path
from concurrent.futures import ProcessPoolExecutor, as_completed
from typing import Callable

from .db_loader   import load_csv
from .tag_stripper import strip_tag, set_tag
from .kmer_filter  import KmerIndex
from .aligner      import batch_align, AlignmentResult
from .store        import ResultsStore


def _format_eta(seconds: float) -> str:
    if seconds < 60:
        return f"{seconds:.0f}s"
    elif seconds < 3600:
        return f"{seconds/60:.1f}min"
    else:
        return f"{seconds/3600:.1f}hr"


def run_pipeline(
    input_path:       str,
    db_path:          str,
    tag_sequence:     str | None = None,
    id_col:           str | None = None,
    seq_col:          str | None = None,
    kmer_threshold:   float = 0.04,
    min_norm_score:   float = 1.0,
    min_aligned_len:  int   = 8,
    matrix:           str   = "blosum62",
    resume:           bool  = True,
    progress_cb:      Callable | None = None,
    workers:          int   = 1,
):
    """
    Full all-vs-all screening pipeline.

    Args:
        input_path:      Path to CSV/TSV antigen database
        db_path:         Path to output SQLite database
        tag_sequence:    His6-ABP tag to strip (uses default if None)
        id_col:          ID column name (auto-detected if None)
        seq_col:         Sequence column name (auto-detected if None)
        kmer_threshold:  Min Jaccard similarity for pre-filter (lower = more sensitive)
        min_norm_score:  Min normalised SW score to record a pair
        min_aligned_len: Min aligned region length (epitope-sized = 8aa)
        matrix:          Substitution matrix: 'blosum62', 'blosum45', 'blosum80'
        resume:          If True, skip queries already in the DB
        progress_cb:     Optional callback(done, total, eta_str) for progress
        workers:         Number of parallel worker processes (experimental)
    """
    start_time = time.time()
    print(f"\n{'='*60}")
    print(f"  Antigen Cross-Reactivity Screener")
    print(f"  Started: {datetime.datetime.now().strftime('%Y-%m-%d %H:%M:%S')}")
    print(f"{'='*60}\n")

    # --- Update tag if provided ---
    if tag_sequence:
        set_tag(tag_sequence)
        print(f"  Tag set to: {tag_sequence}\n")

    # --- Load sequences ---
    print(f"[1/4] Loading sequences from: {input_path}")
    raw_records = load_csv(input_path, id_col=id_col, seq_col=seq_col, verbose=True)
    total = len(raw_records)
    print(f"  Total antigens: {total:,}\n")

    # --- Strip tags and store antigens ---
    print(f"[2/4] Stripping His6-ABP tags...")
    store = ResultsStore(db_path)

    stripped_records = []
    tag_found_count = 0
    antigen_batch = []

    for antigen_id, sequence in raw_records:
        result = strip_tag(sequence)
        stripped = result["stripped"]
        found    = result["tag_found"]
        if found:
            tag_found_count += 1
        antigen_batch.append((antigen_id, sequence, stripped, found))
        stripped_records.append((antigen_id, stripped))

    store.upsert_antigens_batch(antigen_batch)
    print(f"  Tags detected and removed in {tag_found_count:,} / {total:,} sequences\n")

    # --- Build k-mer index ---
    print(f"[3/4] Building k-mer index...")
    index = KmerIndex()
    index.add_batch(stripped_records)
    print(f"  Index built ({len(index):,} sequences)\n")

    # --- All-vs-all SW alignment ---
    print(f"[4/4] Running all-vs-all alignment...")
    print(f"  Parameters: kmer_threshold={kmer_threshold}, "
          f"min_norm_score={min_norm_score}, min_aligned_len={min_aligned_len}, "
          f"matrix={matrix}")

    completed_queries = store.get_completed_queries() if resume else set()
    if completed_queries:
        print(f"  Resuming — {len(completed_queries):,} queries already done\n")

    done = 0
    total_pairs_found = 0
    batch_buffer = []
    BATCH_SIZE = 500   # flush to DB every N results

    for antigen_id, stripped_seq in stripped_records:
        if antigen_id in completed_queries:
            done += 1
            continue

        # Pre-filter with k-mer index
        kmer_candidates = index.query(stripped_seq, threshold=kmer_threshold, exclude_id=antigen_id)

        if kmer_candidates:
            # Resolve sequences for candidates
            candidate_seqs = [
                (cid, index.sequences[cid])
                for cid, _ in kmer_candidates
                if cid in index.sequences
            ]

            # SW alignment scoring
            alignments = batch_align(
                query_seq=stripped_seq,
                query_id=antigen_id,
                candidates=candidate_seqs,
                matrix_name=matrix,
                min_norm_score=min_norm_score,
                min_aligned_len=min_aligned_len,
            )

            for aln in alignments:
                batch_buffer.append(aln.to_dict())
                total_pairs_found += 1

        done += 1

        # Flush to DB periodically
        if len(batch_buffer) >= BATCH_SIZE:
            store.insert_similarities_batch(batch_buffer)
            batch_buffer.clear()

        # Progress reporting
        if done % 100 == 0 or done == total:
            elapsed  = time.time() - start_time
            rate     = done / elapsed if elapsed > 0 else 1
            remaining = (total - done) / rate if rate > 0 else 0
            eta_str  = _format_eta(remaining)
            pct      = 100 * done / total

            print(
                f"  {done:>7,}/{total:,}  ({pct:5.1f}%)  "
                f"pairs found: {total_pairs_found:,}  "
                f"ETA: {eta_str}",
                end="\r",
                flush=True,
            )

            if progress_cb:
                progress_cb(done, total, eta_str)

    # Final flush
    if batch_buffer:
        store.insert_similarities_batch(batch_buffer)

    elapsed_total = time.time() - start_time
    print(f"\n\n{'='*60}")
    print(f"  Run complete!")
    print(f"  Total time:        {_format_eta(elapsed_total)}")
    print(f"  Antigens screened: {total:,}")
    print(f"  Similar pairs found: {total_pairs_found:,}")
    print(f"  Results saved to:  {db_path}")
    print(f"{'='*60}\n")

    store.set_meta("last_run", datetime.datetime.now().isoformat())
    store.set_meta("total_antigens", str(total))
    store.set_meta("total_pairs", str(total_pairs_found))
    store.close()

    return {
        "total_antigens":   total,
        "total_pairs":      total_pairs_found,
        "elapsed_seconds":  elapsed_total,
        "db_path":          db_path,
    }
