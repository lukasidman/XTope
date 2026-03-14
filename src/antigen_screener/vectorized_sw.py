"""
Vectorized Smith-Waterman using NumPy.

Aligns one query against all database sequences simultaneously using
element-wise array operations, replacing the sequential k-mer filter +
per-pair SW pipeline.

Memory budget at 80k sequences (max_len=150):
    db_matrix:  12 MB   (80000, 150) int8
    db_mask:    12 MB   (80000, 150) bool
    DP state:  144 MB   3 × (80000, 150) float32
    Total:    ~170 MB
"""

import time
import datetime
from dataclasses import dataclass
from typing import Callable

import numpy as np

from .sw_fallback import _BLOSUM62_RAW


# ------------------------------------------------------------------ #
#  Amino acid encoding + BLOSUM62 matrix                              #
# ------------------------------------------------------------------ #

_AA_ORDER = "ARNDCQEGHILKMFPSTWYV"
AA_TO_INT: dict[str, int] = {aa: i for i, aa in enumerate(_AA_ORDER)}
_UNKNOWN = len(_AA_ORDER)  # index 20 — for any non-standard residue


def _build_blosum62_array() -> np.ndarray:
    """Parse the raw BLOSUM62 text into a (21, 21) int8 array.

    Row/column 20 is the unknown residue, scoring 0 against everything.
    """
    lines = [line for line in _BLOSUM62_RAW.strip().splitlines() if line.strip()]
    headers = lines[0].split()
    n = len(headers)  # 20

    mat = np.zeros((n + 1, n + 1), dtype=np.int8)
    for row_idx, line in enumerate(lines[1:]):
        parts = line.split()
        for col_idx in range(n):
            mat[row_idx, col_idx] = int(parts[col_idx + 1])

    # Row/col 20 stays zero (unknown residue)
    return mat


BLOSUM62 = _build_blosum62_array()


def encode_sequence(seq: str) -> np.ndarray:
    """Encode an amino acid string as an int8 array of indices into AA_TO_INT."""
    return np.array([AA_TO_INT.get(aa, _UNKNOWN) for aa in seq.upper()], dtype=np.int8)


# ------------------------------------------------------------------ #
#  Encoded database                                                   #
# ------------------------------------------------------------------ #

@dataclass
class EncodedDB:
    """Padded integer-encoded database ready for vectorized SW."""
    matrix: np.ndarray    # (N, max_len) int8
    mask: np.ndarray      # (N, max_len) bool — True where real residue
    ids: list[str]        # row order
    lengths: np.ndarray   # (N,) int32 — real length per sequence
    max_len: int


def encode_database(sequences: dict[str, str]) -> EncodedDB:
    """Convert {id: sequence} to a padded integer matrix.

    Args:
        sequences: Mapping of antigen ID to stripped amino acid sequence.

    Returns:
        EncodedDB with all sequences encoded and zero-padded to uniform length.
    """
    ids = list(sequences.keys())
    n = len(ids)
    if n == 0:
        return EncodedDB(
            matrix=np.empty((0, 0), dtype=np.int8),
            mask=np.empty((0, 0), dtype=bool),
            ids=[],
            lengths=np.empty(0, dtype=np.int32),
            max_len=0,
        )

    seqs = [sequences[sid] for sid in ids]
    lengths = np.array([len(s) for s in seqs], dtype=np.int32)
    max_len = int(lengths.max())

    matrix = np.zeros((n, max_len), dtype=np.int8)
    mask = np.zeros((n, max_len), dtype=bool)

    for i, seq in enumerate(seqs):
        encoded = encode_sequence(seq)
        slen = len(encoded)
        matrix[i, :slen] = encoded
        mask[i, :slen] = True

    return EncodedDB(matrix=matrix, mask=mask, ids=ids, lengths=lengths, max_len=max_len)


# ------------------------------------------------------------------ #
#  Core vectorized SW                                                 #
# ------------------------------------------------------------------ #

def score_one_vs_all(
    query: np.ndarray,
    targets_matrix: np.ndarray,
    targets_mask: np.ndarray,
    gap_open: int = 10,
    gap_extend: int = 1,
) -> np.ndarray:
    """SW-align one query against a set of target sequences simultaneously.

    Uses affine gap penalties with full vectorization. The E matrix (gap in
    target / horizontal gap) has a left-to-right dependency within each row
    that prevents pure element-wise computation. We handle this with a
    column-wise loop over max_len positions — each step is a vectorized
    operation across all N sequences, so the effective Python overhead is
    O(m × max_len) cheap NumPy calls, not O(m × max_len × N).

    The caller controls which sequences are passed in. For upper-triangle
    computation, only sequences with index > query_index need to be included,
    so the arrays shrink as the pipeline progresses — nearly halving total work.

    Args:
        query: int8 array of encoded query residues (length m).
        targets_matrix: (N, max_len) int8 array of encoded target sequences.
        targets_mask: (N, max_len) bool array — True where real residue.
        gap_open: Gap opening penalty (positive value, applied as negative).
        gap_extend: Gap extension penalty (positive value, applied as negative).

    Returns:
        (N,) float32 array of raw SW scores, one per target sequence.
    """
    n, max_len = targets_matrix.shape
    m = len(query)

    if n == 0 or m == 0:
        return np.zeros(n, dtype=np.float32)

    # DP state — only the previous row is kept.
    # Pre-allocate all working arrays once to avoid per-iteration allocation.
    H_prev = np.zeros((n, max_len), dtype=np.float32)
    H_curr = np.zeros((n, max_len), dtype=np.float32)
    F = np.full((n, max_len), -1e6, dtype=np.float32)
    diag = np.empty((n, max_len), dtype=np.float32)
    sub = np.empty((n, max_len), dtype=np.float32)
    E_col = np.empty(n, dtype=np.float32)
    tmp_col = np.empty(n, dtype=np.float32)  # scratch for E_col loop
    best = np.zeros(n, dtype=np.float32)

    go = np.float32(gap_open)
    ge = np.float32(gap_extend)

    for i in range(m):
        # Substitution scores: BLOSUM62[query[i], target[j]] for all (n, max_len)
        sub[:] = BLOSUM62[query[i], targets_matrix]

        # Diagonal: H[i-1, j-1] — shift H_prev right by 1
        diag[:, 0] = 0.0
        diag[:, 1:] = H_prev[:, :-1]

        # Match/mismatch candidate: diag + sub → H_curr (reuse as scratch)
        np.add(diag, sub, out=H_curr)

        # F[i,j] = max(F[i-1,j] - gap_extend, H[i-1,j] - gap_open)
        # Gap in query (vertical gap) — fully parallel across j
        np.subtract(F, ge, out=F)
        np.subtract(H_prev, go, out=diag)  # reuse diag as temp
        np.maximum(F, diag, out=F)

        # H_curr = max(0, match, F) — before horizontal gap correction
        np.maximum(H_curr, F, out=H_curr)
        np.maximum(H_curr, np.float32(0.0), out=H_curr)

        # E[i,j] = max(E[i,j-1] - gap_extend, H[i,j-1] - gap_open)
        # Gap in target (horizontal gap) — left-to-right dependency.
        # Loop over columns; each iteration is vectorized across all N sequences.
        E_col[:] = -1e6
        for j in range(1, max_len):
            np.subtract(E_col, ge, out=E_col)
            np.subtract(H_curr[:, j - 1], go, out=tmp_col)
            np.maximum(E_col, tmp_col, out=E_col)
            np.maximum(H_curr[:, j], E_col, out=H_curr[:, j])

        # Zero out padded positions
        H_curr *= targets_mask

        # Track running max per sequence
        best[:] = np.maximum(best, H_curr.max(axis=1))

        # Swap buffers — H_curr becomes H_prev for next iteration
        H_prev, H_curr = H_curr, H_prev

    return best


# ------------------------------------------------------------------ #
#  Pipeline entry point                                               #
# ------------------------------------------------------------------ #

def run_vectorized_pipeline(
    sequences: dict[str, str],
    store,
    gap_open: int = 10,
    gap_extend: int = 1,
    min_norm_score: float = 1.0,
    resume: bool = True,
    progress_cb: Callable | None = None,
    batch_size: int = 500,
) -> dict:
    """Score all pairs using vectorized SW alignment.

    Replaces the k-mer filter + per-pair SW pipeline with brute-force
    vectorized alignment. Every pair gets a full SW score — no pre-filter,
    no missed pairs.

    Args:
        sequences: {antigen_id: stripped_sequence} for all antigens.
        store: ResultsStore instance to write results into.
        gap_open: Gap opening penalty.
        gap_extend: Gap extension penalty.
        min_norm_score: Minimum normalised score to store a pair.
        resume: If True, skip already-completed query IDs.
        progress_cb: Optional callback(done, total, eta_str).
        batch_size: Flush results to DB every N results.

    Returns:
        Dict with run statistics.
    """
    start_time = time.time()

    # Sort sequences longest-to-shortest. This enables a key optimisation:
    # since we process the upper triangle (query_idx < target_idx), all
    # remaining targets are always shorter or equal to the current query.
    # We can shrink the array width (max_len) after each query, reducing
    # both the column loop iterations and array operation sizes.
    sorted_items = sorted(sequences.items(), key=lambda x: len(x[1]), reverse=True)
    sorted_seqs = {sid: seq for sid, seq in sorted_items}
    ids = list(sorted_seqs.keys())
    n = len(ids)

    print(f"  Encoding {n:,} sequences (sorted longest → shortest)...")
    db = encode_database(sorted_seqs)

    # Populate the antigens table so export_csv can JOIN on it.
    # The vectorized pipeline only receives stripped sequences, so we store
    # the stripped sequence as both 'sequence' and 'stripped' fields.
    antigen_records = [
        (sid, seq, seq, True) for sid, seq in sorted_seqs.items()
    ]
    store.upsert_antigens_batch(antigen_records)

    completed = store.get_completed_queries() if resume else set()
    if completed:
        print(f"  Resuming — {len(completed):,} queries already done")

    # Upper-triangle total: sum of (n-1) + (n-2) + ... + 1 = n*(n-1)/2
    # Each query qi compares against (n - qi - 1) targets.
    # Track progress by comparisons done, not just queries done, so the
    # ETA accounts for later queries being progressively cheaper.
    total_comparisons = n * (n - 1) // 2
    comparisons_done = 0

    total_pairs_found = 0
    batch_buffer: list[dict] = []
    done = 0

    for qi, query_id in enumerate(ids):
        targets_this_query = n - qi - 1

        if query_id in completed:
            done += 1
            comparisons_done += targets_this_query
            continue

        query_encoded = db.matrix[qi, :db.lengths[qi]]

        # Only score against sequences with index > qi (upper triangle).
        # NumPy slicing creates a view — no data is copied.
        target_start = qi + 1

        # Because sequences are sorted longest→shortest, the longest
        # remaining target is at index target_start. Its length gives
        # the effective max_len — all columns beyond that are padding
        # zeros that we can skip entirely.
        if target_start < n:
            effective_max_len = int(db.lengths[target_start])
        else:
            effective_max_len = 0

        targets_matrix = db.matrix[target_start:, :effective_max_len]
        targets_mask = db.mask[target_start:, :effective_max_len]

        raw_scores = score_one_vs_all(
            query_encoded, targets_matrix, targets_mask, gap_open, gap_extend,
        )

        # Collect results above threshold
        query_len = int(db.lengths[qi])
        for ri in range(len(raw_scores)):
            ti = target_start + ri  # original index in db
            target_len = int(db.lengths[ti])
            min_len = min(query_len, target_len)
            if min_len == 0:
                continue

            raw = float(raw_scores[ri])
            norm = raw / min_len

            if norm >= min_norm_score:
                batch_buffer.append({
                    "query_id": query_id,
                    "target_id": ids[ti],
                    "raw_score": int(raw),
                    "normalized_score": round(norm, 4),
                    "query_length": query_len,
                    "target_length": target_len,
                    "aligned_region_len": min_len,  # approximation — no traceback
                })
                total_pairs_found += 1

        done += 1
        comparisons_done += targets_this_query

        # Flush to DB periodically
        if len(batch_buffer) >= batch_size:
            store.insert_similarities_batch(batch_buffer)
            batch_buffer.clear()

        # Progress reporting — based on comparisons, not query count,
        # so ETA reflects that later queries are progressively cheaper.
        if done % 10 == 0 or done == n:
            elapsed = time.time() - start_time
            if comparisons_done > 0:
                time_per_comp = elapsed / comparisons_done
                remaining_comps = total_comparisons - comparisons_done
                remaining_secs = remaining_comps * time_per_comp
            else:
                remaining_secs = 0
            eta = _format_eta(remaining_secs)
            pct = 100 * comparisons_done / total_comparisons if total_comparisons > 0 else 100

            print(
                f"  {done:>7,}/{n:,} queries  ({pct:5.1f}%)  "
                f"pairs found: {total_pairs_found:,}  "
                f"ETA: {eta}",
                end="\r",
                flush=True,
            )

            if progress_cb:
                progress_cb(done, n, eta)

    # Final flush
    if batch_buffer:
        store.insert_similarities_batch(batch_buffer)

    elapsed_total = time.time() - start_time
    print()  # newline after \r progress
    print(
        f"  Complete: {done:,} queries, {comparisons_done:,} comparisons, "
        f"{total_pairs_found:,} similar pairs found in {_format_eta(elapsed_total)}"
    )

    return {
        "total_pairs": total_pairs_found,
        "total_queries": done,
        "total_comparisons": comparisons_done,
        "elapsed_seconds": elapsed_total,
    }


def _format_eta(seconds: float) -> str:
    if seconds < 60:
        return f"{seconds:.0f}s"
    elif seconds < 3600:
        return f"{seconds/60:.1f}min"
    else:
        return f"{seconds/3600:.1f}hr"
