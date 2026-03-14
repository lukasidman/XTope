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
    db: EncodedDB,
    gap_open: int = 10,
    gap_extend: int = 1,
) -> np.ndarray:
    """SW-align one query against all N database sequences simultaneously.

    Uses affine gap penalties with full vectorization. The E matrix (gap in
    target / horizontal gap) has a left-to-right dependency within each row
    that prevents pure element-wise computation. We handle this with a
    column-wise loop over max_len positions — each step is a vectorized
    operation across all N sequences, so the effective Python overhead is
    O(m × max_len) cheap NumPy calls, not O(m × max_len × N).

    For 80k sequences of length ~150, this means ~150 × 150 = 22,500 NumPy
    operations per query — each operating on 80k elements in C. Total wall
    time per query is ~0.5-1.0 sec on a modern CPU.

    Args:
        query: int8 array of encoded query residues (length m).
        db: EncodedDB containing all target sequences.
        gap_open: Gap opening penalty (positive value, applied as negative).
        gap_extend: Gap extension penalty (positive value, applied as negative).

    Returns:
        (N,) float32 array of raw SW scores, one per database sequence.
    """
    n = len(db.ids)
    max_len = db.max_len
    m = len(query)

    if n == 0 or m == 0:
        return np.zeros(n, dtype=np.float32)

    # DP state — only the previous row is kept
    H_prev = np.zeros((n, max_len), dtype=np.float32)
    F = np.full((n, max_len), -1e6, dtype=np.float32)
    best = np.zeros(n, dtype=np.float32)

    go = np.float32(gap_open)
    ge = np.float32(gap_extend)

    for i in range(m):
        # Substitution scores: BLOSUM62[query[i], target[j]] for all (n, max_len)
        sub = BLOSUM62[query[i], db.matrix].astype(np.float32)

        # Diagonal: H[i-1, j-1] — shift H_prev right by 1
        diag = np.empty_like(H_prev)
        diag[:, 0] = 0.0
        diag[:, 1:] = H_prev[:, :-1]

        # Match/mismatch candidate
        match = diag + sub

        # F[i,j] = max(F[i-1,j] - gap_extend, H[i-1,j] - gap_open)
        # Gap in query (vertical gap) — fully parallel across j
        F = np.maximum(F - ge, H_prev - go)

        # H_curr = max(0, match, F) — before horizontal gap correction
        H_curr = np.maximum(np.float32(0.0), np.maximum(match, F))

        # E[i,j] = max(E[i,j-1] - gap_extend, H[i,j-1] - gap_open)
        # Gap in target (horizontal gap) — left-to-right dependency.
        # Loop over columns; each iteration is vectorized across all N sequences.
        E_col = np.full(n, -1e6, dtype=np.float32)
        for j in range(1, max_len):
            E_col = np.maximum(E_col - ge, H_curr[:, j - 1] - go)
            H_curr[:, j] = np.maximum(H_curr[:, j], E_col)

        # Zero out padded positions
        H_curr *= db.mask

        # Track running max per sequence
        best = np.maximum(best, H_curr.max(axis=1))

        H_prev = H_curr

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
    ids = list(sequences.keys())
    n = len(ids)

    print(f"  Encoding {n:,} sequences...")
    db = encode_database(sequences)

    # Build ID → row index mapping for upper-triangle slicing
    id_to_idx = {sid: i for i, sid in enumerate(ids)}

    completed = store.get_completed_queries() if resume else set()
    if completed:
        print(f"  Resuming — {len(completed):,} queries already done")

    total_pairs_found = 0
    batch_buffer: list[dict] = []
    done = 0

    for qi, query_id in enumerate(ids):
        if query_id in completed:
            done += 1
            continue

        query_encoded = db.matrix[qi, :db.lengths[qi]]

        # Score this query against ALL sequences
        raw_scores = score_one_vs_all(query_encoded, db, gap_open, gap_extend)

        # Only keep upper triangle (qi < ti) to avoid duplicate pairs
        query_len = int(db.lengths[qi])
        for ti in range(qi + 1, n):
            target_len = int(db.lengths[ti])
            min_len = min(query_len, target_len)
            if min_len == 0:
                continue

            raw = float(raw_scores[ti])
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

        # Flush to DB periodically
        if len(batch_buffer) >= batch_size:
            store.insert_similarities_batch(batch_buffer)
            batch_buffer.clear()

        # Progress reporting
        if done % 10 == 0 or done == n:
            elapsed = time.time() - start_time
            rate = done / elapsed if elapsed > 0 else 1
            remaining = (n - done) / rate if rate > 0 else 0
            eta = _format_eta(remaining)
            pct = 100 * done / n

            print(
                f"  {done:>7,}/{n:,}  ({pct:5.1f}%)  "
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

    return {
        "total_pairs": total_pairs_found,
        "elapsed_seconds": elapsed_total,
    }


def _format_eta(seconds: float) -> str:
    if seconds < 60:
        return f"{seconds:.0f}s"
    elif seconds < 3600:
        return f"{seconds/60:.1f}min"
    else:
        return f"{seconds/3600:.1f}hr"
