"""
Smith-Waterman alignment scoring module.
Uses a pure-Python implementation for per-pair alignment and traceback.
For bulk scoring, use the batched NumPy/CuPy approach in the pipeline.
"""

from dataclasses import dataclass

from .sw_fallback import sw_score as _sw_fallback
from .evalue import (
    get_karlin_altschul_params,
    compute_bit_score,
    compute_evalue,
)

GAP_OPEN   = 10
GAP_EXTEND = 1

DEFAULT_MATRIX = "blosum62"


def _align_one(query: str, target: str, matrix_name: str = "blosum62") -> tuple:
    """Returns (score, approx_end_pos)."""
    score, end_j = _sw_fallback(query.upper(), target.upper(), GAP_OPEN, GAP_EXTEND)
    return score, end_j


@dataclass
class AlignmentResult:
    query_id:           str
    target_id:          str
    raw_score:          int
    bit_score:          float
    evalue:             float
    query_length:       int
    target_length:      int
    aligned_region_len: int
    query_seq:          str
    target_seq:         str

    def to_dict(self) -> dict:
        return {
            "query_id":           self.query_id,
            "target_id":          self.target_id,
            "raw_score":          self.raw_score,
            "bit_score":          round(self.bit_score, 1),
            "evalue":             self.evalue,
            "query_length":       self.query_length,
            "target_length":      self.target_length,
            "aligned_region_len": self.aligned_region_len,
            "query_seq":          self.query_seq,
            "target_seq":         self.target_seq,
        }


def align_pair(
    query_seq:   str,
    target_seq:  str,
    query_id:    str = "query",
    target_id:   str = "target",
    matrix_name: str = DEFAULT_MATRIX,
    db_size:     int = 1_000_000,
) -> AlignmentResult:
    """Align two sequences and compute E-value and bit-score.

    Args:
        query_seq: Query amino acid sequence.
        target_seq: Target amino acid sequence.
        query_id: Identifier for the query.
        target_id: Identifier for the target.
        matrix_name: Substitution matrix name.
        db_size: Total database size in residues (for E-value computation).

    Returns:
        AlignmentResult with raw score, bit-score, and E-value.
    """
    score, aligned_len = _align_one(query_seq, target_seq, matrix_name)
    lambda_, K, H = get_karlin_altschul_params(GAP_OPEN, GAP_EXTEND, matrix_name)

    bit = compute_bit_score(score, lambda_, K)
    evalue = compute_evalue(score, len(query_seq), db_size, lambda_, K, H)

    return AlignmentResult(
        query_id=query_id, target_id=target_id,
        raw_score=score, bit_score=bit, evalue=evalue,
        query_length=len(query_seq), target_length=len(target_seq),
        aligned_region_len=aligned_len,
        query_seq=query_seq, target_seq=target_seq,
    )


def batch_align(
    query_seq:       str,
    query_id:        str,
    candidates:      list,
    matrix_name:     str   = DEFAULT_MATRIX,
    max_evalue:      float = 0.01,
    min_aligned_len: int   = 8,
    db_size:         int   = 1_000_000,
) -> list:
    """Align a query against multiple candidates, filtering by E-value.

    Args:
        query_seq: Query amino acid sequence.
        query_id: Identifier for the query.
        candidates: List of (target_id, target_seq) tuples.
        matrix_name: Substitution matrix name.
        max_evalue: Maximum E-value to keep a result (lower = stricter).
        min_aligned_len: Minimum aligned region length in residues.
        db_size: Total database size in residues (for E-value computation).

    Returns:
        List of AlignmentResult, sorted by E-value ascending (best first).
    """
    lambda_, K, H = get_karlin_altschul_params(GAP_OPEN, GAP_EXTEND, matrix_name)

    results = []
    for target_id, target_seq in candidates:
        score, aligned_len = _align_one(query_seq, target_seq, matrix_name)

        evalue = compute_evalue(score, len(query_seq), db_size, lambda_, K, H)
        bit = compute_bit_score(score, lambda_, K)

        if evalue <= max_evalue and aligned_len >= min_aligned_len:
            results.append(AlignmentResult(
                query_id=query_id, target_id=target_id,
                raw_score=score, bit_score=bit, evalue=evalue,
                query_length=len(query_seq), target_length=len(target_seq),
                aligned_region_len=aligned_len,
                query_seq=query_seq, target_seq=target_seq,
            ))
    results.sort(key=lambda x: x.evalue)
    return results
