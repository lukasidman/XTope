"""
Smith-Waterman alignment scoring module.
Attempts to use parasail (fast, SIMD-accelerated) and falls back to
the pure Python implementation if parasail is not installed.

Install parasail for ~50x speedup on large databases:
    pip install parasail
"""

from dataclasses import dataclass

try:
    import parasail as _parasail
    _HAS_PARASAIL = True
except ImportError:
    _HAS_PARASAIL = False

from .sw_fallback import sw_score as _sw_fallback

GAP_OPEN   = 10
GAP_EXTEND = 1

DEFAULT_MATRIX = "blosum62"


def _align_one(query: str, target: str, matrix_name: str = "blosum62") -> tuple:
    """Returns (score, approx_end_pos)."""
    if _HAS_PARASAIL:
        matrices = {
            "blosum62": _parasail.blosum62,
            "blosum45": _parasail.blosum45,
            "blosum80": _parasail.blosum80,
        }
        matrix = matrices.get(matrix_name, _parasail.blosum62)
        r = _parasail.sw_trace_striped_16(
            query.upper(), target.upper(), GAP_OPEN, GAP_EXTEND, matrix
        )
        return r.score, max(r.end_query, r.end_ref) + 1
    else:
        score, end_j = _sw_fallback(query.upper(), target.upper(), GAP_OPEN, GAP_EXTEND)
        return score, end_j


@dataclass
class AlignmentResult:
    query_id:           str
    target_id:          str
    raw_score:          int
    normalized_score:   float
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
            "normalized_score":   round(self.normalized_score, 4),
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
) -> AlignmentResult:
    score, aligned_len = _align_one(query_seq, target_seq, matrix_name)
    min_len = min(len(query_seq), len(target_seq))
    normalized = score / min_len if min_len > 0 else 0.0
    return AlignmentResult(
        query_id=query_id, target_id=target_id,
        raw_score=score, normalized_score=normalized,
        query_length=len(query_seq), target_length=len(target_seq),
        aligned_region_len=aligned_len,
        query_seq=query_seq, target_seq=target_seq,
    )


def batch_align(
    query_seq:       str,
    query_id:        str,
    candidates:      list,
    matrix_name:     str   = DEFAULT_MATRIX,
    min_norm_score:  float = 1.0,
    min_aligned_len: int   = 8,
) -> list:
    results = []
    for target_id, target_seq in candidates:
        score, aligned_len = _align_one(query_seq, target_seq, matrix_name)
        min_len = min(len(query_seq), len(target_seq))
        normalized = score / min_len if min_len > 0 else 0.0
        if normalized >= min_norm_score and aligned_len >= min_aligned_len:
            results.append(AlignmentResult(
                query_id=query_id, target_id=target_id,
                raw_score=score, normalized_score=normalized,
                query_length=len(query_seq), target_length=len(target_seq),
                aligned_region_len=aligned_len,
                query_seq=query_seq, target_seq=target_seq,
            ))
    results.sort(key=lambda x: x.normalized_score, reverse=True)
    return results
