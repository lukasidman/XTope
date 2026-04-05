"""
E-value and bit-score computation using Karlin-Altschul statistics.

Converts raw Smith-Waterman scores into statistically meaningful metrics:
- Bit-score: normalised score comparable across different scoring systems
- E-value: expected number of hits this good by chance in the search space

Reference: Karlin & Altschul, "Methods for assessing the statistical significance
of molecular sequence features", PNAS (1990).

Gap penalty convention mapping
------------------------------
Our SW implementation (sw_fallback.py) uses affine gap penalties where:
    gap of length k costs: gap_open + (k-1) * gap_extend

NCBI BLAST uses a different convention:
    gap of length k costs: existence + k * extension

The mapping is: existence = gap_open - gap_extend, extension = gap_extend.
So our default (gap_open=10, gap_extend=1) maps to NCBI (existence=9, extension=1).

The Karlin-Altschul parameters (K, lambda, H) below are indexed by our convention
(gap_open, gap_extend) and internally mapped to the correct NCBI row.
"""

from __future__ import annotations

import math

# Published Karlin-Altschul parameters for BLOSUM62.
# Source: NCBI BLAST C toolkit (blast_stat.c), indexed by NCBI convention
# (existence, extension). We map from our (gap_open, gap_extend) at lookup time.
#
# Each entry: (lambda, K, H)
#   lambda: decay rate of score distribution
#   K:      search space scaling constant
#   H:      relative entropy of the scoring system (bits per aligned position)
_BLOSUM62_PARAMS_NCBI: dict[tuple[int, int], tuple[float, float, float]] = {
    # (existence, extension): (lambda, K, H)
    (11, 1): (0.267, 0.041, 0.140),
    (10, 1): (0.262, 0.040, 0.140),
    (9,  1): (0.253, 0.038, 0.130),
    (8,  1): (0.244, 0.035, 0.120),
    (7,  1): (0.230, 0.032, 0.100),
    (11, 2): (0.199, 0.040, 0.110),
    (10, 2): (0.186, 0.035, 0.090),
    (9,  2): (0.170, 0.028, 0.070),
    (12, 1): (0.270, 0.042, 0.140),
    (13, 1): (0.272, 0.042, 0.140),
}


def get_karlin_altschul_params(
    gap_open: int = 10,
    gap_extend: int = 1,
    matrix: str = "blosum62",
) -> tuple[float, float, float]:
    """Look up Karlin-Altschul parameters for the given scoring system.

    Args:
        gap_open: Gap opening penalty (our convention: first position cost).
        gap_extend: Gap extension penalty (each additional position).
        matrix: Substitution matrix name. Only 'blosum62' is currently supported.

    Returns:
        (lambda_, K, H) tuple of Karlin-Altschul statistical parameters.

    Raises:
        ValueError: If the matrix or gap penalty combination is not supported.
    """
    if matrix.lower() != "blosum62":
        raise ValueError(
            f"Karlin-Altschul parameters not available for matrix '{matrix}'. "
            f"Only 'blosum62' is currently supported."
        )

    # Map our gap convention to NCBI convention
    ncbi_existence = gap_open - gap_extend
    ncbi_extension = gap_extend

    params = _BLOSUM62_PARAMS_NCBI.get((ncbi_existence, ncbi_extension))
    if params is None:
        raise ValueError(
            f"No Karlin-Altschul parameters for gap_open={gap_open}, "
            f"gap_extend={gap_extend} (NCBI existence={ncbi_existence}, "
            f"extension={ncbi_extension}). Available combinations: "
            + ", ".join(f"open={e+x},ext={x}" for e, x in sorted(_BLOSUM62_PARAMS_NCBI))
        )

    return params


def effective_length(seq_len: int, K: float, H: float, db_size: int) -> int:
    """Compute effective sequence length after edge-effect correction.

    The Karlin-Altschul statistics assume infinite sequences. For finite
    sequences, the effective search space is reduced by the expected length
    of a high-scoring segment pair (HSP). This correction is:

        m' = m - (ln(K * m * n) / H)

    where m is the query length and n is the database size in residues.
    The correction is small for typical sequence/database sizes.

    Args:
        seq_len: Actual sequence length.
        K: Karlin-Altschul K parameter.
        H: Relative entropy.
        db_size: Total database size in residues.

    Returns:
        Effective length, clamped to a minimum of 1.
    """
    if seq_len <= 0 or db_size <= 0 or H <= 0:
        return max(seq_len, 1)

    correction = math.log(K * seq_len * db_size) / H
    eff = seq_len - correction
    return max(int(eff), 1)


def compute_bit_score(raw_score: int, lambda_: float, K: float) -> float:
    """Convert a raw SW score to a bit-score.

    The bit-score normalises raw alignment scores into a common scale
    (bits of information) that is independent of the scoring system,
    sequence lengths, and database size. Higher is better.

        S' = (lambda * S - ln(K)) / ln(2)

    Bit-scores are directly comparable across different searches and
    scoring parameters. A bit-score of 50 always represents the same
    statistical significance regardless of the matrix or gap penalties used.

    Args:
        raw_score: Integer raw Smith-Waterman score.
        lambda_: Karlin-Altschul lambda parameter.
        K: Karlin-Altschul K parameter.

    Returns:
        Bit-score (float). Higher values indicate stronger similarity.
    """
    if raw_score <= 0:
        return 0.0
    return (lambda_ * raw_score - math.log(K)) / math.log(2)


def compute_evalue(
    raw_score: int,
    query_len: int,
    db_size: int,
    lambda_: float,
    K: float,
    H: float,
) -> float:
    """Compute E-value for a pairwise alignment in a database context.

    The E-value is the expected number of alignments with score >= S that
    would occur by chance when searching a database of the given size.
    Lower is better: E < 0.001 is very significant, E > 10 is noise.

        E = K * m' * D' * exp(-lambda * S)

    where m' and D' are the effective query length and effective database
    size after edge-effect corrections.

    Args:
        raw_score: Integer raw Smith-Waterman score.
        query_len: Length of the query sequence.
        db_size: Total database size in residues (sum of all sequence lengths).
        lambda_: Karlin-Altschul lambda parameter.
        K: Karlin-Altschul K parameter.
        H: Relative entropy (for effective length correction).

    Returns:
        E-value (float). Lower values indicate stronger significance.
    """
    if raw_score <= 0:
        return float(db_size * query_len)  # worst case

    m_eff = effective_length(query_len, K, H, db_size)
    d_eff = effective_length(db_size, K, H, query_len)

    return K * m_eff * d_eff * math.exp(-lambda_ * raw_score)


def evalue_significance(evalue: float) -> str:
    """Classify an E-value into a human-readable significance tier.

    Based on BLAST conventions, calibrated for large databases of short
    peptide sequences (50-150 aa).

    Args:
        evalue: The E-value to classify.

    Returns:
        Significance label string.
    """
    if evalue < 1e-10:
        return "very high"
    elif evalue < 1e-4:
        return "strong"
    elif evalue < 1e-1:
        return "moderate"
    elif evalue < 1.0:
        return "weak"
    else:
        return "not significant"


def format_evalue(evalue: float) -> str:
    """Format an E-value for display.

    Uses scientific notation for very small values, fixed notation for larger.

    Args:
        evalue: The E-value to format.

    Returns:
        Formatted string.
    """
    if evalue == 0.0:
        return "0.0"
    elif evalue < 1e-100:
        return f"{evalue:.1e}"
    elif evalue < 0.001:
        return f"{evalue:.2e}"
    elif evalue < 10:
        return f"{evalue:.4f}"
    elif evalue < 1000:
        return f"{evalue:.1f}"
    else:
        return f"{evalue:.1e}"
