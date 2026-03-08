"""Physicochemical similarity scoring for aligned antigen regions.

Provides hydrophobicity (Kyte-Doolittle), charge, and isoelectric point
calculations to supplement sequence-based Smith-Waterman scores with
biophysical property comparisons.
"""

from __future__ import annotations

# Kyte-Doolittle hydrophobicity scale
KD_SCALE: dict[str, float] = {
    "A": 1.8,
    "R": -4.5,
    "N": -3.5,
    "D": -3.5,
    "C": 2.5,
    "Q": -3.5,
    "E": -3.5,
    "G": -0.4,
    "H": -3.2,
    "I": 4.5,
    "L": 3.8,
    "K": -3.9,
    "M": 1.9,
    "F": 2.8,
    "P": -1.6,
    "S": -0.8,
    "T": -0.7,
    "W": -0.9,
    "Y": -1.3,
    "V": 4.2,
}


def hydrophobicity_profile(sequence: str, window: int = 5) -> list[float]:
    """Kyte-Doolittle sliding-window hydrophobicity.

    Args:
        sequence: Amino acid sequence (single-letter codes).
        window: Sliding window size.

    Returns:
        List of average hydrophobicity values per window position.
    """
    raise NotImplementedError


def charge_profile(sequence: str, window: int = 5) -> list[float]:
    """Net charge per window (K, R = +1; D, E = -1).

    Args:
        sequence: Amino acid sequence (single-letter codes).
        window: Sliding window size.

    Returns:
        List of net charge values per window position.
    """
    raise NotImplementedError


def isoelectric_point(sequence: str) -> float:
    """Approximate isoelectric point (pI) of a peptide sequence.

    Args:
        sequence: Amino acid sequence (single-letter codes).

    Returns:
        Estimated pI value.
    """
    raise NotImplementedError


def physicochemical_similarity(seq_a: str, seq_b: str) -> dict[str, float]:
    """Compare physicochemical properties of two aligned sequences.

    Both sequences must be the same length (the aligned region).

    Args:
        seq_a: First aligned sequence.
        seq_b: Second aligned sequence.

    Returns:
        Dictionary with keys:
            hydrophobicity_correlation: Pearson correlation (-1 to 1).
            charge_correlation: Pearson correlation (-1 to 1).
            pi_difference: Absolute difference in pI values.
            composite_score: Weighted blend of the above metrics.
    """
    raise NotImplementedError
