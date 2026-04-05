"""Physicochemical similarity scoring for aligned antigen regions.

Provides hydrophobicity (Kyte-Doolittle), charge, and isoelectric point
calculations to supplement sequence-based Smith-Waterman scores with
biophysical property comparisons.

Two modes of use:
  - Post-alignment:  physicochemical_similarity(seq_a, seq_b) on two equal-length
                     aligned regions. Pearson correlation on the residue-level profiles.
  - Global / rescue: sliding_window_correlation(seq_a, seq_b) for pairs where
                     SW produced no useful alignment (CONV-type). Cross-correlates
                     the full hydrophobicity and charge profiles along all offsets
                     and takes the maximum, so shifted or reversed arrangements
                     are still detected.
"""

from __future__ import annotations

import math


# ---------------------------------------------------------------------------
# Physicochemical property scales
# ---------------------------------------------------------------------------

# Kyte-Doolittle hydrophobicity scale (Kyte & Doolittle, 1982)
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

# Net charge at neutral pH: K, R = +1; D, E = -1; H ~ +0.1 (protonation partial)
CHARGE_SCALE: dict[str, float] = {
    "K": 1.0,
    "R": 1.0,
    "D": -1.0,
    "E": -1.0,
    "H": 0.1,  # partially protonated at pH 7
}

# Binary hydrophobicity classification for pattern-level comparison.
# Uses a biologically inclusive definition: aliphatic hydrophobics + aromatics.
# Aromatic residues (F, W, Y) are treated as hydrophobic here because:
#   - They preferentially partition to the hydrophobic core
#   - Antibody epitopes dominated by aromatics are recognised as hydrophobic patches
#   - KD scale places W (-0.9) and Y (-1.3) as slightly polar, which causes
#     systematic under-detection of hydrophobic patterns using those residues.
# This set is used ONLY for binary-pattern cross-correlation; the KD scale is
# still used for continuous post-alignment physicochemical scoring.
HYDROPHOBIC_SET: frozenset[str] = frozenset("ACFILMVWYF")  # all 20 aliphatic + aromatic AAs

# pKa values for isoelectric point calculation (Bjellqvist scale)
# N-terminus, C-terminus, and ionisable side chains
_PKA = {
    "nterm": 8.0,
    "cterm": 3.1,
    "D": 3.9,
    "E": 4.07,
    "H": 6.04,
    "C": 8.14,
    "Y": 10.46,
    "K": 10.54,
    "R": 12.48,
}
_POSITIVE_RESIDUES = {"H", "K", "R"}
_NEGATIVE_RESIDUES = {"D", "E", "C", "Y"}


# ---------------------------------------------------------------------------
# Internal helpers
# ---------------------------------------------------------------------------

def _pearson(xs: list[float], ys: list[float]) -> float:
    """Pearson correlation coefficient between two equal-length lists.

    Returns 0.0 if either list has zero variance (constant profile).
    """
    n = len(xs)
    if n < 2:
        return 0.0

    mx = sum(xs) / n
    my = sum(ys) / n

    num = sum((x - mx) * (y - my) for x, y in zip(xs, ys))
    denom_x = math.sqrt(sum((x - mx) ** 2 for x in xs))
    denom_y = math.sqrt(sum((y - my) ** 2 for y in ys))

    if denom_x == 0.0 or denom_y == 0.0:
        return 0.0
    return num / (denom_x * denom_y)


def _sliding_window_avg(values: list[float], window: int) -> list[float]:
    """Return sliding-window averages over *values* with the given window size.

    Output length = len(values) - window + 1.
    Uses a running sum for O(N) performance.
    """
    n = len(values)
    if n < window:
        return []

    running = sum(values[:window])
    result = [running / window]
    for i in range(window, n):
        running += values[i] - values[i - window]
        result.append(running / window)
    return result


def _per_residue_values(sequence: str, scale: dict[str, float], default: float = 0.0) -> list[float]:
    """Map each residue to a scalar value using *scale*. Unknown residues use *default*."""
    return [scale.get(aa, default) for aa in sequence.upper()]


def _net_charge_at_ph(sequence: str, ph: float) -> float:
    """Compute the net charge of *sequence* at a given pH.

    Uses the Henderson-Hasselbalch equation for each ionisable group.
    """
    seq = sequence.upper()

    # N-terminus: positive below pKa
    charge = 1.0 / (1.0 + 10.0 ** (ph - _PKA["nterm"]))
    # C-terminus: negative above pKa
    charge -= 1.0 / (1.0 + 10.0 ** (_PKA["cterm"] - ph))

    for aa in seq:
        if aa in _POSITIVE_RESIDUES:
            pka = _PKA.get(aa)
            if pka is not None:
                charge += 1.0 / (1.0 + 10.0 ** (ph - pka))
        elif aa in _NEGATIVE_RESIDUES:
            pka = _PKA.get(aa)
            if pka is not None:
                charge -= 1.0 / (1.0 + 10.0 ** (pka - ph))

    return charge


# ---------------------------------------------------------------------------
# Public API
# ---------------------------------------------------------------------------

def binary_hydrophobicity_profile(sequence: str) -> list[float]:
    """Binary hydrophobicity profile: 1.0 if residue is in HYDROPHOBIC_SET, else 0.0.

    Operates at per-residue resolution (no window averaging) to preserve the
    exact H/P pattern for cross-correlation. Use this for detecting convergent
    sequences (CONV-type) that share a hydrophobic/polar pattern but use
    different specific amino acids in each position.

    Args:
        sequence: Amino acid sequence (single-letter codes).

    Returns:
        List of 1.0 (hydrophobic) or 0.0 (polar/charged) per residue.
    """
    return [1.0 if aa.upper() in HYDROPHOBIC_SET else 0.0 for aa in sequence]


def hydrophobicity_profile(sequence: str, window: int = 5) -> list[float]:
    """Kyte-Doolittle sliding-window hydrophobicity.

    Args:
        sequence: Amino acid sequence (single-letter codes).
        window: Sliding window size (default 5).

    Returns:
        List of average hydrophobicity values per window position.
        Length = max(0, len(sequence) - window + 1).
    """
    if not sequence:
        return []
    values = _per_residue_values(sequence, KD_SCALE, default=0.0)
    return _sliding_window_avg(values, window)


def charge_profile(sequence: str, window: int = 5) -> list[float]:
    """Net charge per sliding window (K, R = +1; D, E = -1; H ≈ +0.1).

    Args:
        sequence: Amino acid sequence (single-letter codes).
        window: Sliding window size (default 5).

    Returns:
        List of net charge sum values per window position.
        Length = max(0, len(sequence) - window + 1).
    """
    if not sequence:
        return []
    values = _per_residue_values(sequence, CHARGE_SCALE, default=0.0)
    return _sliding_window_avg(values, window)


def isoelectric_point(sequence: str) -> float:
    """Approximate isoelectric point (pI) of a peptide sequence.

    Uses binary search over pH 0–14, locating the pH at which the
    Henderson-Hasselbalch net charge equals zero. Accurate to ±0.01 pH unit.

    Args:
        sequence: Amino acid sequence (single-letter codes).

    Returns:
        Estimated pI value in the range [0.0, 14.0].
        Returns 7.0 for empty sequences.
    """
    if not sequence:
        return 7.0

    lo, hi = 0.0, 14.0
    for _ in range(50):            # 50 bisection steps → ~14/2^50 ≈ 1e-14 precision
        mid = (lo + hi) / 2.0
        charge = _net_charge_at_ph(sequence, mid)
        if charge > 0:
            lo = mid
        else:
            hi = mid
    return (lo + hi) / 2.0


def max_cross_correlation(profile_a: list[float], profile_b: list[float]) -> float:
    """Maximum Pearson cross-correlation between two property profiles at all offsets.

    Slides the shorter profile along the longer one and computes the Pearson
    correlation at each offset. Returns the maximum value found.

    This is the key primitive for the global / CONV-rescue scan: two sequences
    can share a similar physicochemical envelope at a shifted position and this
    function will still detect it.

    Args:
        profile_a: First sliding-window property profile.
        profile_b: Second sliding-window property profile.

    Returns:
        Maximum Pearson correlation across all valid offsets, in [-1, 1].
        Returns 0.0 if either profile is too short to produce any valid window.
    """
    if not profile_a or not profile_b:
        return 0.0

    # Ensure a is the longer one
    if len(profile_b) > len(profile_a):
        profile_a, profile_b = profile_b, profile_a

    short_len = len(profile_b)
    long_len = len(profile_a)
    n_offsets = long_len - short_len + 1

    if n_offsets <= 0:
        return _pearson(profile_a, profile_b)

    best = -1.0
    for offset in range(n_offsets):
        window = profile_a[offset: offset + short_len]
        r = _pearson(window, profile_b)
        if r > best:
            best = r
    return best


def physicochemical_similarity(seq_a: str, seq_b: str) -> dict[str, float]:
    """Compare physicochemical properties of two aligned sequences.

    Designed for the *post-alignment* mode: both sequences should be the
    same length (the aligned region extracted by SW traceback). Gap characters
    ('-') are stripped before computing profiles so insertions do not distort
    the property averages.

    For the global / CONV-rescue use case (no reliable alignment), call this
    function with the full sequences and ignore hydrophobicity_correlation /
    charge_correlation; use sliding_window_correlation instead.

    Args:
        seq_a: First aligned sequence (may contain '-' gap characters).
        seq_b: Second aligned sequence (may contain '-' gap characters).

    Returns:
        Dictionary with keys:
            hydrophobicity_correlation: Pearson correlation on per-residue
                KD values of the ungapped aligned region (-1 to 1).
            charge_correlation: Pearson correlation on per-residue charge
                values of the ungapped aligned region (-1 to 1).
            pi_difference: |pI(seq_a) - pI(seq_b)| (ungapped sequences).
            hydrophobicity_cross_corr: Maximum sliding-window cross-correlation
                of hydrophobicity profiles (for misaligned / CONV-type pairs).
            charge_cross_corr: Maximum sliding-window cross-correlation of
                charge profiles.
            composite_score: Weighted blend of the above metrics (0 to 1).
                Higher = more physicochemically similar.
    """
    # Strip gap characters for pI and direct correlations
    clean_a = seq_a.replace("-", "").upper()
    clean_b = seq_b.replace("-", "").upper()

    # --- Per-residue profiles on the ungapped aligned region ---
    # For direct (post-alignment) correlation, we compare position-by-position
    # after removing gaps, keeping only matched columns (no '-' in either).
    matched_hydro_a: list[float] = []
    matched_hydro_b: list[float] = []
    matched_charge_a: list[float] = []
    matched_charge_b: list[float] = []

    # Walk aligned columns; include only non-gap paired columns
    len_a = len(seq_a)
    len_b = len(seq_b)
    min_len = min(len_a, len_b)
    for i in range(min_len):
        aa_a = seq_a[i].upper()
        aa_b = seq_b[i].upper()
        if aa_a != "-" and aa_b != "-":
            matched_hydro_a.append(KD_SCALE.get(aa_a, 0.0))
            matched_hydro_b.append(KD_SCALE.get(aa_b, 0.0))
            matched_charge_a.append(CHARGE_SCALE.get(aa_a, 0.0))
            matched_charge_b.append(CHARGE_SCALE.get(aa_b, 0.0))

    hydro_corr = _pearson(matched_hydro_a, matched_hydro_b)
    charge_corr = _pearson(matched_charge_a, matched_charge_b)

    # --- pI difference (global property; always computed on ungapped seqs) ---
    pi_a = isoelectric_point(clean_a)
    pi_b = isoelectric_point(clean_b)
    pi_diff = abs(pi_a - pi_b)

    # --- Sliding-window cross-correlation (for global / CONV-rescue mode) ---
    window = 5
    hydro_profile_a = hydrophobicity_profile(clean_a, window=window)
    hydro_profile_b = hydrophobicity_profile(clean_b, window=window)
    charge_profile_a = charge_profile(clean_a, window=window)
    charge_profile_b = charge_profile(clean_b, window=window)

    hydro_cross = max_cross_correlation(hydro_profile_a, hydro_profile_b)
    charge_cross = max_cross_correlation(charge_profile_a, charge_profile_b)

    # --- Binary H/P pattern cross-correlation ---
    # Compares the hydrophobic/polar pattern at per-residue resolution using the
    # inclusive HYDROPHOBIC_SET (aromatics + aliphatics) rather than KD values.
    # This is the primary signal for convergent CONV-type sequences, which share
    # a H/P pattern but use different amino acids (e.g. I vs L vs V in hydrophobic
    # positions, or W/Y that KD rates as slightly polar but are biologically
    # hydrophobic at an epitope surface).
    bin_a = binary_hydrophobicity_profile(clean_a)
    bin_b = binary_hydrophobicity_profile(clean_b)
    binary_hydro_cross = max_cross_correlation(bin_a, bin_b)

    # --- Composite score ---
    # Weighted blend. Binary H/P cross-correlation is the strongest single
    # predictor of convergent physicochemical similarity and is the key signal
    # for rescuing CONV-type pairs that SW misses.
    # KD hydrophobicity cross-corr adds nuance for sequences where the actual
    # magnitude of hydrophobicity (not just H/P class) is conserved.
    # Charge cross-corr catches charge-patterned pairs.
    # pI difference is a global penalty (normalised: 0 diff → 1.0, 7+ → 0.0).
    # Direct post-alignment correlations are a secondary refinement signal.
    pi_score = max(0.0, 1.0 - pi_diff / 7.0)
    aligned_bonus = (hydro_corr + charge_corr) / 2.0

    w_binary  = 0.40   # binary H/P pattern — primary CONV signal
    w_hydro   = 0.20   # KD sliding-window — magnitude refinement
    w_charge  = 0.20   # charge pattern
    w_pi      = 0.10   # global pI similarity
    w_aligned = 0.10   # post-alignment direct correlations

    composite = (
        w_binary  * max(0.0, binary_hydro_cross)
        + w_hydro   * max(0.0, hydro_cross)
        + w_charge  * max(0.0, charge_cross)
        + w_pi      * pi_score
        + w_aligned * max(0.0, aligned_bonus)
    )

    return {
        "hydrophobicity_correlation": hydro_corr,
        "charge_correlation": charge_corr,
        "pi_difference": pi_diff,
        "hydrophobicity_cross_corr": hydro_cross,
        "charge_cross_corr": charge_cross,
        "binary_hydro_cross_corr": binary_hydro_cross,
        "composite_score": composite,
    }
