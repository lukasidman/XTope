"""Tests for the physicochemical similarity module.

Covers:
  - hydrophobicity_profile
  - charge_profile
  - isoelectric_point
  - max_cross_correlation
  - physicochemical_similarity (post-alignment and global/CONV modes)
"""

import math
import pytest

from xtope.physicochemical import (
    KD_SCALE,
    CHARGE_SCALE,
    HYDROPHOBIC_SET,
    hydrophobicity_profile,
    binary_hydrophobicity_profile,
    charge_profile,
    isoelectric_point,
    max_cross_correlation,
    physicochemical_similarity,
    _pearson,
    _sliding_window_avg,
)


# ---------------------------------------------------------------------------
# _pearson helper
# ---------------------------------------------------------------------------

class TestPearson:
    def test_perfect_positive(self):
        xs = [1.0, 2.0, 3.0, 4.0]
        assert _pearson(xs, xs) == pytest.approx(1.0)

    def test_perfect_negative(self):
        xs = [1.0, 2.0, 3.0]
        ys = [3.0, 2.0, 1.0]
        assert _pearson(xs, ys) == pytest.approx(-1.0)

    def test_zero_correlation(self):
        # Orthogonal-ish vectors
        xs = [1.0, -1.0, 1.0, -1.0]
        ys = [1.0, 1.0, -1.0, -1.0]
        r = _pearson(xs, ys)
        assert abs(r) < 0.05  # near-zero

    def test_constant_returns_zero(self):
        xs = [2.0, 2.0, 2.0]
        ys = [1.0, 2.0, 3.0]
        assert _pearson(xs, ys) == 0.0

    def test_single_element_returns_zero(self):
        assert _pearson([1.0], [1.0]) == 0.0

    def test_empty_returns_zero(self):
        assert _pearson([], []) == 0.0


# ---------------------------------------------------------------------------
# _sliding_window_avg helper
# ---------------------------------------------------------------------------

class TestSlidingWindowAvg:
    def test_window_equals_length(self):
        vals = [1.0, 2.0, 3.0]
        result = _sliding_window_avg(vals, 3)
        assert result == pytest.approx([2.0])

    def test_window_one(self):
        vals = [1.0, 3.0, 5.0]
        assert _sliding_window_avg(vals, 1) == pytest.approx([1.0, 3.0, 5.0])

    def test_standard_window(self):
        vals = [1.0, 2.0, 3.0, 4.0, 5.0]
        result = _sliding_window_avg(vals, 3)
        assert len(result) == 3
        assert result[0] == pytest.approx(2.0)
        assert result[1] == pytest.approx(3.0)
        assert result[2] == pytest.approx(4.0)

    def test_window_larger_than_sequence_returns_empty(self):
        assert _sliding_window_avg([1.0, 2.0], 5) == []

    def test_empty_returns_empty(self):
        assert _sliding_window_avg([], 3) == []


# ---------------------------------------------------------------------------
# binary_hydrophobicity_profile
# ---------------------------------------------------------------------------

class TestBinaryHydrophobicityProfile:
    def test_empty_sequence(self):
        assert binary_hydrophobicity_profile("") == []

    def test_all_hydrophobic(self):
        # A, I, L, V, F, W, Y are all in HYDROPHOBIC_SET
        result = binary_hydrophobicity_profile("AILVFW")
        assert all(v == 1.0 for v in result)

    def test_all_polar(self):
        # D, E, K, R, N, Q, S, T, H, G, P are not hydrophobic
        result = binary_hydrophobicity_profile("DEKRN")
        assert all(v == 0.0 for v in result)

    def test_aromatic_counted_as_hydrophobic(self):
        # W and Y are aromatic; should be 1.0 even though KD scale rates them near 0
        result = binary_hydrophobicity_profile("WY", )
        assert result == [1.0, 1.0]

    def test_mixed(self):
        result = binary_hydrophobicity_profile("AKID")  # A=H, K=P, I=H, D=P
        assert result == [1.0, 0.0, 1.0, 0.0]

    def test_lowercase_handled(self):
        upper = binary_hydrophobicity_profile("AILVFW")
        lower = binary_hydrophobicity_profile("ailvfw")
        assert upper == lower

    def test_length_equals_sequence(self):
        seq = "ACDEFGHIKLMNPQRSTVWY"
        assert len(binary_hydrophobicity_profile(seq)) == len(seq)

    def test_conv_pattern_similarity(self):
        """Two sequences with same H/P pattern but different residues should correlate highly."""
        # Position 0,2,4 = hydrophobic; 1,3,5 = polar
        seq_a = "ISISSV"  # I=H, S=P, I=H, S=P, S=P, V=H  → H P H P P H
        seq_b = "VSVSV"   # V=H, S=P, V=H, S=P, V=H        → H P H P H
        ba = binary_hydrophobicity_profile(seq_a)
        bb = binary_hydrophobicity_profile(seq_b)
        corr = max_cross_correlation(ba, bb)
        assert corr > 0.5


# ---------------------------------------------------------------------------
# hydrophobicity_profile
# ---------------------------------------------------------------------------

class TestHydrophobicityProfile:
    def test_empty_sequence(self):
        assert hydrophobicity_profile("") == []

    def test_single_residue_with_window_gt_1_returns_empty(self):
        assert hydrophobicity_profile("A", window=5) == []

    def test_window_1_returns_per_residue_values(self):
        seq = "AILVG"
        result = hydrophobicity_profile(seq, window=1)
        expected = [KD_SCALE[aa] for aa in seq]
        assert result == pytest.approx(expected)

    def test_known_value_window_5(self):
        # A(1.8) + I(4.5) + L(3.8) + V(4.2) + G(-0.4) / 5 = 13.9/5 = 2.78
        seq = "AILVG"
        result = hydrophobicity_profile(seq, window=5)
        assert len(result) == 1
        assert result[0] == pytest.approx(2.78)

    def test_lowercase_handled(self):
        upper = hydrophobicity_profile("AILVG", window=1)
        lower = hydrophobicity_profile("ailvg", window=1)
        assert upper == pytest.approx(lower)

    def test_unknown_residue_treated_as_zero(self):
        result = hydrophobicity_profile("XA", window=1)
        assert result[0] == pytest.approx(0.0)
        assert result[1] == pytest.approx(KD_SCALE["A"])

    def test_output_length(self):
        seq = "ACDEFGHIKLMNPQRSTVWY"  # all 20 standard AAs
        result = hydrophobicity_profile(seq, window=5)
        assert len(result) == len(seq) - 5 + 1


# ---------------------------------------------------------------------------
# charge_profile
# ---------------------------------------------------------------------------

class TestChargeProfile:
    def test_empty_sequence(self):
        assert charge_profile("") == []

    def test_neutral_residues_give_zero(self):
        result = charge_profile("AAVVLL", window=1)
        assert all(v == pytest.approx(0.0) for v in result)

    def test_positive_residues(self):
        result = charge_profile("KR", window=1)
        assert result[0] == pytest.approx(CHARGE_SCALE["K"])
        assert result[1] == pytest.approx(CHARGE_SCALE["R"])

    def test_negative_residues(self):
        result = charge_profile("DE", window=1)
        assert result[0] == pytest.approx(CHARGE_SCALE["D"])
        assert result[1] == pytest.approx(CHARGE_SCALE["E"])

    def test_window_5_known_value(self):
        # K(+1) R(+1) D(-1) E(-1) A(0) → sum = 0, avg = 0
        result = charge_profile("KRDEA", window=5)
        assert len(result) == 1
        assert result[0] == pytest.approx(0.0)

    def test_output_length(self):
        seq = "KRDEKRDEAA"
        result = charge_profile(seq, window=3)
        assert len(result) == len(seq) - 3 + 1


# ---------------------------------------------------------------------------
# isoelectric_point
# ---------------------------------------------------------------------------

class TestIsoelectricPoint:
    def test_empty_returns_7(self):
        assert isoelectric_point("") == pytest.approx(7.0)

    def test_acidic_peptide_low_pi(self):
        # EEEDD — very acidic, pI should be < 4
        pi = isoelectric_point("EEEEDD")
        assert pi < 4.5

    def test_basic_peptide_high_pi(self):
        # KKKKR — very basic, pI should be > 10
        pi = isoelectric_point("KKKKR")
        assert pi > 10.0

    def test_neutral_peptide_near_7(self):
        # Mix of equal charge; pI should be in a reasonable midrange
        pi = isoelectric_point("AAVVLLIIMM")  # all hydrophobic, no charges
        # Only termini contribute — expect near 5–7 range
        assert 3.0 < pi < 9.0

    def test_result_in_valid_range(self):
        for seq in ["ACDEFGHIKLMNPQRSTVWY", "KKKKEEEE", "MHHHHHHGSSG"]:
            pi = isoelectric_point(seq)
            assert 0.0 <= pi <= 14.0

    def test_known_approximate(self):
        # Polyarginine — pI should be very high (near pKa of R = 12.48)
        pi = isoelectric_point("RRRR")
        assert pi > 11.0

    def test_lowercase_input(self):
        pi_upper = isoelectric_point("KRDEA")
        pi_lower = isoelectric_point("krdea")
        assert pi_upper == pytest.approx(pi_lower, abs=0.01)


# ---------------------------------------------------------------------------
# max_cross_correlation
# ---------------------------------------------------------------------------

class TestMaxCrossCorrelation:
    def test_identical_profiles(self):
        p = [1.0, 2.0, 3.0, 2.0, 1.0]
        assert max_cross_correlation(p, p) == pytest.approx(1.0)

    def test_empty_profiles_return_zero(self):
        assert max_cross_correlation([], [1.0, 2.0]) == 0.0
        assert max_cross_correlation([1.0, 2.0], []) == 0.0

    def test_shifted_match_detected(self):
        # Profile b is the same as a[2:] — a shifted match
        long_p = [0.0, 0.0, 1.0, 2.0, 3.0, 2.0, 1.0]
        short_p = [1.0, 2.0, 3.0, 2.0, 1.0]
        corr = max_cross_correlation(long_p, short_p)
        assert corr == pytest.approx(1.0)

    def test_anti_correlated_profiles(self):
        p = [1.0, 2.0, 3.0]
        q = [-1.0, -2.0, -3.0]
        corr = max_cross_correlation(p, q)
        assert corr == pytest.approx(-1.0)

    def test_commutative(self):
        p = [1.0, 2.0, 3.0, 4.0, 5.0]
        q = [2.0, 3.0, 4.0]
        assert max_cross_correlation(p, q) == pytest.approx(
            max_cross_correlation(q, p)
        )

    def test_equal_length_profiles(self):
        p = [1.0, -1.0, 1.0, -1.0]
        q = [1.0, -1.0, 1.0, -1.0]
        assert max_cross_correlation(p, q) == pytest.approx(1.0)


# ---------------------------------------------------------------------------
# physicochemical_similarity
# ---------------------------------------------------------------------------

class TestPhysicochemicalSimilarity:
    def test_identical_sequences(self):
        seq = "ACDEFGHIKLMNPQRSTVWY"
        result = physicochemical_similarity(seq, seq)
        assert result["hydrophobicity_correlation"] == pytest.approx(1.0)
        assert result["charge_correlation"] == pytest.approx(1.0, abs=0.01)
        assert result["pi_difference"] == pytest.approx(0.0)
        assert result["composite_score"] > 0.8

    def test_result_keys_present(self):
        result = physicochemical_similarity("ACDE", "KLMN")
        expected_keys = {
            "hydrophobicity_correlation",
            "charge_correlation",
            "pi_difference",
            "hydrophobicity_cross_corr",
            "charge_cross_corr",
            "binary_hydro_cross_corr",
            "composite_score",
        }
        assert set(result.keys()) == expected_keys

    def test_composite_score_in_range(self):
        result = physicochemical_similarity("ACDEFGHIK", "LMNPQRSTVW")
        assert 0.0 <= result["composite_score"] <= 1.0

    def test_gap_characters_stripped(self):
        # Gaps should be ignored — result should match gapless comparison
        seq_a = "ACD-EF"
        seq_b = "ACDEF"
        result_with_gap = physicochemical_similarity(seq_a, seq_b)
        result_clean = physicochemical_similarity("ACDEF", "ACDEF")
        # pI should be same for the clean ungapped sequences
        assert result_with_gap["pi_difference"] == pytest.approx(
            result_clean["pi_difference"], abs=0.1
        )

    def test_very_different_sequences_low_composite(self):
        # Highly hydrophobic vs highly charged — should give low composite
        hydrophobic = "IIIIIVVVVLLLL"
        charged = "KKKKEEEERRRR"
        result = physicochemical_similarity(hydrophobic, charged)
        assert result["composite_score"] < 0.7

    def test_conv_family_detected(self):
        """CONV family: same hydrophobic/polar pattern, different residues.

        Both sequences follow HHPPHHPP... pattern but use different amino acids.
        Hydrophobic cross-correlation should be high even if direct correlation
        of per-residue values is imperfect.
        """
        # H = hydrophobic (I/L/V), P = polar (S/T/N)
        conv_a = "IISSIISSTTIISSII"   # H=I or L, P=S or T
        conv_b = "LLTTLLTTSSLLTTLL"
        result = physicochemical_similarity(conv_a, conv_b)
        assert result["hydrophobicity_cross_corr"] > 0.6

    def test_empty_sequences(self):
        result = physicochemical_similarity("", "")
        assert result["pi_difference"] == pytest.approx(0.0)
        assert result["composite_score"] >= 0.0

    def test_single_residue(self):
        result = physicochemical_similarity("A", "V")
        # Both hydrophobic — should be reasonably similar in cross-corr
        assert isinstance(result["composite_score"], float)

    def test_pi_difference_acidic_vs_basic(self):
        # EEEE is very acidic, RRRR is very basic
        result = physicochemical_similarity("EEEEEEEE", "RRRRRRRR")
        assert result["pi_difference"] > 6.0
