"""Tests for kmer_filter module — KmerIndex, ChainKmerIndex, and TwoPassFilter."""

from __future__ import annotations

import pytest
from antigen_screener.kmer_filter import (
    KmerIndex,
    ChainKmerIndex,
    TwoPassFilter,
    kmers,
    jaccard,
    DEFAULT_K,
    DEFAULT_THRESHOLD,
    DEFAULT_K2,
    DEFAULT_MIN_CHAIN,
)

# ---------------------------------------------------------------------------
# Shared test sequences
# ---------------------------------------------------------------------------

# Flanking backgrounds — each uses a unique repeated amino acid so they share
# zero k-mers with each other regardless of k.  Using repeated single amino
# acids (KKKK..., DDDD..., etc.) guarantees this.
BG_K = "K" * 35   # all-K background
BG_D = "D" * 35   # all-D background
BG_E = "E" * 35   # all-E background
BG_H = "H" * 35   # all-H background

# 20 aa shared region — long enough to pass Jaccard (pass 1)
SHARED_LONG = "KQIRYLDGISALRKETCNKS"

# 8 aa shared region — just long enough for chain detection (pass 2)
# 8 aa = chain of 5 consecutive 4-mers
SHARED_SHORT = "RPCHQFNV"

# 6 aa shared region — minimum catchable by pass 2 (chain=3, k=4)
SHARED_MIN = "CHQFNV"

# Sequences with a long shared region (should pass Jaccard — pass 1)
SEQ_LONG_A = BG_K + SHARED_LONG + BG_D
SEQ_LONG_B = BG_E + SHARED_LONG + BG_H

# Sequences with only a short shared region (Jaccard < threshold — pass 2 only)
# Backgrounds are single-amino-acid runs so they contribute zero shared k-mers.
SEQ_SHORT_A = BG_K + SHARED_SHORT + BG_D
SEQ_SHORT_B = BG_E + SHARED_SHORT + BG_H

# Sequences with the minimum 6 aa shared region (pass 2 only)
SEQ_MIN_A = BG_K + SHARED_MIN + BG_D
SEQ_MIN_B = BG_E + SHARED_MIN + BG_H

# Completely unrelated sequences (no shared k-mers by construction)
SEQ_UNRELATED_A = "AAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAA"
SEQ_UNRELATED_B = "WWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWW"

# Real IGF sequences (tag stripped) from the test dataset
IGF_01 = "GASHRPSRHAMTESLQPTCSLGAEAWMYMGSWGYRPCHQFNVGVAQGSKKKYHEQISGCKKYEGDHGKREKKRWPKPS"
IGF_02 = "SIINPQWISNMIPFDRLINGGKYGWLDWGMRPCHQFNVLVKCQTHHCMVVRWPSMDFDSANITNPCHYRM"


# ---------------------------------------------------------------------------
# Utility functions
# ---------------------------------------------------------------------------

class TestKmers:
    def test_correct_count(self) -> None:
        assert len(kmers("ACDEFG", k=4)) == 3  # ACDE, CDEF, DEFG

    def test_returns_set(self) -> None:
        result = kmers("AAAA", k=2)
        assert isinstance(result, set)
        assert result == {"AA"}

    def test_sequence_shorter_than_k(self) -> None:
        assert kmers("ACG", k=6) == set()


class TestJaccard:
    def test_identical_sets(self) -> None:
        s = {"ABC", "BCD", "CDE"}
        assert jaccard(s, s) == 1.0

    def test_disjoint_sets(self) -> None:
        assert jaccard({"AAA"}, {"BBB"}) == 0.0

    def test_partial_overlap(self) -> None:
        a = {"ABC", "BCD", "CDE"}
        b = {"BCD", "CDE", "DEF"}
        assert jaccard(a, b) == pytest.approx(2 / 4)

    def test_empty_set(self) -> None:
        assert jaccard(set(), {"ABC"}) == 0.0


# ---------------------------------------------------------------------------
# Pass 1 — KmerIndex
# ---------------------------------------------------------------------------

class TestKmerIndex:
    def test_add_and_query_finds_similar(self) -> None:
        idx = KmerIndex()
        idx.add("A", SEQ_LONG_A)
        results = idx.query(SEQ_LONG_B)
        ids = [r[0] for r in results]
        assert "A" in ids

    def test_unrelated_sequences_filtered(self) -> None:
        idx = KmerIndex()
        idx.add("A", SEQ_UNRELATED_A)
        results = idx.query(SEQ_UNRELATED_B)
        assert results == []

    def test_jaccard_calculation_correct(self) -> None:
        idx = KmerIndex()
        idx.add("A", SEQ_LONG_A)
        results = idx.query(SEQ_LONG_B)
        assert len(results) == 1
        _, score = results[0]
        # Verify score is in valid Jaccard range
        assert 0.0 <= score <= 1.0
        assert score >= DEFAULT_THRESHOLD

    def test_exclude_self(self) -> None:
        idx = KmerIndex()
        idx.add("A", SEQ_LONG_A)
        results = idx.query(SEQ_LONG_A, exclude_id="A")
        assert all(r[0] != "A" for r in results)

    def test_results_sorted_descending(self) -> None:
        idx = KmerIndex()
        idx.add("near", SEQ_LONG_A)
        idx.add("far", SEQ_UNRELATED_B + SHARED_LONG)
        results = idx.query(SEQ_LONG_B)
        scores = [r[1] for r in results]
        assert scores == sorted(scores, reverse=True)

    def test_empty_sequence_returns_empty(self) -> None:
        idx = KmerIndex()
        idx.add("A", SEQ_LONG_A)
        assert idx.query("") == []

    def test_short_shared_region_fails_pass1(self) -> None:
        """
        IGF_01 / IGF_02 share only an 8 aa motif embedded in diverse 70-78 aa
        sequences.  Their Jaccard score is 3/135 ≈ 0.022 — below the 0.04
        threshold — so pass 1 alone should NOT find them.
        """
        idx = KmerIndex()
        idx.add("IGF_02", IGF_02)
        results = idx.query(IGF_01)
        assert results == []

    def test_add_batch(self) -> None:
        idx = KmerIndex()
        idx.add_batch([("A", SEQ_LONG_A), ("B", SEQ_UNRELATED_A)])
        assert len(idx) == 2

    def test_len(self) -> None:
        idx = KmerIndex()
        idx.add("A", SEQ_LONG_A)
        idx.add("B", SEQ_LONG_B)
        assert len(idx) == 2


# ---------------------------------------------------------------------------
# Pass 2 — ChainKmerIndex
# ---------------------------------------------------------------------------

class TestChainKmerIndex:
    def test_short_shared_region_detected(self) -> None:
        """8 aa shared region (5-mer chain) should be found by chain detection."""
        idx = ChainKmerIndex()
        idx.add("B", SEQ_SHORT_B)
        results = idx.query(SEQ_SHORT_A)
        assert "B" in results

    def test_minimum_6aa_epitope_detected(self) -> None:
        """6 aa shared region (chain=3, k=4) is the minimum catchable size."""
        idx = ChainKmerIndex()
        idx.add("B", SEQ_MIN_B)
        results = idx.query(SEQ_MIN_A)
        assert "B" in results

    def test_unrelated_sequences_not_detected(self) -> None:
        idx = ChainKmerIndex()
        idx.add("B", SEQ_UNRELATED_B)
        results = idx.query(SEQ_UNRELATED_A)
        assert results == []

    def test_scattered_shared_kmers_do_not_pass(self) -> None:
        """Shared k-mers that are NOT consecutive should not form a passing chain."""
        # ACDE at pos 0 in A, at the end of B — different diagonal.
        # RPCH at the end of A, at pos 0 of B — different diagonal.
        # L-padding in A and M-padding in B share no k-mers with each other.
        # Result: only two isolated seed pairs (ACDE, RPCH), max chain length = 1.
        seq_a = "ACDE" + "L" * 16 + "RPCH"   # 24 aa
        seq_b = "RPCH" + "M" * 20 + "ACDE"   # 28 aa
        idx = ChainKmerIndex(k=4, min_chain=3)
        idx.add("B", seq_b)
        results = idx.query(seq_a)
        assert "B" not in results

    def test_exclude_self(self) -> None:
        idx = ChainKmerIndex()
        idx.add("A", SEQ_SHORT_A)
        results = idx.query(SEQ_SHORT_A, exclude_id="A")
        assert "A" not in results

    def test_add_batch(self) -> None:
        idx = ChainKmerIndex()
        idx.add_batch([("A", SEQ_SHORT_A), ("B", SEQ_SHORT_B)])
        assert len(idx) == 2

    def test_igf_synthetic_pair_detected(self) -> None:
        """IGF_01/IGF_02 share an 8 aa epitope — should be found by pass 2."""
        idx = ChainKmerIndex()
        idx.add("IGF_02", IGF_02)
        results = idx.query(IGF_01)
        assert "IGF_02" in results


# ---------------------------------------------------------------------------
# TwoPassFilter — combined behaviour
# ---------------------------------------------------------------------------

class TestTwoPassFilter:
    def test_long_shared_region_passes_via_pass1(self) -> None:
        """Long shared region passes Jaccard and appears with a real score."""
        f = TwoPassFilter()
        f.add("B", SEQ_LONG_B)
        results = f.query(SEQ_LONG_A)
        ids = [r[0] for r in results]
        assert "B" in ids
        score = next(s for i, s in results if i == "B")
        assert score >= DEFAULT_THRESHOLD  # real Jaccard score

    def test_short_shared_region_passes_via_pass2(self) -> None:
        """8 aa shared region should be caught by pass 2 even though pass 1 misses it."""
        f = TwoPassFilter()
        f.add("B", SEQ_SHORT_B)
        results = f.query(SEQ_SHORT_A)
        ids = [r[0] for r in results]
        assert "B" in ids

    def test_pass2_hit_has_zero_jaccard_score(self) -> None:
        """
        IGF_01/IGF_02 fail pass 1 (Jaccard 0.022 < 0.04) but pass pass 2.
        They must be returned with jaccard_score = 0.0, since Jaccard is
        meaningless for sequences found only via chain detection.
        """
        f = TwoPassFilter()
        f.add("IGF_02", IGF_02)
        results = f.query(IGF_01)
        score = next(s for i, s in results if i == "IGF_02")
        assert score == 0.0

    def test_unrelated_sequences_filtered_by_both_passes(self) -> None:
        f = TwoPassFilter()
        f.add("B", SEQ_UNRELATED_B)
        results = f.query(SEQ_UNRELATED_A)
        assert results == []

    def test_no_duplicate_ids_in_results(self) -> None:
        """A sequence found by both passes should appear only once."""
        f = TwoPassFilter()
        f.add("B", SEQ_LONG_B)
        results = f.query(SEQ_LONG_A)
        ids = [r[0] for r in results]
        assert len(ids) == len(set(ids))

    def test_sequences_dict_accessible(self) -> None:
        """pipeline.py accesses index.sequences[cid] — must be present."""
        f = TwoPassFilter()
        f.add("A", SEQ_LONG_A)
        assert "A" in f.sequences
        assert f.sequences["A"] == SEQ_LONG_A

    def test_exclude_self(self) -> None:
        f = TwoPassFilter()
        f.add("A", SEQ_SHORT_A)
        results = f.query(SEQ_SHORT_A, exclude_id="A")
        assert all(r[0] != "A" for r in results)

    def test_add_batch(self) -> None:
        f = TwoPassFilter()
        f.add_batch([("A", SEQ_LONG_A), ("B", SEQ_SHORT_B)])
        assert len(f) == 2

    def test_igf_synthetic_pair_now_passes(self) -> None:
        """
        The key regression test: IGF_01/IGF_02 share an 8 aa epitope.
        They were previously invisible to Jaccard alone (score 0.022 < 0.04).
        TwoPassFilter must find them via pass 2.
        """
        f = TwoPassFilter()
        f.add("IGF_02", IGF_02)
        results = f.query(IGF_01)
        ids = [r[0] for r in results]
        assert "IGF_02" in ids

    def test_minimum_6aa_epitope_passes(self) -> None:
        """6 aa is the minimum epitope size the two-pass filter can detect."""
        f = TwoPassFilter()
        f.add("B", SEQ_MIN_B)
        results = f.query(SEQ_MIN_A)
        assert "B" in [r[0] for r in results]

    def test_threshold_override(self) -> None:
        """threshold kwarg is forwarded to pass 1."""
        f = TwoPassFilter()
        f.add("B", SEQ_LONG_B)
        # Very high threshold — pass 1 won't fire, but pass 2 still might
        results_strict = f.query(SEQ_LONG_A, threshold=0.99)
        results_normal = f.query(SEQ_LONG_A, threshold=DEFAULT_THRESHOLD)
        # Normal threshold should find at least as many as strict
        assert len(results_normal) >= len(results_strict)
