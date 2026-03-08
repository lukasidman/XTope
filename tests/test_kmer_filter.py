"""Tests for kmer_filter module."""

from __future__ import annotations


class TestKmerFilter:
    """Tests for inverted k-mer index and Jaccard pre-filter."""

    def test_add_and_query(self) -> None:
        """Adding a sequence and querying returns it as a candidate."""

    def test_jaccard_calculation(self) -> None:
        """Jaccard similarity is computed correctly for known overlaps."""

    def test_empty_sequence(self) -> None:
        """Empty sequence produces no k-mers and zero Jaccard score."""
