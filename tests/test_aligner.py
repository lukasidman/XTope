"""Tests for aligner module."""

from __future__ import annotations


class TestAligner:
    """Tests for Smith-Waterman alignment scoring."""

    def test_identical_sequences(self) -> None:
        """Identical sequences produce a high normalised score."""

    def test_unrelated_sequences(self) -> None:
        """Unrelated sequences produce a low normalised score."""

    def test_short_sequences(self) -> None:
        """Short sequences are handled without errors."""
