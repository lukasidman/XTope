"""Tests for tag_stripper module."""

from __future__ import annotations


class TestTagStripper:
    """Tests for His6-ABP tag detection and removal."""

    def test_exact_match(self) -> None:
        """Tag is removed when it matches exactly at the N-terminus."""

    def test_sw_match(self) -> None:
        """Tag is detected via Smith-Waterman when exact match fails."""

    def test_no_tag_found(self) -> None:
        """Sequence is returned unchanged when no tag is present."""

    def test_very_short_sequence(self) -> None:
        """Edge case: sequence shorter than the tag."""
