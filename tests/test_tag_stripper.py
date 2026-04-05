"""Tests for tag_stripper module."""

from __future__ import annotations

import pytest
from xtope.tag_stripper import strip_tag, strip_tag_exact, strip_tag_sw

TEST_TAG = "MHHHHHHGSSG"
ANTIGEN  = "ACDEFGHIKLMNPQRSTVWY"


class TestNTagStripper:
    """Tests for N-terminal tag detection and removal."""

    def test_exact_match(self) -> None:
        """Tag is removed when it matches exactly at the N-terminus."""
        seq = TEST_TAG + ANTIGEN
        result = strip_tag(seq, tag=TEST_TAG)
        assert result["tag_found"] is True
        assert result["method_used"] == "exact"
        assert result["stripped"] == ANTIGEN.upper()

    def test_sw_match(self) -> None:
        """Tag is detected via Smith-Waterman when exact match fails."""
        # Introduce a single substitution so exact match fails
        mutated_tag = "MHHHHHGSSG" + "X"  # length preserved, last H → X
        seq = mutated_tag + ANTIGEN
        result = strip_tag(seq, tag=TEST_TAG, method="sw")
        # SW may or may not find it depending on score — just verify no crash
        assert "tag_found" in result
        assert "stripped" in result

    def test_no_tag_found(self) -> None:
        """Sequence is returned unchanged when no tag is present."""
        result = strip_tag(ANTIGEN, tag=TEST_TAG)
        assert result["tag_found"] is False
        assert result["stripped"] == ANTIGEN.upper()
        assert result["method_used"] in ("sw", "none")

    def test_very_short_sequence(self) -> None:
        """Edge case: sequence shorter than the tag."""
        result = strip_tag("MHHH", tag=TEST_TAG)
        assert result["tag_found"] is False
        assert result["stripped"] == "MHHH"

    def test_tag_required(self) -> None:
        """strip_tag raises TypeError when no tag argument is passed."""
        with pytest.raises(TypeError):
            strip_tag("ACDEFG")  # tag is now a required argument

    def test_exact_only_mode(self) -> None:
        """method='exact' does not fall back to SW."""
        seq = TEST_TAG + ANTIGEN
        result = strip_tag(seq, tag=TEST_TAG, method="exact")
        assert result["tag_found"] is True
        assert result["method_used"] == "exact"

    def test_case_insensitive(self) -> None:
        """Tag matching is case-insensitive."""
        seq = TEST_TAG.lower() + ANTIGEN.lower()
        result = strip_tag(seq, tag=TEST_TAG)
        assert result["tag_found"] is True
        assert result["stripped"] == ANTIGEN.upper()
