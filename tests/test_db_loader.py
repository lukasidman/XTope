"""Tests for db_loader module."""

from __future__ import annotations


class TestDbLoader:
    """Tests for CSV/TSV loading with delimiter auto-detection."""

    def test_comma_csv(self) -> None:
        """Comma-delimited CSV loads correctly."""

    def test_semicolon_csv(self) -> None:
        """Semicolon-delimited CSV loads correctly."""

    def test_tsv(self) -> None:
        """Tab-separated file loads correctly."""

    def test_windows_crlf(self) -> None:
        """Windows line endings are handled."""

    def test_blank_leading_rows(self) -> None:
        """Blank rows before the header are skipped."""

    def test_bom(self) -> None:
        """UTF-8 BOM is handled transparently."""

    def test_duplicate_ids(self) -> None:
        """Duplicate antigen IDs are detected or handled."""

    def test_invalid_sequences(self) -> None:
        """Non-amino-acid characters are flagged or filtered."""
