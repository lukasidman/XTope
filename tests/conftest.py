"""Shared fixtures for antigen_screener tests."""

from __future__ import annotations

import pytest


@pytest.fixture
def sample_sequence() -> str:
    """A short antigen sequence with the His6-ABP tag for testing."""
    return "MHHHHHHGSSGVKQTLNFDLLKLAGDVESNPGPAGSK"


@pytest.fixture
def sample_sequence_no_tag() -> str:
    """A short antigen sequence without the tag."""
    return "VKQTLNFDLLKLAGDVESNPGPAGSK"


@pytest.fixture
def tmp_db(tmp_path) -> str:
    """Path to a temporary SQLite database."""
    return str(tmp_path / "test.db")
