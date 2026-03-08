"""Tests for pipeline module."""

from __future__ import annotations


class TestPipeline:
    """Integration tests for the full screening pipeline."""

    def test_full_pipeline_small_dataset(self) -> None:
        """End-to-end pipeline run on a small synthetic dataset."""
