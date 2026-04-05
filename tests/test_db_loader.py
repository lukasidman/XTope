"""Tests for db_loader module."""

from __future__ import annotations

import pytest
from pathlib import Path

from xtope.db_loader import (
    load_csv,
    load_fasta,
    load_sequences,
    sniff_delimiter,
    clean_sequence,
    is_valid_sequence,
    detect_columns,
    iter_csv,
)

# --- Helpers ---

SEQ_A = "MHHHHHHGSSGVKQTLNFDLLKLAGDVESNPGPAGSK"
SEQ_B = "MHHHHHHGSSGACDEFGHIKLMNPQRSTVWYACDEFGH"


def _write(tmp_path: Path, name: str, content: str) -> Path:
    p = tmp_path / name
    p.write_text(content, encoding="utf-8")
    return p


def _write_bytes(tmp_path: Path, name: str, data: bytes) -> Path:
    p = tmp_path / name
    p.write_bytes(data)
    return p


# --- CSV delimiter detection ---


class TestSniffDelimiter:
    def test_comma(self, tmp_path: Path) -> None:
        f = _write(tmp_path, "c.csv", "name,sequence\nAG1,ACDEFGHIKLMNPQRSTVWY\n")
        assert sniff_delimiter(f) == ","

    def test_semicolon(self, tmp_path: Path) -> None:
        f = _write(tmp_path, "s.csv", ";name;sequence\nAG1;ACDEFGHIKLMNPQRSTVWY\n")
        assert sniff_delimiter(f) == ";"

    def test_tab(self, tmp_path: Path) -> None:
        f = _write(tmp_path, "t.tsv", "name\tsequence\nAG1\tACDEFGHIKLMNPQRSTVWY\n")
        assert sniff_delimiter(f) == "\t"

    def test_skips_blank_leading_lines(self, tmp_path: Path) -> None:
        f = _write(tmp_path, "b.csv", "\n\n;name;sequence\nAG1;ACDEFGHIKLMNPQRSTVWY\n")
        assert sniff_delimiter(f) == ";"

    def test_empty_file(self, tmp_path: Path) -> None:
        f = _write(tmp_path, "e.csv", "")
        # Should not crash, returns default
        assert sniff_delimiter(f) == ","


# --- CSV loading ---


class TestLoadCsvComma:
    def test_basic_comma(self, tmp_path: Path) -> None:
        f = _write(tmp_path, "c.csv", f"name,sequence\nAG1,{SEQ_A}\nAG2,{SEQ_B}\n")
        records = load_csv(f, verbose=False)
        assert len(records) == 2
        assert records[0][0] == "AG1"
        assert records[1][0] == "AG2"


class TestLoadCsvSemicolon:
    def test_semicolon_with_leading_empty_col(self, tmp_path: Path) -> None:
        """The real-world format: ;name;sequence with empty first column."""
        content = f";name;sequence\n;AG1;{SEQ_A}\n;AG2;{SEQ_B}\n"
        f = _write(tmp_path, "s.csv", content)
        records = load_csv(f, verbose=False)
        assert len(records) == 2
        assert records[0][0] == "AG1"

    def test_semicolon_auto_detected(self, tmp_path: Path) -> None:
        """Semicolon delimiter is detected even with .csv extension."""
        content = f"name;sequence\nAG1;{SEQ_A}\n"
        f = _write(tmp_path, "data.csv", content)
        records = load_csv(f, verbose=False)
        assert len(records) == 1
        assert records[0][0] == "AG1"


class TestLoadCsvTsv:
    def test_tab_separated(self, tmp_path: Path) -> None:
        content = f"name\tsequence\nAG1\t{SEQ_A}\n"
        f = _write(tmp_path, "t.tsv", content)
        records = load_csv(f, verbose=False)
        assert len(records) == 1


class TestLoadCsvEdgeCases:
    def test_windows_crlf(self, tmp_path: Path) -> None:
        content = f"name,sequence\r\nAG1,{SEQ_A}\r\nAG2,{SEQ_B}\r\n"
        f = _write(tmp_path, "w.csv", content)
        records = load_csv(f, verbose=False)
        assert len(records) == 2

    def test_blank_leading_rows(self, tmp_path: Path) -> None:
        content = f"\n\n\nname,sequence\nAG1,{SEQ_A}\n"
        f = _write(tmp_path, "b.csv", content)
        records = load_csv(f, verbose=False)
        assert len(records) == 1
        assert records[0][0] == "AG1"

    def test_bom(self, tmp_path: Path) -> None:
        content = f"name,sequence\nAG1,{SEQ_A}\n"
        f = _write_bytes(tmp_path, "bom.csv", content.encode("utf-8-sig"))
        records = load_csv(f, verbose=False)
        assert len(records) == 1

    def test_duplicate_ids_renamed(self, tmp_path: Path) -> None:
        content = f"name,sequence\nAG1,{SEQ_A}\nAG1,{SEQ_B}\n"
        f = _write(tmp_path, "d.csv", content)
        records = load_csv(f, verbose=False)
        assert len(records) == 2
        ids = [r[0] for r in records]
        assert len(set(ids)) == 2  # all unique

    def test_invalid_sequences_skipped(self, tmp_path: Path) -> None:
        content = "name,sequence\nAG1,SHORT\nAG2," + SEQ_A + "\n"
        f = _write(tmp_path, "inv.csv", content)
        records = load_csv(f, verbose=False)
        assert len(records) == 1
        assert records[0][0] == "AG2"

    def test_blank_rows_in_body_skipped(self, tmp_path: Path) -> None:
        content = f"name,sequence\nAG1,{SEQ_A}\n\n\nAG2,{SEQ_B}\n"
        f = _write(tmp_path, "gaps.csv", content)
        records = load_csv(f, verbose=False)
        assert len(records) == 2

    def test_file_not_found(self) -> None:
        with pytest.raises(FileNotFoundError):
            load_csv("/nonexistent/path.csv", verbose=False)

    def test_explicit_delimiter_overrides_sniffer(self, tmp_path: Path) -> None:
        content = f"name\tsequence\nAG1\t{SEQ_A}\n"
        f = _write(tmp_path, "t.csv", content)
        records = load_csv(f, delimiter="\t", verbose=False)
        assert len(records) == 1

    def test_empty_file_returns_empty(self, tmp_path: Path) -> None:
        f = _write(tmp_path, "empty.csv", "")
        records = load_csv(f, verbose=False)
        assert records == []

    def test_sanity_warning_on_no_valid_seqs(self, tmp_path: Path, capsys) -> None:
        """Non-empty file that yields 0 valid records triggers stderr warning."""
        content = "name,sequence\nAG1,XX\n"
        f = _write(tmp_path, "bad.csv", content)
        records = load_csv(f, verbose=False)
        assert records == []
        captured = capsys.readouterr()
        assert "no valid sequences" in captured.err.lower()


# --- FASTA loading ---


class TestLoadFasta:
    def test_basic_fasta(self, tmp_path: Path) -> None:
        content = f">AG1 some description\n{SEQ_A}\n>AG2\n{SEQ_B}\n"
        f = _write(tmp_path, "seqs.fasta", content)
        records = load_fasta(f, verbose=False)
        assert len(records) == 2
        assert records[0][0] == "AG1"
        assert records[1][0] == "AG2"

    def test_multiline_sequence(self, tmp_path: Path) -> None:
        half = len(SEQ_A) // 2
        content = f">AG1\n{SEQ_A[:half]}\n{SEQ_A[half:]}\n"
        f = _write(tmp_path, "multi.fasta", content)
        records = load_fasta(f, verbose=False)
        assert len(records) == 1
        assert records[0][1] == clean_sequence(SEQ_A)

    def test_skips_short_sequences(self, tmp_path: Path) -> None:
        content = ">SHORT\nACD\n>GOOD\n" + SEQ_A + "\n"
        f = _write(tmp_path, "short.fasta", content)
        records = load_fasta(f, verbose=False, min_len=8)
        assert len(records) == 1
        assert records[0][0] == "GOOD"

    def test_duplicate_ids(self, tmp_path: Path) -> None:
        content = f">AG1\n{SEQ_A}\n>AG1\n{SEQ_B}\n"
        f = _write(tmp_path, "dup.fasta", content)
        records = load_fasta(f, verbose=False)
        assert len(records) == 2
        ids = [r[0] for r in records]
        assert len(set(ids)) == 2

    def test_empty_fasta(self, tmp_path: Path) -> None:
        f = _write(tmp_path, "empty.fasta", "")
        records = load_fasta(f, verbose=False)
        assert records == []

    def test_fa_extension(self, tmp_path: Path) -> None:
        content = f">AG1\n{SEQ_A}\n"
        f = _write(tmp_path, "seqs.fa", content)
        records = load_fasta(f, verbose=False)
        assert len(records) == 1

    def test_fasta_with_bom(self, tmp_path: Path) -> None:
        content = f">AG1\n{SEQ_A}\n"
        f = _write_bytes(tmp_path, "bom.fasta", content.encode("utf-8-sig"))
        records = load_fasta(f, verbose=False)
        assert len(records) == 1

    def test_max_len(self, tmp_path: Path) -> None:
        content = f">AG1\n{SEQ_A}\n"
        f = _write(tmp_path, "long.fasta", content)
        records = load_fasta(f, verbose=False, max_len=10)
        assert records == []

    def test_file_not_found(self) -> None:
        with pytest.raises(FileNotFoundError):
            load_fasta("/nonexistent/seqs.fasta", verbose=False)


# --- Unified loader ---


class TestLoadSequences:
    def test_dispatches_to_csv(self, tmp_path: Path) -> None:
        f = _write(tmp_path, "data.csv", f"name,sequence\nAG1,{SEQ_A}\n")
        records = load_sequences(f, verbose=False)
        assert len(records) == 1

    def test_dispatches_to_fasta(self, tmp_path: Path) -> None:
        f = _write(tmp_path, "data.fasta", f">AG1\n{SEQ_A}\n")
        records = load_sequences(f, verbose=False)
        assert len(records) == 1

    def test_dispatches_fa_extension(self, tmp_path: Path) -> None:
        f = _write(tmp_path, "data.fa", f">AG1\n{SEQ_A}\n")
        records = load_sequences(f, verbose=False)
        assert len(records) == 1

    def test_dispatches_faa_extension(self, tmp_path: Path) -> None:
        f = _write(tmp_path, "data.faa", f">AG1\n{SEQ_A}\n")
        records = load_sequences(f, verbose=False)
        assert len(records) == 1


# --- Iterator ---


class TestIterCsv:
    def test_basic(self, tmp_path: Path) -> None:
        f = _write(tmp_path, "i.csv", f"name,sequence\nAG1,{SEQ_A}\nAG2,{SEQ_B}\n")
        results = list(iter_csv(f))
        assert len(results) == 2

    def test_semicolon(self, tmp_path: Path) -> None:
        f = _write(tmp_path, "i.csv", f"name;sequence\nAG1;{SEQ_A}\n")
        results = list(iter_csv(f))
        assert len(results) == 1

    def test_blank_leading_rows(self, tmp_path: Path) -> None:
        f = _write(tmp_path, "i.csv", f"\n\nname,sequence\nAG1,{SEQ_A}\n")
        results = list(iter_csv(f))
        assert len(results) == 1


# --- Utility functions ---


class TestCleanSequence:
    def test_strips_non_aa(self) -> None:
        assert clean_sequence("  acD-E*FG 123 ") == "ACDEFG"

    def test_uppercase(self) -> None:
        assert clean_sequence("acdefgh") == "ACDEFGH"


class TestIsValidSequence:
    def test_valid(self) -> None:
        assert is_valid_sequence("ACDEFGHIKLMNPQ", min_len=8)

    def test_too_short(self) -> None:
        assert not is_valid_sequence("ACD", min_len=8)

    def test_invalid_char(self) -> None:
        assert not is_valid_sequence("ACDEFGHX", min_len=8)


class TestDetectColumns:
    def test_standard_headers(self) -> None:
        id_idx, seq_idx = detect_columns(["name", "sequence"])
        assert id_idx == 0
        assert seq_idx == 1

    def test_with_empty_first_col(self) -> None:
        id_idx, seq_idx = detect_columns(["", "name", "sequence"])
        assert id_idx == 1
        assert seq_idx == 2

    def test_two_col_fallback(self) -> None:
        id_idx, seq_idx = detect_columns(["foo", "bar"])
        assert id_idx == 0
        assert seq_idx == 1

    def test_unrecognized_raises(self) -> None:
        with pytest.raises(ValueError, match="Could not detect"):
            detect_columns(["foo", "bar", "baz"])
