"""
Database loader module.
Reads antigen CSV/TSV files, validates sequences, and prepares records
for indexing. Flexible column detection handles varied input formats.
"""

import csv
import re
import sys
from pathlib import Path
from typing import Iterator

# Valid single-letter amino acid codes
VALID_AA = set("ACDEFGHIKLMNPQRSTVWY")

# Common column name variants (case-insensitive)
ID_COLS  = {"id", "antigen_id", "name", "accession", "gene", "protein_id", "entry"}
SEQ_COLS = {"sequence", "seq", "aa_sequence", "protein_sequence", "peptide", "amino_acid_sequence"}


def detect_columns(headers: list[str]) -> tuple[int, int]:
    """
    Auto-detect ID and sequence columns from header names.
    Returns (id_col_index, seq_col_index).
    Raises ValueError if columns can't be determined.
    """
    headers_lower = [h.strip().lower() for h in headers]

    id_idx = None
    seq_idx = None

    for i, h in enumerate(headers_lower):
        if h in ID_COLS and id_idx is None:
            id_idx = i
        if h in SEQ_COLS and seq_idx is None:
            seq_idx = i

    # Fallback: if only 2 columns, assume first=id, second=sequence
    if id_idx is None and seq_idx is None and len(headers) == 2:
        id_idx, seq_idx = 0, 1

    if id_idx is None or seq_idx is None:
        raise ValueError(
            f"Could not detect ID and sequence columns from headers: {headers}\n"
            f"Expected one of {ID_COLS} for ID and {SEQ_COLS} for sequence.\n"
            f"Use --id-col and --seq-col to specify column names manually."
        )

    return id_idx, seq_idx


def clean_sequence(seq: str) -> str:
    """Strip whitespace, convert to uppercase, remove non-AA characters."""
    seq = seq.strip().upper()
    seq = re.sub(r"[^ACDEFGHIKLMNPQRSTVWY]", "", seq)
    return seq


def is_valid_sequence(seq: str, min_len: int = 8) -> bool:
    return len(seq) >= min_len and all(c in VALID_AA for c in seq)


def load_csv(
    filepath: str | Path,
    id_col: str | None = None,
    seq_col: str | None = None,
    delimiter: str | None = None,
    min_len: int = 8,
    max_len: int | None = None,
    verbose: bool = True,
) -> list[tuple[str, str]]:
    """
    Load antigen sequences from a CSV or TSV file.

    Args:
        filepath:  Path to the input file
        id_col:    Column name for antigen IDs (auto-detected if None)
        seq_col:   Column name for sequences (auto-detected if None)
        delimiter: ',' or '\\t' — auto-detected from extension if None
        min_len:   Minimum sequence length to include
        max_len:   Maximum sequence length (None = no limit)
        verbose:   Print loading summary

    Returns:
        List of (antigen_id, sequence) tuples, deduplicated by ID.
    """
    filepath = Path(filepath)
    if not filepath.exists():
        raise FileNotFoundError(f"Input file not found: {filepath}")

    # Auto-detect delimiter
    if delimiter is None:
        delimiter = "\t" if filepath.suffix.lower() in {".tsv", ".tab"} else ","

    records = []
    seen_ids: set[str] = set()
    skipped_invalid = 0
    skipped_duplicate = 0

    with open(filepath, newline="", encoding="utf-8-sig") as f:
        reader = csv.reader(f, delimiter=delimiter)
        headers = next(reader)

        # Resolve column indices
        if id_col and seq_col:
            headers_lower = [h.strip().lower() for h in headers]
            try:
                id_idx  = headers_lower.index(id_col.lower())
                seq_idx = headers_lower.index(seq_col.lower())
            except ValueError as e:
                raise ValueError(f"Specified column not found: {e}")
        else:
            id_idx, seq_idx = detect_columns(headers)

        if verbose:
            print(f"  Detected columns — ID: '{headers[id_idx]}', Sequence: '{headers[seq_idx]}'")

        for row_num, row in enumerate(reader, start=2):
            if len(row) <= max(id_idx, seq_idx):
                continue  # malformed row

            antigen_id = row[id_idx].strip()
            raw_seq    = row[seq_idx]

            if not antigen_id:
                antigen_id = f"antigen_{row_num}"

            seq = clean_sequence(raw_seq)

            if not is_valid_sequence(seq, min_len):
                skipped_invalid += 1
                continue

            if max_len and len(seq) > max_len:
                skipped_invalid += 1
                continue

            if antigen_id in seen_ids:
                # Make ID unique by appending row number
                antigen_id = f"{antigen_id}_{row_num}"
                skipped_duplicate += 1

            seen_ids.add(antigen_id)
            records.append((antigen_id, seq))

    if verbose:
        print(f"  Loaded {len(records):,} sequences")
        if skipped_invalid:
            print(f"  Skipped {skipped_invalid:,} invalid/too-short sequences")
        if skipped_duplicate:
            print(f"  Renamed {skipped_duplicate:,} duplicate IDs")

    return records


def iter_csv(
    filepath: str | Path,
    id_col: str | None = None,
    seq_col: str | None = None,
    delimiter: str | None = None,
) -> Iterator[tuple[str, str]]:
    """Memory-efficient iterator version of load_csv (no filtering)."""
    filepath = Path(filepath)
    if delimiter is None:
        delimiter = "\t" if filepath.suffix.lower() in {".tsv", ".tab"} else ","

    with open(filepath, newline="", encoding="utf-8-sig") as f:
        reader = csv.reader(f, delimiter=delimiter)
        headers = next(reader)

        if id_col and seq_col:
            headers_lower = [h.strip().lower() for h in headers]
            id_idx  = headers_lower.index(id_col.lower())
            seq_idx = headers_lower.index(seq_col.lower())
        else:
            id_idx, seq_idx = detect_columns(headers)

        for row in reader:
            if len(row) <= max(id_idx, seq_idx):
                continue
            yield row[id_idx].strip(), clean_sequence(row[seq_idx])
