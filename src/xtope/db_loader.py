"""
Database loader module.
Reads antigen CSV/TSV/FASTA files, validates sequences, and prepares records
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

# Delimiters that csv.Sniffer is allowed to pick
ACCEPTED_DELIMITERS = {",", ";", "\t"}


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


def sniff_delimiter(filepath: Path) -> str:
    """
    Detect the delimiter of a CSV/TSV file using csv.Sniffer.

    Reads up to 8192 bytes (skipping blank leading lines) and sniffs
    the delimiter.  Falls back to comma if Sniffer can't decide or
    picks something outside the accepted set (comma, semicolon, tab).
    """
    with open(filepath, newline="", encoding="utf-8-sig") as f:
        sample_lines: list[str] = []
        chars_read = 0
        for line in f:
            if chars_read > 8192:
                break
            # skip blank leading lines when building the sample
            if not sample_lines and line.strip() == "":
                continue
            sample_lines.append(line)
            chars_read += len(line)

    sample = "".join(sample_lines)
    if not sample.strip():
        return ","  # empty file, delimiter doesn't matter

    try:
        dialect = csv.Sniffer().sniff(sample, delimiters=",;\t")
        if dialect.delimiter in ACCEPTED_DELIMITERS:
            return dialect.delimiter
    except csv.Error:
        pass

    # Heuristic fallback: count occurrences in the first non-blank line
    first_line = sample_lines[0] if sample_lines else ""
    counts = {d: first_line.count(d) for d in ACCEPTED_DELIMITERS}
    best = max(counts, key=counts.get)  # type: ignore[arg-type]
    if counts[best] > 0:
        return best

    return ","


def _skip_blank_leading_rows(f) -> str | None:
    """
    Advance a file object past any blank leading lines.
    Returns the first non-blank line, or None if the file is empty.
    """
    for line in f:
        if line.strip():
            return line
    return None


def clean_sequence(seq: str) -> str:
    """Strip whitespace, convert to uppercase, remove non-AA characters."""
    seq = seq.strip().upper()
    seq = re.sub(r"[^ACDEFGHIKLMNPQRSTVWY]", "", seq)
    return seq


def is_valid_sequence(seq: str, min_len: int = 8) -> bool:
    """Check that a sequence has only valid AA characters and meets minimum length."""
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
        filepath:  Path to the input file.
        id_col:    Column name for antigen IDs (auto-detected if None).
        seq_col:   Column name for sequences (auto-detected if None).
        delimiter: ',', ';', or '\\t' — auto-detected with csv.Sniffer if None.
        min_len:   Minimum sequence length to include.
        max_len:   Maximum sequence length (None = no limit).
        verbose:   Print loading summary.

    Returns:
        List of (antigen_id, sequence) tuples, deduplicated by ID.
    """
    filepath = Path(filepath)
    if not filepath.exists():
        raise FileNotFoundError(f"Input file not found: {filepath}")

    # Auto-detect delimiter using Sniffer (replaces extension-only check)
    if delimiter is None:
        delimiter = sniff_delimiter(filepath)

    if verbose:
        delim_name = {",": "comma", ";": "semicolon", "\t": "tab"}.get(delimiter, repr(delimiter))
        print(f"  Detected delimiter: {delim_name}")

    records: list[tuple[str, str]] = []
    seen_ids: set[str] = set()
    skipped_invalid = 0
    skipped_duplicate = 0

    with open(filepath, newline="", encoding="utf-8-sig") as f:
        # Skip blank leading rows before the header
        first_line = _skip_blank_leading_rows(f)
        if first_line is None:
            if verbose:
                print("  Warning: file is empty")
            return []

        # Parse the first non-blank line as the header
        # We need to re-wrap remaining lines with the header prepended
        import itertools
        remaining = itertools.chain([first_line], f)
        reader = csv.reader(remaining, delimiter=delimiter)
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
            print(f"  Detected columns — ID: '{headers[id_idx].strip()}', Sequence: '{headers[seq_idx].strip()}'")

        for row_num, row in enumerate(reader, start=2):
            # Skip blank rows in the data body
            if not any(cell.strip() for cell in row):
                continue

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

    # Sanity check: warn if no records were loaded from a non-empty file
    if not records and filepath.stat().st_size > 0:
        sys.stderr.write(
            f"Warning: no valid sequences loaded from {filepath}. "
            f"Check that the file format, delimiter, and column names are correct.\n"
        )

    return records


def load_fasta(
    filepath: str | Path,
    min_len: int = 8,
    max_len: int | None = None,
    verbose: bool = True,
) -> list[tuple[str, str]]:
    """
    Load antigen sequences from a FASTA file.

    Parses standard FASTA format where header lines start with '>'
    and the ID is the first whitespace-delimited token after '>'.

    Args:
        filepath: Path to the .fasta / .fa / .faa file.
        min_len:  Minimum sequence length to include.
        max_len:  Maximum sequence length (None = no limit).
        verbose:  Print loading summary.

    Returns:
        List of (antigen_id, sequence) tuples.
    """
    filepath = Path(filepath)
    if not filepath.exists():
        raise FileNotFoundError(f"Input file not found: {filepath}")

    records: list[tuple[str, str]] = []
    seen_ids: set[str] = set()
    skipped_invalid = 0
    skipped_duplicate = 0

    current_id: str | None = None
    current_seq_parts: list[str] = []

    def _flush() -> None:
        nonlocal current_id, current_seq_parts, skipped_invalid, skipped_duplicate
        if current_id is None:
            return
        seq = clean_sequence("".join(current_seq_parts))
        current_seq_parts = []

        if not is_valid_sequence(seq, min_len):
            skipped_invalid += 1
            current_id = None
            return

        if max_len and len(seq) > max_len:
            skipped_invalid += 1
            current_id = None
            return

        aid = current_id
        if aid in seen_ids:
            skipped_duplicate += 1
            aid = f"{aid}_dup{skipped_duplicate}"
        seen_ids.add(aid)
        records.append((aid, seq))
        current_id = None

    with open(filepath, encoding="utf-8-sig") as f:
        for line in f:
            line = line.rstrip("\n\r")
            if line.startswith(">"):
                # Flush previous record
                _flush()
                # Parse header: >ID description...
                header = line[1:].strip()
                current_id = header.split()[0] if header else None
                current_seq_parts = []
            elif current_id is not None:
                # Sequence continuation line
                current_seq_parts.append(line.strip())

        # Flush the last record
        _flush()

    if verbose:
        print(f"  Loaded {len(records):,} sequences from FASTA")
        if skipped_invalid:
            print(f"  Skipped {skipped_invalid:,} invalid/too-short sequences")
        if skipped_duplicate:
            print(f"  Renamed {skipped_duplicate:,} duplicate IDs")

    if not records and filepath.stat().st_size > 0:
        sys.stderr.write(
            f"Warning: no valid sequences loaded from {filepath}. "
            f"Check that the file is in valid FASTA format.\n"
        )

    return records


# File extensions recognized as FASTA
FASTA_EXTENSIONS = {".fasta", ".fa", ".faa", ".fas", ".fna"}


def load_sequences(
    filepath: str | Path,
    id_col: str | None = None,
    seq_col: str | None = None,
    delimiter: str | None = None,
    min_len: int = 8,
    max_len: int | None = None,
    verbose: bool = True,
) -> list[tuple[str, str]]:
    """
    Load antigen sequences, auto-detecting file format.

    Dispatches to load_fasta() for FASTA files (.fasta, .fa, .faa, .fas, .fna)
    and load_csv() for everything else (CSV, TSV, etc.).

    Args:
        filepath:  Path to the input file.
        id_col:    Column name for antigen IDs (CSV only, auto-detected if None).
        seq_col:   Column name for sequences (CSV only, auto-detected if None).
        delimiter: Delimiter override (CSV only, auto-detected if None).
        min_len:   Minimum sequence length to include.
        max_len:   Maximum sequence length (None = no limit).
        verbose:   Print loading summary.

    Returns:
        List of (antigen_id, sequence) tuples.
    """
    filepath = Path(filepath)
    if filepath.suffix.lower() in FASTA_EXTENSIONS:
        return load_fasta(filepath, min_len=min_len, max_len=max_len, verbose=verbose)
    return load_csv(
        filepath, id_col=id_col, seq_col=seq_col,
        delimiter=delimiter, min_len=min_len, max_len=max_len, verbose=verbose,
    )


def iter_csv(
    filepath: str | Path,
    id_col: str | None = None,
    seq_col: str | None = None,
    delimiter: str | None = None,
) -> Iterator[tuple[str, str]]:
    """Memory-efficient iterator version of load_csv (no filtering)."""
    filepath = Path(filepath)
    if delimiter is None:
        delimiter = sniff_delimiter(filepath)

    with open(filepath, newline="", encoding="utf-8-sig") as f:
        first_line = _skip_blank_leading_rows(f)
        if first_line is None:
            return

        import itertools
        remaining = itertools.chain([first_line], f)
        reader = csv.reader(remaining, delimiter=delimiter)
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
