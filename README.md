# Antigen Cross-Reactivity Screener

Find antigens in your database that share enough sequence similarity with
any other antigen to potentially be recognised by the same antibody.

Designed for purification troubleshooting: when you're struggling to produce
a functional antibody against a target, find alternative antigens from your
existing database that might act as suitable capture column substitutes.

---

## How it works

```
All antigens (CSV / TSV / FASTA)
       │
       ▼
Strip His6-ABP N-terminal tag (exact match + Smith-Waterman fallback)
       │
       ▼
┌──────────────────────────────────────────────────┐
│  Choose one of two backends:                     │
│                                                  │
│  "kmer" (default)          "vectorized"          │
│  Build k-mer index  →      Batched NumPy         │
│  Jaccard pre-filter →      all-vs-all SW         │
│  Per-pair SW                (no pre-filter)       │
└──────────────────────────────────────────────────┘
       │
       ▼
Score with E-values & bit-scores (Karlin-Altschul statistics)
       │
       ▼
Store significant pairs in SQLite (resumable)
       │
       ▼
Query anytime — instant lookup from precomputed results
```

---

## Installation

Requires **Python 3.10+** and **numpy**.

```bash
# Editable install (recommended for development)
pip install -e .

# Or install core dependencies only
pip install -r requirements.txt

# Optional: Excel export support
pip install -r requirements-optional.txt

# Optional: development tools (pytest, ruff, mypy)
pip install -e ".[dev]"
```

---

## Quick start

### 1. Set your His6-ABP tag

Pass the tag via the `--tag` CLI flag. The default is `MHHHHHHGSSG`.

```bash
python -m antigen_screener run --input antigens.csv --db results.db --tag MHHHHHHGSSG
```

### 2. Prepare your input file

The tool accepts **CSV**, **TSV**, and **FASTA** files. Delimiters (comma, semicolon, tab) are auto-detected using `csv.Sniffer` — no configuration needed.

Your file needs at minimum two columns — one for IDs and one for sequences.
Column names are auto-detected. Supported names:

- **ID column:** `id`, `antigen_id`, `name`, `accession`, `gene`, `protein_id`, `entry`
- **Sequence column:** `sequence`, `seq`, `aa_sequence`, `protein_sequence`, `peptide`

If your column names differ, use `--id-col` and `--seq-col` flags.

Example CSV:
```
antigen_id,sequence
AG_001,MHHHHHHGSSGARKLVTPQIYWDG...
AG_002,MHHHHHHGSSGPVALKEQRTMWDNF...
```

Example semicolon-delimited (auto-detected):
```
;name;sequence
AG_001;MHHHHHHGSSGARKLVTPQIYWDG...
AG_002;MHHHHHHGSSGPVALKEQRTMWDNF...
```

Example FASTA:
```
>AG_001
MHHHHHHGSSGARKLVTPQIYWDG...
>AG_002
MHHHHHHGSSGPVALKEQRTMWDNF...
```

### 3. Run the precomputation

**Default backend (k-mer pre-filter + per-pair SW):**

```bash
python -m antigen_screener run --input antigens.csv --db results.db
```

**Vectorized backend (NumPy batched all-vs-all, no pre-filter — more sensitive, slower on large datasets):**

```bash
python -m antigen_screener run --input antigens.csv --db results.db --backend vectorized
```

With custom column names:
```bash
python -m antigen_screener run \
  --input antigens.tsv \
  --db results.db \
  --id-col gene_name \
  --seq-col aa_sequence
```

This will run for hours/days on 80,000 sequences — that's expected and fine.
**It is resumable**: if interrupted, re-run the same command and it picks up
where it left off. Use `--no-resume` to start fresh.

### 4. Query results

Look up precomputed similar antigens for a known ID:
```bash
python -m antigen_screener query --db results.db --id AG_001
```

Run a live alignment for a new sequence not in the database:
```bash
python -m antigen_screener query --db results.db --seq MHHHHHHGSSGAKLTPVQIYWDG
```

Control the number of results and E-value cutoff:
```bash
python -m antigen_screener query --db results.db --id AG_001 --top-n 50 --max-evalue 0.001
```

### 5. Export to CSV

```bash
python -m antigen_screener export --db results.db --output results.csv
```

Filter to only high-confidence pairs:
```bash
python -m antigen_screener export --db results.db --output results.csv --max-evalue 0.001
```

### 6. Check database stats

```bash
python -m antigen_screener stats --db results.db
```

### 7. Score interpretation guide

```bash
python -m antigen_screener help-scores
```

---

## Commands reference

| Command | Purpose |
|---------|---------|
| `run` | Run all-vs-all precomputation across the full database |
| `query` | Look up precomputed results by ID, or run live alignment for a new sequence |
| `export` | Export similarity pairs to CSV |
| `stats` | Show database statistics (antigen count, pairs, last run) |
| `help-scores` | Display E-value and bit-score interpretation guide |

### `run` flags

| Flag | Default | Description |
|------|---------|-------------|
| `--input` | *(required)* | Path to CSV, TSV, or FASTA antigen database |
| `--db` | *(required)* | Path to output SQLite database |
| `--tag` | `MHHHHHHGSSG` | His6-ABP tag sequence to strip |
| `--id-col` | *(auto-detect)* | Column name for antigen IDs |
| `--seq-col` | *(auto-detect)* | Column name for sequences |
| `--kmer-threshold` | `0.04` | Jaccard threshold for k-mer pre-filter (kmer backend only) |
| `--max-evalue` | `0.01` | Maximum E-value to record a pair |
| `--min-aligned` | `8` | Minimum aligned region length in amino acids |
| `--matrix` | `blosum62` | Substitution matrix (`blosum62`, `blosum45`, `blosum80`) |
| `--no-resume` | off | Start fresh even if the database has partial results |
| `--backend` | `kmer` | `kmer` (k-mer filter + per-pair SW) or `vectorized` (NumPy batched all-vs-all) |

### `query` flags

| Flag | Default | Description |
|------|---------|-------------|
| `--db` | *(required)* | Path to SQLite database |
| `--id` | — | Antigen ID to look up precomputed results |
| `--seq` | — | Raw sequence to run live alignment (mutually exclusive with `--id`) |
| `--tag` | `MHHHHHHGSSG` | Tag to strip from `--seq` in live mode |
| `--max-evalue` | `0.01` | Maximum E-value threshold |
| `--min-aligned` | `8` | Minimum aligned region length |
| `--top-n` | `25` | Number of results to show |

### `export` flags

| Flag | Default | Description |
|------|---------|-------------|
| `--db` | *(required)* | Path to SQLite database |
| `--output` | *(db name).csv* | Output CSV path |
| `--max-evalue` | `10.0` | Maximum E-value to include |

---

## Scoring

This tool uses **E-values** and **bit-scores** based on the Karlin-Altschul statistical framework (the same approach used by BLAST and Diamond). These replace the earlier normalised score, which suffered from length bias on short sequences.

**E-value** — the number of alignments this good you'd expect by chance in a database of this size. Lower is better.

| E-value | Significance | Interpretation |
|---------|-------------|----------------|
| < 1e-10 | Very high | Almost certainly share an epitope |
| 1e-10 – 1e-4 | Strong | Likely cross-reactive |
| 1e-4 – 0.1 | Moderate | Worth investigating |
| 0.1 – 1.0 | Weak | Possibly coincidental |
| > 1.0 | Not significant | Probably random |

**Bit-score** — normalised alignment score independent of sequence length and database size. Higher is better.

| Bit-score | Interpretation |
|-----------|----------------|
| > 50 | Very strong similarity |
| 30 – 50 | Strong similarity |
| 20 – 30 | Moderate similarity |
| < 20 | Weak or no similarity |

Run `python -m antigen_screener help-scores` for the full guide including Karlin-Altschul parameters.

---

## Backends

The tool offers two alignment backends, selectable via `--backend`:

**`kmer`** (default) — Builds a k-mer inverted index (k=6) and uses Jaccard similarity as a pre-filter to eliminate >99% of sequence pairs before running Smith-Waterman alignment on the remaining candidates. Fast on large datasets but may miss pairs sharing only short epitopes (below the Jaccard threshold).

**`vectorized`** — Runs batched Smith-Waterman alignment on all pairs simultaneously using NumPy array operations. No pre-filter, so every pair is scored — more sensitive, but slower on very large datasets. Recommended for smaller databases or when sensitivity matters more than speed.

---

## Testing

```bash
# Install dev dependencies
pip install -e ".[dev]"

# Run all tests
pytest tests/ -v

# Run with coverage
pytest tests/ --cov=antigen_screener --cov-report=term-missing

# Type check
mypy src/antigen_screener/

# Lint
ruff check src/ tests/
```

Generate synthetic test data and run a quick pipeline test:
```bash
python -m antigen_screener run --input antigens.csv --db test.db --tag MHHHHHHGSSG
python -m antigen_screener query --db test.db --id AG1
python -m antigen_screener stats --db test.db
```

---

## Output CSV columns

| Column | Description |
|--------|-------------|
| `query_id` | First antigen in the pair |
| `target_id` | Second antigen in the pair |
| `evalue` | E-value (expected number of chance hits this good) |
| `bit_score` | Bit-score (normalised, length-independent) |
| `raw_score` | Raw Smith-Waterman score |
| `aligned_region_len` | Length of the locally aligned region (amino acids) |
| `query_length` | Length of query sequence after tag removal |
| `target_length` | Length of target sequence after tag removal |
| `query_stripped` | Query sequence after tag removal |
| `target_stripped` | Target sequence after tag removal |

---

## Project structure

```
antigen-screener/
├── README.md
├── CLAUDE.md                        # Development instructions
├── pyproject.toml                   # Package config (setuptools, Python ≥3.10)
├── requirements.txt                 # Core: numpy
├── requirements-optional.txt        # Optional: openpyxl (Excel export)
├── requirements-dev.txt             # Dev: pytest, pytest-cov, mypy, ruff
├── src/
│   └── antigen_screener/
│       ├── __init__.py              # Package version
│       ├── __main__.py              # CLI: run, query, export, stats, help-scores
│       ├── pipeline.py              # All-vs-all orchestrator with resume/checkpoint
│       ├── tag_stripper.py          # His6-ABP tag detection and removal
│       ├── kmer_filter.py           # K-mer inverted index, Jaccard pre-filter
│       ├── aligner.py               # Smith-Waterman alignment scoring
│       ├── sw_fallback.py           # Pure-Python SW implementation
│       ├── vectorized_sw.py         # Batched NumPy all-vs-all SW pipeline
│       ├── evalue.py                # E-value & bit-score (Karlin-Altschul)
│       ├── physicochemical.py       # Hydrophobicity, charge, pI scoring
│       ├── db_loader.py             # CSV / TSV / FASTA loader (auto-detect)
│       ├── store.py                 # SQLite results store (WAL mode)
│       └── generate_test_data.py    # Synthetic test data generator
└── tests/
    ├── conftest.py                  # Shared pytest fixtures
    ├── test_tag_stripper.py
    ├── test_kmer_filter.py
    ├── test_aligner.py
    ├── test_db_loader.py
    ├── test_store.py
    └── test_pipeline.py
```

---

## Dependencies

| Dependency | Required | Purpose |
|------------|----------|---------|
| `numpy` | Yes | Array operations, vectorized SW backend |
| `openpyxl` | No | Excel (.xlsx) export |
| `pytest` | Dev only | Test runner |
| `pytest-cov` | Dev only | Coverage reporting |
| `mypy` | Dev only | Type checking |
| `ruff` | Dev only | Linting |
