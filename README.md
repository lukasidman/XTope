```
  ██╗  ██╗████████╗ ██████╗ ██████╗ ███████╗
  ╚██╗██╔╝╚══██╔══╝██╔═══██╗██╔══██╗██╔════╝
   ╚███╔╝    ██║   ██║   ██║██████╔╝█████╗
   ██╔██╗    ██║   ██║   ██║██╔═══╝ ██╔══╝
  ██╔╝ ██╗   ██║   ╚██████╔╝██║     ███████╗
  ╚═╝  ╚═╝   ╚═╝    ╚═════╝ ╚═╝     ╚══════╝

  Cross-Reactivity Antigen Screener
```

Find antigens in your database that share enough sequence similarity with
any other antigen to potentially be recognised by the same antibody.

Designed for purification troubleshooting: when you're struggling to produce
a functional antibody against a target, find alternative antigens from your
existing database that might act as suitable capture column substitutes.

---

## How it works

```
1. Prepare your antigen database (CSV / TSV / FASTA)
       │
       ▼
2. (Optional) Strip N-terminal purification tag — opt-in via --tag
       │
       ▼
3. Run all-vs-all alignment
   ┌──────────────────────────────────────────────────┐
   │  vectorized (default)      kmer                  │
   │  RA diagonal pre-filter →  Jaccard pre-filter →  │
   │  Batched NumPy SW          Per-pair SW            │
   └──────────────────────────────────────────────────┘
       │
       ▼
4. Score with E-values & bit-scores (Karlin-Altschul statistics)
       │
       ▼
5. Store significant pairs in SQLite (resumable — safe to interrupt)
       │
       ▼
6. Query anytime — instant lookup from precomputed results
       │
       ▼
7. Export to CSV for downstream analysis
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

### 1. Prepare your input file

The tool accepts **CSV**, **TSV**, and **FASTA** files. Delimiters (comma, semicolon, tab) are auto-detected using `csv.Sniffer` — no configuration needed.

Your file needs at minimum two columns — one for IDs and one for sequences.
Column names are auto-detected. Supported names:

- **ID column:** `id`, `antigen_id`, `name`, `accession`, `gene`, `protein_id`, `entry`
- **Sequence column:** `sequence`, `seq`, `aa_sequence`, `protein_sequence`, `peptide`

If your column names differ, use `--id-col` and `--seq-col` flags.

Example CSV:
```
antigen_id,sequence
AG_001,ARKLVTPQIYWDG...
AG_002,PVALKEQRTMWDNF...
```

Example semicolon-delimited (auto-detected):
```
;name;sequence
AG_001;ARKLVTPQIYWDG...
AG_002;PVALKEQRTMWDNF...
```

Example FASTA:
```
>AG_001
ARKLVTPQIYWDG...
>AG_002
PVALKEQRTMWDNF...
```

### 2. N-terminal tag stripping (optional)

N-tag stripping is **disabled by default** — sequences are used exactly as provided.

If your sequences carry an N-terminal purification tag that you want removed before
alignment, pass it via `--tag`:

```bash
python -m xtope run --input antigens.csv --db results.db --tag MHHHHHHGSSG
```

The tool will try exact string matching first, then fall back to Smith-Waterman
alignment within the first 60 residues if the exact match fails. At startup it
reports how many sequences the tag was found in.

If you don't pass `--tag`, the tool logs `N-tag stripping: disabled` and proceeds
with full sequences — appropriate when sequences are already tag-free.

### 3. Run the precomputation

**Default (vectorized NumPy + RA diagonal pre-filter):**

```bash
python -m xtope run --input antigens.csv --db results.db
```

This runs batched NumPy Smith-Waterman alignment with a reduced-alphabet diagonal pre-filter.
The pre-filter requires contiguous local similarity before a pair reaches SW — this is biologically
appropriate for epitope-level cross-reactivity (antibodies bind contiguous regions, not gapped matches).

**With N-tag stripping:**

```bash
python -m xtope run --input antigens.csv --db results.db --tag MHHHHHHGSSG
```

**Exhaustive mode — no pre-filter, every pair scored:**

```bash
python -m xtope run --input antigens.csv --db results.db --no-prefilter
```

Scores all pairs including those with only gapped or scattered similarity. Slower on large datasets;
useful when you want to be certain nothing is missed.

**Quick screen — k-mer filter + per-pair SW:**

```bash
python -m xtope run --input antigens.csv --db results.db --backend kmer
```

Fastest option for very large datasets. Uses a two-pass k-mer index (Jaccard + chain detection)
to eliminate most pairs before alignment. May miss some matches the vectorized backend finds.

With custom column names:
```bash
python -m xtope run \
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
python -m xtope query --db results.db --id AG_001
```

Run a live alignment for a new sequence not in the database:
```bash
python -m xtope query --db results.db --seq AKLTPVQIYWDG
```

Strip an N-tag from the query sequence in live mode:
```bash
python -m xtope query --db results.db --seq MHHHHHHGSSGAKLTPVQIYWDG --tag MHHHHHHGSSG
```

Control the number of results and E-value cutoff:
```bash
python -m xtope query --db results.db --id AG_001 --top-n 50 --max-evalue 0.001
```

### 5. Export to CSV

```bash
python -m xtope export --db results.db --output results.csv
```

Filter to only high-confidence pairs:
```bash
python -m xtope export --db results.db --output results.csv --max-evalue 0.001
```

### 6. Check database stats

```bash
python -m xtope stats --db results.db
```

### 7. Score interpretation guide

```bash
python -m xtope help-scores
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
| `--tag` | *(disabled)* | N-terminal tag sequence to strip before alignment (opt-in) |
| `--id-col` | *(auto-detect)* | Column name for antigen IDs |
| `--seq-col` | *(auto-detect)* | Column name for sequences |
| `--kmer-threshold` | `0.04` | Jaccard threshold for k-mer pre-filter (kmer backend only) |
| `--max-evalue` | `0.01` | Maximum E-value to record a pair |
| `--min-aligned` | `8` | Minimum aligned region length in amino acids |
| `--matrix` | `blosum62` | Substitution matrix (`blosum62`, `blosum45`, `blosum80`) |
| `--no-resume` | off | Start fresh even if the database has partial results |
| `--backend` | `vectorized` | `vectorized` (default, NumPy batched SW) or `kmer` (fast k-mer filter + per-pair SW) |
| `--no-prefilter` | off | Disable the RA diagonal pre-filter; score every pair exhaustively (vectorized only) |

### `query` flags

| Flag | Default | Description |
|------|---------|-------------|
| `--db` | *(required)* | Path to SQLite database |
| `--id` | — | Antigen ID to look up precomputed results |
| `--seq` | — | Raw sequence to run live alignment (mutually exclusive with `--id`) |
| `--tag` | *(disabled)* | N-terminal tag to strip from `--seq` in live mode (opt-in) |
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

Run `python -m xtope help-scores` for the full guide including Karlin-Altschul parameters.

---

## Backends

The tool offers two alignment backends, selectable via `--backend`:

**`vectorized`** (default) — Runs batched Smith-Waterman alignment using NumPy array operations. By default uses a reduced-alphabet diagonal pre-filter to skip pairs with no contiguous local similarity before they reach SW — biologically appropriate for antibody cross-reactivity, where epitopes are contiguous regions. Use `--no-prefilter` to disable the pre-filter and score every pair exhaustively.

**`kmer`** — Builds a k-mer inverted index (k=6) and uses a two-pass filter (Jaccard similarity + diagonal chain detection) to eliminate most pairs before running per-pair Smith-Waterman. Fastest option for very large datasets. Use `--backend kmer` for a quick screen when full sensitivity is less important than speed.

---

## Test data

A set of 2,500 synthetic test sequences is included at `data/test_antigens_2500.csv`. These contain no N-terminal tag and are ready to use directly:

```bash
python -m xtope run --input data/test_antigens_2500.csv --db test.db
python -m xtope query --db test.db --id FAM00_VAR00
python -m xtope stats --db test.db
```

The dataset includes 10 similarity families (15 variants each), 5 short shared-epitope groups, a physicochemically convergent family, mosaic domain-shuffled sequences, and ~2,300 random background sequences across 50–150 aa.

---

## Testing

```bash
# Install dev dependencies
pip install -e ".[dev]"

# Run all tests
pytest tests/ -v

# Run with coverage
pytest tests/ --cov=xtope --cov-report=term-missing

# Type check
mypy src/xtope/

# Lint
ruff check src/ tests/
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
| `query_length` | Length of query sequence (after tag removal if applicable) |
| `target_length` | Length of target sequence (after tag removal if applicable) |
| `query_stripped` | Query sequence used for alignment |
| `target_stripped` | Target sequence used for alignment |

---

## Project structure

```
xtope/
├── README.md
├── DEVELOPMENT.md                   # Development instructions and design decisions
├── pyproject.toml                   # Package config (setuptools, Python ≥3.10)
├── requirements.txt                 # Core: numpy
├── requirements-optional.txt        # Optional: openpyxl (Excel export)
├── requirements-dev.txt             # Dev: pytest, pytest-cov, mypy, ruff
├── data/
│   └── test_antigens_2500.csv       # 2,500 synthetic test sequences (no N-tag)
├── src/
│   └── xtope/
│       ├── __init__.py              # Package version
│       ├── __main__.py              # CLI: run, query, export, stats, help-scores
│       ├── pipeline.py              # All-vs-all orchestrator with resume/checkpoint
│       ├── tag_stripper.py          # N-tag detection and removal (opt-in)
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
