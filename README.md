# Antigen Cross-Reactivity Screener

Find antigens in your database that share enough sequence similarity with
any other antigen to potentially be recognised by the same antibody.

Designed for purification troubleshooting: when you're struggling to produce
a functional antibody against a target, find alternative antigens from your
existing database that might act as suitable capture column substitutes.

---

## How it works

```
All antigens (CSV)
       │
       ▼
Strip His6-ABP N-terminal tag (Smith-Waterman, robust to minor variations)
       │
       ▼
Build k-mer inverted index (k=6, Jaccard similarity pre-filter)
       │
       ▼
For each antigen: k-mer filter → ~handful of candidates
       │
       ▼
Smith-Waterman local alignment (BLOSUM62, parasail SIMD)
       │
       ▼
Store all pairs above threshold in SQLite (resumable)
       │
       ▼
Query anytime — instant lookup from precomputed results
```

---

## Installation

```bash
pip install -r requirements.txt
```

Requires Python 3.11+.

---

## Quick start

### 1. Set your His6-ABP tag

Open `antigen_screener/tag_stripper.py` and update `DEFAULT_TAG`:

```python
DEFAULT_TAG = "MHHHHHHGSSG"   # ← replace with your real His6-ABP sequence
```

### 2. Prepare your CSV

Your file needs at minimum two columns — one for IDs and one for sequences.
Column names are auto-detected. Supported names:

- **ID column:** `id`, `antigen_id`, `name`, `accession`, `gene`, `protein_id`, `entry`
- **Sequence column:** `sequence`, `seq`, `aa_sequence`, `protein_sequence`, `peptide`

If your column names differ, use `--id-col` and `--seq-col` flags.

Example format:
```
antigen_id,sequence
AG_001,MHHHHHHGSSGARKLVTPQIYWDG...
AG_002,MHHHHHHGSSGPVALKEQRTMWDNF...
```

### 3. Run the precomputation

```bash
python -m antigen_screener run --input antigens.csv --db results.db
```

With a custom tag:
```bash
python -m antigen_screener run \
  --input antigens.csv \
  --db results.db \
  --tag MHHHHHHGSSG
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
where it left off.

### 4. Query results

Look up precomputed similar antigens for a known ID:
```bash
python -m antigen_screener query --db results.db --id AG_001
```

Run a live alignment for a new sequence:
```bash
python -m antigen_screener query --db results.db --seq MHHHHHHGSSGAKLTPVQIYWDG
```

### 5. Export to CSV

```bash
python -m antigen_screener export --db results.db --output results.csv
```

Filter to only high-confidence pairs:
```bash
python -m antigen_screener export --db results.db --output results.csv --min-score 2.0
```

### 6. Check database stats

```bash
python -m antigen_screener stats --db results.db
```

---

## Tuning parameters

| Parameter | Default | Effect |
|-----------|---------|--------|
| `--kmer-threshold` | 0.04 | Lower = more sensitive pre-filter, slower. Raise to 0.08+ for speed. |
| `--min-score` | 1.0 | Normalised SW score threshold. Lower = more results. Try 0.5–2.0. |
| `--min-aligned` | 8 | Minimum shared region in amino acids. 8aa ≈ one epitope. |
| `--matrix` | blosum62 | Use `blosum45` for short/distant homologs (<30% identity). |

### Normalised score guide

The score is `raw SW score / min(query_length, target_length)`.

| Score | Approximate meaning |
|-------|---------------------|
| > 3.0 | Very high similarity — almost certainly same epitope |
| 1.5–3.0 | Strong similarity — likely cross-reactive |
| 1.0–1.5 | Moderate — worth investigating |
| < 1.0 | Weak — probably coincidental |

---

## Testing

Generate synthetic test data (500 antigens, 5 planted similarity families):
```bash
python generate_test_data.py
python -m antigen_screener run --input data/test_antigens.csv --db data/test.db --tag MHHHHHHGSSG
python -m antigen_screener query --db data/test.db --id FAM00_VAR00
```

You should see other `FAM00_*` variants at the top of the results.

---

## Output CSV columns

| Column | Description |
|--------|-------------|
| `query_id` | First antigen in the pair |
| `target_id` | Second antigen in the pair |
| `normalized_score` | SW score / min sequence length (main ranking metric) |
| `raw_score` | Raw Smith-Waterman score |
| `aligned_region_len` | Length of the locally aligned region (aa) |
| `query_length` | Length of query stripped sequence |
| `target_length` | Length of target stripped sequence |
| `query_stripped` | Query sequence after tag removal |
| `target_stripped` | Target sequence after tag removal |

---

## File structure

```
antigen_screener/
├── antigen_screener/
│   ├── __init__.py        # Package
│   ├── __main__.py        # CLI entry point
│   ├── tag_stripper.py    # His6-ABP tag detection and removal
│   ├── kmer_filter.py     # K-mer pre-filter (Jaccard similarity)
│   ├── aligner.py         # Smith-Waterman alignment (parasail)
│   ├── db_loader.py       # CSV/TSV loading and validation
│   ├── store.py           # SQLite persistence layer
│   └── pipeline.py        # Full pipeline orchestrator
├── generate_test_data.py  # Synthetic test data generator
├── requirements.txt
└── README.md
```
