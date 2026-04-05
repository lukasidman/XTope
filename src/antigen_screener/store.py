"""
SQLite results store.
Persists all-vs-all precomputed similarity results for fast future queries.
Schema is append-friendly so partial runs can be resumed.
"""

import sqlite3
from pathlib import Path
from contextlib import contextmanager
from typing import Iterator


SCHEMA = """
CREATE TABLE IF NOT EXISTS antigens (
    id          TEXT PRIMARY KEY,
    sequence    TEXT NOT NULL,
    stripped    TEXT NOT NULL,
    tag_found   INTEGER NOT NULL,   -- 0 or 1
    length      INTEGER NOT NULL
);

CREATE TABLE IF NOT EXISTS similarities (
    query_id            TEXT NOT NULL,
    target_id           TEXT NOT NULL,
    raw_score           INTEGER NOT NULL,
    bit_score           REAL NOT NULL,
    evalue              REAL NOT NULL,
    query_length        INTEGER NOT NULL,
    target_length       INTEGER NOT NULL,
    aligned_region_len  INTEGER NOT NULL,
    PRIMARY KEY (query_id, target_id)
);

CREATE INDEX IF NOT EXISTS idx_sim_query  ON similarities(query_id);
CREATE INDEX IF NOT EXISTS idx_sim_target ON similarities(target_id);
CREATE INDEX IF NOT EXISTS idx_sim_evalue ON similarities(evalue ASC);
CREATE INDEX IF NOT EXISTS idx_sim_bitscore ON similarities(bit_score DESC);

CREATE TABLE IF NOT EXISTS run_metadata (
    key     TEXT PRIMARY KEY,
    value   TEXT
);
"""


class ResultsStore:
    def __init__(self, db_path: str | Path):
        self.db_path = Path(db_path)
        self.conn = sqlite3.connect(str(self.db_path), check_same_thread=False)
        self.conn.row_factory = sqlite3.Row
        self.conn.execute("PRAGMA journal_mode=WAL")
        self.conn.execute("PRAGMA synchronous=NORMAL")
        self.conn.execute("PRAGMA cache_size=-64000")   # 64MB cache
        self._init_schema()

    def _init_schema(self):
        self.conn.executescript(SCHEMA)
        # Migration: add columns for existing DBs created before E-value support
        self._migrate_add_columns()
        self.conn.commit()

    def _migrate_add_columns(self):
        """Add new columns to existing databases, handling schema evolution."""
        # Contiguous match columns (from vectorized backend)
        for col, coltype in [
            ("consec_match_len", "INTEGER DEFAULT 0"),
            ("consec_query_start", "INTEGER DEFAULT -1"),
            ("consec_query_end", "INTEGER DEFAULT -1"),
            ("consec_target_start", "INTEGER DEFAULT -1"),
            ("consec_target_end", "INTEGER DEFAULT -1"),
        ]:
            try:
                self.conn.execute(
                    f"ALTER TABLE similarities ADD COLUMN {col} {coltype}"
                )
            except sqlite3.OperationalError:
                pass  # column already exists

        # E-value migration: add bit_score and evalue to old DBs that only have
        # normalized_score. We add them with defaults so existing rows get values.
        for col, coltype in [
            ("bit_score", "REAL DEFAULT 0.0"),
            ("evalue", "REAL DEFAULT 999.0"),
        ]:
            try:
                self.conn.execute(
                    f"ALTER TABLE similarities ADD COLUMN {col} {coltype}"
                )
            except sqlite3.OperationalError:
                pass  # column already exists

        # Physicochemical similarity columns (Track 2 scoring).
        # NULL = not yet scored (pairs from pre-physicochemical runs).
        for col, coltype in [
            ("hydrophobicity_corr",      "REAL"),   # post-alignment KD correlation
            ("charge_corr",              "REAL"),   # post-alignment charge correlation
            ("pi_diff",                  "REAL"),   # |pI_a - pI_b|
            ("hydrophobicity_cross_corr","REAL"),   # KD sliding-window cross-corr
            ("charge_cross_corr",        "REAL"),   # charge sliding-window cross-corr
            ("binary_hydro_cross_corr",  "REAL"),   # binary H/P pattern cross-corr (primary CONV signal)
            ("pc_composite_score",       "REAL"),   # weighted composite (0–1)
            ("detection_source",         "TEXT"),   # 'sequence', 'physicochemical', or 'both'
        ]:
            try:
                self.conn.execute(
                    f"ALTER TABLE similarities ADD COLUMN {col} {coltype}"
                )
            except sqlite3.OperationalError:
                pass  # column already exists

    # ------------------------------------------------------------------ #
    #  Antigen records                                                     #
    # ------------------------------------------------------------------ #

    def upsert_antigen(self, antigen_id: str, sequence: str, stripped: str, tag_found: bool):
        self.conn.execute(
            """INSERT OR REPLACE INTO antigens (id, sequence, stripped, tag_found, length)
               VALUES (?, ?, ?, ?, ?)""",
            (antigen_id, sequence, stripped, int(tag_found), len(stripped)),
        )

    def upsert_antigens_batch(self, records: list[tuple]):
        """records: list of (id, sequence, stripped, tag_found)"""
        self.conn.executemany(
            """INSERT OR REPLACE INTO antigens (id, sequence, stripped, tag_found, length)
               VALUES (?, ?, ?, ?, ?)""",
            [(r[0], r[1], r[2], int(r[3]), len(r[2])) for r in records],
        )
        self.conn.commit()

    def get_all_stripped(self) -> list[tuple[str, str]]:
        """Return all (id, stripped_sequence) pairs."""
        rows = self.conn.execute("SELECT id, stripped FROM antigens ORDER BY id").fetchall()
        return [(r["id"], r["stripped"]) for r in rows]

    def antigen_count(self) -> int:
        return self.conn.execute("SELECT COUNT(*) FROM antigens").fetchone()[0]

    def total_residues(self) -> int:
        """Return total residue count across all stripped sequences."""
        row = self.conn.execute("SELECT COALESCE(SUM(length), 0) FROM antigens").fetchone()
        return row[0]

    # ------------------------------------------------------------------ #
    #  Similarity results                                                  #
    # ------------------------------------------------------------------ #

    def insert_similarities_batch(self, results: list[dict]):
        """Bulk insert similarity results. Ignores duplicates.

        Handles both old-style dicts (without contiguous match / physicochemical
        fields) and new-style dicts for backwards compatibility.
        """
        self.conn.executemany(
            """INSERT OR IGNORE INTO similarities
               (query_id, target_id, raw_score, bit_score, evalue,
                query_length, target_length, aligned_region_len,
                consec_match_len, consec_query_start, consec_query_end,
                consec_target_start, consec_target_end,
                hydrophobicity_corr, charge_corr, pi_diff,
                hydrophobicity_cross_corr, charge_cross_corr,
                binary_hydro_cross_corr, pc_composite_score, detection_source)
               VALUES (:query_id, :target_id, :raw_score, :bit_score, :evalue,
                       :query_length, :target_length, :aligned_region_len,
                       :consec_match_len, :consec_query_start, :consec_query_end,
                       :consec_target_start, :consec_target_end,
                       :hydrophobicity_corr, :charge_corr, :pi_diff,
                       :hydrophobicity_cross_corr, :charge_cross_corr,
                       :binary_hydro_cross_corr, :pc_composite_score,
                       :detection_source)""",
            [
                {
                    "consec_match_len":        r.get("consec_match_len", 0),
                    "consec_query_start":       r.get("consec_query_start", -1),
                    "consec_query_end":         r.get("consec_query_end", -1),
                    "consec_target_start":      r.get("consec_target_start", -1),
                    "consec_target_end":        r.get("consec_target_end", -1),
                    "bit_score":                r.get("bit_score", 0.0),
                    "evalue":                   r.get("evalue", 999.0),
                    "hydrophobicity_corr":      r.get("hydrophobicity_corr"),
                    "charge_corr":              r.get("charge_corr"),
                    "pi_diff":                  r.get("pi_diff"),
                    "hydrophobicity_cross_corr":r.get("hydrophobicity_cross_corr"),
                    "charge_cross_corr":        r.get("charge_cross_corr"),
                    "binary_hydro_cross_corr":  r.get("binary_hydro_cross_corr"),
                    "pc_composite_score":       r.get("pc_composite_score"),
                    "detection_source":         r.get("detection_source"),
                    **r,
                }
                for r in results
            ],
        )
        self.conn.commit()

    def query_similar(
        self,
        antigen_id: str,
        max_evalue: float = 0.01,
        top_n: int = 50,
    ) -> list[dict]:
        """
        Retrieve precomputed similar antigens for a given ID.
        Searches both directions (query->target and target->query), deduplicates.
        Results sorted by E-value ascending (best first).
        """
        rows = self.conn.execute(
            """SELECT
                CASE WHEN query_id = ? THEN target_id ELSE query_id END AS partner_id,
                raw_score, bit_score, evalue, query_length, target_length,
                aligned_region_len,
                hydrophobicity_corr, charge_corr, pi_diff,
                hydrophobicity_cross_corr, charge_cross_corr,
                binary_hydro_cross_corr, pc_composite_score, detection_source
               FROM similarities
               WHERE (query_id = ? OR target_id = ?)
                 AND evalue <= ?
               ORDER BY evalue ASC""",
            (antigen_id, antigen_id, antigen_id, max_evalue),
        ).fetchall()

        # Deduplicate by partner_id, keeping best (lowest) E-value
        seen = {}
        for r in rows:
            pid = r["partner_id"]
            if pid not in seen or r["evalue"] < seen[pid]["evalue"]:
                seen[pid] = dict(r)

        results = sorted(seen.values(), key=lambda x: x["evalue"])
        return results[:top_n]

    def similarity_count(self) -> int:
        return self.conn.execute("SELECT COUNT(*) FROM similarities").fetchone()[0]

    def get_completed_queries(self) -> set[str]:
        """Return set of query IDs that have already been processed (for resume support)."""
        rows = self.conn.execute(
            "SELECT DISTINCT query_id FROM similarities"
        ).fetchall()
        return {r[0] for r in rows}

    # ------------------------------------------------------------------ #
    #  Metadata                                                            #
    # ------------------------------------------------------------------ #

    def set_meta(self, key: str, value: str):
        self.conn.execute(
            "INSERT OR REPLACE INTO run_metadata (key, value) VALUES (?, ?)", (key, value)
        )
        self.conn.commit()

    def get_meta(self, key: str) -> str | None:
        row = self.conn.execute(
            "SELECT value FROM run_metadata WHERE key = ?", (key,)
        ).fetchone()
        return row[0] if row else None

    # ------------------------------------------------------------------ #
    #  Export                                                              #
    # ------------------------------------------------------------------ #

    def export_csv(self, output_path: str | Path, max_evalue: float = 10.0):
        """Export all similarity results to CSV.

        Includes pairs that passed either the SW E-value threshold or the
        contiguous match threshold (or both).
        """
        import csv
        rows = self.conn.execute(
            """SELECT s.*, a1.sequence as query_original, a2.sequence as target_original,
                      a1.stripped as query_stripped, a2.stripped as target_stripped
               FROM similarities s
               JOIN antigens a1 ON s.query_id  = a1.id
               JOIN antigens a2 ON s.target_id = a2.id
               WHERE s.evalue <= ?
               ORDER BY s.evalue ASC""",
            (max_evalue,),
        ).fetchall()

        with open(output_path, "w", newline="") as f:
            writer = csv.writer(f)
            writer.writerow([
                "query_id", "target_id", "evalue", "bit_score", "raw_score",
                "aligned_region_len", "query_length", "target_length",
                "consec_match_len", "consec_query_start", "consec_query_end",
                "consec_target_start", "consec_target_end",
                # Physicochemical Track 2 columns (NULL if not yet scored)
                "hydrophobicity_corr", "charge_corr", "pi_diff",
                "hydrophobicity_cross_corr", "charge_cross_corr",
                "binary_hydro_cross_corr", "pc_composite_score",
                "detection_source",
                "query_stripped", "target_stripped",
            ])
            for r in rows:
                writer.writerow([
                    r["query_id"], r["target_id"], r["evalue"],
                    r["bit_score"], r["raw_score"], r["aligned_region_len"],
                    r["query_length"], r["target_length"],
                    r["consec_match_len"], r["consec_query_start"],
                    r["consec_query_end"], r["consec_target_start"],
                    r["consec_target_end"],
                    r["hydrophobicity_corr"], r["charge_corr"], r["pi_diff"],
                    r["hydrophobicity_cross_corr"], r["charge_cross_corr"],
                    r["binary_hydro_cross_corr"], r["pc_composite_score"],
                    r["detection_source"],
                    r["query_stripped"], r["target_stripped"],
                ])

        return len(rows)

    def close(self):
        self.conn.close()

    def __enter__(self):
        return self

    def __exit__(self, *args):
        self.close()
