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
    normalized_score    REAL NOT NULL,
    query_length        INTEGER NOT NULL,
    target_length       INTEGER NOT NULL,
    aligned_region_len  INTEGER NOT NULL,
    PRIMARY KEY (query_id, target_id)
);

CREATE INDEX IF NOT EXISTS idx_sim_query  ON similarities(query_id);
CREATE INDEX IF NOT EXISTS idx_sim_target ON similarities(target_id);
CREATE INDEX IF NOT EXISTS idx_sim_score  ON similarities(normalized_score DESC);

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
        self.conn.commit()

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

    # ------------------------------------------------------------------ #
    #  Similarity results                                                  #
    # ------------------------------------------------------------------ #

    def insert_similarities_batch(self, results: list[dict]):
        """Bulk insert similarity results. Ignores duplicates."""
        self.conn.executemany(
            """INSERT OR IGNORE INTO similarities
               (query_id, target_id, raw_score, normalized_score,
                query_length, target_length, aligned_region_len)
               VALUES (:query_id, :target_id, :raw_score, :normalized_score,
                       :query_length, :target_length, :aligned_region_len)""",
            results,
        )
        self.conn.commit()

    def query_similar(
        self,
        antigen_id: str,
        min_score: float = 1.0,
        top_n: int = 50,
    ) -> list[dict]:
        """
        Retrieve precomputed similar antigens for a given ID.
        Searches both directions (query→target and target→query), deduplicates.
        """
        rows = self.conn.execute(
            """SELECT
                CASE WHEN query_id = ? THEN target_id ELSE query_id END AS partner_id,
                raw_score, normalized_score, query_length, target_length, aligned_region_len
               FROM similarities
               WHERE (query_id = ? OR target_id = ?)
                 AND normalized_score >= ?
               ORDER BY normalized_score DESC""",
            (antigen_id, antigen_id, antigen_id, min_score),
        ).fetchall()

        # Deduplicate by partner_id, keeping best score
        seen = {}
        for r in rows:
            pid = r["partner_id"]
            if pid not in seen or r["normalized_score"] > seen[pid]["normalized_score"]:
                seen[pid] = dict(r)

        results = sorted(seen.values(), key=lambda x: x["normalized_score"], reverse=True)
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

    def export_csv(self, output_path: str | Path, min_score: float = 0.0):
        """Export all similarity results to CSV."""
        import csv
        rows = self.conn.execute(
            """SELECT s.*, a1.sequence as query_original, a2.sequence as target_original,
                      a1.stripped as query_stripped, a2.stripped as target_stripped
               FROM similarities s
               JOIN antigens a1 ON s.query_id  = a1.id
               JOIN antigens a2 ON s.target_id = a2.id
               WHERE s.normalized_score >= ?
               ORDER BY s.normalized_score DESC""",
            (min_score,),
        ).fetchall()

        with open(output_path, "w", newline="") as f:
            writer = csv.writer(f)
            writer.writerow([
                "query_id", "target_id", "normalized_score", "raw_score",
                "aligned_region_len", "query_length", "target_length",
                "query_stripped", "target_stripped",
            ])
            for r in rows:
                writer.writerow([
                    r["query_id"], r["target_id"], r["normalized_score"],
                    r["raw_score"], r["aligned_region_len"],
                    r["query_length"], r["target_length"],
                    r["query_stripped"], r["target_stripped"],
                ])

        return len(rows)

    def close(self):
        self.conn.close()

    def __enter__(self):
        return self

    def __exit__(self, *args):
        self.close()
