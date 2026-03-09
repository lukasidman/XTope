"""
K-mer pre-filter module.
Rapidly reduces the comparison space using Jaccard similarity on amino acid k-mers,
with an optional second pass that detects short consecutive shared epitopes that
Jaccard alone would miss.

Two-pass strategy
-----------------
Pass 1 — KmerIndex (k=6, Jaccard ≥ 0.04)
    Fast.  Catches shared regions of ~15 aa or longer.

Pass 2 — ChainKmerIndex (k=4, consecutive chain ≥ 3 k-mers)
    Slower but still O(shared_4mers).  Catches shared regions of 6 aa or longer
    (chain of 3 consecutive 4-mers = 3 + 4 - 1 = 6 contiguous amino acids).
    Only considers pairs not already found by Pass 1.

Use TwoPassFilter as a drop-in replacement for KmerIndex to enable both passes.
"""

from collections import defaultdict
from typing import Iterator
import numpy as np


DEFAULT_K = 6          # k-mer length for pass 1
DEFAULT_THRESHOLD = 0.04  # Minimum Jaccard similarity to pass pass-1 filter
DEFAULT_K2 = 4         # k-mer length for pass 2 (chain detection)
DEFAULT_MIN_CHAIN = 3  # Minimum consecutive k-mers to pass pass-2 filter
                       # (3 × k=4 → 6 aa shared region)


def kmers(sequence: str, k: int = DEFAULT_K) -> set[str]:
    """Return the set of all k-mers in a sequence."""
    return {sequence[i:i+k] for i in range(len(sequence) - k + 1)}


def jaccard(set_a: set, set_b: set) -> float:
    """Jaccard similarity between two k-mer sets."""
    if not set_a or not set_b:
        return 0.0
    intersection = len(set_a & set_b)
    union = len(set_a | set_b)
    return intersection / union if union > 0 else 0.0


class KmerIndex:
    """
    Inverted index for fast k-mer based pre-filtering.
    Maps each k-mer → list of antigen IDs containing it.
    Allows O(|query_kmers| * avg_postings) lookup instead of O(n) full scan.
    """

    def __init__(self, k: int = DEFAULT_K):
        self.k = k
        self.index: dict[str, list[str]] = defaultdict(list)
        self.kmer_sets: dict[str, set[str]] = {}   # id → kmer set
        self.sequences: dict[str, str] = {}         # id → stripped sequence

    def add(self, antigen_id: str, sequence: str):
        """Index a single sequence."""
        ks = kmers(sequence, self.k)
        self.kmer_sets[antigen_id] = ks
        self.sequences[antigen_id] = sequence
        for kmer in ks:
            self.index[kmer].append(antigen_id)

    def add_batch(self, records: list[tuple[str, str]]):
        """Bulk index a list of (id, sequence) tuples."""
        for antigen_id, sequence in records:
            self.add(antigen_id, sequence)

    def query(
        self,
        sequence: str,
        threshold: float = DEFAULT_THRESHOLD,
        exclude_id: str | None = None,
    ) -> list[tuple[str, float]]:
        """
        Find all indexed sequences with Jaccard similarity >= threshold.

        Returns list of (antigen_id, jaccard_score) sorted descending.
        """
        query_kmers = kmers(sequence, self.k)
        if not query_kmers:
            return []

        # Count shared k-mers via inverted index (much faster than full scan)
        candidate_counts: dict[str, int] = defaultdict(int)
        for kmer in query_kmers:
            for aid in self.index.get(kmer, []):
                if aid != exclude_id:
                    candidate_counts[aid] += 1

        # Compute Jaccard only for candidates with at least 1 shared k-mer
        results = []
        for aid, shared_count in candidate_counts.items():
            target_kmers = self.kmer_sets[aid]
            union = len(query_kmers) + len(target_kmers) - shared_count
            score = shared_count / union if union > 0 else 0.0
            if score >= threshold:
                results.append((aid, score))

        results.sort(key=lambda x: x[1], reverse=True)
        return results

    def __len__(self):
        return len(self.sequences)


# ---------------------------------------------------------------------------
# Pass 2 — consecutive chain detection
# ---------------------------------------------------------------------------

class ChainKmerIndex:
    """
    k=4 inverted index that detects short consecutive shared epitopes.

    Unlike Jaccard (which ignores position), this index checks whether shared
    k-mers form a consecutive chain — i.e., k-mer at position i in the query
    is matched by a k-mer at position j in the target, AND the k-mer at i+1
    is matched at j+1, and so on.  A chain of length N covers N + k - 1
    contiguous amino acids, so with k=4 and min_chain=3 the shortest passing
    shared region is 6 aa.

    Used as Pass 2 in TwoPassFilter to catch short epitopes that Jaccard misses.
    """

    def __init__(self, k: int = DEFAULT_K2, min_chain: int = DEFAULT_MIN_CHAIN):
        self.k = k
        self.min_chain = min_chain
        # kmer → list of (antigen_id, position_in_sequence)
        self.index: dict[str, list[tuple[str, int]]] = defaultdict(list)
        self.sequences: dict[str, str] = {}

    def add(self, antigen_id: str, sequence: str) -> None:
        """Index a single sequence, storing per-position k-mer entries."""
        self.sequences[antigen_id] = sequence
        for i in range(len(sequence) - self.k + 1):
            self.index[sequence[i:i + self.k]].append((antigen_id, i))

    def add_batch(self, records: list[tuple[str, str]]) -> None:
        """Bulk index a list of (id, sequence) tuples."""
        for antigen_id, sequence in records:
            self.add(antigen_id, sequence)

    def _max_chain(self, seed_pairs: set[tuple[int, int]]) -> int:
        """
        Given a set of (query_pos, target_pos) seed pairs, return the length
        of the longest consecutive diagonal run.

        Two pairs (i, j) and (i+1, j+1) are consecutive — they correspond to
        adjacent overlapping k-mers on the same diagonal, meaning the sequences
        share a contiguous region of length chain + k - 1.
        """
        max_len = 0
        for qp, tp in seed_pairs:
            # Only start a chain from the beginning of a run
            if (qp - 1, tp - 1) in seed_pairs:
                continue
            length = 1
            while (qp + length, tp + length) in seed_pairs:
                length += 1
            if length > max_len:
                max_len = length
        return max_len

    def query(
        self,
        sequence: str,
        exclude_id: str | None = None,
    ) -> list[str]:
        """
        Return IDs of indexed sequences that share a consecutive k-mer chain
        of length >= min_chain with the query sequence.

        Args:
            sequence:   Query sequence (tag-stripped).
            exclude_id: Antigen ID to skip (typically the query itself).

        Returns:
            List of passing antigen IDs (unordered).
        """
        k = self.k

        # Build query position map: kmer → positions in query sequence
        query_pos: dict[str, list[int]] = defaultdict(list)
        for i in range(len(sequence) - k + 1):
            query_pos[sequence[i:i + k]].append(i)

        # Accumulate (query_pos, target_pos) seed pairs per candidate
        candidate_seeds: dict[str, set[tuple[int, int]]] = defaultdict(set)
        for kmer, qpos_list in query_pos.items():
            for (target_id, tpos) in self.index.get(kmer, []):
                if target_id == exclude_id:
                    continue
                for qpos in qpos_list:
                    candidate_seeds[target_id].add((qpos, tpos))

        # Filter candidates whose seed pairs form a long enough chain
        passing = []
        for target_id, seeds in candidate_seeds.items():
            if self._max_chain(seeds) >= self.min_chain:
                passing.append(target_id)

        return passing

    def __len__(self) -> int:
        return len(self.sequences)


# ---------------------------------------------------------------------------
# Two-pass filter — drop-in replacement for KmerIndex
# ---------------------------------------------------------------------------

class TwoPassFilter:
    """
    Two-pass k-mer pre-filter combining Jaccard (pass 1) and consecutive chain
    detection (pass 2).  Drop-in replacement for KmerIndex.

    Pass 1 (KmerIndex, k=6, Jaccard ≥ threshold)
        Catches shared regions of ~15 aa or longer quickly.

    Pass 2 (ChainKmerIndex, k=4, chain ≥ min_chain)
        Catches short epitopes (≥ 6 aa) missed by Jaccard alone.
        Only runs for candidates not already found by Pass 1.

    The return type of query() is identical to KmerIndex.query():
        list[tuple[str, float]]
    Pass-2 hits are returned with jaccard_score = 0.0 (they have no Jaccard
    score by definition; SW alignment provides the authoritative score).
    """

    def __init__(
        self,
        k1: int = DEFAULT_K,
        threshold: float = DEFAULT_THRESHOLD,
        k2: int = DEFAULT_K2,
        min_chain: int = DEFAULT_MIN_CHAIN,
    ):
        self.threshold = threshold
        self.pass1 = KmerIndex(k=k1)
        self.pass2 = ChainKmerIndex(k=k2, min_chain=min_chain)
        # Mirror of KmerIndex.sequences for pipeline compatibility
        self.sequences: dict[str, str] = self.pass1.sequences

    def add(self, antigen_id: str, sequence: str) -> None:
        """Index a sequence in both passes."""
        self.pass1.add(antigen_id, sequence)
        self.pass2.add(antigen_id, sequence)

    def add_batch(self, records: list[tuple[str, str]]) -> None:
        """Bulk index a list of (id, sequence) tuples."""
        for antigen_id, sequence in records:
            self.add(antigen_id, sequence)

    def query(
        self,
        sequence: str,
        threshold: float | None = None,
        exclude_id: str | None = None,
    ) -> list[tuple[str, float]]:
        """
        Find candidates passing either the Jaccard filter (pass 1) or the
        consecutive chain filter (pass 2).

        Args:
            sequence:   Query sequence (tag-stripped).
            threshold:  Jaccard threshold override for pass 1.
            exclude_id: Antigen ID to exclude (typically the query itself).

        Returns:
            list of (antigen_id, jaccard_score) sorted by score descending.
            Pass-2 hits that did not reach pass 1 carry jaccard_score = 0.0.
        """
        t = threshold if threshold is not None else self.threshold

        # Pass 1 — Jaccard
        pass1_results = self.pass1.query(sequence, threshold=t, exclude_id=exclude_id)
        pass1_ids = {aid for aid, _ in pass1_results}

        # Pass 2 — chain detection (only for IDs not already in pass 1)
        pass2_ids = self.pass2.query(sequence, exclude_id=exclude_id)
        pass2_new = [aid for aid in pass2_ids if aid not in pass1_ids]

        # Merge: pass 1 results keep their Jaccard score; pass 2 additions get 0.0
        results: list[tuple[str, float]] = list(pass1_results)
        results.extend((aid, 0.0) for aid in pass2_new)

        return results

    def __len__(self) -> int:
        return len(self.pass1)
