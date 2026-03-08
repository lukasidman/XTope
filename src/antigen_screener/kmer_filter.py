"""
K-mer pre-filter module.
Rapidly reduces the comparison space using Jaccard similarity on amino acid k-mers.
This is a speed filter only — threshold loosely to avoid missing true positives.
"""

from collections import defaultdict
from typing import Iterator
import numpy as np


DEFAULT_K = 6          # k-mer length (6-8 recommended for short peptides)
DEFAULT_THRESHOLD = 0.04  # Minimum Jaccard similarity to pass filter


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
