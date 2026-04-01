"""
Reduced-alphabet diagonal chain filter for pre-screening before batched NumPy SW.

Goal: cheaply detect pairs that share a ~10aa local region (even with conservative
substitutions) so we can skip obviously unrelated pairs in the expensive SW stage.
This is NOT trying to be a sensitive final filter — it's a generous screen to
eliminate ~60% of background pairs that share no local similarity whatsoever.

Strategy:
  1. Convert sequences to a reduced alphabet (~11 BLOSUM-derived groups)
  2. Build an inverted index of (reduced_kmer, position) tuples
  3. For each query, find target sequences sharing k-mers on the same diagonal
     (diagonal = query_pos - target_pos), indicating a contiguous local match
  4. Pairs with a diagonal chain of >= min_chain reduced k-mers pass the filter

Why reduced alphabet works here:
  A 10aa conserved epitope with 2-3 conservative substitutions (I→L, D→E, K→R)
  still produces ~5 consecutive reduced-alphabet 5-mers on the same diagonal.
  With 11 groups and k=5, random match probability per position is (1/11)^5 ≈ 0.0006%,
  so random diagonal chains of length >=2 are rare — good selectivity.

Why this is different from the abandoned ReducedAlphabetChainIndex:
  That attempt used 6 groups (too coarse → flooding) and was benchmarked as a
  sensitivity-optimizing filter (trying to catch everything Jaccard missed). This
  uses 11 groups (more selective) and optimizes for *cheap rejection* of obvious
  non-matches before batched NumPy SW, where a false positive just means a sequence
  stays in the batch (low marginal cost) rather than triggering an expensive per-pair
  SW call (high marginal cost).
"""

from collections import defaultdict
from typing import Iterator


# ---------------------------------------------------------------------------
# Reduced alphabet: 11 BLOSUM62-derived similarity groups
# ---------------------------------------------------------------------------
# Groups are based on BLOSUM62 positive-score clusters. Amino acids within a
# group have positive or zero BLOSUM62 scores for most pairwise comparisons,
# meaning conservative substitutions within a group are common in evolution.

REDUCED_ALPHABET = {
    # Aliphatic hydrophobic (BLOSUM62: I↔L=2, I↔V=3, L↔V=1, M↔I=1, M↔L=2)
    'I': 'h', 'L': 'h', 'V': 'h', 'M': 'h',
    # Small hydrophobic + alanine (A↔G=0 borderline, but A is more hydrophobic)
    'A': 'a',
    # Aromatic (F↔Y=3, F↔W=1, W↔Y=2)
    'F': 'r', 'W': 'r', 'Y': 'r',
    # Small polar (S↔T=1, S↔N=1)
    'S': 's', 'T': 's',
    # Acidic (D↔E=2)
    'D': 'd', 'E': 'd',
    # Amide (N↔Q=0 borderline, but chemically similar)
    'N': 'n', 'Q': 'n',
    # Basic (K↔R=2)
    'K': 'b', 'R': 'b',
    # Histidine (unique: aromatic + basic, pH-dependent ionisation)
    'H': 'i',
    # Cysteine (unique: disulfide bonding)
    'C': 'c',
    # Glycine (unique: backbone flexibility, smallest side chain)
    'G': 'g',
    # Proline (unique: helix breaker, ring structure)
    'P': 'p',
}

# Number of distinct groups (for selectivity calculations)
N_GROUPS = len(set(REDUCED_ALPHABET.values()))  # should be 11

DEFAULT_K = 5          # k-mer length in reduced alphabet space
DEFAULT_MIN_CHAIN = 2  # minimum consecutive k-mers on same diagonal
                       # chain=2 with k=5 → shared region of 6aa in reduced space
                       # generous enough to catch 8-10aa epitopes with 2-3 substitutions


def reduce_sequence(sequence: str) -> str:
    """Convert an amino acid sequence to reduced alphabet.

    Unknown characters (X, B, Z, etc.) are mapped to '?' and will not match
    any indexed k-mer, which is the desired behaviour (unknown = no evidence).

    Args:
        sequence: Amino acid sequence (uppercase).

    Returns:
        Reduced-alphabet string of same length.
    """
    return ''.join(REDUCED_ALPHABET.get(aa, '?') for aa in sequence.upper())


def reduced_kmers_with_pos(sequence: str, k: int = DEFAULT_K) -> list[tuple[str, int]]:
    """Extract (reduced_kmer, position) tuples from a sequence.

    Args:
        sequence: Original amino acid sequence.
        k: k-mer length in reduced alphabet space.

    Returns:
        List of (reduced_kmer, position) tuples.
    """
    reduced = reduce_sequence(sequence)
    result = []
    for i in range(len(reduced) - k + 1):
        kmer = reduced[i:i + k]
        if '?' not in kmer:  # skip k-mers containing unknown residues
            result.append((kmer, i))
    return result


class RADiagonalFilter:
    """
    Reduced-Alphabet Diagonal chain filter.

    Indexes sequences in reduced-alphabet space and detects pairs sharing
    consecutive k-mers on the same diagonal (= contiguous local match in
    reduced space). Designed as a generous pre-filter for batched NumPy SW.

    Usage:
        filt = RADiagonalFilter(k=5, min_chain=2)
        filt.add_batch([(id1, seq1), (id2, seq2), ...])

        # For each query, get set of IDs that pass the filter
        passing_ids = filt.query(query_seq, exclude_id=query_id)
    """

    def __init__(self, k: int = DEFAULT_K, min_chain: int = DEFAULT_MIN_CHAIN):
        self.k = k
        self.min_chain = min_chain
        # Inverted index: reduced_kmer → list of (antigen_id, position)
        self.index: dict[str, list[tuple[str, int]]] = defaultdict(list)
        self.sequences: dict[str, str] = {}  # id → original sequence
        self._reduced: dict[str, str] = {}   # id → reduced sequence (for debugging)

    def add(self, antigen_id: str, sequence: str) -> None:
        """Index a single sequence."""
        self.sequences[antigen_id] = sequence
        reduced = reduce_sequence(sequence)
        self._reduced[antigen_id] = reduced
        for kmer, pos in reduced_kmers_with_pos(sequence, self.k):
            self.index[kmer].append((antigen_id, pos))

    def add_batch(self, records: list[tuple[str, str]]) -> None:
        """Bulk index a list of (id, sequence) tuples."""
        for antigen_id, sequence in records:
            self.add(antigen_id, sequence)

    def _find_diagonal_chains(
        self,
        seed_pairs: dict[str, list[tuple[int, int]]],
    ) -> dict[str, int]:
        """Find the longest diagonal chain per candidate.

        A diagonal chain is a run of consecutive (query_pos, target_pos) pairs
        where both positions increment by 1 each step — meaning the reduced
        k-mers are adjacent and on the same alignment diagonal.

        Args:
            seed_pairs: target_id → list of (query_pos, target_pos) seed hits.

        Returns:
            target_id → longest chain length.
        """
        result = {}
        for target_id, pairs in seed_pairs.items():
            if len(pairs) < self.min_chain:
                # Can't possibly form a long enough chain
                result[target_id] = len(pairs) if pairs else 0
                continue

            # Group by diagonal (query_pos - target_pos)
            diagonals: dict[int, list[int]] = defaultdict(list)
            for qp, tp in pairs:
                diag = qp - tp
                diagonals[diag].append(qp)

            max_chain = 0
            for diag, positions in diagonals.items():
                if len(positions) < self.min_chain:
                    if len(positions) > max_chain:
                        max_chain = len(positions)
                    continue

                # Sort positions and find longest consecutive run
                positions.sort()
                chain = 1
                best = 1
                for i in range(1, len(positions)):
                    if positions[i] == positions[i - 1] + 1:
                        chain += 1
                        if chain > best:
                            best = chain
                    else:
                        chain = 1
                if best > max_chain:
                    max_chain = best

            result[target_id] = max_chain
        return result

    def query(
        self,
        sequence: str,
        exclude_id: str | None = None,
    ) -> set[str]:
        """Find IDs of sequences that share a local region with the query.

        Args:
            sequence: Query amino acid sequence (tag-stripped).
            exclude_id: ID to skip (typically the query itself).

        Returns:
            Set of antigen IDs passing the diagonal chain filter.
        """
        query_kmer_pos = reduced_kmers_with_pos(sequence, self.k)
        if not query_kmer_pos:
            return set()

        # Collect seed hits: target_id → [(query_pos, target_pos), ...]
        seed_pairs: dict[str, list[tuple[int, int]]] = defaultdict(list)
        for kmer, qpos in query_kmer_pos:
            for target_id, tpos in self.index.get(kmer, []):
                if target_id != exclude_id:
                    seed_pairs[target_id].append((qpos, tpos))

        # Find diagonal chains and filter
        chain_lengths = self._find_diagonal_chains(seed_pairs)
        return {
            tid for tid, chain_len in chain_lengths.items()
            if chain_len >= self.min_chain
        }

    def query_with_details(
        self,
        sequence: str,
        exclude_id: str | None = None,
    ) -> list[tuple[str, int, int]]:
        """Like query() but returns (target_id, chain_length, diagonal) for analysis.

        Args:
            sequence: Query amino acid sequence (tag-stripped).
            exclude_id: ID to skip.

        Returns:
            List of (target_id, max_chain_length, best_diagonal) tuples,
            sorted by chain length descending. Includes ALL targets with
            at least 1 seed hit (not just those passing min_chain).
        """
        query_kmer_pos = reduced_kmers_with_pos(sequence, self.k)
        if not query_kmer_pos:
            return []

        seed_pairs: dict[str, list[tuple[int, int]]] = defaultdict(list)
        for kmer, qpos in query_kmer_pos:
            for target_id, tpos in self.index.get(kmer, []):
                if target_id != exclude_id:
                    seed_pairs[target_id].append((qpos, tpos))

        results = []
        for target_id, pairs in seed_pairs.items():
            # Find best diagonal
            diagonals: dict[int, list[int]] = defaultdict(list)
            for qp, tp in pairs:
                diag = qp - tp
                diagonals[diag].append(qp)

            best_chain = 0
            best_diag = 0
            for diag, positions in diagonals.items():
                positions.sort()
                chain = 1
                run_best = 1
                for i in range(1, len(positions)):
                    if positions[i] == positions[i - 1] + 1:
                        chain += 1
                        if chain > run_best:
                            run_best = chain
                    else:
                        chain = 1
                if run_best > best_chain:
                    best_chain = run_best
                    best_diag = diag

            results.append((target_id, best_chain, best_diag))

        results.sort(key=lambda x: x[1], reverse=True)
        return results

    def __len__(self) -> int:
        return len(self.sequences)
