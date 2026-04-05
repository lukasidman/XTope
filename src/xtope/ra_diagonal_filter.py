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

Short-region rescue (BLOSUM62 scoring):
  Pairs that fail the chain threshold but have at least one seed hit on a diagonal
  get a second chance. The original (full-alphabet) residues along the best diagonal
  region are scored with BLOSUM62. If the mean per-position score is high enough
  (default >= 3.0), the pair is force-included in the SW batch. This catches short
  (5-8 aa) regions of near-identity that produce only 1 reduced-alphabet k-mer hit
  — too short for a chain, but biologically meaningful.

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

from .sw_fallback import BLOSUM62 as _BLOSUM62_DICT


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

# ---------------------------------------------------------------------------
# BLOSUM62 rescue for short high-similarity regions
# ---------------------------------------------------------------------------
# When a pair fails the chain threshold (chain < min_chain) but has at least
# one seed hit, we score the original residues along the best diagonal using
# BLOSUM62. If the mean per-position score is high enough, the pair is rescued.
#
# Threshold rationale:
#   BLOSUM62 identity scores range from 4 (A, V) to 11 (W), with most at 4-6.
#   A mean of 3.0 per position requires near-identity across the region:
#     - 5 identical residues averaging 4.5 → mean 4.5 (passes easily)
#     - 4 identical + 1 conservative sub (score ~1) → mean ~3.4 (passes)
#     - 3 identical + 2 mismatches (score ~-2) → mean ~1.4 (fails)
#   This is deliberately strict — we only rescue pairs where the shared
#   region is genuinely high-identity, not merely conserved in properties
#   (the reduced alphabet already handles that).
DEFAULT_RESCUE_THRESHOLD = 3.0  # min mean BLOSUM62 score per position
DEFAULT_RESCUE_MIN_REGION = 5   # min region length (aa) to attempt rescue


def _score_diagonal_region(
    query_seq: str,
    target_seq: str,
    query_start: int,
    target_start: int,
    length: int,
) -> float:
    """Score a diagonal region between two sequences using BLOSUM62.

    Extracts `length` residues from each sequence starting at the given
    positions and computes the mean BLOSUM62 score per position.

    Args:
        query_seq: Full query amino acid sequence.
        target_seq: Full target amino acid sequence.
        query_start: Start position in query.
        target_start: Start position in target.
        length: Number of positions to score.

    Returns:
        Mean BLOSUM62 score per position, or -999.0 if the region
        extends beyond either sequence.
    """
    if (query_start < 0 or target_start < 0
            or query_start + length > len(query_seq)
            or target_start + length > len(target_seq)):
        return -999.0

    total = 0.0
    for offset in range(length):
        q_aa = query_seq[query_start + offset].upper()
        t_aa = target_seq[target_start + offset].upper()
        total += _BLOSUM62_DICT.get((q_aa, t_aa), -4)
    return total / length


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

    def __init__(
        self,
        k: int = DEFAULT_K,
        min_chain: int = DEFAULT_MIN_CHAIN,
        rescue_threshold: float = DEFAULT_RESCUE_THRESHOLD,
        rescue_min_region: int = DEFAULT_RESCUE_MIN_REGION,
    ):
        self.k = k
        self.min_chain = min_chain
        self.rescue_threshold = rescue_threshold
        self.rescue_min_region = rescue_min_region
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
    ) -> dict[str, tuple[int, int, int]]:
        """Find the longest diagonal chain per candidate.

        A diagonal chain is a run of consecutive (query_pos, target_pos) pairs
        where both positions increment by 1 each step — meaning the reduced
        k-mers are adjacent and on the same alignment diagonal.

        Args:
            seed_pairs: target_id → list of (query_pos, target_pos) seed hits.

        Returns:
            target_id → (chain_length, best_diagonal, chain_start_qpos).
            best_diagonal is query_pos - target_pos for the diagonal with the
            longest chain. chain_start_qpos is the query position where the
            best chain begins.
        """
        result: dict[str, tuple[int, int, int]] = {}
        for target_id, pairs in seed_pairs.items():
            # Group by diagonal (query_pos - target_pos)
            diagonals: dict[int, list[int]] = defaultdict(list)
            for qp, tp in pairs:
                diag = qp - tp
                diagonals[diag].append(qp)

            max_chain = 0
            best_diag = 0
            best_start_qpos = 0
            for diag, positions in diagonals.items():
                positions.sort()
                chain = 1
                chain_start = positions[0]
                run_best = 1
                run_best_start = positions[0]
                for i in range(1, len(positions)):
                    if positions[i] == positions[i - 1] + 1:
                        chain += 1
                        if chain > run_best:
                            run_best = chain
                            run_best_start = chain_start
                    else:
                        chain = 1
                        chain_start = positions[i]
                if run_best > max_chain:
                    max_chain = run_best
                    best_diag = diag
                    best_start_qpos = run_best_start

            result[target_id] = (max_chain, best_diag, best_start_qpos)
        return result

    def query(
        self,
        sequence: str,
        exclude_id: str | None = None,
    ) -> set[str]:
        """Find IDs of sequences that share a local region with the query.

        Includes both normal chain-threshold passes and BLOSUM62-rescued pairs.
        Rescued pairs are those with a short chain (below min_chain) but whose
        original residues along the best diagonal score highly in BLOSUM62.

        Args:
            sequence: Query amino acid sequence (tag-stripped).
            exclude_id: ID to skip (typically the query itself).

        Returns:
            Set of antigen IDs passing the filter (chain or rescue).
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
        chain_info = self._find_diagonal_chains(seed_pairs)
        passing: set[str] = set()

        for tid, (chain_len, best_diag, chain_start_qpos) in chain_info.items():
            if chain_len >= self.min_chain:
                passing.add(tid)
            else:
                # Rescue attempt: score original residues along best diagonal
                # The seed region covers k positions starting at each seed's
                # query position. For a chain of length c, the region is
                # chain_start_qpos .. chain_start_qpos + c + k - 1
                region_len = chain_len + self.k - 1
                if region_len < self.rescue_min_region:
                    continue

                target_start = chain_start_qpos - best_diag
                mean_score = _score_diagonal_region(
                    sequence,
                    self.sequences[tid],
                    chain_start_qpos,
                    target_start,
                    region_len,
                )
                if mean_score >= self.rescue_threshold:
                    passing.add(tid)

        return passing

    def query_with_details(
        self,
        sequence: str,
        exclude_id: str | None = None,
    ) -> list[dict]:
        """Like query() but returns detailed info including rescue scoring.

        Args:
            sequence: Query amino acid sequence (tag-stripped).
            exclude_id: ID to skip.

        Returns:
            List of dicts with keys:
                target_id: str
                chain_length: int
                diagonal: int
                chain_start_qpos: int
                region_len: int
                blosum62_mean: float  (mean per-position BLOSUM62 score, or None)
                passed_chain: bool    (met chain threshold)
                rescued: bool         (failed chain but passed BLOSUM62 rescue)
            Sorted by chain_length descending, then blosum62_mean descending.
            Includes ALL targets with at least 1 seed hit.
        """
        query_kmer_pos = reduced_kmers_with_pos(sequence, self.k)
        if not query_kmer_pos:
            return []

        seed_pairs: dict[str, list[tuple[int, int]]] = defaultdict(list)
        for kmer, qpos in query_kmer_pos:
            for target_id, tpos in self.index.get(kmer, []):
                if target_id != exclude_id:
                    seed_pairs[target_id].append((qpos, tpos))

        chain_info = self._find_diagonal_chains(seed_pairs)
        results = []
        for target_id, (chain_len, best_diag, chain_start_qpos) in chain_info.items():
            region_len = chain_len + self.k - 1
            passed_chain = chain_len >= self.min_chain

            # Score original residues regardless of chain pass/fail
            blosum62_mean = None
            rescued = False
            if region_len >= self.rescue_min_region:
                target_start = chain_start_qpos - best_diag
                blosum62_mean = _score_diagonal_region(
                    sequence,
                    self.sequences[target_id],
                    chain_start_qpos,
                    target_start,
                    region_len,
                )
                if not passed_chain and blosum62_mean >= self.rescue_threshold:
                    rescued = True

            results.append({
                "target_id": target_id,
                "chain_length": chain_len,
                "diagonal": best_diag,
                "chain_start_qpos": chain_start_qpos,
                "region_len": region_len,
                "blosum62_mean": blosum62_mean,
                "passed_chain": passed_chain,
                "rescued": rescued,
            })

        results.sort(key=lambda x: (x["chain_length"], x["blosum62_mean"] or -999),
                      reverse=True)
        return results

    def __len__(self) -> int:
        return len(self.sequences)
