"""
Pure Python Smith-Waterman local alignment.
Used as a fallback when parasail is not installed.
Slower than parasail (~50x) but functionally identical.
For large databases (100,000+ antigens), strongly recommend installing parasail.
"""

import numpy as np
from typing import NamedTuple

# BLOSUM62 as a flat dict (subset of commonly needed pairs)
# Full matrix loaded from data below
_BLOSUM62_RAW = """
   A  R  N  D  C  Q  E  G  H  I  L  K  M  F  P  S  T  W  Y  V
A  4 -1 -2 -2  0 -1 -1  0 -2 -1 -1 -1 -1 -2 -1  1  0 -3 -2  0
R -1  5  0 -2 -3  1  0 -2  0 -3 -2  2 -1 -3 -2 -1 -1 -3 -2 -3
N -2  0  6  1 -3  0  0  0  1 -3 -3  0 -2 -3 -2  1  0 -4 -2 -3
D -2 -2  1  6 -3  0  2 -1 -1 -3 -4 -1 -3 -3 -1  0 -1 -4 -3 -3
C  0 -3 -3 -3  9 -3 -4 -3 -3 -1 -1 -3 -1 -2 -3 -1 -1 -2 -2 -1
Q -1  1  0  0 -3  5  2 -2  0 -3 -2  1  0 -3 -1  0 -1 -2 -1 -2
E -1  0  0  2 -4  2  5 -2  0 -3 -3  1 -2 -3 -1  0 -1 -3 -2 -2
G  0 -2  0 -1 -3 -2 -2  6 -2 -4 -4 -2 -3 -3 -2  0 -2 -2 -3 -3
H -2  0  1 -1 -3  0  0 -2  8 -3 -3 -1 -2 -1 -2 -1 -2 -2  2 -3
I -1 -3 -3 -3 -1 -3 -3 -4 -3  4  2 -3  1  0 -3 -2 -1 -3 -1  3
L -1 -2 -3 -4 -1 -2 -3 -4 -3  2  4 -2  2  0 -3 -2 -1 -2 -1  1
K -1  2  0 -1 -3  1  1 -2 -1 -3 -2  5 -1 -3 -1  0 -1 -3 -2 -2
M -1 -1 -2 -3 -1  0 -2 -3 -2  1  2 -1  5  0 -2 -1 -1 -1 -1  1
F -2 -3 -3 -3 -2 -3 -3 -3 -1  0  0 -3  0  6 -4 -2 -2  1  3 -1
P -1 -2 -2 -1 -3 -1 -1 -2 -2 -3 -3 -1 -2 -4  7 -1 -1 -4 -3 -2
S  1 -1  1  0 -1  0  0  0 -1 -2 -2  0 -1 -2 -1  4  1 -3 -2 -2
T  0 -1  0 -1 -1 -1 -1 -2 -2 -1 -1 -1 -1 -2 -1  1  5 -2 -2  0
W -3 -3 -4 -4 -2 -2 -3 -2 -2 -3 -2 -3 -1  1 -4 -3 -2 11  2 -3
Y -2 -2 -2 -3 -2 -1 -2 -3  2 -1 -1 -2 -1  3 -3 -2 -2  2  7 -1
V  0 -3 -3 -3 -1 -2 -2 -3 -3  3  1 -2  1 -1 -2 -2  0 -3 -1  4
"""

def _parse_blosum62():
    lines = [l for l in _BLOSUM62_RAW.strip().splitlines() if l.strip()]
    headers = lines[0].split()
    matrix = {}
    for line in lines[1:]:
        parts = line.split()
        aa1 = parts[0]
        for j, aa2 in enumerate(headers):
            matrix[(aa1, aa2)] = int(parts[j+1])
    return matrix

BLOSUM62 = _parse_blosum62()


def sw_score(seq1: str, seq2: str, gap_open: int = 10, gap_extend: int = 1) -> tuple[int, int]:
    """
    Smith-Waterman local alignment.
    Returns (score, end_pos_in_seq2) — approximates parasail's interface.
    Uses affine gap penalties.
    """
    m, n = len(seq1), len(seq2)
    if m == 0 or n == 0:
        return 0, 0

    # H = match scores, E = gap in seq2, F = gap in seq1
    NEG_INF = -10**6
    H = np.zeros((m+1, n+1), dtype=np.int32)
    E = np.full((m+1, n+1), NEG_INF, dtype=np.int32)
    F = np.full((m+1, n+1), NEG_INF, dtype=np.int32)

    best_score = 0
    best_j = 0

    for i in range(1, m+1):
        for j in range(1, n+1):
            aa1 = seq1[i-1]
            aa2 = seq2[j-1]
            sub = BLOSUM62.get((aa1, aa2), BLOSUM62.get((aa2, aa1), -4))

            E[i][j] = max(E[i][j-1] - gap_extend, H[i][j-1] - gap_open)
            F[i][j] = max(F[i-1][j] - gap_extend, H[i-1][j] - gap_open)
            H[i][j] = max(0, H[i-1][j-1] + sub, E[i][j], F[i][j])

            if H[i][j] > best_score:
                best_score = int(H[i][j])
                best_j = j

    return best_score, best_j
