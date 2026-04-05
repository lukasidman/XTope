"""
N-tag stripping module.
Detects and removes an N-terminal purification tag from antigen sequences.
Uses exact string matching first, then Smith-Waterman local alignment as fallback.

N-tag stripping is opt-in — sequences are used as-is unless a tag sequence is
explicitly provided via the --tag CLI flag or the strip_tag() function.
"""

from .sw_fallback import sw_score as _sw_fallback

SEARCH_WINDOW = 60   # aa from N-terminus to search
MIN_TAG_SCORE = 20   # minimum SW score to confirm tag match
GAP_OPEN      = 10
GAP_EXTEND    = 1


def strip_tag_exact(sequence: str, tag: str) -> tuple[str, bool]:
    """Remove tag from sequence by exact string match.

    Args:
        sequence: Amino acid sequence (may contain tag at N-terminus).
        tag:      Tag sequence to search for and remove.

    Returns:
        Tuple of (stripped_sequence, tag_found).
    """
    seq_upper = sequence.upper()
    tag_upper = tag.upper()
    idx = seq_upper.find(tag_upper)
    if idx != -1:
        return sequence[idx + len(tag):], True
    return sequence, False


def strip_tag_sw(sequence: str, tag: str) -> tuple[str, bool, float]:
    """Remove tag from sequence using Smith-Waterman alignment as fallback.

    Searches only the first SEARCH_WINDOW residues from the N-terminus.

    Args:
        sequence: Amino acid sequence.
        tag:      Tag sequence to locate and remove.

    Returns:
        Tuple of (stripped_sequence, tag_found, sw_score).
    """
    window = sequence[:SEARCH_WINDOW].upper()
    tag_upper = tag.upper()

    score, end_in_window = _sw_fallback(tag_upper, window, GAP_OPEN, GAP_EXTEND)

    if score < MIN_TAG_SCORE:
        return sequence, False, score

    stripped = sequence[end_in_window:]
    return stripped, True, score


def strip_tag(sequence: str, tag: str, method: str = "sw") -> dict:
    """Strip an N-terminal tag from a sequence.

    Tries exact match first; falls back to Smith-Waterman if method='sw'.

    Args:
        sequence: Amino acid sequence to process.
        tag:      Tag sequence to detect and remove.
        method:   'sw' (default) uses SW fallback after exact match fails;
                  'exact' skips SW fallback.

    Returns:
        Dict with keys:
            original    - input sequence (uppercased)
            stripped    - sequence after tag removal (or original if not found)
            tag_found   - True if tag was detected and removed
            method_used - 'exact', 'sw', or 'none'
            score       - SW score if SW was used, else None
    """
    sequence = sequence.strip().upper()

    stripped_exact, found_exact = strip_tag_exact(sequence, tag)
    if found_exact:
        return {"original": sequence, "stripped": stripped_exact,
                "tag_found": True, "method_used": "exact", "score": None}

    if method == "sw":
        stripped_sw, found_sw, score = strip_tag_sw(sequence, tag)
        return {"original": sequence, "stripped": stripped_sw,
                "tag_found": found_sw, "method_used": "sw" if found_sw else "none",
                "score": score}

    return {"original": sequence, "stripped": sequence,
            "tag_found": False, "method_used": "none", "score": None}
