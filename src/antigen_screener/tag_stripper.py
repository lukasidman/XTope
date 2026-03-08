"""
Tag stripping module.
Detects and removes the His6-ABP tag from the N-terminus of antigen sequences.
Uses Smith-Waterman local alignment for robust detection (falls back to
pure Python if parasail is not installed).
"""

try:
    import parasail as _parasail
    _HAS_PARASAIL = True
except ImportError:
    _HAS_PARASAIL = False

from .sw_fallback import sw_score as _sw_fallback

# ── Update this to your real His6-ABP tag sequence ──────────────────────────
DEFAULT_TAG = "MHHHHHHGSSG"

SEARCH_WINDOW = 60   # aa from N-terminus to search
MIN_TAG_SCORE = 20   # minimum SW score to confirm tag match
GAP_OPEN      = 10
GAP_EXTEND    = 1


def strip_tag_exact(sequence: str, tag: str = DEFAULT_TAG):
    seq_upper = sequence.upper()
    tag_upper = tag.upper()
    idx = seq_upper.find(tag_upper)
    if idx != -1:
        return sequence[idx + len(tag):], True
    return sequence, False


def strip_tag_sw(sequence: str, tag: str = DEFAULT_TAG):
    window = sequence[:SEARCH_WINDOW].upper()
    tag_upper = tag.upper()

    if _HAS_PARASAIL:
        r = _parasail.sw_trace_striped_16(
            tag_upper, window, GAP_OPEN, GAP_EXTEND, _parasail.blosum62
        )
        score = r.score
        end_in_window = r.end_query + 1
    else:
        score, end_in_window = _sw_fallback(tag_upper, window, GAP_OPEN, GAP_EXTEND)

    if score < MIN_TAG_SCORE:
        return sequence, False, score

    stripped = sequence[end_in_window:]
    return stripped, True, score


def strip_tag(sequence: str, tag: str = DEFAULT_TAG, method: str = "sw") -> dict:
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


def set_tag(new_tag: str):
    global DEFAULT_TAG
    DEFAULT_TAG = new_tag.strip().upper()
