"""Cross-reactivity antigen screener.

Identifies cross-reactive antigens in a database of short antigen sequences
(50-150 amino acids) by combining k-mer pre-filtering with Smith-Waterman
alignment and optional physicochemical similarity scoring.
"""

__version__ = "0.1.0"

import sys

if sys.version_info >= (3, 13):
    import warnings

    warnings.warn(
        "Python 3.13+ detected. The 'parasail' library is not available for this version. "
        "Install Python 3.12 for best performance. The pure-Python fallback will be used.",
        UserWarning,
        stacklevel=2,
    )
