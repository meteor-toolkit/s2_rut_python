"""
Internal helper to make vendored third-party S2-RUT code importable.

This module adjusts sys.path to include the vendored subtree.
Do not modify without updating the subtree integration.
"""

import os
import sys

_BASE_DIR = os.path.abspath(os.path.dirname(os.path.dirname(__file__)))
_VENDOR_DIR = os.path.join(
    _BASE_DIR,
    "third-party",
    "Source",
)

if _VENDOR_DIR not in sys.path:
    sys.path.insert(0, _VENDOR_DIR)
