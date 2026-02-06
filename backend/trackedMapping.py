"""
Backward-compatible re-export from the tracked-mapping package.

Install:  pip install -e ~/Python/TrackedDict

Author: Eremey Valetov
"""

from tracked_mapping import TrackedDict, TrackedList

__all__ = ["TrackedDict", "TrackedList"]
