"""
Backward-compatible re-export from the tracked-mapping package.

Install:  pip install -e ~/Python/TrackedDict

Author: Eremey Valetov
"""

from tracked_dict import TrackedDict, TrackedList

__all__ = ["TrackedDict", "TrackedList"]
