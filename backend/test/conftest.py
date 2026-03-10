"""
Shared pytest configuration and fixtures for FELsim tests.

Author: Eremey Valetov
"""

import sys
from pathlib import Path

import pytest

# Ensure backend/ is importable for all tests
_BACKEND = Path(__file__).resolve().parent.parent
if str(_BACKEND) not in sys.path:
    sys.path.insert(0, str(_BACKEND))


def pytest_configure(config):
    config.addinivalue_line("markers", "visual: tests requiring interactive display")
