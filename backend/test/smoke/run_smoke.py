#!/usr/bin/env python3
"""
Dynamic CLI runner for FELsim smoke tests.

Displays a spinner while each test runs, then shows pass/fail/skip status
with ANSI colours. Returns exit code 0 if all passed, 1 if any failed.

Usage:
    python run_smoke.py

Author: Eremey Valetov
"""

import importlib
import sys
import time
import threading
import itertools

# ---------------------------------------------------------------------------
# ANSI helpers
# ---------------------------------------------------------------------------

_IS_TTY = sys.stdout.isatty()

def _ansi(code):
    return f"\033[{code}m" if _IS_TTY else ""

GREEN  = _ansi("32")
RED    = _ansi("31")
YELLOW = _ansi("33")
BOLD   = _ansi("1")
RESET  = _ansi("0")

SPINNER_CHARS = "⠋⠙⠹⠸⠼⠴⠦⠧⠇⠏"


# ---------------------------------------------------------------------------
# Test discovery
# ---------------------------------------------------------------------------

def _collect_tests():
    """Import the test module and return [(name, func), ...]."""
    mod = importlib.import_module("test_smoke")
    tests = []
    for name in sorted(dir(mod)):
        if name.startswith("test_"):
            tests.append((name, getattr(mod, name)))
    return tests


# ---------------------------------------------------------------------------
# Runner with live spinner
# ---------------------------------------------------------------------------

class _SkipTest(Exception):
    """Raised when a test is skipped (mirrors pytest.skip)."""

def _run_test(func):
    """Run a single test function. Returns (status, message).
    status: 'pass', 'fail', or 'skip'
    """
    try:
        func()
        return "pass", ""
    except _SkipTest as e:
        return "skip", str(e)
    except Exception as e:
        # pytest.skip raises Skipped, which is a subclass of Exception
        ename = type(e).__name__
        if ename in ("Skipped", "skip"):
            return "skip", str(e)
        return "fail", f"{ename}: {e}"


def _spinner_loop(label, stop_event):
    """Background thread that animates a spinner on the current line."""
    for ch in itertools.cycle(SPINNER_CHARS):
        if stop_event.is_set():
            break
        sys.stdout.write(f"\r{ch} Running {label}...")
        sys.stdout.flush()
        time.sleep(0.08)


def _run_all(tests):
    passed = failed = skipped = 0
    results = []

    for name, func in tests:
        if _IS_TTY:
            stop = threading.Event()
            spinner = threading.Thread(target=_spinner_loop, args=(name, stop), daemon=True)
            spinner.start()

        status, msg = _run_test(func)

        if _IS_TTY:
            stop.set()
            spinner.join()
            sys.stdout.write("\r\033[K")  # clear spinner line
            sys.stdout.flush()

        if status == "pass":
            passed += 1
            print(f"  {GREEN}✓{RESET} {name}")
        elif status == "skip":
            skipped += 1
            reason = f": {msg}" if msg else ""
            print(f"  {YELLOW}-{RESET} {name}{YELLOW}{reason}{RESET}")
        else:
            failed += 1
            print(f"  {RED}✗{RESET} {name}: {RED}{msg}{RESET}")

        results.append((name, status, msg))

    # Summary
    print()
    parts = [f"{GREEN}{passed} passed{RESET}"]
    if failed:
        parts.append(f"{RED}{failed} failed{RESET}")
    if skipped:
        parts.append(f"{YELLOW}{skipped} skipped{RESET}")
    print(f"  {BOLD}{'  '.join(parts)}{RESET}")

    return failed == 0


# ---------------------------------------------------------------------------
# Entry point
# ---------------------------------------------------------------------------

def main():
    # Ensure we can import test_smoke and backend modules
    smoke_dir = str(__import__("pathlib").Path(__file__).resolve().parent)
    backend_dir = str(__import__("pathlib").Path(__file__).resolve().parent.parent.parent)
    for p in (smoke_dir, backend_dir):
        if p not in sys.path:
            sys.path.insert(0, p)

    print(f"\n  {BOLD}FELsim Smoke Tests{RESET}\n")

    tests = _collect_tests()
    if not tests:
        print("  No tests found.")
        sys.exit(1)

    success = _run_all(tests)
    print()
    sys.exit(0 if success else 1)


if __name__ == "__main__":
    main()
