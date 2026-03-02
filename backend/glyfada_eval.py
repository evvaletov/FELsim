#!/usr/bin/env python3
"""Glyfada DH evaluator for FELsim beam optimization.

Supports two modes:
  1. Persistent worker (fast): glyfada_worker.py imports this module and calls
     evaluate(params, timeout) directly — no subprocess overhead.
  2. Subprocess (legacy): glyfada spawns
     ``python model.py --timeout T '{"var1": val1, ...}'``

Loads a pickled objective function from _objective.pkl in the working
directory, evaluates it with the parameter values.  Glyfada maximizes,
so MSE is negated.

Author: Eremey Valetov
"""

import json
import math
import sys
import os

# Cached state for persistent worker mode (loaded once, reused)
_state = None
_emit_constraints = None


def _load_state():
    """Load and cache the pickled objective function."""
    global _state, _emit_constraints
    if _state is None:
        pkl_path = os.path.join(os.path.dirname(os.path.abspath(__file__)), "_objective.pkl")
        try:
            import cloudpickle as pickle
        except ImportError:
            import pickle
        with open(pkl_path, "rb") as f:
            _state = pickle.load(f)

        flag_path = os.path.join(os.path.dirname(os.path.abspath(__file__)), ".emit_constraints")
        _emit_constraints = os.path.exists(flag_path)
    return _state


def evaluate(params, timeout=1200):
    """Persistent worker fast path — called directly, no subprocess.

    Parameters
    ----------
    params : dict
        Variable names → values, e.g. {"Ic": 4.2, "I": 1.5, ...}
    timeout : float
        Timeout in seconds (unused here; objective is fast).

    Returns
    -------
    dict
        {"objective_1": -mse, ...} with optional constraint values.
    """
    state = _load_state()
    variable_vals = [float(params[name]) for name in state["variable_names"]]

    try:
        mse = float(state["objective_func"](variable_vals))
    except Exception:
        mse = 1e6

    if not math.isfinite(mse):
        mse = 1e6

    result = {"objective_1": -mse}
    if _emit_constraints:
        result["constraint_stable"] = 0.0 if (math.isfinite(mse) and mse < 1e4) else 1.0
    return result


def main():
    """Subprocess CLI entry point (legacy DH protocol)."""
    args = list(sys.argv[1:])

    if "--timeout" in args:
        idx = args.index("--timeout")
        args.pop(idx)
        if idx < len(args):
            args.pop(idx)

    if not args:
        print("F no JSON argument", file=sys.stderr)
        sys.exit(1)

    params = json.loads(args[0])
    result = evaluate(params)
    print(json.dumps(result))


if __name__ == "__main__":
    main()
