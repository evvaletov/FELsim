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
_logger_cache = None


def _get_logger():
    """Lazily obtain a logger (avoid import cost on every evaluation)."""
    global _logger_cache
    if _logger_cache is None:
        import logging
        _logger_cache = logging.getLogger(__name__)
    return _logger_cache


def _load_state():
    """Load and cache the pickled objective function."""
    global _state
    if _state is None:
        pkl_path = os.path.join(os.path.dirname(os.path.abspath(__file__)), "_objective.pkl")
        import cloudpickle
        with open(pkl_path, "rb") as f:
            _state = cloudpickle.load(f)
    return _state


def _format_objectives(raw, n_objectives):
    """Normalize objective function return value to a dict.

    Accepts:
      - float/int → single objective
      - list/tuple → positional objectives
      - dict → values extracted in insertion order, keys replaced with objective_N

    Each value is negated (glyfada maximizes). Non-finite values after negation
    are replaced with -1e6 (penalty in glyfada's maximization space).
    """
    if isinstance(raw, dict):
        result = {f"objective_{i+1}": -float(v) for i, (_, v) in enumerate(raw.items())}
    elif isinstance(raw, (list, tuple)):
        result = {f"objective_{i+1}": -float(v) for i, v in enumerate(raw)}
    else:
        result = {"objective_1": -float(raw)}

    # Sanitize non-finite values (NaN, Inf) to a penalty
    for k, v in result.items():
        if not math.isfinite(v):
            result[k] = -1e6

    n_obj_keys = sum(1 for k in result if k.startswith("objective_"))
    if n_obj_keys != n_objectives:
        _logger = _get_logger()
        _logger.warning(
            f"Objective count mismatch: expected {n_objectives}, got {n_obj_keys} "
            f"from return type {type(raw).__name__}"
        )

    return result


def evaluate(params, timeout=1200):
    """Persistent worker fast path — called directly, no subprocess.

    Parameters
    ----------
    params : dict
        Variable names -> values, e.g. {"Ic": 4.2, "I": 1.5, ...}
    timeout : float
        Timeout in seconds (unused here; objective is fast).

    Returns
    -------
    dict
        {"objective_1": -mse, ...} with optional constraint values.
    """
    state = _load_state()
    variable_names = state["variable_names"]
    n_objectives = state.get("n_objectives", 1)
    variable_vals = [float(params[name]) for name in variable_names]

    raw = state["objective_func"](variable_vals)
    result = _format_objectives(raw, n_objectives)

    # Binary feasibility constraint placeholder.
    # The objective function already encodes stability penalties (large MSE
    # for unstable solutions). This exposes a binary feasibility signal to
    # glyfada's constraint mechanism: feasible (0) if the un-negated first
    # objective is finite and below the limit, infeasible (1) otherwise.
    # A single constraint value is emitted regardless of how many constraint
    # specs are provided — all specs check the same underlying condition.
    constraints = state.get("constraints")
    if constraints:
        obj_val = -result.get("objective_1", float("nan"))  # un-negate
        limit = constraints[0].get("limit", 1e4)
        feasible = math.isfinite(obj_val) and obj_val < limit
        for i in range(1, len(constraints) + 1):
            result[f"constraint_{i}"] = 0.0 if feasible else 1.0

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
