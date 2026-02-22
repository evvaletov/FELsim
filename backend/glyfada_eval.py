#!/usr/bin/env python3
"""Glyfada DH evaluator for FELsim beam optimization.

Called by glyfada as:
    python model.py --timeout T '{"var1": val1, "var2": val2, ...}'

Loads a pickled objective function from _objective.pkl in the working
directory, evaluates it with the parameter values, and prints the result
as JSON to stdout.  Glyfada maximizes, so MSE is negated.

Author: Eremey Valetov
"""

import json
import sys
import os


def main():
    args = list(sys.argv[1:])

    # Consume --timeout argument (glyfada always passes it)
    if "--timeout" in args:
        idx = args.index("--timeout")
        args.pop(idx)
        if idx < len(args):
            args.pop(idx)

    if not args:
        print("F no JSON argument", file=sys.stderr)
        sys.exit(1)

    params = json.loads(args[0])

    # Load pickled objective
    pkl_path = os.path.join(os.path.dirname(os.path.abspath(__file__)), "_objective.pkl")
    try:
        import cloudpickle as pickle
    except ImportError:
        import pickle

    with open(pkl_path, "rb") as f:
        state = pickle.load(f)

    objective_func = state["objective_func"]
    variable_names = state["variable_names"]

    # Build variable values in the correct order
    variable_vals = [float(params[name]) for name in variable_names]

    try:
        mse = float(objective_func(variable_vals))
    except Exception as e:
        print(f"F evaluation failed: {e}", file=sys.stderr)
        sys.exit(1)

    # Glyfada maximizes — negate MSE
    print(json.dumps({"objective_1": -mse}))


if __name__ == "__main__":
    main()
