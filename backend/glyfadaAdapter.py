"""Adapter for the glyfada evolutionary optimizer.

Delegates to glyfada's Python API (``glyfada.optimize()``) for binary
discovery, config construction, MPI launching, and result parsing.
The adapter handles FELsim-specific concerns: pickling the objective
function, copying the evaluator script, and translating results into
the format expected by beamOptimizer.

Author: Eremey Valetov
"""

import os
import shutil
import sys
import tempfile
import types

import numpy as np

from loggingConfig import get_logger_with_fallback

EVAL_SCRIPT = os.path.join(os.path.dirname(os.path.abspath(__file__)), "glyfada_eval.py")


class GlyfadaOptimizer:
    """Run glyfada evolutionary optimization on a pickled objective function.

    Parameters
    ----------
    objective_func : callable
        Objective function ``f(variable_vals) -> float`` (MSE to minimize).
    variable_names : list[str]
        Ordered variable names.
    bounds : list[tuple | dict]
        Per-variable bounds. Each entry is either a ``(min, max)`` tuple
        (continuous) or a dict in glyfada format::

            {"min": 0, "max": 10, "type": "integer"}

    default_values : list[float]
        Starting / default values for each variable.
    pop_size : int
        Population size per generation.
    max_gen : int
        Maximum number of generations.
    sigma : float
        Gaussian mutation std dev as fraction of parameter range.
    n_processes : int
        Number of MPI ranks.
    timeout_minutes : float
        Per-evaluation timeout.
    algorithm : str
        Glyfada algorithm: ``"auto"`` selects ULS for single-objective.
    n_objectives : int
        Number of objectives (1 = single-objective, >1 = multi-objective).
    constraints : list[dict] | None
        Glyfada constraint specs, e.g.
        ``[{"name": "stability", "type": "<=", "limit": 1e4}]``.
    seed_from_results : str | None
        Path to previous Pareto CSV for warm starting.
    callback : callable | None
        Progress callback ``(gen, hypervolume) -> None``.
    debug : bool | None
        Debug flag.
    extra_config : dict | None
        Additional glyfada config keys forwarded to the Python API.
    """

    def __init__(self, objective_func, variable_names, bounds, default_values,
                 pop_size=50, max_gen=100, sigma=0.05, n_processes=1,
                 timeout_minutes=2, algorithm="auto", n_objectives=1,
                 constraints=None, seed_from_results=None, callback=None,
                 debug=None, extra_config=None):
        self.objective_func = objective_func
        self.variable_names = variable_names
        self.bounds = bounds
        self.default_values = default_values
        self.pop_size = pop_size
        self.max_gen = max_gen
        self.sigma = sigma
        self.n_processes = n_processes
        self.timeout_minutes = timeout_minutes
        self.algorithm = algorithm
        self.n_objectives = n_objectives
        self.constraints = constraints
        self.seed_from_results = seed_from_results
        self.callback = callback
        self.extra_config = extra_config or {}
        self.logger, self.debug = get_logger_with_fallback(__name__, debug)

    def _pickle_objective(self, work_dir):
        """Serialize the objective function, variable names, and metadata."""
        import cloudpickle

        state = {
            "objective_func": self.objective_func,
            "variable_names": self.variable_names,
            "n_objectives": self.n_objectives,
        }
        if self.constraints:
            state["constraints"] = self.constraints
        pkl_path = os.path.join(work_dir, "_objective.pkl")
        with open(pkl_path, "wb") as f:
            cloudpickle.dump(state, f)
        return pkl_path

    def _build_parameters(self):
        """Build the glyfada parameter list from bounds + variable names."""
        parameters = []
        for name, bound, default in zip(
            self.variable_names, self.bounds, self.default_values
        ):
            if isinstance(bound, dict):
                p = dict(bound)
                p.setdefault("name", name)
                p.setdefault("default", default)
                parameters.append(p)
            else:
                lo, hi = bound
                parameters.append({
                    "name": name,
                    "min": lo,
                    "max": hi,
                    "default": default,
                })
        return parameters

    def _build_source_command(self):
        """Environment setup so the worker subprocess finds the right Python."""
        python_bin_dir = os.path.dirname(sys.executable)
        backend_dir = os.path.dirname(os.path.abspath(__file__))
        return (
            f"export PATH={python_bin_dir}:$PATH && "
            f"export PYTHONPATH={backend_dir}:$PYTHONPATH"
        )

    def optimize(self):
        """Run glyfada optimization.

        Returns
        -------
        result : SimpleNamespace
            For single-objective: scipy.optimize.OptimizeResult-compatible
            with ``x``, ``fun``, ``success``, ``message``, ``nfev``.
            For multi-objective: includes ``pareto_front``, ``objectives``,
            ``parameters``, ``summary``.
            Always includes ``work_dir``.
        """
        try:
            import glyfada
        except ImportError:
            raise ImportError(
                "glyfada Python package not installed. "
                "Install from ~/ML/paradiseo/glyfada/python: pip install -e ."
            )

        work_dir = tempfile.mkdtemp(prefix="glyfada_felsim_")
        self.logger.info(f"Glyfada working directory: {work_dir}")

        try:
            self._pickle_objective(work_dir)
            shutil.copy2(EVAL_SCRIPT, os.path.join(work_dir, "model.py"))

            parameters = self._build_parameters()

            # Assemble kwargs for glyfada.optimize()
            api_kwargs = {
                "parameters": parameters,
                "n_objectives": self.n_objectives,
                "evaluator": "dh",
                "program_file": "model.py",
                "workdir": work_dir,
                "algorithm": self.extra_config.get("algorithm", self.algorithm),
                "pop_size": self.pop_size,
                "max_gen": self.max_gen,
                "n_ranks": self.n_processes,
                "timeout_minutes": self.timeout_minutes,
                "sigma": self.sigma,
                "persistent_worker": True,
                "source_command": self._build_source_command(),
                "print_all_results": True,
                "omp_num_threads": max(1, (os.cpu_count() or 4) // self.n_processes),
            }

            if self.callback:
                api_kwargs["callback"] = self.callback
            if self.constraints:
                api_kwargs["constraints"] = self.constraints
            if self.seed_from_results:
                api_kwargs["seed_from_results"] = self.seed_from_results

            # Forward extra config (algorithm already handled above)
            for k, v in self.extra_config.items():
                if k != "algorithm" and k not in api_kwargs:
                    api_kwargs[k] = v

            if self.debug:
                import json
                safe = {k: v for k, v in api_kwargs.items() if k != "callback"}
                self.logger.debug(f"glyfada.optimize() kwargs: {json.dumps(safe, indent=2, default=str)}")

            glyfada_result = glyfada.optimize(**api_kwargs)
            return self._wrap_result(glyfada_result, work_dir)

        except Exception:
            self.logger.error(f"Glyfada failed. Working directory preserved: {work_dir}")
            raise

    def _wrap_result(self, glyfada_result, work_dir):
        """Convert glyfada.OptimizeResult to FELsim-expected format."""
        result = types.SimpleNamespace()
        result.pareto_front = glyfada_result.pareto_front
        result.objectives = glyfada_result.objectives
        result.parameters = glyfada_result.parameters
        result.summary = glyfada_result.summary
        result.work_dir = work_dir
        result.method = "glyfada"
        result.success = True
        result.message = "Glyfada optimization completed"
        result.nit = self.max_gen
        # Try to read actual nfev from glyfada summary; fall back to estimate.
        # With MPI multistart, evaluations are distributed across ranks, not multiplied.
        nfev_from_summary = getattr(glyfada_result, 'summary', {}).get('total_evaluations', None)
        result.nfev = nfev_from_summary if nfev_from_summary is not None else self.pop_size * self.max_gen

        if self.n_objectives == 1 and len(glyfada_result.objectives) > 0:
            # Glyfada maximizes (negated MSE), so objectives are negative.
            # Find the row whose objective is closest to zero (smallest MSE).
            obj_col = glyfada_result.objectives[:, 0]
            best_idx = int(np.argmax(obj_col))  # max of negated MSE = best
            best_params = glyfada_result.parameters[best_idx]

            result.x = best_params
            result.fun = float(-obj_col[best_idx])  # un-negate to get MSE
            result.best_params = dict(zip(self.variable_names, best_params))
            result.best_objective = result.fun

            self.logger.info(
                f"Glyfada optimization complete: MSE = {result.fun:.6e}, "
                f"params = {result.best_params}"
            )
        elif self.n_objectives == 1:
            self.logger.warning("Pareto front is empty (0 rows) — optimization likely failed")
            result.x = None
            result.fun = None
            result.success = False
            result.message = "Glyfada optimization returned empty Pareto front"
        else:
            result.x = None
            result.fun = None

        return result
