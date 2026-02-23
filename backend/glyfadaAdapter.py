"""Adapter for the glyfada evolutionary optimizer.

Wraps the glyfada C++ MPI binary (~/ML/paradiseo/glyfada/build/optimiser)
using its DH (generic) evaluator protocol.  For each fitness evaluation,
glyfada spawns ``python model.py --timeout T '<json_params>'`` and reads
the result from stdout.

Author: Eremey Valetov
"""

import json
import os
import shutil
import subprocess
import tempfile
import glob as globmod
import numpy as np
from loggingConfig import get_logger_with_fallback

GLYFADA_BINARY = os.path.expanduser(
    "~/ML/paradiseo/glyfada/build/optimiser"
)


class GlyfadaOptimizer:
    """Run glyfada evolutionary optimization on a pickled objective function.

    Parameters
    ----------
    objective_func : callable
        Objective function ``f(variable_vals) -> float`` (MSE to minimize).
    variable_names : list[str]
        Ordered variable names.
    bounds : list[tuple]
        (min, max) for each variable.
    default_values : list[float]
        Starting / default values for each variable.
    pop_size : int
        Population size per generation.
    max_gen : int
        Maximum number of generations.
    sigma : float
        Gaussian mutation std dev as fraction of parameter range.
    n_processes : int
        Number of MPI processes (multistart ranks).
    timeout_minutes : float
        Per-evaluation timeout.
    algorithm : str
        Glyfada algorithm: ``"auto"`` selects ULS for single-objective.
    debug : bool or None
        Debug flag.
    """

    def __init__(self, objective_func, variable_names, bounds, default_values,
                 pop_size=50, max_gen=100, sigma=0.05, n_processes=1,
                 timeout_minutes=2, algorithm="auto", debug=None):
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
        self.logger, self.debug = get_logger_with_fallback(__name__, debug)

    def _build_config(self, work_dir):
        """Build glyfada JSON configuration."""
        parameters = []
        for name, (lo, hi), default in zip(
            self.variable_names, self.bounds, self.default_values
        ):
            parameters.append({
                "name": name,
                "type": "continuous",
                "min_value": lo,
                "max_value": hi,
                "default_value": default
            })

        # Make subprocess Python match the running interpreter
        import sys
        python_bin_dir = os.path.dirname(sys.executable)
        backend_dir = os.path.dirname(os.path.abspath(__file__))
        source_cmd = (
            f"export PATH={python_bin_dir}:$PATH && "
            f"export PYTHONPATH={backend_dir}"
        )

        config = {
            "mode": "multistart",
            "evaluator": "dh",
            "algorithm": self.algorithm,
            "n_objectives": 1,
            "program_file": "model.py",
            "config_file": "parameters.json",
            "program_directory": work_dir,
            "source_command": source_cmd,
            "timeout_minutes": self.timeout_minutes,
            "pop_size": self.pop_size,
            "max_gen": self.max_gen,
            "run_limit_type": "max_gen",
            "sigma": self.sigma,
            "p_cross": 0.25,
            "p_mut": 0.35,
            "p_change": 1.0,
            "eta_c": 30.0,
            "mut_epsilon": 0.01,
            "print_all_results": True,
            "timein_seconds": 0,
            "omp_num_threads": os.cpu_count() or 4,
            "tournament_size": min(15, self.pop_size),
            "parameters": parameters,
        }
        return config

    def _pickle_objective(self, work_dir):
        """Serialize the objective function and variable names."""
        try:
            import cloudpickle as pickle
        except ImportError:
            import pickle

        state = {
            "objective_func": self.objective_func,
            "variable_names": self.variable_names,
        }
        pkl_path = os.path.join(work_dir, "_objective.pkl")
        with open(pkl_path, "wb") as f:
            pickle.dump(state, f)
        return pkl_path

    def _parse_results(self, work_dir):
        """Parse the best solution from glyfada CSV output.

        Returns (best_objective, best_params_dict) or raises on failure.
        """
        # Look for all_evaluated_solutions or pareto_frontier CSVs
        patterns = [
            os.path.join(work_dir, "all_evaluated_solutions_*.csv"),
            os.path.join(work_dir, "pareto_frontier_*.csv"),
        ]

        csv_files = []
        for pattern in patterns:
            csv_files.extend(globmod.glob(pattern))

        if not csv_files:
            raise RuntimeError(
                f"No glyfada output CSV found in {work_dir}. "
                f"Check glyfada stderr for errors."
            )

        best_obj = -np.inf
        best_params = None

        for csv_file in csv_files:
            with open(csv_file, 'r') as f:
                for line in f:
                    line = line.strip()
                    if not line or line.startswith('#'):
                        continue
                    values = [float(v) for v in line.split(',')]
                    # Format: objective_1, param_1, param_2, ...
                    obj_val = values[0]
                    param_vals = values[1:]
                    if obj_val > best_obj:
                        best_obj = obj_val
                        best_params = param_vals

        if best_params is None:
            raise RuntimeError("No valid solutions found in glyfada CSV output")

        # Negate back (glyfada maximizes, we negated MSE)
        best_mse = -best_obj
        params_dict = dict(zip(self.variable_names, best_params))

        return best_mse, params_dict

    def optimize(self):
        """Run glyfada optimization.

        Returns
        -------
        result : object
            scipy.optimize.OptimizeResult-compatible object with attributes:
            x, fun, success, message, nfev (estimated).
        """
        if not os.path.isfile(GLYFADA_BINARY):
            raise FileNotFoundError(
                f"Glyfada binary not found at {GLYFADA_BINARY}. "
                f"Build glyfada first: cd ~/ML/paradiseo/glyfada/build && make"
            )

        work_dir = tempfile.mkdtemp(prefix="glyfada_felsim_")
        self.logger.info(f"Glyfada working directory: {work_dir}")

        try:
            # Copy evaluator script
            eval_script = os.path.join(
                os.path.dirname(os.path.abspath(__file__)), "glyfada_eval.py"
            )
            shutil.copy(eval_script, os.path.join(work_dir, "model.py"))

            # Pickle objective
            self._pickle_objective(work_dir)

            # Write config
            config = self._build_config(work_dir)
            config_path = os.path.join(work_dir, "parameters.json")
            with open(config_path, 'w') as f:
                json.dump(config, f, indent=2)

            if self.debug:
                self.logger.debug(f"Glyfada config: {json.dumps(config, indent=2)}")

            # Run glyfada
            cmd = [
                "mpirun", "-np", str(self.n_processes),
                GLYFADA_BINARY, f"--config={config_path}"
            ]
            self.logger.info(f"Running: {' '.join(cmd)}")

            proc = subprocess.run(
                cmd, cwd=work_dir,
                capture_output=True, text=True,
                timeout=self.timeout_minutes * 60 * (self.max_gen + 10)
            )

            if proc.returncode != 0:
                self.logger.error(f"Glyfada stderr:\n{proc.stderr}")
                raise RuntimeError(
                    f"Glyfada exited with code {proc.returncode}: {proc.stderr[:500]}"
                )

            if self.debug and proc.stdout:
                self.logger.debug(f"Glyfada stdout:\n{proc.stdout[:2000]}")

            # Parse results
            best_mse, params_dict = self._parse_results(work_dir)
            x_best = [params_dict[name] for name in self.variable_names]

            self.logger.info(
                f"Glyfada optimization complete: MSE = {best_mse:.6e}, "
                f"params = {params_dict}"
            )

            # Return scipy-compatible result
            from types import SimpleNamespace
            result = SimpleNamespace(
                x=np.array(x_best),
                fun=best_mse,
                success=True,
                message="Glyfada optimization completed",
                nfev=self.pop_size * self.max_gen * self.n_processes,
                nit=self.max_gen,
                method="glyfada",
                work_dir=work_dir
            )
            return result

        except Exception:
            self.logger.error(f"Glyfada failed. Working directory preserved: {work_dir}")
            raise
