"""
Regression tests for beamOptimizer convergence.

Verifies that the optimizer drives MSE toward zero on a simple FODO cell
problem with known solution.

Author: Eremey Valetov
"""
import numpy as np
import pytest

from beamline import driftLattice, qpfLattice, qpdLattice
from beamOptimizer import beamOptimizer
from ebeam import beam

KE = 45.0  # MeV


def _make_fodo_beamline():
    """Drift-QF-Drift-QD-Drift at 45 MeV."""
    elements = [
        driftLattice(0.3),
        qpfLattice(0.1, 2.0),
        driftLattice(0.4),
        qpdLattice(0.1, 2.0),
        driftLattice(0.3),
    ]
    for el in elements:
        el.setE(KE)
    return elements


def _make_particles(n=500):
    ebeam = beam()
    std_dev = [1.0, 0.1, 1.0, 0.1, 1.0, 0.5]
    return np.array(ebeam.gen_6d_gaussian(0, std_dev, n))


class TestOptimizerConvergence:
    def test_nelder_mead_reduces_mse(self):
        """Optimizer should reduce MSE from its initial value."""
        beamline = _make_fodo_beamline()
        particles = _make_particles()
        ebeam_obj = beam()

        opt = beamOptimizer(beamline, particles)

        # Optimize QF current to hit a target x-envelope at the end
        segmentVar = {
            1: ["current_QF", "current", lambda x: x],
        }
        objectives = {
            4: [{"measure": ["x", "envelope"],
                 "goal": 1.5, "weight": 1.0}],
        }
        startPoint = {
            "current_QF": {"start": 1.0, "bounds": (0.1, 10.0)},
        }

        result = opt.calc(
            method='Nelder-Mead',
            segmentVar=segmentVar,
            startPoint=startPoint,
            objectives=objectives,
        )

        assert result.fun < opt.plotMSE[0], \
            f"Final MSE {result.fun} should be less than initial {opt.plotMSE[0]}"
        assert len(opt.plotIterate) > 1, "Optimizer should run multiple iterations"

    def test_zero_goals_returns_inf(self):
        """Optimizer with no objectives should return inf."""
        beamline = _make_fodo_beamline()
        particles = _make_particles(50)

        opt = beamOptimizer(beamline, particles)
        opt.segmentVar = {}
        opt.objectives = {}
        opt.variablesToOptimize = []
        opt.trackGoals = {}
        opt.trackVariables = []
        opt.plotMSE = []
        opt.plotIterate = []
        opt.iterationTrack = 0

        result = opt._optiSpeed([])
        assert result == np.inf

    def test_invalid_index_raises(self):
        """Out-of-bounds segment index should raise IndexError."""
        beamline = _make_fodo_beamline()
        particles = _make_particles(50)

        opt = beamOptimizer(beamline, particles)

        segmentVar = {
            99: ["current", "current", lambda x: x],
        }
        objectives = {4: [{"measure": ["x", "envelope"], "goal": 1.0, "weight": 1.0}]}
        startPoint = {"current": {"start": 1.0, "bounds": (0.1, 10.0)}}

        with pytest.raises(IndexError):
            opt.calc('Nelder-Mead', segmentVar, startPoint, objectives)

    def test_invalid_objective_index_raises(self):
        """Out-of-bounds objective index should raise TypeError."""
        beamline = _make_fodo_beamline()
        particles = _make_particles(50)

        opt = beamOptimizer(beamline, particles)

        segmentVar = {1: ["current_QF", "current", lambda x: x]}
        objectives = {99: [{"measure": ["x", "envelope"], "goal": 1.0, "weight": 1.0}]}
        startPoint = {"current_QF": {"start": 1.0, "bounds": (0.1, 10.0)}}

        with pytest.raises(TypeError):
            opt.calc('Nelder-Mead', segmentVar, startPoint, objectives)
