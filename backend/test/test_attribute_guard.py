"""I6.2b: Attribute typo guard for beamline element classes.

Validates that all setattr targets in the codebase reference legitimate
element attributes.  Catches silent bugs where a misspelled attribute
(e.g. 'currnet' instead of 'current') creates a new attribute that the
physics routines never read.

Author: Eremey Valetov
"""

import ast
import sys
from pathlib import Path

import pytest

_BACKEND = Path(__file__).resolve().parent.parent
if str(_BACKEND) not in sys.path:
    sys.path.insert(0, str(_BACKEND))

from beamline import lattice, driftLattice, qpfLattice, qpdLattice, dipole, dipole_wedge

ELEMENT_CLASSES = [driftLattice, qpfLattice, qpdLattice, dipole, dipole_wedge]


def _collect_instance_attrs(cls):
    """Instantiate element class and return all attribute names."""
    constructors = {
        driftLattice: lambda: driftLattice(0.1),
        qpfLattice: lambda: qpfLattice(current=1.0, length=0.1),
        qpdLattice: lambda: qpdLattice(current=1.0, length=0.1),
        dipole: lambda: dipole(length=0.1, angle=1.5),
        dipole_wedge: lambda: dipole_wedge(length=0.01, angle=1.0),
    }
    obj = constructors[cls]()
    return set(vars(obj).keys())


def _collect_class_attrs(cls):
    """Collect class-level attributes (not inherited from object)."""
    attrs = set()
    for klass in cls.__mro__:
        if klass is object:
            break
        for k, v in vars(klass).items():
            if not k.startswith('_') and not callable(v):
                attrs.add(k)
    return attrs


def _valid_attributes(cls):
    """Full set of legitimate attribute names for an element class."""
    return _collect_instance_attrs(cls) | _collect_class_attrs(cls)


# Build the complete whitelist: union of all element class attributes
ALL_VALID_ATTRS = set()
for cls in ELEMENT_CLASSES:
    ALL_VALID_ATTRS |= _valid_attributes(cls)


class TestElementAttributes:
    """Verify element class attribute inventories."""

    @pytest.mark.parametrize("cls", ELEMENT_CLASSES, ids=lambda c: c.__name__)
    def test_base_attributes_present(self, cls):
        """Every element must have the core lattice attributes."""
        attrs = _valid_attributes(cls)
        required = {'name', 'E', 'E0', 'Q', 'M', 'C', 'f', 'gamma', 'beta',
                     'length', 'color', 'chromatic', 'aperture_x', 'aperture_y'}
        missing = required - attrs
        assert not missing, f"{cls.__name__} missing base attributes: {missing}"

    def test_quad_has_current_and_G(self):
        for cls in (qpfLattice, qpdLattice):
            attrs = _valid_attributes(cls)
            assert 'current' in attrs
            assert 'G' in attrs

    def test_dipole_has_angle(self):
        attrs = _valid_attributes(dipole)
        assert 'angle' in attrs

    def test_dipole_wedge_has_wedge_attrs(self):
        attrs = _valid_attributes(dipole_wedge)
        for attr in ('angle', 'dipole_length', 'dipole_angle', 'pole_gap'):
            assert attr in attrs, f"dipole_wedge missing '{attr}'"


class TestSetAttrTargets:
    """Scan source code for setattr() calls and validate attribute names."""

    def _find_setattr_calls(self, filepath):
        """Parse a Python file and return (line, attr_name) for setattr calls
        where the attribute name is a string literal."""
        source = filepath.read_text()
        try:
            tree = ast.parse(source)
        except SyntaxError:
            return []

        results = []
        for node in ast.walk(tree):
            if isinstance(node, ast.Call) and isinstance(node.func, ast.Name):
                if node.func.id == 'setattr' and len(node.args) >= 2:
                    attr_arg = node.args[1]
                    if isinstance(attr_arg, ast.Constant) and isinstance(attr_arg.value, str):
                        results.append((node.lineno, attr_arg.value))
        return results

    def test_optimizer_setattr_uses_valid_attrs(self):
        """beamOptimizer.py setattr targets should be valid element attributes."""
        path = _BACKEND / 'beamOptimizer.py'
        if not path.exists():
            pytest.skip("beamOptimizer.py not found")
        calls = self._find_setattr_calls(path)
        for line, attr in calls:
            assert attr in ALL_VALID_ATTRS, (
                f"beamOptimizer.py:{line}: setattr target '{attr}' "
                f"is not a known element attribute"
            )

    def test_felapi_setattr_uses_valid_attrs(self):
        """felAPI.py setattr targets should be valid element attributes."""
        path = _BACKEND / 'felAPI.py'
        if not path.exists():
            pytest.skip("felAPI.py not found")
        calls = self._find_setattr_calls(path)
        for line, attr in calls:
            assert attr in ALL_VALID_ATTRS, (
                f"felAPI.py:{line}: setattr target '{attr}' "
                f"is not a known element attribute"
            )

    def test_algebraic_setattr_uses_valid_attrs(self):
        """AlgebraicOptimization.py setattr targets should be valid element attributes."""
        path = _BACKEND / 'AlgebraicOptimization.py'
        if not path.exists():
            pytest.skip("AlgebraicOptimization.py not found")
        calls = self._find_setattr_calls(path)
        for line, attr in calls:
            assert attr in ALL_VALID_ATTRS, (
                f"AlgebraicOptimization.py:{line}: setattr target '{attr}' "
                f"is not a known element attribute"
            )

    def test_all_backend_setattr_targets(self):
        """Scan all backend .py files for setattr with literal attribute names."""
        unknown = []
        for pyfile in sorted(_BACKEND.glob('*.py')):
            calls = self._find_setattr_calls(pyfile)
            for line, attr in calls:
                if attr not in ALL_VALID_ATTRS and not attr.startswith('_'):
                    unknown.append(f"{pyfile.name}:{line}: '{attr}'")
        if unknown:
            msg = "setattr targets not in element attribute whitelist:\n"
            msg += "\n".join(f"  {u}" for u in unknown)
            pytest.fail(msg)


class TestAttributeConsistency:
    """Cross-check that element classes don't silently shadow attributes."""

    def test_no_typo_variants(self):
        """Check for suspicious near-duplicates in attribute names."""
        import difflib
        for cls in ELEMENT_CLASSES:
            attrs = sorted(_valid_attributes(cls))
            for i, a in enumerate(attrs):
                close = difflib.get_close_matches(a, attrs[:i] + attrs[i+1:],
                                                  n=1, cutoff=0.85)
                if close:
                    # Known legitimate pairs
                    legit_pairs = {
                        frozenset({'E', 'E0'}),
                        frozenset({'M', 'M_AMU'}),
                        frozenset({'startPos', 'endPos'}),
                        frozenset({'aperture_x', 'aperture_y'}),
                        frozenset({'dipole_angle', 'dipole_length'}),
                    }
                    if frozenset({a, close[0]}) not in legit_pairs:
                        pytest.fail(
                            f"{cls.__name__}: '{a}' suspiciously similar "
                            f"to '{close[0]}'"
                        )
