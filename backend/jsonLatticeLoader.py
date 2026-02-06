"""
Load beamline lattice configurations from JSON files.

Produces two output formats:
  - parse_beamline()  → list of dicts (BeamlineBuilder-compatible)
  - create_beamline() → list of beamline.py class instances

Uses tracked_mapping.TrackedDict for unconsumed-data reporting and
optionally validates against the JSON Schema in var/lattice_schema_v1.json.

Author: Eremey Valetov
"""

import json
import os
from pathlib import Path

from tracked_mapping import TrackedDict
from beamline import driftLattice, qpfLattice, qpdLattice, dipole, dipole_wedge
from loggingConfig import get_logger_with_fallback

SUPPORTED_FORMAT_VERSION = 1

# Map JSON element types to internal short names used by BeamlineBuilder / adapters.
_TYPE_ALIASES = {
    "QUADRUPOLE": None,  # resolved via polarity
    "QPF": "QPF",
    "QPD": "QPD",
    "DIPOLE": "DPH",
    "DPH": "DPH",
    "DIPOLE_WEDGE": "DPW",
    "DPW": "DPW",
    "SOLENOID": "SOL",
    "SOL": "SOL",
    "RF_CAVITY": "RFC",
    "RFC": "RFC",
    "SEXTUPOLE": "SXT",
    "SXT": "SXT",
    "UNDULATOR": "UND",
    "UND": "UND",
    "BPM": "BPM",
    "OTR": "OTR",
    "CORRECTOR_V": "STV",
    "STV": "STV",
    "CORRECTOR_H": "STH",
    "STH": "STH",
    "SPECTROMETER": "SPC",
    "SPC": "SPC",
    "XRS": "XRS",
    "BSW": "BSW",
    "DRIFT": "DRIFT",
}


class JsonLatticeLoader:
    """Load a FELsim JSON lattice file into beamline representations."""

    def __init__(self, file_path, validate_schema=True, debug=None):
        self.file_path = str(file_path)
        self.logger, self.debug = get_logger_with_fallback(__name__, debug)

        with open(self.file_path) as f:
            raw = json.load(f)

        if validate_schema:
            self._validate_schema(raw)

        self._tracked = TrackedDict(raw)
        self._beamline = self._tracked["beamline"]

        fv = self._beamline["metadata"]["format_version"]
        if fv != SUPPORTED_FORMAT_VERSION:
            raise ValueError(
                f"Unsupported format_version {fv} (expected {SUPPORTED_FORMAT_VERSION})"
            )

    # ------------------------------------------------------------------
    # Public API
    # ------------------------------------------------------------------

    def parse_beamline(self):
        """Return list of dicts compatible with BeamlineBuilder / adapters.

        Dict keys: type, length, current, angle, wedge_angle, gap_wedge,
        pole_gap, enge_fct, z_start, z_end.
        Drifts are auto-inserted between elements.
        """
        elements = self._positioned_elements()
        result = []
        prev_z_end = 0.0

        for elem in elements:
            z_start = elem["s_start_m"]
            z_end = elem["s_end_m"]

            if z_start > prev_z_end:
                result.append({"type": "DRIFT", "length": z_start - prev_z_end})

            result.append(self._element_to_dict(elem))
            prev_z_end = z_end

        self._report_unaccessed()
        return result

    def create_beamline(self):
        """Return list of beamline.py class instances (driftLattice, qpfLattice, etc.).

        Mirrors the output of ExcelElements.create_beamline().
        """
        elements = self._positioned_elements()
        result = []
        prev_z_end = 0.0

        for elem in elements:
            z_start = elem["s_start_m"]
            z_end = elem["s_end_m"]

            if z_start > prev_z_end:
                result.append(driftLattice(z_start - prev_z_end))

            obj = self._element_to_object(elem)
            if obj is not None:
                result.append(obj)
            prev_z_end = z_end

        self._report_unaccessed()
        return result

    @property
    def metadata(self):
        return self._beamline["metadata"]

    @property
    def beam_parameters(self):
        return self._beamline["beam_parameters"]

    # ------------------------------------------------------------------
    # Internal helpers
    # ------------------------------------------------------------------

    def _positioned_elements(self):
        """Return elements with valid s_start/s_end, sorted by s_start."""
        # Read top-level sections we understand; mark others as accessed.
        meta = self._beamline["metadata"]
        _ = meta["name"], meta["version"], meta["reference_energy_mev"], meta["particle_type"]
        meta.mark_accessed("description", "author", "date")

        bp = self._beamline["beam_parameters"]
        _ = bp["particle"]["type"], bp["particle"]["kinetic_energy_mev"]
        _ = bp["particle"]["mass_mev"], bp["particle"]["charge_e"]
        _ = bp["rf_frequency_hz"]

        self._beamline.mark_accessed(
            "lattice_structure", "global_settings", "simulator_specific"
        )

        raw_elements = self._beamline["elements"]
        positioned = []
        for elem in raw_elements:
            s_start = elem["s_start_m"]
            s_end = elem["s_end_m"]
            if s_start is not None and s_end is not None and s_end > s_start:
                positioned.append(elem)
            elif s_start is not None and s_end is not None and s_start == s_end:
                # Zero-length element (BPM, corrector, etc.) — include but skip drift
                positioned.append(elem)
            else:
                # Example/placeholder element with null positions — skip
                self.logger.debug(f"Skipping element {elem['name']!r} (no valid position)")
                elem.mark_all_accessed()
        positioned.sort(key=lambda e: e["s_start_m"].raw if hasattr(e["s_start_m"], 'raw') else e["s_start_m"])
        return positioned

    def _resolve_type(self, elem):
        """Resolve JSON element type + polarity to internal short name."""
        json_type = elem["type"]
        if json_type == "QUADRUPOLE":
            polarity = elem["polarity"]
            return "QPF" if polarity == "focusing" else "QPD"
        short = _TYPE_ALIASES.get(json_type, json_type)
        return short

    def _element_to_dict(self, elem):
        """Convert a tracked element to a BeamlineBuilder-compatible dict."""
        internal_type = self._resolve_type(elem)
        length = elem["length_m"]
        z_start = elem["s_start_m"]
        z_end = elem["s_end_m"]
        params = elem["parameters"]

        current = params.get("current_a", 0) or 0
        angle = 0.0
        wedge_angle = 0.0
        gap_wedge = 0.0
        pole_gap = 0.0
        enge_fct = ""

        if internal_type == "DPH":
            angle = params.get("bending_angle_deg", 0) or 0
            pole_gap = params.get("pole_gap_m", 0) or 0
            params.mark_accessed("dipole_length_m")
        elif internal_type == "DPW":
            angle = params.get("dipole_angle_deg", 0) or 0
            wedge_angle = params.get("wedge_angle_deg", 0) or 0
            gap_wedge = length
            pole_gap = params.get("pole_gap_m", 0) or 0
            params.mark_accessed("dipole_length_m")
            enge_fct = self._get_enge(elem)

        # Mark remaining standard fields as consumed
        elem.mark_accessed("name", "aperture_m", "optimization", "fringe_fields", "metadata")
        params.mark_all_accessed()

        return {
            "type": internal_type,
            "length": length,
            "current": current,
            "angle": angle,
            "wedge_angle": wedge_angle,
            "gap_wedge": gap_wedge,
            "pole_gap": pole_gap,
            "enge_fct": enge_fct,
            "z_start": z_start,
            "z_end": z_end,
        }

    def _element_to_object(self, elem):
        """Convert a tracked element to a beamline.py class instance."""
        internal_type = self._resolve_type(elem)
        length = elem["length_m"]
        z_start = elem["s_start_m"]
        z_end = elem["s_end_m"]
        params = elem["parameters"]

        elem.mark_accessed("name", "aperture_m", "optimization", "fringe_fields", "metadata")

        if internal_type == "DRIFT":
            params.mark_all_accessed()
            if length > 0:
                return driftLattice(length)
            return None

        elif internal_type == "QPF":
            current = params.get("current_a", 0) or 0
            params.mark_all_accessed()
            return qpfLattice(current=current, length=length)

        elif internal_type == "QPD":
            current = params.get("current_a", 0) or 0
            params.mark_all_accessed()
            return qpdLattice(current=current, length=length)

        elif internal_type == "DPH":
            angle = params.get("bending_angle_deg", 0) or 0
            dipole_length = params.get("dipole_length_m", length) or length
            params.mark_all_accessed()
            return dipole(length=dipole_length, angle=angle)

        elif internal_type == "DPW":
            wedge_angle = params.get("wedge_angle_deg", 0) or 0
            dipole_angle = params.get("dipole_angle_deg", 0) or 0
            dipole_length = params.get("dipole_length_m", 0) or 0
            pole_gap = params.get("pole_gap_m", 0) or 0
            params.mark_all_accessed()
            enge_fct = self._get_enge(elem)
            return dipole_wedge(
                length=length, angle=wedge_angle,
                dipole_length=dipole_length, dipole_angle=dipole_angle,
                pole_gap=pole_gap, enge_fct=enge_fct,
            )

        else:
            # Diagnostic / passive elements (BPM, OTR, correctors, undulator, etc.)
            # treated as zero-length or drift
            params.mark_all_accessed()
            if length > 0:
                return driftLattice(length)
            return None

    def _get_enge(self, elem):
        """Extract Enge coefficients from a tracked element."""
        if "fringe_fields" not in elem:
            return []
        ff = elem["fringe_fields"]
        coeffs = ff.get("enge_coefficients")
        ff.mark_all_accessed()
        if coeffs is None:
            return []
        if hasattr(coeffs, "raw"):
            return list(coeffs.raw)
        return list(coeffs) if coeffs else []

    def _validate_schema(self, raw):
        """Validate raw JSON against the lattice schema (if jsonschema is installed)."""
        try:
            import jsonschema
        except ImportError:
            self.logger.debug("jsonschema not installed; skipping schema validation")
            return

        schema_path = Path(self.file_path).resolve().parent.parent / "var" / "lattice_schema_v1.json"
        if not schema_path.exists():
            # Try relative to backend/
            schema_path = Path(__file__).resolve().parent.parent / "var" / "lattice_schema_v1.json"
        if not schema_path.exists():
            self.logger.debug(f"Schema file not found at {schema_path}; skipping validation")
            return

        with open(schema_path) as f:
            schema = json.load(f)

        try:
            jsonschema.validate(raw, schema)
            self.logger.debug("JSON schema validation passed")
        except jsonschema.ValidationError as e:
            raise ValueError(f"JSON schema validation failed: {e.message}") from e

    def _report_unaccessed(self):
        """Log any unconsumed data paths."""
        unaccessed = self._tracked.unaccessed()
        if unaccessed:
            self.logger.info(f"JSON lattice: {len(unaccessed)} unhandled field(s):")
            for path in unaccessed:
                self.logger.info(f"  {path}")
