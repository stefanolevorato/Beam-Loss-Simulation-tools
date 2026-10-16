"""Geometry model shared by the historical section editor.

The simulator snapshot shipped with the recovery is self-contained, while the
standalone editor snapshot imports this module.  This file keeps that editor
API available: sequential sections, material/index presets, independent
detector active areas, and JSON round trips.  Quantities ending in ``_m`` are
SI metres.
"""

from __future__ import annotations

import json
import math
from copy import deepcopy
from dataclasses import dataclass, field, replace
from pathlib import Path
from typing import Any, Mapping, Optional, Sequence, Tuple

import numpy as np


SECTION_SHAPES = (
    "rectangle",
    "taper",
    "segment",
    "tapered",
    "segmented",
    "chamfered",
    "polygonal",
)


def fused_silica_index(wavelength_m: float | np.ndarray) -> float | np.ndarray:
    """Return the fused-silica Sellmeier index for wavelength in metres."""

    value = np.asarray(wavelength_m, dtype=float)
    um2 = (value * 1e6) ** 2
    n2 = (
        1.0
        + 0.6961663 * um2 / (um2 - 0.0684043**2)
        + 0.4079426 * um2 / (um2 - 0.1162414**2)
        + 0.8974794 * um2 / (um2 - 9.896161**2)
    )
    result = np.sqrt(n2)
    return float(result) if result.ndim == 0 else result


@dataclass(frozen=True)
class MaterialSpec:
    """Optical material preset with an optional custom refractive index."""

    preset: str = "fused_silica"
    custom_n: Optional[float] = None

    @property
    def refractive_index(self) -> float:
        return self.effective_index(450e-9)

    def index(self, wavelength_m: float | np.ndarray) -> float | np.ndarray:
        key = str(self.preset).strip().lower().replace("-", "_")
        value = np.asarray(wavelength_m, dtype=float)
        if key in {"fused_silica", "quartz", "silica"}:
            result = fused_silica_index(value)
        elif key in {"air", "vacuum"}:
            result = np.ones_like(value)
        elif key in {"water"}:
            result = np.full_like(value, 1.333)
        elif key in {"pmma", "acrylic"}:
            result = np.full_like(value, 1.49)
        elif key == "custom":
            if self.custom_n is None or not math.isfinite(float(self.custom_n)):
                raise ValueError("custom material requires a finite custom_n")
            result = np.full_like(value, float(self.custom_n))
        else:
            # Preserve unknown historical presets as custom-index materials.
            if self.custom_n is None:
                raise ValueError(f"unknown material preset: {self.preset}")
            result = np.full_like(value, float(self.custom_n))
        return float(result) if np.asarray(result).ndim == 0 else result

    def effective_index(self, wavelength_m: float = 450e-9) -> float:
        return float(self.index(wavelength_m))

    def to_dict(self) -> dict[str, Any]:
        return {
            "preset": self.preset,
            "custom_n": self.custom_n,
            "refractive_index": self.refractive_index,
        }

    @classmethod
    def from_dict(cls, value: Mapping[str, Any] | None) -> "MaterialSpec":
        value = value if isinstance(value, Mapping) else {}
        preset = str(value.get("preset", value.get("name", "fused_silica")))
        custom = value.get("custom_n", value.get("refractive_index"))
        if preset.strip().lower() not in {"custom"} and custom is not None:
            # Older JSON files used refractive_index only for custom materials.
            custom = custom if preset.strip().lower() in {"unknown"} else None
        return cls(preset, None if custom in (None, "") else float(custom))


@dataclass(frozen=True)
class ActiveDetectorArea:
    """Independent active area on a terminal x face."""

    enabled: bool = True
    mode: str = "match_exit_surface"
    width_m: Optional[float] = None
    height_m: Optional[float] = None
    center_y_m: float = 0.0
    center_z_m: float = 0.0
    shape: str = "rectangle"
    radius_m: Optional[float] = None
    semi_axis_a_m: Optional[float] = None
    semi_axis_b_m: Optional[float] = None

    def __post_init__(self) -> None:
        if self.mode == "disabled":
            object.__setattr__(self, "enabled", False)

    def resolved(self, face_width: float, face_height: float) -> Tuple[float, float, float, float]:
        if self.mode == "match_exit_surface" or self.shape == "full_exit_surface":
            return face_width, face_height, self.center_y_m, self.center_z_m
        width = face_width if self.width_m is None else float(self.width_m)
        height = face_height if self.height_m is None else float(self.height_m)
        if self.shape == "circle":
            radius = self.radius_m if self.radius_m is not None else min(width, height) / 2.0
            width = height = 2.0 * radius
        elif self.shape == "ellipse":
            width = 2.0 * (self.semi_axis_a_m if self.semi_axis_a_m is not None else width / 2.0)
            height = 2.0 * (self.semi_axis_b_m if self.semi_axis_b_m is not None else height / 2.0)
        return width, height, self.center_y_m, self.center_z_m

    def contains(self, y: float, z: float, face_width: float, face_height: float) -> bool:
        if not self.enabled or self.mode == "disabled":
            return False
        width, height, cy, cz = self.resolved(face_width, face_height)
        if self.shape == "circle":
            radius = self.radius_m if self.radius_m is not None else min(width, height) / 2.0
            return (y - cy) ** 2 + (z - cz) ** 2 <= radius**2
        if self.shape == "ellipse":
            a = self.semi_axis_a_m if self.semi_axis_a_m is not None else width / 2.0
            b = self.semi_axis_b_m if self.semi_axis_b_m is not None else height / 2.0
            return ((y - cy) / max(a, 1e-30)) ** 2 + ((z - cz) / max(b, 1e-30)) ** 2 <= 1.0
        return abs(y - cy) <= width / 2.0 and abs(z - cz) <= height / 2.0

    def to_dict(self) -> dict[str, Any]:
        return {
            "enabled": bool(self.enabled and self.mode != "disabled"),
            "mode": self.mode,
            "shape": self.shape,
            "width_m": self.width_m,
            "height_m": self.height_m,
            "radius_m": self.radius_m,
            "semi_axis_a_m": self.semi_axis_a_m,
            "semi_axis_b_m": self.semi_axis_b_m,
            "center_y_m": self.center_y_m,
            "center_z_m": self.center_z_m,
        }

    @classmethod
    def from_dict(cls, value: Mapping[str, Any] | None) -> "ActiveDetectorArea":
        value = value if isinstance(value, Mapping) else {}
        enabled = bool(value.get("enabled", True))
        mode = str(value.get("mode", "match_exit_surface"))
        if not enabled:
            mode = "disabled"
        return cls(
            enabled=enabled,
            mode=mode,
            shape=str(value.get("shape", "rectangle")),
            width_m=None if value.get("width_m") in (None, "") else float(value["width_m"]),
            height_m=None if value.get("height_m") in (None, "") else float(value["height_m"]),
            radius_m=None if value.get("radius_m") in (None, "") else float(value["radius_m"]),
            semi_axis_a_m=None if value.get("semi_axis_a_m") in (None, "") else float(value["semi_axis_a_m"]),
            semi_axis_b_m=None if value.get("semi_axis_b_m") in (None, "") else float(value["semi_axis_b_m"]),
            center_y_m=float(value.get("center_y_m", 0.0)),
            center_z_m=float(value.get("center_z_m", 0.0)),
        )


@dataclass(frozen=True)
class GeometrySection:
    length_m: float = 0.100
    width_start_m: float = 0.020
    width_end_m: Optional[float] = None
    thickness_start_m: float = 0.010
    thickness_end_m: Optional[float] = None
    material: MaterialSpec = field(default_factory=MaterialSpec)
    name: str = "section_1"
    shape: str = "rectangle"

    def __post_init__(self) -> None:
        if self.shape not in SECTION_SHAPES:
            raise ValueError(f"unsupported section shape: {self.shape}")
        if self.length_m <= 0 or self.width_start_m <= 0 or self.thickness_start_m <= 0:
            raise ValueError("section length and dimensions must be positive")

    def dimensions_at(self, fraction: float) -> Tuple[float, float]:
        fraction = min(1.0, max(0.0, float(fraction)))
        width_end = self.width_start_m if self.width_end_m is None else self.width_end_m
        thickness_end = self.thickness_start_m if self.thickness_end_m is None else self.thickness_end_m
        return (
            (1.0 - fraction) * self.width_start_m + fraction * width_end,
            (1.0 - fraction) * self.thickness_start_m + fraction * thickness_end,
        )

    def to_dict(self) -> dict[str, Any]:
        return {
            "name": self.name,
            "shape": self.shape,
            "length_m": self.length_m,
            "width_start_m": self.width_start_m,
            "width_end_m": self.width_end_m,
            "thickness_start_m": self.thickness_start_m,
            "thickness_end_m": self.thickness_end_m,
            "material": self.material.to_dict(),
        }

    @classmethod
    def from_dict(cls, value: Mapping[str, Any] | None) -> "GeometrySection":
        value = value if isinstance(value, Mapping) else {}
        width = float(value.get("width_start_m", value.get("width_y_m", 0.020)))
        thickness = float(value.get("thickness_start_m", value.get("thickness_z_m", 0.010)))
        return cls(
            length_m=float(value.get("length_m", 0.100)),
            width_start_m=width,
            width_end_m=None if value.get("width_end_m") in (None, "") else float(value["width_end_m"]),
            thickness_start_m=thickness,
            thickness_end_m=None if value.get("thickness_end_m") in (None, "") else float(value["thickness_end_m"]),
            material=MaterialSpec.from_dict(value.get("material")),
            name=str(value.get("name", "section_1")),
            shape=str(value.get("shape", value.get("type", "rectangle"))).lower(),
        )


@dataclass(frozen=True)
class GeometryModel:
    length_x_m: float = 0.100
    width_y_m: float = 0.020
    thickness_z_m: float = 0.010
    sections: Tuple[GeometrySection, ...] = field(default_factory=tuple)
    detector_plus_x: ActiveDetectorArea = field(default_factory=ActiveDetectorArea)
    detector_minus_x: ActiveDetectorArea = field(default_factory=ActiveDetectorArea)

    def __post_init__(self) -> None:
        sections = tuple(self.sections) if self.sections else (GeometrySection(length_m=self.length_x_m, width_start_m=self.width_y_m, thickness_start_m=self.thickness_z_m),)
        object.__setattr__(self, "sections", sections)
        if self.length_x_m <= 0 or self.width_y_m <= 0 or self.thickness_z_m <= 0:
            raise ValueError("geometry dimensions must be positive")

    @property
    def half_x(self) -> float:
        return self.length_x_m / 2.0

    @property
    def half_y(self) -> float:
        return self.width_y_m / 2.0

    @property
    def half_z(self) -> float:
        return self.thickness_z_m / 2.0

    @property
    def x_edges(self) -> Tuple[float, ...]:
        total = sum(section.length_m for section in self.sections)
        scale = self.length_x_m / total if total else 1.0
        edges = [-self.half_x]
        for section in self.sections:
            edges.append(edges[-1] + section.length_m * scale)
        edges[-1] = self.half_x
        return tuple(edges)

    def section_index_at_x(self, x: float) -> int:
        return min(len(self.sections) - 1, max(0, int(np.searchsorted(self.x_edges[1:-1], float(x), side="right"))))

    def dimensions_at_x(self, x: float) -> Tuple[float, float]:
        index = self.section_index_at_x(x)
        x0, x1 = self.x_edges[index:index + 2]
        return self.sections[index].dimensions_at((float(x) - x0) / max(x1 - x0, 1e-30))

    def material_at_x(self, x: float) -> MaterialSpec:
        return self.sections[self.section_index_at_x(x)].material

    def face_dimensions(self, face: str) -> Tuple[float, float]:
        return self.dimensions_at_x(self.half_x if face == "plus_x" else -self.half_x)

    def active_area(self, face: str) -> ActiveDetectorArea:
        return self.detector_plus_x if face == "plus_x" else self.detector_minus_x

    def contains(self, point: Sequence[float], tolerance: float = 1e-12) -> bool:
        x, y, z = point
        width, thickness = self.dimensions_at_x(float(x))
        return -self.half_x - tolerance <= x <= self.half_x + tolerance and abs(y) <= width / 2 + tolerance and abs(z) <= thickness / 2 + tolerance

    def to_dict(self) -> dict[str, Any]:
        return {
            "length_x_m": self.length_x_m,
            "width_y_m": self.width_y_m,
            "thickness_z_m": self.thickness_z_m,
            "sections": [section.to_dict() for section in self.sections],
            "detector_plus_x": self.detector_plus_x.to_dict(),
            "detector_minus_x": self.detector_minus_x.to_dict(),
            "plus_x_active_area": self.detector_plus_x.to_dict(),
            "minus_x_active_area": self.detector_minus_x.to_dict(),
        }

    @classmethod
    def from_dict(cls, value: Mapping[str, Any] | None) -> "GeometryModel":
        value = value if isinstance(value, Mapping) else {}
        sections = tuple(GeometrySection.from_dict(item) for item in value.get("sections", ()))
        plus = value.get("detector_plus_x", value.get("plus_x_active_area", {}))
        minus = value.get("detector_minus_x", value.get("minus_x_active_area", {}))
        return cls(
            length_x_m=float(value.get("length_x_m", 0.100)),
            width_y_m=float(value.get("width_y_m", 0.020)),
            thickness_z_m=float(value.get("thickness_z_m", 0.010)),
            sections=sections,
            detector_plus_x=ActiveDetectorArea.from_dict(plus),
            detector_minus_x=ActiveDetectorArea.from_dict(minus),
        )

    @classmethod
    def from_json(cls, path: str | Path) -> "GeometryModel":
        return cls.from_dict(json.loads(Path(path).read_text(encoding="utf-8")))

    def to_json(self, path: str | Path) -> Path:
        destination = Path(path)
        destination.write_text(json.dumps(self.to_dict(), indent=2) + "\n", encoding="utf-8")
        return destination


BarGeometry = GeometryModel


# ---------------------------------------------------------------------------
# Shared assembly model used by the v13 GUI
# ---------------------------------------------------------------------------

ASSEMBLY_SECTION_SHAPES = (
    "box",
    "tapered_rectangular_prism",
    "cylinder",
    "elliptical_cylinder",
    "truncated_cone",
)
ASSEMBLY_AXES = ("X", "Y", "Z")

# These are deliberately shared by the model, GUI controls, render scene, and
# optical adapter.  Keeping the user-facing wording here prevents the old
# ``+x/-x sensitive area`` terminology from leaking back into the assembly
# editors.
TERMINAL_SURFACE_LABELS = {
    "open": "Open",
    "mirror": "Mirror",
    "detector": "Detector/Photosensor coupling",
}
TERMINAL_SURFACE_CHOICES = tuple(TERMINAL_SURFACE_LABELS.values())

# Versioned, detector-specific envelopes are deliberately kept beside the
# canonical model.  Both GUI tabs use these helpers, so a file cannot be
# accidentally interpreted as the other detector assembly.
ASSEMBLY_GEOMETRY_FORMAT = "quartz-assembly-geometry"
ASSEMBLY_GEOMETRY_VERSION = 27
ASSEMBLY_GEOMETRY_TYPES = ("cherenkov", "scintillator")

# v32 is the complete local detector-assembly envelope.  The older v27
# envelope remains intentionally available for compatibility with existing
# typed geometry files; v32 adds the photosensor/readout state and makes the
# complete-assembly intent explicit instead of treating the sensor as a
# separate profile.
COMPLETE_DETECTOR_ASSEMBLY_GEOMETRY_FORMAT = "complete-detector-assembly-geometry"
COMPLETE_DETECTOR_ASSEMBLY_GEOMETRY_VERSION = 32
COMPLETE_DETECTOR_ASSEMBLY_GEOMETRY_MARKER = "complete-detector-assembly-geometry"

# v62 adds explicit shape-specific readout dimensions while retaining the
# v61 ``dimensions_mm.width/height`` envelope.  The retained dimensions are
# always full transverse extents; circle/ellipse controls are radii/semiaxes.
SENSOR_READOUT_GEOMETRY_VERSION = 2

# Display-only Mylar styling.  These values never enter the transport model;
# the explicit epsilon is applied only to copied overlay polygons to prevent
# z-fighting with the canonical physical faces.
MYLAR_DISPLAY_EPSILON_MM = 0.05
MYLAR_DISPLAY_COLOR = "#f2c14e"
MYLAR_DISPLAY_EDGE = "#b26a00"
MYLAR_DISPLAY_HATCH = "gold-diagonal-stripes"


def resolve_sensor_readout_dimensions_mm(
    shape: str, *, width_mm: float | None = None, height_mm: float | None = None,
    radius_mm: float | None = None, semi_axis_a_mm: float | None = None,
    semi_axis_b_mm: float | None = None,
) -> dict[str, float | str]:
    """Normalize one readout to full extents plus explicit shape dimensions.

    v61 width/height values were full transverse dimensions.  Consequently a
    legacy circular 2 mm x 2 mm readout migrates to radius 1 mm, not 2 mm.
    Ellipse width/height similarly migrate to semiaxes width/2 and height/2.
    """

    key = str(shape).strip().lower()
    aliases = {"rectangle": "rectangular", "circle": "circular", "ellipse": "elliptical"}
    key = aliases.get(key, key)
    if key not in {"rectangular", "circular", "elliptical"}:
        raise ValueError("sensor/readout shape must be rectangular, circular, or elliptical")

    def positive(value: float | None, label: str) -> float:
        if value is None:
            raise ValueError(f"{label} is required for {key} readout")
        result = float(value)
        if not math.isfinite(result) or result <= 0.0:
            raise ValueError(f"{label} must be finite and positive")
        return result

    if key == "circular":
        radius = positive(radius_mm if radius_mm is not None else (
            min(float(width_mm), float(height_mm)) / 2.0
            if width_mm is not None and height_mm is not None else None
        ), "readout radius")
        return {"shape": key, "width_mm": 2.0 * radius, "height_mm": 2.0 * radius,
                "radius_mm": radius, "semi_axis_a_mm": radius, "semi_axis_b_mm": radius}
    if key == "elliptical":
        a = positive(semi_axis_a_mm if semi_axis_a_mm is not None else (
            float(width_mm) / 2.0 if width_mm is not None else None
        ), "readout semi-axis A")
        b = positive(semi_axis_b_mm if semi_axis_b_mm is not None else (
            float(height_mm) / 2.0 if height_mm is not None else None
        ), "readout semi-axis B")
        return {"shape": key, "width_mm": 2.0 * a, "height_mm": 2.0 * b,
                "radius_mm": min(a, b), "semi_axis_a_mm": a, "semi_axis_b_mm": b}
    return {"shape": key, "width_mm": positive(width_mm, "readout width"),
            "height_mm": positive(height_mm, "readout height"),
            "radius_mm": min(float(width_mm), float(height_mm)) / 2.0,
            "semi_axis_a_mm": float(width_mm) / 2.0, "semi_axis_b_mm": float(height_mm) / 2.0}


def _is_legacy_guide_terminal(value: Any) -> bool:
    text = str(value or "").strip().lower().replace("_", " ")
    return text in {"guide", "light guide", "light-guide coupling", "optical/light guide coupling", "optical/light-guide coupling"}


def normalize_terminal_surface(value: Any, default: str = "open") -> str:
    """Return the canonical terminal state used by every assembly consumer."""

    text = str(value if value is not None else default).strip().lower()
    aliases = {
        "open": "open", "none": "open", "escape": "open",
        "mirror": "mirror", "reflector": "mirror",
        "detector": "detector", "photosensor": "detector",
        "detector_photosensor": "detector", "detector/photosensor": "detector",
        "detector/photosensor coupling": "detector",
    }
    # v24 exposed the guide as a fourth terminal mode.  It actually meant
    # “read out through the configured guide”, now derived from geometry.
    if _is_legacy_guide_terminal(value):
        return "detector"
    result = aliases.get(text, text)
    if result not in TERMINAL_SURFACE_LABELS:
        raise ValueError(f"unsupported terminal surface state: {value}")
    return result


def terminal_surface_label(value: Any) -> str:
    return TERMINAL_SURFACE_LABELS[normalize_terminal_surface(value)]


def terminal_surfaces_from_payload(
    payload: Mapping[str, Any], *, default_plus: str = "detector", default_minus: str = "open",
    require_explicit_endpoints: bool = False,
) -> tuple[str, str]:
    """Normalize terminal states from current and legacy geometry envelopes.

    The scalar ``terminal_plus_surface``/``terminal_minus_surface`` fields
    are canonical.  ``terminal_surfaces`` is their structured vXX form and
    ``opposite_end_mode`` is the legacy Run-control alias for -X.  When both
    a scalar and structured field are present they must agree semantically;
    stale legacy values are ignored whenever either modern representation
    supplies -X.
    """

    if not isinstance(payload, Mapping):
        raise ValueError("assembly geometry JSON must be an object")

    nested: Mapping[str, Any] = {}
    if "terminal_surfaces" in payload:
        value = payload["terminal_surfaces"]
        if not isinstance(value, Mapping):
            raise ValueError("terminal_surfaces must be an object")
        nested = value
        allowed_terminal_keys = {"+X", "-X", "plus_x", "minus_x", "plus-x", "minus-x", "plus", "minus"}
        unknown_terminal_keys = set(nested) - allowed_terminal_keys
        if unknown_terminal_keys:
            raise ValueError(f"terminal_surfaces contains unsupported keys: {', '.join(sorted(map(str, unknown_terminal_keys)))}")
        if not nested:
            raise ValueError("terminal_surfaces must specify at least one endpoint")

    def nested_face(keys: tuple[str, ...], label: str) -> tuple[bool, str | None]:
        found = [(key, nested[key]) for key in keys if key in nested]
        if not found:
            return False, None
        normalized = []
        for key, value in found:
            if value is None or not str(value).strip():
                raise ValueError(f"terminal_surfaces {key} must specify a terminal state")
            try:
                normalized.append((key, normalize_terminal_surface(value)))
            except (TypeError, ValueError) as exc:
                raise ValueError(f"terminal_surfaces {key} is invalid: {exc}") from exc
        if len({state for _key, state in normalized}) != 1:
            raise ValueError(f"terminal_surfaces contains conflicting {label} values")
        return True, normalized[0][1]

    plus_nested_present, plus_nested = nested_face(("+X", "plus_x", "plus-x", "plus"), "+X")
    minus_nested_present, minus_nested = nested_face(("-X", "minus_x", "minus-x", "minus"), "-X")

    def scalar(key: str) -> tuple[bool, str | None]:
        if key not in payload:
            return False, None
        if payload[key] is None or not str(payload[key]).strip():
            raise ValueError(f"{key} must specify Open, Mirror, or Detector/Photosensor coupling")
        try:
            return True, normalize_terminal_surface(payload[key])
        except (TypeError, ValueError) as exc:
            raise ValueError(f"{key} is invalid: {exc}") from exc

    plus_scalar_present, plus_scalar = scalar("terminal_plus_surface")
    minus_scalar_present, minus_scalar = scalar("terminal_minus_surface")
    if plus_scalar_present and plus_nested_present and plus_scalar != plus_nested:
        raise ValueError("terminal_surfaces +X disagrees with terminal_plus_surface")
    if minus_scalar_present and minus_nested_present and minus_scalar != minus_nested:
        raise ValueError("terminal_surfaces -X disagrees with terminal_minus_surface")

    if require_explicit_endpoints:
        if not (minus_scalar_present or minus_nested_present or "opposite_end_mode" in payload):
            raise ValueError("geometry JSON is missing a -X terminal surface")
        if not (plus_scalar_present or plus_nested_present) and "opposite_end_mode" not in payload:
            raise ValueError("geometry JSON is missing a +X terminal surface")

    if plus_scalar_present:
        plus = plus_scalar
    elif plus_nested_present:
        plus = plus_nested
    else:
        plus = normalize_terminal_surface(default_plus)

    if minus_scalar_present:
        minus = minus_scalar
    elif minus_nested_present:
        minus = minus_nested
    elif "opposite_end_mode" in payload:
        if payload["opposite_end_mode"] is None or not str(payload["opposite_end_mode"]).strip():
            raise ValueError("opposite_end_mode must specify the legacy -X terminal state")
        try:
            minus = normalize_terminal_surface(payload["opposite_end_mode"])
        except (TypeError, ValueError) as exc:
            raise ValueError(f"opposite_end_mode is invalid: {exc}") from exc
    else:
        minus = normalize_terminal_surface(default_minus)
    return plus, minus


@dataclass
class AssemblySection:
    """One ordered, editable section shared by both detector assemblies.

    The GUI stores this model in millimetres.  ``to_simulation_geometry`` is
    intentionally lazy-imported so this geometry module remains useful to the
    standalone editor without introducing a simulator import cycle.
    """

    name: str = "Section 1"
    shape: str = "box"
    material: str = "fused silica"
    length_x_mm: float = 100.0
    width_y_mm: float = 20.0
    thickness_z_mm: float = 10.0
    width_y_start_mm: float = 20.0
    width_y_end_mm: float = 20.0
    thickness_z_start_mm: float = 10.0
    thickness_z_end_mm: float = 10.0
    axis: str = "Z"
    radius_mm: float = 10.0
    length_mm: float = 100.0
    semi_axis_a_mm: float = 10.0
    semi_axis_b_mm: float = 5.0
    radius_start_mm: float = 10.0
    radius_end_mm: float = 5.0
    center_x_mm: float = 0.0
    center_y_mm: float = 0.0
    center_z_mm: float = 0.0

    def validate(self) -> None:
        if not str(self.name).strip():
            raise ValueError("geometry section name cannot be empty")
        self.shape = str(self.shape).strip().lower()
        self.axis = str(self.axis).strip().upper()
        if self.shape not in ASSEMBLY_SECTION_SHAPES:
            raise ValueError(f"unsupported geometry shape: {self.shape}")
        if self.axis not in ASSEMBLY_AXES:
            raise ValueError("section axis must be X, Y, or Z")
        for key in ("center_x_mm", "center_y_mm", "center_z_mm"):
            if not math.isfinite(float(getattr(self, key))):
                raise ValueError(f"{key} must be finite")
        positive = {
            "box": ("length_x_mm", "width_y_mm", "thickness_z_mm"),
            "tapered_rectangular_prism": ("length_x_mm", "width_y_start_mm", "width_y_end_mm", "thickness_z_start_mm", "thickness_z_end_mm"),
            "cylinder": ("radius_mm", "length_mm"),
            "elliptical_cylinder": ("semi_axis_a_mm", "semi_axis_b_mm", "length_mm"),
            "truncated_cone": ("radius_start_mm", "radius_end_mm", "length_mm"),
        }[self.shape]
        for key in positive:
            value = float(getattr(self, key))
            if not math.isfinite(value) or value <= 0.0:
                raise ValueError(f"{key} must be positive")

    def to_dict(self) -> dict[str, Any]:
        self.validate()
        return {
            key: getattr(self, key)
            for key in (
                "name", "shape", "material", "length_x_mm", "width_y_mm", "thickness_z_mm",
                "width_y_start_mm", "width_y_end_mm", "thickness_z_start_mm", "thickness_z_end_mm",
                "axis", "radius_mm", "length_mm", "semi_axis_a_mm", "semi_axis_b_mm",
                "radius_start_mm", "radius_end_mm", "center_x_mm", "center_y_mm", "center_z_mm",
            )
        }

    @classmethod
    def from_dict(cls, payload: Mapping[str, Any] | None, index: int = 0) -> "AssemblySection":
        raw = payload if isinstance(payload, Mapping) else {}
        def number(key: str, default: float) -> float:
            try:
                value = float(raw.get(key, default))
                return value if math.isfinite(value) else default
            except (TypeError, ValueError):
                return default
        shape = str(raw.get("shape", "box")).strip().lower()
        aliases = {"rectangular": "box", "rectangular_prism": "box", "tapered box": "tapered_rectangular_prism"}
        shape = aliases.get(shape, shape)
        section = cls(
            name=str(raw.get("name", f"Section {index + 1}")), shape=shape,
            material=str(raw.get("material", "fused silica")),
            length_x_mm=number("length_x_mm", 100.0), width_y_mm=number("width_y_mm", 20.0),
            thickness_z_mm=number("thickness_z_mm", 10.0),
            width_y_start_mm=number("width_y_start_mm", number("width_y_mm", 20.0)),
            width_y_end_mm=number("width_y_end_mm", number("width_y_mm", 20.0)),
            thickness_z_start_mm=number("thickness_z_start_mm", number("thickness_z_mm", 10.0)),
            thickness_z_end_mm=number("thickness_z_end_mm", number("thickness_z_mm", 10.0)),
            axis=str(raw.get("axis", "Z")).upper(), radius_mm=number("radius_mm", 10.0),
            length_mm=number("length_mm", 100.0), semi_axis_a_mm=number("semi_axis_a_mm", 10.0),
            semi_axis_b_mm=number("semi_axis_b_mm", 5.0), radius_start_mm=number("radius_start_mm", 10.0),
            radius_end_mm=number("radius_end_mm", 5.0), center_x_mm=number("center_x_mm", 0.0),
            center_y_mm=number("center_y_mm", 0.0), center_z_mm=number("center_z_mm", 0.0),
        )
        section.validate()
        return section

    def to_simulation_geometry(self):
        import quartz_cherenkov_sim as simulator
        self.validate()
        material = {
            "fused silica": "fused_silica", "quartz": "fused_silica", "air": "air",
            "scintillator": "scintillator", "plastic scintillator": "scintillator",
            "light guide": "light_guide", "light_guide": "light_guide", "custom": "custom",
        }.get(str(self.material).strip().lower())
        if material is None:
            raise ValueError(f"unsupported section material: {self.material}")
        return simulator.SolidGeometry(
            shape=self.shape, material=material,
            center_x_m=self.center_x_mm * 1e-3, center_y_m=self.center_y_mm * 1e-3, center_z_m=self.center_z_mm * 1e-3,
            length_x_m=self.length_x_mm * 1e-3, width_y_m=self.width_y_mm * 1e-3, thickness_z_m=self.thickness_z_mm * 1e-3,
            width_y_start_m=self.width_y_start_mm * 1e-3, width_y_end_m=self.width_y_end_mm * 1e-3,
            thickness_z_start_m=self.thickness_z_start_mm * 1e-3, thickness_z_end_m=self.thickness_z_end_mm * 1e-3,
            axis=self.axis, radius_m=self.radius_mm * 1e-3, length_m=self.length_mm * 1e-3,
            semi_axis_a_m=self.semi_axis_a_mm * 1e-3, semi_axis_b_m=self.semi_axis_b_mm * 1e-3,
            radius_start_m=self.radius_start_mm * 1e-3, radius_end_m=self.radius_end_mm * 1e-3,
        )


def replace_simulation_terminal_metadata(solid: Any, terminal_metadata: Sequence[Mapping[str, Any]],
                                         sensor_components: Sequence[Mapping[str, Any]] = ()) -> Any:
    """Attach canonical terminal metadata without changing local dimensions."""

    # The simulator dataclasses deliberately keep this field passive: optical
    # transport reads the same state through OpticalConfig, while renderers
    # and diagnostics can inspect the exact model that produced the solid.
    return replace(
        solid,
        terminal_interfaces=tuple(dict(item) for item in terminal_metadata),
        sensor_components=tuple(dict(item) for item in sensor_components),
    )


@dataclass
class AssemblyGeometryModel:
    """Complete local assembly description shared by local and global views."""

    name: str
    sections: list[AssemblySection] = field(default_factory=list)
    guide_sections: list[AssemblySection] = field(default_factory=list)
    guide_enabled: bool = False
    terminal_plus_surface: str = "detector_photosensor"
    terminal_minus_surface: str = "open"
    mylar_enabled: bool = True
    mylar_reflectivity: float = 0.95
    impact_x_mm: float = 0.0
    impact_y_mm: float = 0.0

    @staticmethod
    def _section_is_valid(section: AssemblySection) -> bool:
        """Return whether a configured optional section can be built.

        Optional guide rows are deliberately ignored when malformed.  This
        lets an in-progress/empty editor row behave as an absent guide rather
        than leaking a placeholder into every downstream consumer.
        """

        try:
            section.validate()
            material = str(section.material).strip().lower().replace("_", " ")
            return material in {"light guide", "fused silica", "quartz", "scintillator", "air", "custom"}
        except (TypeError, ValueError):
            return False

    @property
    def valid_guide_sections(self) -> tuple[AssemblySection, ...]:
        """Configured guide sections that are safe to append to active geometry."""

        if not self.guide_enabled:
            return ()
        return tuple(section for section in self.guide_sections if self._section_is_valid(section))

    @property
    def has_light_guide(self) -> bool:
        """Whether the explicit enable-and-configured rule is satisfied."""

        return bool(self.guide_enabled and self.valid_guide_sections)

    def guide_attached_at_terminal(self, face: str) -> bool:
        """Whether a valid enabled guide actually meets the selected X face."""

        if not self.has_light_guide or not self.sections:
            return False
        key = str(face).strip().lower().replace("-", "_")
        plus = key in {"plus", "plus_x", "+x", "x+"}
        active_start = _assembly_section_bounds_mm(self.sections[0])[0]
        active_end = _assembly_section_bounds_mm(self.sections[-1])[1]
        guide_start = _assembly_section_bounds_mm(self.valid_guide_sections[0])[0]
        guide_end = _assembly_section_bounds_mm(self.valid_guide_sections[-1])[1]
        return abs((guide_start if plus else guide_end) - (active_end if plus else active_start)) < 1e-9

    def guide_attachment_faces(self) -> tuple[str, ...]:
        """Return terminal faces physically joined to the active geometry."""

        return tuple(face for face in ("minus_x", "plus_x") if self.guide_attached_at_terminal(face))

    def active_sections(self) -> tuple[AssemblySection, ...]:
        return tuple(self.sections) + self.valid_guide_sections

    def terminal_state(self, face: str) -> str:
        """Return the canonical state for a real YZ terminal face."""

        key = str(face).strip().lower().replace("-", "_")
        if key in {"plus", "plus_x", "+x", "x+"}:
            requested = normalize_terminal_surface(self.terminal_plus_surface, "detector")
        elif key in {"minus", "minus_x", "-x", "x-"}:
            requested = normalize_terminal_surface(self.terminal_minus_surface, "open")
        else:
            raise ValueError("face must identify plus_x or minus_x")
        # Keep an explicitly requested guide state even when its guide is
        # missing.  That distinction is required by the scene contract:
        # invalid guide coupling must not silently become a photosensor.
        return requested

    def terminal_coupling_status(self, face: str) -> dict[str, Any]:
        """Describe the derived readout meaning of one terminal surface."""

        state = self.terminal_state(face)
        if state == "detector":
            guide_present = self.has_light_guide
            guide_attached = self.guide_attached_at_terminal(face)
            valid = not guide_present or guide_attached
            return {
                "state": state, "valid": valid, "guide_present": guide_present,
                "guide_attached": guide_attached,
                "label": TERMINAL_SURFACE_LABELS[state] if valid else "Detector/Photosensor coupling (guide output is at the other terminal)",
                "sensor_at_active_terminal": not guide_present and valid,
                "sensor_at_guide_output": guide_present and guide_attached and valid,
                "derived_mode": "light_guide" if guide_present and guide_attached and valid else "direct",
            }
        return {
            "state": state, "valid": True, "guide_present": self.has_light_guide,
            "guide_attached": self.guide_attached_at_terminal(face) if self.has_light_guide else False,
            "label": TERMINAL_SURFACE_LABELS[state],
            "sensor_at_active_terminal": False,
            "sensor_at_guide_output": False,
            "derived_mode": "none",
        }

    def readout_status(self) -> dict[str, Any]:
        """Return the single canonical readout decision for this assembly."""

        selected = tuple(face for face in ("minus_x", "plus_x") if self.terminal_state(face) == "detector")
        if len(selected) == 0:
            return {
                "valid": False, "code": "guide_present_no_readout" if self.has_light_guide else "no_readout",
                "message": "Guide present but no photosensor readout terminal selected." if self.has_light_guide else "No Detector/Photosensor coupling terminal selected.",
                "face": None, "mode": "none", "optical_path": (),
            }
        if len(selected) > 1:
            return {
                "valid": False, "code": "multiple_readout_terminals",
                "message": "Multiple photosensor readout terminals selected; choose one Detector/Photosensor coupling terminal.",
                "face": None, "mode": "none", "optical_path": (),
            }
        face = selected[0]
        coupling = self.terminal_coupling_status(face)
        if not coupling["valid"]:
            return {
                "valid": False, "code": "readout_not_at_guide_output",
                "message": "A valid guide is present, but the selected photosensor terminal is not its output.",
                "face": face, "mode": "none", "optical_path": (),
            }
        mode = str(coupling["derived_mode"])
        active_label = "Active scintillator" if str(self.name).strip().lower().startswith("scint") else "Active quartz radiator"
        path = (active_label, "Complete light guide", "Photosensor") if mode == "light_guide" else (active_label, "Photosensor")
        return {
            "valid": True, "code": "ok", "message": "Valid photosensor readout configuration.",
            "face": face, "mode": mode, "optical_path": path,
        }

    def terminal_dimensions_mm(self, face: str) -> dict[str, float | str]:
        """Return the actual terminal Y/Z dimensions, including taper ends."""

        sections = self.active_sections()
        if not sections:
            raise ValueError(f"{self.name} requires at least one active geometry section")
        plus = str(face).strip().lower().replace("-", "_") in {"plus", "plus_x", "+x", "x+"}
        section = sections[-1] if plus else sections[0]
        end = "end" if plus else "start"
        if section.shape == "box":
            width, height = section.width_y_mm, section.thickness_z_mm
        elif section.shape == "tapered_rectangular_prism":
            width = getattr(section, f"width_y_{end}_mm")
            height = getattr(section, f"thickness_z_{end}_mm")
        elif section.shape == "cylinder":
            width = height = 2.0 * section.radius_mm
        elif section.shape == "elliptical_cylinder":
            width, height = 2.0 * section.semi_axis_a_mm, 2.0 * section.semi_axis_b_mm
        elif section.shape == "truncated_cone":
            radius = section.radius_end_mm if plus else section.radius_start_mm
            width = height = 2.0 * radius
        else:  # validate() provides the detailed error for callers
            raise ValueError(f"unsupported terminal shape: {section.shape}")
        return {
            "face": "plus_x" if plus else "minus_x",
            "width_y_mm": float(width), "thickness_z_mm": float(height),
            "center_y_mm": float(section.center_y_mm), "center_z_mm": float(section.center_z_mm),
            "section": section.name,
        }

    def readout_terminal_dimensions_mm(self, face: str) -> dict[str, float | str]:
        """Return the physical face where the selected readout window couples.

        A direct photosensor couples to the active terminal.  With an enabled
        guide, the photosensor couples to the guide's outside output instead.
        Keeping this distinction in the canonical model prevents run
        validation from comparing a readout against an unrelated legacy
        portal size.
        """

        key = str(face).strip().lower().replace("-", "_")
        plus = key in {"plus", "plus_x", "+x", "x+"}
        if key not in {"plus", "plus_x", "+x", "x+", "minus", "minus_x", "-x", "x-"}:
            raise ValueError("face must identify plus_x or minus_x")
        if self.has_light_guide and self.guide_attached_at_terminal(face):
            sections = self.valid_guide_sections
            section = sections[-1] if plus else sections[0]
            end = "end" if plus else "start"
            if section.shape == "box":
                width, height = section.width_y_mm, section.thickness_z_mm
            elif section.shape == "tapered_rectangular_prism":
                width = getattr(section, f"width_y_{end}_mm")
                height = getattr(section, f"thickness_z_{end}_mm")
            elif section.shape == "cylinder":
                width = height = 2.0 * section.radius_mm
            elif section.shape == "elliptical_cylinder":
                width, height = 2.0 * section.semi_axis_a_mm, 2.0 * section.semi_axis_b_mm
            else:
                radius = section.radius_end_mm if plus else section.radius_start_mm
                width = height = 2.0 * radius
            return {
                "face": "plus_x" if plus else "minus_x", "width_y_mm": float(width),
                "thickness_z_mm": float(height), "center_y_mm": float(section.center_y_mm),
                "center_z_mm": float(section.center_z_mm), "section": section.name,
                "target": "guide output" if plus else "guide input/output terminal",
            }
        return {**self.terminal_dimensions_mm(face), "target": "active terminal face"}

    def validate_complete_assembly(
        self, *, sensor_shape: str = "rectangular", sensor_width_mm: float = 20.0,
        sensor_height_mm: float = 10.0, sensor_depth_mm: float = 2.0,
        sensor_radius_mm: float | None = None, sensor_semi_axis_a_mm: float | None = None,
        sensor_semi_axis_b_mm: float | None = None,
    ) -> dict[str, Any]:
        """Validate one complete assembly, including its real readout window.

        This is the run-time authority for the current editor model.  It
        validates the active/guide transport geometry and terminal routing, and
        checks the detector-specific photosensor window against the actual
        active terminal or guide-output dimensions.
        """

        self.validate()
        transport = self.to_simulation_geometry(
            sensor_shape=sensor_shape, sensor_width_mm=sensor_width_mm,
            sensor_height_mm=sensor_height_mm, sensor_depth_mm=sensor_depth_mm,
            sensor_radius_mm=sensor_radius_mm, sensor_semi_axis_a_mm=sensor_semi_axis_a_mm,
            sensor_semi_axis_b_mm=sensor_semi_axis_b_mm,
        )
        transport.validate()
        readout = self.readout_status()
        if not readout["valid"]:
            raise ValueError(f"{self.name} assembly has no valid photosensor readout: {readout['message']}")
        face = str(readout["face"])
        terminal = self.readout_terminal_dimensions_mm(face)
        shape = str(sensor_shape).strip().lower()
        if shape in {"circle", "ellipse"}:
            shape = "circular" if shape == "circle" else "elliptical"
        if shape not in {"rectangular", "circular", "elliptical"}:
            raise ValueError(f"{self.name} photosensor shape is unsupported: {sensor_shape}")
        dimensions = resolve_sensor_readout_dimensions_mm(
            sensor_shape, width_mm=sensor_width_mm, height_mm=sensor_height_mm,
            radius_mm=sensor_radius_mm, semi_axis_a_mm=sensor_semi_axis_a_mm,
            semi_axis_b_mm=sensor_semi_axis_b_mm,
        )
        width, height, depth = float(dimensions["width_mm"]), float(dimensions["height_mm"]), float(sensor_depth_mm)
        if not math.isfinite(depth) or depth <= 0.0:
            raise ValueError(f"{self.name} photosensor depth must be finite and positive")
        actual_width = float(terminal["width_y_mm"])
        actual_height = float(terminal["thickness_z_mm"])
        if shape == "circular":
            radius = float(dimensions["radius_mm"])
            boundary = [(radius * math.cos(2.0 * math.pi * i / 64), radius * math.sin(2.0 * math.pi * i / 64)) for i in range(64)]
        elif shape == "elliptical":
            a, b = float(dimensions["semi_axis_a_mm"]), float(dimensions["semi_axis_b_mm"])
            boundary = [(a * math.cos(2.0 * math.pi * i / 64), b * math.sin(2.0 * math.pi * i / 64)) for i in range(64)]
        else:
            boundary = [(-width / 2.0, -height / 2.0), (width / 2.0, -height / 2.0), (width / 2.0, height / 2.0), (-width / 2.0, height / 2.0)]
        terminal_profile = str(self._terminal_profile(face))
        if terminal_profile == "ellipse":
            outside = [((y / (actual_width / 2.0)) ** 2 + (z / (actual_height / 2.0)) ** 2) > 1.0 + 1e-9 for y, z in boundary]
        else:
            outside = [abs(y) > actual_width / 2.0 + 1e-9 or abs(z) > actual_height / 2.0 + 1e-9 for y, z in boundary]
        if any(outside):
            raise ValueError(
                f"{self.name} {face} photosensor window ({width:g} × {height:g} mm) does not fit "
                f"the actual {face} readout terminal of section '{terminal['section']}' "
                f"({actual_width:g} × {actual_height:g} mm)"
            )
        return {"assembly": self.name, "readout": readout, "terminal": terminal,
                "sensor": {"shape": shape, "width_mm": width, "height_mm": height, "depth_mm": depth,
                           "radius_mm": float(dimensions["radius_mm"]),
                           "semi_axis_a_mm": float(dimensions["semi_axis_a_mm"]),
                           "semi_axis_b_mm": float(dimensions["semi_axis_b_mm"])} }

    def _terminal_profile(self, face: str) -> str:
        """Return the actual terminal cross-section kind for readout checks."""

        key = str(face).strip().lower().replace("-", "_")
        plus = key in {"plus", "plus_x", "+x", "x+"}
        sections = self.valid_guide_sections if self.has_light_guide and self.guide_attached_at_terminal(face) else self.sections
        section = sections[-1] if plus else sections[0]
        return "rectangle" if section.shape in {"box", "tapered_rectangular_prism"} else "ellipse"

    def validate(self) -> None:
        if not self.sections:
            raise ValueError(f"{self.name} requires at least one active geometry section")
        for section in self.sections:
            section.validate()
        for section in self.valid_guide_sections:
            section.validate()
        if not 0.0 <= float(self.mylar_reflectivity) <= 1.0:
            raise ValueError("Mylar reflectivity must be between 0 and 1")
        for key in ("impact_x_mm", "impact_y_mm"):
            if not math.isfinite(float(getattr(self, key))):
                raise ValueError(f"{key} must be finite")

    def to_simulation_geometry(self, *, sensor_shape: str = "rectangular", sensor_width_mm: float = 20.0,
                               sensor_height_mm: float = 10.0, sensor_depth_mm: float = 2.0,
                               sensor_radius_mm: float | None = None, sensor_semi_axis_a_mm: float | None = None,
                               sensor_semi_axis_b_mm: float | None = None):
        import quartz_cherenkov_sim as simulator
        self.validate()
        if self.guide_enabled and not self.guide_sections:
            raise ValueError(f"{self.name} guide transport is enabled but no guide sections are configured")
        if self.guide_enabled and len(self.valid_guide_sections) != len(self.guide_sections):
            invalid = next(index for index, section in enumerate(self.guide_sections) if not self._section_is_valid(section))
            raise ValueError(
                f"{self.name} guide section {invalid + 1} has unsupported shape, axis, material, or dimensions; "
                "the run was rejected before transport"
            )
        terminal_metadata = []
        for face in ("minus_x", "plus_x"):
            requested_state = self.terminal_state(face)
            coupling_status = self.terminal_coupling_status(face)
            # The transport adapter understands only Open/Mirror/Detector.
            # Invalid readout selections remain visible in metadata but cannot
            # create a detector portal.
            readout = self.readout_status()
            effective_state = requested_state
            if requested_state == "detector" and (not readout["valid"] or readout["face"] != face):
                effective_state = "open"
            terminal_metadata.append({
                "face": face, "axis_face": "YZ", "state": effective_state,
                "requested_state": requested_state, "label": coupling_status["label"],
                "valid": coupling_status["valid"], "dimensions_mm": self.terminal_dimensions_mm(face),
            })
        terminal_metadata = tuple(terminal_metadata)
        # Carry the same passive canonical sensor transform into transport.
        # Terminal state remains the only photon-interaction control.
        sensor_components = []
        for face in ("minus_x", "plus_x"):
            if self.terminal_state(face) != "detector":
                continue
            transform = canonical_sensor_transform_mm(
                self, face=face, width_mm=sensor_width_mm, height_mm=sensor_height_mm,
                depth_mm=sensor_depth_mm, shape=sensor_shape, radius_mm=sensor_radius_mm,
                semi_axis_a_mm=sensor_semi_axis_a_mm, semi_axis_b_mm=sensor_semi_axis_b_mm,
            )
            if transform is not None:
                sensor_components.append({
                    "face": transform["face"], "state": "detector",
                    "coupling_position_mm": list(transform["coupling_position_mm"]),
                    "position_mm": list(transform["position_mm"]),
                    "normal_axis": transform["normal_axis"], "normal_sign": transform["normal_sign"],
                    "target": transform["target"], "target_section": transform["target_section"],
                    "solid": transform["solid"], "status": self.readout_status()["code"],
                    # The transport aperture is the same finite window drawn
                    # by the canonical sensor component.  v50 omitted this
                    # handoff, leaving _active_area_values() with zero-sized
                    # dimensions after a complete assembly conversion.
                    "shape": transform["solid"]["shape"],
                    "dimensions_mm": dict(transform["solid"]["dimensions_mm"]),
                    "readout_radius_mm": float(transform["solid"]["readout_radius_mm"]),
                    "readout_semiaxes_mm": dict(transform["solid"]["readout_semiaxes_mm"]),
                    "window_faces_mm": transform["solid"]["coupling_face_mm"],
                })
        # The editor preserves the semantic order “active sections, then
        # guide sections”.  Transport, however, traverses physical -X→+X
        # order; sorting here also makes a valid guide attached at -X a real
        # compound instead of a false overlap.
        transport_sections = tuple(sorted(self.active_sections(), key=lambda section: _assembly_section_bounds_mm(section)[0]))
        solids = tuple(section.to_simulation_geometry() for section in transport_sections)
        if len(solids) == 1:
            return replace_simulation_terminal_metadata(solids[0], terminal_metadata, sensor_components)
        return simulator.CompoundGeometry(
            solids, allow_interface_mismatch=True,
            wrapped_external=bool(self.mylar_enabled), mylar_reflectivity=float(self.mylar_reflectivity),
            terminal_interfaces=terminal_metadata,
            sensor_components=tuple(sensor_components),
        )

    def to_dict(self) -> dict[str, Any]:
        self.validate()
        return {
            "model": "AssemblyGeometryModel", "version": ASSEMBLY_GEOMETRY_VERSION, "name": self.name,
            "sections": [section.to_dict() for section in self.sections],
            # Preserve the configured ordered guide list.  Serialization is a
            # committed operation and therefore must not silently discard a
            # configured row; malformed optional rows fail here.
            "guide_sections": [section.to_dict() for section in self.guide_sections],
            "guide_enabled": bool(self.guide_enabled),
            "guide_present": bool(self.has_light_guide),
            "terminal_plus_surface": terminal_surface_label(self.terminal_plus_surface),
            "terminal_minus_surface": terminal_surface_label(self.terminal_minus_surface),
            "mylar_enabled": bool(self.mylar_enabled), "mylar_reflectivity": float(self.mylar_reflectivity),
            "impact_x_mm": float(self.impact_x_mm), "impact_y_mm": float(self.impact_y_mm),
        }

    @classmethod
    def from_dict(cls, payload: Mapping[str, Any], name: str | None = None) -> "AssemblyGeometryModel":
        sections = [AssemblySection.from_dict(item, i) for i, item in enumerate(payload.get("sections", ())) if isinstance(item, Mapping)]
        guides = [AssemblySection.from_dict(item, i) for i, item in enumerate(payload.get("guide_sections", ())) if isinstance(item, Mapping)]
        guide_enabled = str(payload.get("guide_enabled", False)).strip().lower() in {"1", "true", "yes", "on", "enabled"}
        plus_surface, minus_surface = terminal_surfaces_from_payload(payload)
        legacy_guide_present = guide_enabled and bool(guides)
        if legacy_guide_present and _is_legacy_guide_terminal(payload.get("terminal_plus_surface")):
            plus_surface = "Detector/Photosensor coupling"
        if legacy_guide_present and _is_legacy_guide_terminal(payload.get("terminal_minus_surface")):
            minus_surface = "Detector/Photosensor coupling"
        model = cls(
            name=str(name or payload.get("name", "Assembly")), sections=sections, guide_sections=guides,
            guide_enabled=guide_enabled,
            terminal_plus_surface=plus_surface,
            terminal_minus_surface=minus_surface,
            mylar_enabled=str(payload.get("mylar_enabled", True)).strip().lower() in {"1", "true", "yes", "on", "enabled"},
            mylar_reflectivity=float(payload.get("mylar_reflectivity", 0.95)),
            impact_x_mm=float(payload.get("impact_x_mm", 0.0)), impact_y_mm=float(payload.get("impact_y_mm", 0.0)),
        )
        model.validate()
        return model

    def complete_scene(self, *, sensor_label: str = "Photosensor", sensor_coupling: str = "direct",
                       sensor_automatic_position: bool = True, sensor_shape: str = "rectangular",
                       sensor_width_mm: float | None = None, sensor_height_mm: float | None = None,
                       sensor_radius_mm: float | None = None, sensor_semi_axis_a_mm: float | None = None,
                       sensor_semi_axis_b_mm: float | None = None) -> dict[str, Any]:
        """Build the canonical complete local assembly scene.

        The returned dictionary is intentionally renderer-neutral.  It is the
        same object family consumed by the Tk preview, headless tests, global
        scene metadata, and the optical geometry adapter.  ``views`` contains
        synchronized physical-coordinate payloads for XZ, XY, and YZ.
        """

        return build_complete_assembly_scene(
            self, sensor_label=sensor_label, sensor_coupling=sensor_coupling,
            sensor_automatic_position=sensor_automatic_position, sensor_shape=sensor_shape,
            sensor_width_mm=sensor_width_mm, sensor_height_mm=sensor_height_mm,
            sensor_radius_mm=sensor_radius_mm, sensor_semi_axis_a_mm=sensor_semi_axis_a_mm,
            sensor_semi_axis_b_mm=sensor_semi_axis_b_mm,
        )


def _strict_serialized_section(payload: Any, index: int, role: str) -> AssemblySection:
    """Parse one section without the forgiving editor defaults.

    ``AssemblySection.from_dict`` intentionally supports old/in-progress GUI
    rows.  A file load has a stronger contract: malformed numeric input and
    missing shape-specific taper dimensions must be rejected before commit.
    """

    if not isinstance(payload, Mapping):
        raise ValueError(f"{role} section {index + 1} must be a JSON object")
    required = {"name", "shape", "material"}
    missing = sorted(required - set(payload))
    if missing:
        raise ValueError(f"{role} section {index + 1} is missing: {', '.join(missing)}")
    shape = str(payload["shape"]).strip().lower()
    aliases = {"rectangular": "box", "rectangular_prism": "box", "tapered box": "tapered_rectangular_prism"}
    shape = aliases.get(shape, shape)
    if shape not in ASSEMBLY_SECTION_SHAPES:
        raise ValueError(f"{role} section {index + 1} has unsupported shape: {shape}")
    if role == "guide":
        material = str(payload["material"]).strip().lower().replace("_", " ")
        if material not in {"light guide", "fused silica", "quartz", "scintillator", "air", "custom"}:
            raise ValueError(f"guide section {index + 1} has unsupported material: {payload['material']}")
    shape_required = {
        "box": ("length_x_mm", "width_y_mm", "thickness_z_mm"),
        "tapered_rectangular_prism": ("length_x_mm", "width_y_start_mm", "width_y_end_mm", "thickness_z_start_mm", "thickness_z_end_mm"),
        "cylinder": ("radius_mm", "length_mm"),
        "elliptical_cylinder": ("semi_axis_a_mm", "semi_axis_b_mm", "length_mm"),
        "truncated_cone": ("radius_start_mm", "radius_end_mm", "length_mm"),
    }[shape]
    missing = [key for key in shape_required if key not in payload]
    if missing:
        raise ValueError(f"{role} section {index + 1} is missing dimensions: {', '.join(missing)}")
    numeric_fields = (
        "length_x_mm", "width_y_mm", "thickness_z_mm", "width_y_start_mm", "width_y_end_mm",
        "thickness_z_start_mm", "thickness_z_end_mm", "radius_mm", "length_mm", "semi_axis_a_mm",
        "semi_axis_b_mm", "radius_start_mm", "radius_end_mm", "center_x_mm", "center_y_mm", "center_z_mm",
    )
    clean = dict(payload)
    for key in numeric_fields:
        if key not in payload:
            continue
        try:
            value = float(payload[key])
        except (TypeError, ValueError) as exc:
            raise ValueError(f"{role} section {index + 1} {key} must be numeric") from exc
        if not math.isfinite(value):
            raise ValueError(f"{role} section {index + 1} {key} must be finite")
        clean[key] = value
    section = AssemblySection.from_dict(clean, index)
    section.validate()
    return section


def assembly_geometry_to_dict(model: AssemblyGeometryModel, assembly_type: str) -> dict[str, Any]:
    """Return the versioned, detector-specific canonical assembly payload."""

    kind = str(assembly_type).strip().lower()
    if kind not in ASSEMBLY_GEOMETRY_TYPES:
        raise ValueError(f"unsupported assembly type: {assembly_type}")
    model.validate()
    # Validate every configured guide section, including disabled/inactive
    # rows, before writing.  This makes save deterministic and lossless.
    for index, section in enumerate(model.guide_sections):
        _strict_serialized_section(section.to_dict(), index, "guide")
    payload = model.to_dict()
    payload.update({
        "format": ASSEMBLY_GEOMETRY_FORMAT,
        "schema_version": ASSEMBLY_GEOMETRY_VERSION,
        "assembly_type": kind,
        "detector_type": kind,
        "type": kind,
        "terminal_surfaces": {
            "+X": terminal_surface_label(model.terminal_plus_surface),
            "-X": terminal_surface_label(model.terminal_minus_surface),
        },
        # v27 canonical terminals are authoritative; keep the old detector
        # control alias for clients that only understand the -X Run mode.
        "opposite_end_mode": normalize_terminal_surface(model.terminal_minus_surface),
        "metadata": {
            "coordinate_system": "assembly-local millimetres",
            "section_order": "active sections followed by ordered guide sections",
            "terminal_faces": ["minus_x", "plus_x"],
            "impact_coordinates": ["impact_x_mm", "impact_y_mm"],
            "impact_z": "derived; no independent impact Z is serialized",
        },
    })
    return payload


def assembly_geometry_from_dict(payload: Mapping[str, Any], expected_type: str, *, name: str | None = None) -> AssemblyGeometryModel:
    """Strictly parse and validate a detector-specific assembly payload.

    The function has no side effects.  Callers can therefore parse a file,
    build all derived metadata, and only then replace their live model.
    Unmarked v26 payloads remain readable when their name identifies the same
    detector; a clearly opposite detector name is rejected.
    """

    kind = str(expected_type).strip().lower()
    if kind not in ASSEMBLY_GEOMETRY_TYPES:
        raise ValueError(f"unsupported expected assembly type: {expected_type}")
    if not isinstance(payload, Mapping):
        raise ValueError("assembly geometry JSON must be an object")
    markers = [payload[key] for key in ("assembly_type", "detector_type", "type") if key in payload]
    if any(str(marker).strip().lower() != kind for marker in markers):
        raise ValueError(f"geometry file is for {markers[0] if markers else 'an incompatible detector'}, not {kind}")
    marker = markers[0] if markers else None
    if payload.get("format") is not None and payload.get("format") != ASSEMBLY_GEOMETRY_FORMAT:
        raise ValueError("unsupported geometry JSON format")
    version = payload.get("schema_version", payload.get("version"))
    if version is not None:
        try:
            version_number = int(version)
        except (TypeError, ValueError) as exc:
            raise ValueError("geometry JSON schema version must be an integer") from exc
        if version_number > ASSEMBLY_GEOMETRY_VERSION:
            raise ValueError(f"geometry JSON schema version {version_number} is newer than supported version {ASSEMBLY_GEOMETRY_VERSION}")
    file_name = str(payload.get("name", "")).strip().lower()
    if marker is None and file_name:
        opposite = "scintillator" if kind == "cherenkov" else "cherenkov"
        if opposite in file_name:
            raise ValueError(f"untyped geometry file is identified as {opposite}, not {kind}")
    required = ("sections", "guide_sections", "guide_enabled", "mylar_enabled", "mylar_reflectivity", "impact_x_mm", "impact_y_mm")
    missing = [key for key in required if key not in payload]
    if missing:
        raise ValueError(f"assembly geometry JSON is missing: {', '.join(missing)}")
    if "impact_z_mm" in payload:
        raise ValueError("independent impact Z is not supported")
    if not isinstance(payload["sections"], (list, tuple)) or not payload["sections"]:
        raise ValueError("assembly geometry JSON requires at least one active section")
    if not isinstance(payload["guide_sections"], (list, tuple)):
        raise ValueError("guide_sections must be an ordered JSON array")
    sections = [_strict_serialized_section(item, index, "active") for index, item in enumerate(payload["sections"])]
    guides = [_strict_serialized_section(item, index, "guide") for index, item in enumerate(payload["guide_sections"])]
    def strict_bool(value: Any, key: str) -> bool:
        if isinstance(value, bool):
            return value
        text = str(value).strip().lower()
        if text in {"1", "true", "yes", "on", "enabled"}:
            return True
        if text in {"0", "false", "no", "off", "disabled"}:
            return False
        raise ValueError(f"{key} must be boolean")
    try:
        guide_enabled = strict_bool(payload["guide_enabled"], "guide_enabled")
        mylar_enabled = strict_bool(payload["mylar_enabled"], "mylar_enabled")
        mylar_reflectivity = float(payload["mylar_reflectivity"])
        impact_x = float(payload["impact_x_mm"])
        impact_y = float(payload["impact_y_mm"])
    except (TypeError, ValueError) as exc:
        raise ValueError("assembly-level numeric or boolean field is invalid") from exc
    if not all(math.isfinite(value) for value in (mylar_reflectivity, impact_x, impact_y)):
        raise ValueError("assembly-level numeric fields must be finite")
    try:
        plus, minus = terminal_surfaces_from_payload(payload, require_explicit_endpoints=True)
    except ValueError as exc:
        raise ValueError(f"invalid terminal surface state: {exc}") from exc
    model = AssemblyGeometryModel(
        name=str(name or payload.get("name", kind.title())), sections=sections, guide_sections=guides,
        guide_enabled=guide_enabled, terminal_plus_surface=plus, terminal_minus_surface=minus,
        mylar_enabled=mylar_enabled, mylar_reflectivity=mylar_reflectivity,
        impact_x_mm=impact_x, impact_y_mm=impact_y,
    )
    model.validate()
    return model


def _strict_complete_sensor_settings(payload: Any) -> dict[str, Any]:
    """Validate the readout portion of a complete v32 assembly file.

    Sensor placement is retained as diagnostic metadata only.  The canonical
    terminal/guide transform is recomputed by the caller after this validation
    succeeds, so a hand-edited placement cannot detach a sensor from the
    assembly.
    """

    if not isinstance(payload, Mapping):
        raise ValueError("photosensor/readout must be an object")
    required = {
        "present", "shape", "dimensions_mm", "quantum_efficiency",
        "acceptance_deg", "optical_window_coupling_index",
        "collection_efficiency", "timing_readout", "automatic_position",
    }
    missing = sorted(required - set(payload))
    if missing:
        raise ValueError(f"photosensor/readout is missing: {', '.join(missing)}")

    def strict_bool(value: Any, key: str) -> bool:
        if isinstance(value, bool):
            return value
        text = str(value).strip().lower()
        if text in {"1", "true", "yes", "on", "enabled"}:
            return True
        if text in {"0", "false", "no", "off", "disabled"}:
            return False
        raise ValueError(f"{key} must be boolean")

    def number(value: Any, key: str, *, positive: bool = False, fraction: bool = False) -> float:
        try:
            result = float(value)
        except (TypeError, ValueError) as exc:
            raise ValueError(f"{key} must be numeric") from exc
        if not math.isfinite(result) or (positive and result <= 0.0) or (fraction and not 0.0 <= result <= 1.0):
            if fraction:
                raise ValueError(f"{key} must be between 0 and 1")
            if positive:
                raise ValueError(f"{key} must be positive")
            raise ValueError(f"{key} must be finite")
        return result

    shape = str(payload["shape"]).strip().lower()
    if shape not in {"rectangular", "circular", "elliptical"}:
        raise ValueError("photosensor/readout shape must be rectangular, circular, or elliptical")
    dimensions = payload["dimensions_mm"]
    if not isinstance(dimensions, Mapping):
        raise ValueError("photosensor/readout dimensions_mm must be an object")
    dimensions_clean = {
        "width": number(dimensions.get("width"), "photosensor width", positive=True),
        "height": number(dimensions.get("height"), "photosensor height", positive=True),
        "depth": number(dimensions.get("depth"), "photosensor depth", positive=True),
    }
    # v61 and earlier carried only full width/height.  Additive v62 keys take
    # precedence; absent keys are a lossless one-time migration from those
    # full extents (circle diameter -> radius; ellipse extents -> semiaxes).
    radius_raw = payload.get("readout_radius_mm")
    axes_raw = payload.get("readout_semiaxes_mm")
    radius = None if radius_raw in (None, "") else number(radius_raw, "readout_radius_mm", positive=True)
    axis_a = axis_b = None
    if axes_raw is not None:
        if not isinstance(axes_raw, Mapping):
            raise ValueError("readout_semiaxes_mm must be an object")
        axis_a = number(axes_raw.get("a"), "readout_semiaxes_mm.a", positive=True)
        axis_b = number(axes_raw.get("b"), "readout_semiaxes_mm.b", positive=True)
    normalized = resolve_sensor_readout_dimensions_mm(
        shape, width_mm=dimensions_clean["width"], height_mm=dimensions_clean["height"],
        radius_mm=radius, semi_axis_a_mm=axis_a, semi_axis_b_mm=axis_b,
    )
    dimensions_clean["width"] = float(normalized["width_mm"])
    dimensions_clean["height"] = float(normalized["height_mm"])
    timing = payload["timing_readout"]
    if not isinstance(timing, Mapping):
        raise ValueError("timing_readout must be an object")
    timing_clean = {
        "sigma_ns": number(timing.get("sigma_ns", 0.0), "timing_readout.sigma_ns"),
        "mode": str(timing.get("mode", "configured")),
    }
    if timing_clean["sigma_ns"] < 0.0:
        raise ValueError("timing_readout.sigma_ns must be non-negative")
    result = dict(payload)
    result["present"] = strict_bool(payload["present"], "photosensor.present")
    result["automatic_position"] = strict_bool(payload["automatic_position"], "photosensor.automatic_position")
    result["shape"] = shape
    result["dimensions_mm"] = dimensions_clean
    result["readout_radius_mm"] = float(normalized["radius_mm"])
    result["readout_semiaxes_mm"] = {"a": float(normalized["semi_axis_a_mm"]), "b": float(normalized["semi_axis_b_mm"])}
    result["quantum_efficiency"] = number(payload["quantum_efficiency"], "quantum_efficiency", fraction=True)
    result["acceptance_deg"] = number(payload["acceptance_deg"], "acceptance_deg", positive=True)
    if result["acceptance_deg"] > 180.0:
        raise ValueError("acceptance_deg must not exceed 180")
    result["optical_window_coupling_index"] = number(payload["optical_window_coupling_index"], "optical_window_coupling_index", positive=True)
    result["collection_efficiency"] = number(payload["collection_efficiency"], "collection_efficiency", fraction=True)
    result["timing_readout"] = timing_clean
    if "position_mm" in payload:
        position = payload["position_mm"]
        if not isinstance(position, (list, tuple)) or len(position) != 3:
            raise ValueError("photosensor position_mm must contain three values")
        result["position_mm"] = [number(value, f"position_mm[{index}]") for index, value in enumerate(position)]
    if "coupled_terminal" in payload and payload["coupled_terminal"] is not None:
        if str(payload["coupled_terminal"]) not in {"plus_x", "minus_x"}:
            raise ValueError("coupled_terminal must be plus_x, minus_x, or null")
    # Keep detector-specific response fields lossless while validating the
    # common numeric values above.  Unknown keys are diagnostic extensions.
    return result


def complete_detector_assembly_geometry_to_dict(
    model: AssemblyGeometryModel,
    assembly_type: str,
    photosensor: Mapping[str, Any],
    *,
    local_settings: Mapping[str, Any] | None = None,
    derived_diagnostics: Mapping[str, Any] | None = None,
) -> dict[str, Any]:
    """Serialize one complete detector assembly, including readout state.

    ``photosensor`` is the detector-specific readout state from the GUI.  The
    optional derived sections are explicitly diagnostic and are never needed
    to reconstruct the canonical geometry.
    """

    kind = str(assembly_type).strip().lower()
    if kind not in ASSEMBLY_GEOMETRY_TYPES:
        raise ValueError(f"unsupported assembly type: {assembly_type}")
    model.validate()
    for index, section in enumerate(model.guide_sections):
        _strict_serialized_section(section.to_dict(), index, "guide")
    sensor = _strict_complete_sensor_settings(photosensor)
    payload = {
        "format": COMPLETE_DETECTOR_ASSEMBLY_GEOMETRY_FORMAT,
        "type_marker": COMPLETE_DETECTOR_ASSEMBLY_GEOMETRY_MARKER,
        "schema_version": COMPLETE_DETECTOR_ASSEMBLY_GEOMETRY_VERSION,
        "detector_type": kind,
        "assembly_type": kind,
        "name": model.name,
        "active_sections": [dict(section.to_dict(), role="active", order=index) for index, section in enumerate(model.sections)],
        # ``sections`` is retained as the canonical v27-compatible field name.
        "sections": [section.to_dict() for section in model.sections],
        "light_guide": {
            "enabled": bool(model.guide_enabled),
            "present": bool(model.has_light_guide),
            "sections": [dict(section.to_dict(), role="guide", order=index) for index, section in enumerate(model.guide_sections)],
        },
        "guide_enabled": bool(model.guide_enabled),
        "guide_sections": [section.to_dict() for section in model.guide_sections],
        "terminal_surfaces": {
            "+X": terminal_surface_label(model.terminal_plus_surface),
            "-X": terminal_surface_label(model.terminal_minus_surface),
        },
        "terminal_plus_surface": terminal_surface_label(model.terminal_plus_surface),
        "terminal_minus_surface": terminal_surface_label(model.terminal_minus_surface),
        # Retain the Run-tab spelling as a -X compatibility alias.  Endpoint
        # fields and the structured block remain the geometry authority.
        "opposite_end_mode": normalize_terminal_surface(model.terminal_minus_surface),
        "mylar": {
            "enabled": bool(model.mylar_enabled),
            "reflectivity": float(model.mylar_reflectivity),
            "overlay_surface_count": len(derive_mylar_overlay_surfaces(model)),
        },
        "mylar_enabled": bool(model.mylar_enabled),
        "mylar_reflectivity": float(model.mylar_reflectivity),
        "impact_xy_mm": {"x": float(model.impact_x_mm), "y": float(model.impact_y_mm)},
        "impact_x_mm": float(model.impact_x_mm),
        "impact_y_mm": float(model.impact_y_mm),
        "photosensor_readout": sensor,
        "local_assembly_settings": dict(local_settings or {
            "coordinate_system": "assembly-local millimetres",
            "global_placement_offsets_excluded": True,
            "sensor_transform_source": "canonical terminal/guide output; recompute on load",
        }),
        "derived_diagnostics": dict(derived_diagnostics or {}),
        "metadata": {
            "section_order": "ordered active sections, then ordered optional guide sections",
            "terminal_faces": ["+X", "-X"],
            "impact_coordinates": ["impact_x_mm", "impact_y_mm"],
            "impact_z": "derived; no independent impact Z is serialized",
        },
    }
    return payload


def complete_detector_assembly_geometry_from_dict(
    payload: Mapping[str, Any], expected_type: str, *, name: str | None = None,
) -> tuple[AssemblyGeometryModel, dict[str, Any]]:
    """Strictly parse a complete v32 envelope without mutating live state."""

    kind = str(expected_type).strip().lower()
    if kind not in ASSEMBLY_GEOMETRY_TYPES:
        raise ValueError(f"unsupported expected assembly type: {expected_type}")
    if not isinstance(payload, Mapping):
        raise ValueError("complete detector assembly JSON must be an object")
    if payload.get("format") != COMPLETE_DETECTOR_ASSEMBLY_GEOMETRY_FORMAT or payload.get("type_marker") != COMPLETE_DETECTOR_ASSEMBLY_GEOMETRY_MARKER:
        raise ValueError("not a complete-detector-assembly-geometry JSON file")
    try:
        version = int(payload.get("schema_version"))
    except (TypeError, ValueError) as exc:
        raise ValueError("complete geometry schema_version must be an integer") from exc
    if version > COMPLETE_DETECTOR_ASSEMBLY_GEOMETRY_VERSION:
        raise ValueError(f"complete geometry schema version {version} is newer than supported version {COMPLETE_DETECTOR_ASSEMBLY_GEOMETRY_VERSION}")
    for key in ("detector_type", "assembly_type"):
        if str(payload.get(key, "")).strip().lower() != kind:
            raise ValueError(f"complete geometry file is not for {kind}")
    required = ("sections", "guide_sections", "guide_enabled", "mylar_enabled", "mylar_reflectivity", "impact_x_mm", "impact_y_mm", "photosensor_readout")
    missing = [key for key in required if key not in payload]
    if missing:
        raise ValueError(f"complete geometry JSON is missing: {', '.join(missing)}")
    if "impact_z_mm" in payload:
        raise ValueError("independent impact Z is not supported")
    if not isinstance(payload["sections"], (list, tuple)) or not payload["sections"]:
        raise ValueError("complete geometry JSON requires at least one active section")
    if not isinstance(payload["guide_sections"], (list, tuple)):
        raise ValueError("guide_sections must be an ordered JSON array")
    sections = [_strict_serialized_section(item, index, "active") for index, item in enumerate(payload["sections"])]
    guides = [_strict_serialized_section(item, index, "guide") for index, item in enumerate(payload["guide_sections"])]

    def strict_bool(value: Any, key: str) -> bool:
        if isinstance(value, bool):
            return value
        text = str(value).strip().lower()
        if text in {"1", "true", "yes", "on", "enabled"}:
            return True
        if text in {"0", "false", "no", "off", "disabled"}:
            return False
        raise ValueError(f"{key} must be boolean")

    try:
        guide_enabled = strict_bool(payload["guide_enabled"], "guide_enabled")
        mylar_enabled = strict_bool(payload["mylar_enabled"], "mylar_enabled")
        mylar_reflectivity = float(payload["mylar_reflectivity"])
        impact_x = float(payload["impact_x_mm"])
        impact_y = float(payload["impact_y_mm"])
    except (TypeError, ValueError) as exc:
        raise ValueError("complete geometry assembly-level field is invalid") from exc
    if not all(math.isfinite(value) for value in (mylar_reflectivity, impact_x, impact_y)):
        raise ValueError("complete geometry assembly-level numeric fields must be finite")
    if "active_sections" in payload:
        if not isinstance(payload["active_sections"], (list, tuple)) or len(payload["active_sections"]) != len(sections):
            raise ValueError("active_sections must match the ordered sections array")
        for index, item in enumerate(payload["active_sections"]):
            if not isinstance(item, Mapping) or str(item.get("role", "active")) != "active" or int(item.get("order", index)) != index:
                raise ValueError("active_sections role/order metadata is invalid")
    if "light_guide" in payload:
        guide_block = payload["light_guide"]
        if not isinstance(guide_block, Mapping):
            raise ValueError("light_guide must be an object")
        if strict_bool(guide_block.get("enabled"), "light_guide.enabled") != guide_enabled:
            raise ValueError("light_guide.enabled disagrees with guide_enabled")
        if not isinstance(guide_block.get("sections"), (list, tuple)) or len(guide_block["sections"]) != len(guides):
            raise ValueError("light_guide.sections must match guide_sections")
        present = strict_bool(guide_block.get("present", False), "light_guide.present")
        if present != bool(guide_enabled and guides):
            raise ValueError("light_guide.present is inconsistent with enabled valid sections")
        for index, item in enumerate(guide_block["sections"]):
            if not isinstance(item, Mapping) or str(item.get("role", "guide")) != "guide" or int(item.get("order", index)) != index:
                raise ValueError("light_guide section role/order metadata is invalid")
    if "mylar" in payload:
        mylar_block = payload["mylar"]
        if not isinstance(mylar_block, Mapping) or strict_bool(mylar_block.get("enabled"), "mylar.enabled") != mylar_enabled or abs(float(mylar_block.get("reflectivity")) - mylar_reflectivity) > 1e-12:
            raise ValueError("mylar block disagrees with Mylar fields")
    if "impact_xy_mm" in payload:
        impact_block = payload["impact_xy_mm"]
        if not isinstance(impact_block, Mapping) or float(impact_block.get("x")) != impact_x or float(impact_block.get("y")) != impact_y:
            raise ValueError("impact_xy_mm disagrees with impact_x_mm/impact_y_mm")
    if "local_assembly_settings" in payload:
        local_settings = payload["local_assembly_settings"]
        if not isinstance(local_settings, Mapping):
            raise ValueError("local_assembly_settings must be an object")
        if "global_placement_offsets_excluded" in local_settings and not strict_bool(local_settings["global_placement_offsets_excluded"], "global_placement_offsets_excluded"):
            raise ValueError("global placement offsets cannot be loaded as local geometry")
    # Validate cross-field agreement after normalization, so older vXX aliases
    # such as ``plus_x`` and lower-case values retain meaning while conflicting
    # geometry is rejected before commit.
    plus, minus = terminal_surfaces_from_payload(payload, require_explicit_endpoints=True)
    model = AssemblyGeometryModel(
        name=str(name or payload.get("name", kind.title())), sections=sections, guide_sections=guides,
        guide_enabled=guide_enabled, terminal_plus_surface=plus, terminal_minus_surface=minus,
        mylar_enabled=mylar_enabled, mylar_reflectivity=mylar_reflectivity,
        impact_x_mm=impact_x, impact_y_mm=impact_y,
    )
    model.validate()
    sensor = _strict_complete_sensor_settings(payload["photosensor_readout"])
    # Validate the derived scene now, while this function remains side-effect
    # free.  The returned placement is deliberately not used as authority.
    model.complete_scene(
        sensor_label="Photosensor" if kind == "cherenkov" else "Scintillator readout",
        sensor_shape=sensor["shape"], sensor_width_mm=sensor["dimensions_mm"]["width"],
        sensor_height_mm=sensor["dimensions_mm"]["height"],
    )
    return model, sensor


# Short aliases make the v32 API discoverable alongside the v27 helpers.
complete_assembly_geometry_to_dict = complete_detector_assembly_geometry_to_dict
complete_assembly_geometry_from_dict = complete_detector_assembly_geometry_from_dict


def _assembly_section_bounds_mm(section: AssemblySection) -> tuple[float, float]:
    if section.shape in {"box", "tapered_rectangular_prism"}:
        length = section.length_x_mm
    else:
        length = section.length_mm if section.axis == "X" else section.length_x_mm
    return section.center_x_mm - length / 2.0, section.center_x_mm + length / 2.0


def _assembly_terminal_x_mm(model: AssemblyGeometryModel, face: str) -> float:
    sections = model.active_sections()
    low, _ = _assembly_section_bounds_mm(sections[0])
    _, high = _assembly_section_bounds_mm(sections[-1])
    key = str(face).strip().lower().replace("-", "_")
    return float(high if key in {"plus", "plus_x", "+x", "x+"} else low)


def canonical_readout_frame_mm(model: AssemblyGeometryModel, face: str) -> dict[str, Any]:
    """Return the actual coupling/output face used by the readout solid.

    A direct sensor couples to the selected active terminal.  When a valid
    guide is attached at that terminal, the coupling plane is the guide's
    outside output instead.  Keeping this small frame helper beside the
    sensor transform prevents a marker, sensor body, and transport portal
    from quietly choosing three different planes.
    """

    key = str(face).strip().lower().replace("-", "_")
    if key not in {"minus", "minus_x", "-x", "x-", "plus", "plus_x", "+x", "x+"}:
        raise ValueError("face must identify plus_x or minus_x")
    plus = key in {"plus", "plus_x", "+x", "x+"}
    readout = model.readout_status()
    selected = str(readout.get("face") or "") == ("plus_x" if plus else "minus_x")
    guide_at_selected_terminal = bool(model.has_light_guide and model.guide_attached_at_terminal(face))
    use_guide_output = bool(selected and readout.get("valid") and guide_at_selected_terminal)
    target_sections = list(model.valid_guide_sections if use_guide_output else model.active_sections())
    if not target_sections:
        raise ValueError(f"{model.name} has no geometry at readout face {face}")
    target_section = target_sections[-1] if plus else target_sections[0]
    low_x, high_x = _assembly_section_bounds_mm(target_section)
    dimensions = model.readout_terminal_dimensions_mm(face) if use_guide_output or selected else model.terminal_dimensions_mm(face)
    return {
        "face": "plus_x" if plus else "minus_x",
        "x_mm": float(high_x if plus else low_x),
        "center_y_mm": float(dimensions["center_y_mm"]),
        "center_z_mm": float(dimensions["center_z_mm"]),
        "width_y_mm": float(dimensions["width_y_mm"]),
        "thickness_z_mm": float(dimensions["thickness_z_mm"]),
        "target": "guide output" if use_guide_output and plus else (
            "guide input/output terminal" if use_guide_output else "active terminal face"
        ),
        "target_section": target_section.name,
        "guide_output": use_guide_output,
        "selected_readout": selected,
    }


def _assembly_section_dimensions(section: AssemblySection, end: str) -> tuple[float, float]:
    if section.shape == "box":
        return section.width_y_mm, section.thickness_z_mm
    if section.shape == "tapered_rectangular_prism":
        return getattr(section, f"width_y_{end}_mm"), getattr(section, f"thickness_z_{end}_mm")
    if section.shape == "cylinder":
        return 2.0 * section.radius_mm, 2.0 * section.radius_mm
    if section.shape == "elliptical_cylinder":
        return 2.0 * section.semi_axis_a_mm, 2.0 * section.semi_axis_b_mm
    radius = section.radius_end_mm if end == "end" else section.radius_start_mm
    return 2.0 * radius, 2.0 * radius


def _ellipse_points(center_y: float, center_z: float, width: float, height: float, count: int = 32) -> list[list[float]]:
    return [
        [center_y + width / 2.0 * math.cos(2.0 * math.pi * i / count),
         center_z + height / 2.0 * math.sin(2.0 * math.pi * i / count)]
        for i in range(count + 1)
    ]


def canonical_section_faces_mm(section: AssemblySection, samples: int = 32) -> list[list[list[float]]]:
    """Return renderer faces from the same solid used by optical transport.

    This is intentionally a model-layer helper.  Renderers may project or
    shade these polygons, but they do not manufacture a second simplified
    geometry model for Mylar or for the physical volume.
    """

    if samples < 4:
        raise ValueError("surface face samples must be at least 4")
    solid = section.to_simulation_geometry()
    rings = [
        [[float(value) * 1e3 for value in point] for point in ring]
        for ring in solid.surface_rings(samples=int(samples))
    ]
    if len(rings) != 2 or any(len(ring) < 3 for ring in rings):
        raise ValueError("solid surface rings must contain two polygonal caps")
    first, second = rings
    faces: list[list[list[float]]] = [first, list(reversed(second))]
    count = min(len(first), len(second))
    for index in range(count):
        next_index = (index + 1) % count
        faces.append([first[index], first[next_index], second[next_index], second[index]])
    return faces


def _face_is_x_cap(face: Sequence[Sequence[float]]) -> bool:
    values = [float(point[0]) for point in face]
    tolerance = max(1e-8, max(abs(value) for value in values) * 1e-9)
    return max(values) - min(values) <= tolerance


def _faces_share_external_interface(first: Sequence[Sequence[float]], second: Sequence[Sequence[float]]) -> bool:
    """Return whether two X-cap polygons occupy the same physical plane.

    The transverse overlap test handles rectangular, tapered, and sampled
    rotational sections without depending on face ordering.
    """

    if not (_face_is_x_cap(first) and _face_is_x_cap(second)):
        return False
    x_first = sum(float(point[0]) for point in first) / len(first)
    x_second = sum(float(point[0]) for point in second) / len(second)
    if abs(x_first - x_second) > max(1e-8, abs(x_first) * 1e-9, abs(x_second) * 1e-9):
        return False
    first_y = [float(point[1]) for point in first]
    first_z = [float(point[2]) for point in first]
    second_y = [float(point[1]) for point in second]
    second_z = [float(point[2]) for point in second]
    overlap_y = min(max(first_y), max(second_y)) - max(min(first_y), min(second_y))
    overlap_z = min(max(first_z), max(second_z)) - max(min(first_z), min(second_z))
    return overlap_y > 1e-8 and overlap_z > 1e-8


def _display_offset_face(face: Sequence[Sequence[float]], section: AssemblySection,
                         epsilon_mm: float = MYLAR_DISPLAY_EPSILON_MM) -> list[list[float]]:
    """Copy one eligible physical face with a display-only outward offset."""

    points = np.asarray(face, dtype=float)
    center = np.asarray((section.center_x_mm, section.center_y_mm, section.center_z_mm), dtype=float)
    direction = np.mean(points, axis=0) - center
    norm = float(np.linalg.norm(direction))
    if norm <= 1e-10 and len(points) >= 3:
        direction = np.cross(points[1] - points[0], points[2] - points[0])
        norm = float(np.linalg.norm(direction))
    if norm <= 1e-10:
        direction = np.asarray((0.0, 0.0, 1.0))
        norm = 1.0
    offset = direction / norm * float(epsilon_mm)
    return [[float(value) for value in point + offset] for point in points]


def derive_mylar_overlay_surfaces(model: AssemblyGeometryModel, *, samples: int = 32) -> list[dict[str, Any]]:
    """Derive explicit display-only Mylar surfaces from canonical section faces.

    Eligible faces are external side faces plus genuinely open external
    terminal caps.  Shared active/guide caps and any Mirror or Detector
    terminal cap are excluded.  Sensor bodies/windows are not section faces
    and therefore can never receive a Mylar overlay here.
    """

    model.validate()
    sections = list(model.active_sections())
    face_sets = [canonical_section_faces_mm(section, samples=samples) for section in sections]
    excluded: dict[tuple[int, int], str] = {}

    # Shared caps are optical interfaces, irrespective of terminal labels.
    for first_index, first_faces in enumerate(face_sets):
        for second_index in range(first_index + 1, len(face_sets)):
            for first_face_index, first_face in enumerate(first_faces):
                if not _face_is_x_cap(first_face):
                    continue
                for second_face_index, second_face in enumerate(face_sets[second_index]):
                    if _faces_share_external_interface(first_face, second_face):
                        excluded[(first_index, first_face_index)] = "active/light-guide optical interface"
                        excluded[(second_index, second_face_index)] = "active/light-guide optical interface"

    low_x = min(_assembly_section_bounds_mm(section)[0] for section in sections)
    high_x = max(_assembly_section_bounds_mm(section)[1] for section in sections)
    for section_index, section in enumerate(sections):
        section_low, section_high = _assembly_section_bounds_mm(section)
        for face_index, face in enumerate(face_sets[section_index]):
            if not _face_is_x_cap(face):
                continue
            x_value = sum(float(point[0]) for point in face) / len(face)
            at_low = abs(x_value - low_x) <= max(1e-8, abs(low_x) * 1e-9)
            at_high = abs(x_value - high_x) <= max(1e-8, abs(high_x) * 1e-9)
            terminal_face = "minus_x" if at_low else "plus_x" if at_high else None
            if terminal_face and model.terminal_state(terminal_face) != "open":
                excluded[(section_index, face_index)] = (
                    "mirror terminal face" if model.terminal_state(terminal_face) == "mirror"
                    else "terminal optical/readout interface"
                )

    overlays: list[dict[str, Any]] = []
    for section_index, (section, faces) in enumerate(zip(sections, face_sets)):
        eligible = [index for index in range(len(faces)) if (section_index, index) not in excluded]
        if not model.mylar_enabled:
            eligible = []
        for face_index in eligible:
            source_face = faces[face_index]
            overlays.append({
                "name": f"Mylar • {section.name} • face {face_index}",
                "kind": "mylar_overlay",
                "label": f"Mylar ({model.mylar_reflectivity * 100.0:g}%) • {section.name}",
                "color": MYLAR_DISPLAY_COLOR,
                "edge_color": MYLAR_DISPLAY_EDGE,
                "hatch": MYLAR_DISPLAY_HATCH,
                "stripe_pattern": "diagonal-gold",
                "section_name": section.name,
                "section_index": section_index,
                "source_face_index": face_index,
                "source_faces_mm": [deepcopy(source_face)],
                "faces_mm": [_display_offset_face(source_face, section)],
                "display_offset_mm": MYLAR_DISPLAY_EPSILON_MM,
                "display_only": True,
                "reflectivity": float(model.mylar_reflectivity),
                "excluded_surface_reasons": [],
            })
    return overlays


def sensor_solid_geometry_mm(*, center_mm: Sequence[float], width_mm: float,
                             height_mm: float, depth_mm: float,
                             shape: str = "rectangular", normal_axis: str = "X",
                             radius_mm: float | None = None,
                             semi_axis_a_mm: float | None = None,
                             semi_axis_b_mm: float | None = None,
                             normal_sign: float = 1.0) -> dict[str, Any]:
    """Return one physical sensor body and its coupling window in millimetres.

    Rectangular sensors use a six-faced box. Circular and elliptical sensors
    use a sampled elliptical cylinder. The cap nearest the coupled surface is
    returned as ``coupling_face_mm`` and is reused as the visible window.
    """

    dimensions = resolve_sensor_readout_dimensions_mm(
        shape, width_mm=width_mm, height_mm=height_mm, radius_mm=radius_mm,
        semi_axis_a_mm=semi_axis_a_mm, semi_axis_b_mm=semi_axis_b_mm,
    )
    width = float(dimensions["width_mm"])
    height = float(dimensions["height_mm"])
    depth = float(depth_mm)
    if min(width, height, depth) <= 0.0 or not all(math.isfinite(value) for value in (width, height, depth)):
        raise ValueError("sensor body dimensions must be finite and positive")
    axis = str(normal_axis).strip().upper()
    if axis not in {"X", "Y", "Z"}:
        raise ValueError("sensor body normal_axis must be X, Y, or Z")
    sign = 1.0 if float(normal_sign) >= 0.0 else -1.0
    center = tuple(float(value) for value in center_mm)
    if len(center) != 3 or not all(math.isfinite(value) for value in center):
        raise ValueError("sensor body center must contain three finite values")
    axis_index = {"X": 0, "Y": 1, "Z": 2}[axis]
    near_position = center[axis_index] - sign * depth / 2.0
    far_position = center[axis_index] + sign * depth / 2.0
    shape_key = str(shape).strip().lower()
    if shape_key in {"circular", "circle"}:
        ring_count = 32
        transverse_a = transverse_b = 2.0 * float(dimensions["radius_mm"])
        canonical_shape = "circular"
    elif shape_key in {"elliptical", "ellipse"}:
        ring_count = 32
        transverse_a = 2.0 * float(dimensions["semi_axis_a_mm"])
        transverse_b = 2.0 * float(dimensions["semi_axis_b_mm"])
        canonical_shape = "elliptical"
    else:
        ring_count = 4
        transverse_a, transverse_b = width, height
        canonical_shape = "rectangular"

    def ring(position: float) -> list[list[float]]:
        if ring_count == 4:
            transverse = [(-1.0, -1.0), (1.0, -1.0), (1.0, 1.0), (-1.0, 1.0)]
        else:
            transverse = [
                (math.cos(2.0 * math.pi * index / ring_count),
                 math.sin(2.0 * math.pi * index / ring_count))
                for index in range(ring_count)
            ]
        result = []
        for first, second in transverse:
            u, v = first * transverse_a / 2.0, second * transverse_b / 2.0
            if axis == "X":
                result.append([position, center[1] + u, center[2] + v])
            elif axis == "Y":
                result.append([center[0] + u, position, center[2] + v])
            else:
                result.append([center[0] + u, center[1] + v, position])
        return result

    near = ring(near_position)
    far = ring(far_position)
    faces = [near, far]
    for index in range(len(near)):
        next_index = (index + 1) % len(near)
        faces.append([near[index], near[next_index], far[next_index], far[index]])
    return {
        "shape": canonical_shape,
        "faces_mm": faces,
        "coupling_face_mm": [near],
        "outer_face_mm": [far],
        "center_mm": list(center),
        "dimensions_mm": {"width": width, "height": height, "depth": depth},
        "readout_radius_mm": float(dimensions["radius_mm"]),
        "readout_semiaxes_mm": {"a": float(dimensions["semi_axis_a_mm"]), "b": float(dimensions["semi_axis_b_mm"])},
        "normal_axis": axis,
        "normal_sign": sign,
    }


def canonical_sensor_transform_mm(model: AssemblyGeometryModel, *, face: str,
                                  width_mm: float, height_mm: float,
                                  depth_mm: float, shape: str = "rectangular",
                                  radius_mm: float | None = None,
                                  semi_axis_a_mm: float | None = None,
                                  semi_axis_b_mm: float | None = None) -> dict[str, Any] | None:
    """Return the one terminal/output transform used by all sensor renderers.

    The coupling point is a real X terminal of the complete ordered assembly.
    When a valid enabled guide exists, the selected terminal is resolved to
    that guide's outside output (the last section at +X or first at -X).  The
    sensor body is then built outward from that face; its centre is never used
    as the coupling plane.
    """

    state = model.terminal_state(face)
    readout = model.readout_status()
    # A configured detector terminal always receives a physical sensor.  An
    # ambiguous or invalid readout remains visible and is reported in the
    # status fields instead of causing the component to disappear.
    if state != "detector":
        return None
    key = str(face).strip().lower().replace("-", "_")
    plus = key in {"plus", "plus_x", "+x", "x+"}
    direction = 1.0 if plus else -1.0
    frame = canonical_readout_frame_mm(model, face)
    target_section = next(
        section for section in (model.valid_guide_sections if frame["guide_output"] else model.active_sections())
        if section.name == frame["target_section"]
    )
    target_name = str(frame["target"])
    coupling_x = float(frame["x_mm"])
    terminal_width, terminal_height = float(frame["width_y_mm"]), float(frame["thickness_z_mm"])
    center_y = float(frame["center_y_mm"])
    center_z = float(frame["center_z_mm"])
    body_center = [coupling_x + direction * float(depth_mm) / 2.0, center_y, center_z]
    solid = sensor_solid_geometry_mm(
        center_mm=body_center, width_mm=width_mm, height_mm=height_mm,
        depth_mm=depth_mm, shape=shape, radius_mm=radius_mm,
        semi_axis_a_mm=semi_axis_a_mm, semi_axis_b_mm=semi_axis_b_mm,
        normal_axis="X", normal_sign=direction,
    )
    points = [point for face_points in solid["faces_mm"] for point in face_points]
    bounds = [min(point[i] for point in points) for i in range(3)] + [max(point[i] for point in points) for i in range(3)]
    coupling = [float(coupling_x), center_y, center_z]
    return {
        "face": "plus_x" if plus else "minus_x", "state": state,
        "target": target_name, "target_section": target_section.name,
        "coupling_position_mm": coupling, "position_mm": list(body_center),
        "normal_axis": "X", "normal_sign": direction,
        "normal_vector_mm": [direction, 0.0, 0.0],
        "display_offset_mm": [0.0, 0.0, 0.0],
        "terminal_dimensions_mm": {"width_y_mm": float(terminal_width), "thickness_z_mm": float(terminal_height)},
        "solid": solid, "bounds_mm": bounds,
        "readout_valid": bool(readout["valid"]), "readout_status": readout["code"],
    }


def _orthographic_axis_payload(view_name: str, components: Sequence[Mapping[str, Any]]) -> dict[str, Any]:
    """Build visible, physical-coordinate axes for one complete-scene view.

    The local renderer is canvas based, so axis limits and labels cannot be
    delegated to Matplotlib metadata. This payload is the renderer-neutral
    source for the actual line, arrow, label, and origin marker primitives.
    Its camera bounds are derived from the view's real component polygons.
    """

    if view_name == "xz":
        axis_names = ("X", "Z")
        direction_label = "+Z beam • +X longitudinal"
        beam_label = "+Z beam"
        terminal_label = ""
    elif view_name == "xy":
        axis_names = ("X", "Y")
        direction_label = "+Z beam normal • +X longitudinal"
        beam_label = "+Z beam normal"
        terminal_label = ""
    elif view_name == "yz":
        axis_names = ("Y", "Z")
        direction_label = "+X → +X terminal • −X → −X terminal"
        beam_label = "+Z beam normal"
        terminal_label = "+X / −X terminal direction"
    else:
        raise ValueError(f"unsupported orthographic view: {view_name}")

    points = [
        point for item in components for point in item.get("polygon", ())
        if len(point) >= 2
    ]
    if not points:
        raise ValueError(f"{view_name} view has no physical component bounds")
    min_a = min(float(point[0]) for point in points)
    max_a = max(float(point[0]) for point in points)
    min_b = min(float(point[1]) for point in points)
    max_b = max(float(point[1]) for point in points)
    span_a = max(max_a - min_a, 1.0)
    span_b = max(max_b - min_b, 1.0)
    pad_a = max(span_a * 0.10, 1.0)
    pad_b = max(span_b * 0.10, 1.0)
    # Put the triad in a lower-left data-space corner. The anchor and both
    # endpoints are computed from current scene bounds, never a fixed sketch.
    origin = [min_a + 0.12 * span_a, min_b + 0.13 * span_b]
    axis_length_a = max(0.23 * span_a, 0.5)
    axis_length_b = max(0.23 * span_b, 0.5)
    primitives = [
        {"kind": "line", "axis": axis_names[0], "start": list(origin),
         "end": [origin[0] + axis_length_a, origin[1]], "label": axis_names[0],
         "color": "#c62828"},
        {"kind": "line", "axis": axis_names[1], "start": list(origin),
         "end": [origin[0], origin[1] + axis_length_b], "label": axis_names[1],
         "color": "#1565c0"},
        {"kind": "origin", "point": list(origin), "label": "O",
         "color": "#263238"},
    ]
    return {
        "coordinate_system": "physical local XYZ (mm)",
        "geometry_bounds_mm": [min_a, min_b, max_a, max_b],
        "bounds_mm": [min_a - pad_a, min_b - pad_b, max_a + pad_a, max_b + pad_b],
        "origin_mm": list(origin),
        "primitives": primitives,
        "arrowheads": True,
        "labels": list(axis_names),
        "direction": direction_label,
        "beam": "+Z",
        "beam_label": beam_label,
        "terminal_direction": terminal_label,
        "legend": "O = local view origin • arrows = positive axes",
    }


def build_complete_assembly_scene(model: AssemblyGeometryModel, *, sensor_label: str = "Photosensor",
                                  sensor_coupling: str = "direct", sensor_automatic_position: bool = True,
                                  sensor_shape: str = "rectangular", sensor_width_mm: float | None = None,
                                  sensor_height_mm: float | None = None,
                                  sensor_radius_mm: float | None = None,
                                  sensor_semi_axis_a_mm: float | None = None,
                                  sensor_semi_axis_b_mm: float | None = None) -> dict[str, Any]:
    """Create geometry, terminal interfaces, markers, and three view payloads.

    State-dependent items are real scene components.  Consequently a state
    change changes the component list and interface metadata even when a
    caller chooses not to use the Tk renderer.
    """

    model.validate()
    # ``sensor_coupling`` is retained as a call-compatible legacy argument;
    # v25 never trusts it.  The terminal selection and actual guide sections
    # are the only source of optical routing truth.
    readout = model.readout_status()
    route = str(readout["mode"])
    derived_component = "Complete light guide" if route == "light_guide" else (
        "Active scintillator" if str(model.name).strip().lower().startswith("scint") else "Active quartz radiator"
    )
    sections = list(model.active_sections())
    guide_ids = {id(section) for section in model.valid_guide_sections}
    colors = {"fused silica": "#80b1d3", "quartz": "#80b1d3", "light guide": "#f0a35b",
              "scintillator": "#66c2a5", "air": "#d9d9d9", "custom": "#c5b0d5"}
    components: list[dict[str, Any]] = []
    overlay_components = derive_mylar_overlay_surfaces(model)
    overlays_by_section: dict[int, list[dict[str, Any]]] = {}
    for overlay in overlay_components:
        overlays_by_section.setdefault(int(overlay["section_index"]), []).append(overlay)
    xz: list[dict[str, Any]] = []
    xy: list[dict[str, Any]] = []
    low_x = min(_assembly_section_bounds_mm(section)[0] for section in sections)
    high_x = max(_assembly_section_bounds_mm(section)[1] for section in sections)
    for index, section in enumerate(sections):
        x0, x1 = _assembly_section_bounds_mm(section)
        width0, height0 = _assembly_section_dimensions(section, "start")
        width1, height1 = _assembly_section_dimensions(section, "end")
        cy, cz = section.center_y_mm, section.center_z_mm
        canonical_faces = canonical_section_faces_mm(section)
        section_overlays = overlays_by_section.get(index, [])
        mylar_face_indices = sorted(int(item["source_face_index"]) for item in section_overlays)
        component = {
            "name": section.name, "kind": "guide" if id(section) in guide_ids else "geometry",
            "material": section.material, "color": colors.get(section.material.lower(), "#c5b0d5"),
            "bounds_mm": [x0, cy - width0 / 2, cz - height0 / 2, x1, cy + width1 / 2, cz + height1 / 2],
            "faces_mm": canonical_faces,
            "mylar_faces": ["eligible_external_faces"] if mylar_face_indices else [],
            "mylar_face_indices": mylar_face_indices,
            "mylar_overlay_faces_mm": [face for item in section_overlays for face in item["faces_mm"]],
            "mylar_source_faces_mm": [face for item in section_overlays for face in item["source_faces_mm"]],
        }
        components.append(component)
        xz_item = {"name": section.name, "kind": component["kind"], "color": component["color"],
                   "polygon": [[x0, cz - height0 / 2], [x1, cz - height1 / 2],
                               [x1, cz + height1 / 2], [x0, cz + height0 / 2]]}
        xy_item = {"name": section.name, "kind": component["kind"], "color": component["color"],
                   "polygon": [[x0, cy - width0 / 2], [x1, cy - width1 / 2],
                               [x1, cy + width1 / 2], [x0, cy + width0 / 2]]}
        xz.append(xz_item); xy.append(xy_item)

    terminal_interfaces: list[dict[str, Any]] = []
    yz: list[dict[str, Any]] = []
    invalid_couplings: list[dict[str, Any]] = []
    # The end view carries the same section inventory as the longitudinal
    # views, using each section's terminal Y/Z profile.  This keeps all three
    # synchronized payloads truthful when an optional guide is appended.
    for section in sections:
        width, height = _assembly_section_dimensions(section, "end")
        cy, cz = section.center_y_mm, section.center_z_mm
        yz.append({
            "name": section.name,
            "kind": "guide" if id(section) in guide_ids else "geometry",
            "color": colors.get(section.material.lower(), "#c5b0d5"),
            "polygon": [[cy - width / 2, cz - height / 2], [cy + width / 2, cz - height / 2],
                        [cy + width / 2, cz + height / 2], [cy - width / 2, cz + height / 2]],
        })
    for face, end in (("minus_x", "start"), ("plus_x", "end")):
        state = model.terminal_state(face)
        coupling_status = model.terminal_coupling_status(face)
        is_readout_face = bool(readout["valid"] and readout["face"] == face)
        dimensions = model.terminal_dimensions_mm(face)
        x = low_x if face == "minus_x" else high_x
        width, height = dimensions["width_y_mm"], dimensions["thickness_z_mm"]
        cy, cz = dimensions["center_y_mm"], dimensions["center_z_mm"]
        terminal = {
            "face": face, "axis_face": "YZ", "state": state,
            "label": coupling_status["label"], "valid": coupling_status["valid"],
            "x_mm": x,
            "width_y_mm": width, "thickness_z_mm": height,
            "center_y_mm": cy, "center_z_mm": cz,
            "mylar_excluded": bool(state != "open" or model.mylar_enabled),
            "interface_name": f"{state}_{face}_interface",
        }
        if state == "detector" and not coupling_status["valid"]:
            invalid_couplings.append({"face": face, "label": coupling_status["label"], "reason": "readout_not_at_guide_output"})
        terminal_interfaces.append(terminal)
        yz.append({"name": f"{face}_terminal", "kind": "terminal", "state": state,
                   "label": terminal["label"], "color": "#1b9e77" if face == "plus_x" else "#984ea3",
                   "linestyle": "-" if face == "plus_x" else "--",
                   "polygon": [[cy - width / 2, cz - height / 2], [cy + width / 2, cz - height / 2],
                               [cy + width / 2, cz + height / 2], [cy - width / 2, cz + height / 2]],
                   "dimensions_mm": [width, height]})
        if state == "mirror":
            marker = {"name": f"terminal_mirror_{face}", "kind": "mirror", "face": face,
                      "color": "#d62728", "label": f"Mirror • {face} YZ",
                      "x_mm": x, "polygon": yz[-1]["polygon"]}
            components.append(marker); xz.append({**marker, "polygon": [[x, cz - height / 2], [x, cz + height / 2]]})
            xy.append({**marker, "polygon": [[x, cy - width / 2], [x, cy + width / 2]]})
        elif state == "detector":
            # The purple interface and the sensor must share the same
            # physical plane.  For a guided readout that plane is outside the
            # complete guide, not the active section's terminal cap.
            frame = canonical_readout_frame_mm(model, face)
            target = str(frame["target"])
            readout_x = float(frame["x_mm"]) if is_readout_face else x
            readout_width = float(frame["width_y_mm"]) if is_readout_face else float(width)
            readout_height = float(frame["thickness_z_mm"]) if is_readout_face else float(height)
            readout_cy = float(frame["center_y_mm"]) if is_readout_face else float(cy)
            readout_cz = float(frame["center_z_mm"]) if is_readout_face else float(cz)
            readout_polygon = [
                [readout_cy - readout_width / 2, readout_cz - readout_height / 2],
                [readout_cy + readout_width / 2, readout_cz - readout_height / 2],
                [readout_cy + readout_width / 2, readout_cz + readout_height / 2],
                [readout_cy - readout_width / 2, readout_cz + readout_height / 2],
            ]
            marker = {"name": f"readout_{face}", "kind": "readout_interface", "face": face,
                      "color": "#6a1b9a", "label": f"{coupling_status['label']} • {face} YZ",
                      "x_mm": readout_x, "polygon": readout_polygon, "coupled_target": target,
                      "target_section": frame["target_section"], "guide_output": bool(frame["guide_output"]),
                      "valid": coupling_status["valid"], "readout_valid": bool(readout["valid"]),
                      "status": readout["code"] if not readout["valid"] else "ok"}
            components.append(marker); yz.append(marker)
            direction = 1.0 if face == "plus_x" else -1.0
            width = float(sensor_width_mm if sensor_width_mm is not None else dimensions["width_y_mm"])
            height = float(sensor_height_mm if sensor_height_mm is not None else dimensions["thickness_z_mm"])
            normalized_sensor = resolve_sensor_readout_dimensions_mm(
                sensor_shape, width_mm=width, height_mm=height, radius_mm=sensor_radius_mm,
                semi_axis_a_mm=sensor_semi_axis_a_mm, semi_axis_b_mm=sensor_semi_axis_b_mm,
            )
            width = float(normalized_sensor["width_mm"])
            height = float(normalized_sensor["height_mm"])
            depth = max(1.5, min(4.0, min(width, height) * 0.20))
            transform = canonical_sensor_transform_mm(
                model, face=face, width_mm=width, height_mm=height,
                depth_mm=depth, shape=sensor_shape,
                radius_mm=float(normalized_sensor["radius_mm"]),
                semi_axis_a_mm=float(normalized_sensor["semi_axis_a_mm"]),
                semi_axis_b_mm=float(normalized_sensor["semi_axis_b_mm"]),
            )
            if transform is None:
                continue
            sensor_solid = transform["solid"]
            coupling_position = transform["coupling_position_mm"]
            body_center = transform["position_mm"]
            if sensor_solid["shape"] == "circular":
                dimension_label = f"circular r={sensor_solid['readout_radius_mm']:g} mm"
            elif sensor_solid["shape"] == "elliptical":
                axes = sensor_solid["readout_semiaxes_mm"]
                dimension_label = f"elliptical a={axes['a']:g}, b={axes['b']:g} mm"
            else:
                dims = sensor_solid["dimensions_mm"]
                dimension_label = f"rectangular {dims['width']:g} × {dims['height']:g} mm"
            x_bounds = (transform["bounds_mm"][0], transform["bounds_mm"][3])
            y_bounds = (transform["bounds_mm"][1], transform["bounds_mm"][4])
            z_bounds = (transform["bounds_mm"][2], transform["bounds_mm"][5])
            sensor = {"name": f"sensor_{face}", "kind": "sensor", "face": face,
                      "color": "#9b59b6", "label": f"{sensor_label} • {face} ({dimension_label}; {transform['target']})",
                      "position_mm": list(body_center), "coupling_position_mm": coupling_position,
                      "coupling": route,
                      "routing": route, "derived_component": derived_component,
                      "automatic_position": bool(sensor_automatic_position),
                      "shape": sensor_solid["shape"], "dimensions_mm": sensor_solid["dimensions_mm"],
                      "readout_radius_mm": sensor_solid["readout_radius_mm"],
                      "readout_semiaxes_mm": sensor_solid["readout_semiaxes_mm"],
                      "normal_axis": "X", "normal_sign": direction,
                      "normal_vector_mm": list(transform.get("normal_vector_mm", [direction, 0.0, 0.0])),
                      "display_offset_mm": transform["display_offset_mm"], "terminal_x_mm": float(coupling_position[0]),
                      "coupled_terminal": transform["face"], "coupled_target": transform["target"],
                      "readout_valid": bool(transform.get("readout_valid", readout["valid"])),
                      "readout_status": transform.get("readout_status", readout["code"]),
                      "target_section": transform["target_section"],
                      "solid_faces_mm": sensor_solid["faces_mm"], "faces_mm": sensor_solid["faces_mm"],
                      "coupling_face_mm": sensor_solid["coupling_face_mm"],
                      "window_faces_mm": sensor_solid["coupling_face_mm"],
                      "bounds_mm": [x_bounds[0], y_bounds[0], z_bounds[0], x_bounds[1], y_bounds[1], z_bounds[1]],
                      "transform": {"position_mm": list(body_center), "coupling_position_mm": coupling_position,
                                    "normal_axis": "X", "normal_sign": direction,
                                    "normal_vector_mm": list(transform.get("normal_vector_mm", [direction, 0.0, 0.0])),
                                    "coupled_terminal": transform["face"], "coupled_target": transform["target"],
                                    "target_section": transform["target_section"],
                                    "display_offset_mm": transform["display_offset_mm"]}}
            components.append(sensor)
            sensor_x = body_center[0]
            xz.append({"name": sensor["name"], "kind": "sensor", "color": sensor["color"],
                       "label": sensor["label"], "solid_component": sensor["name"],
                       "polygon": [[x_bounds[0], z_bounds[0]], [x_bounds[1], z_bounds[0]],
                                   [x_bounds[1], z_bounds[1]], [x_bounds[0], z_bounds[1]]]})
            xy.append({"name": sensor["name"], "kind": "sensor", "color": sensor["color"],
                       "label": sensor["label"], "solid_component": sensor["name"],
                       "polygon": [[x_bounds[0], y_bounds[0]], [x_bounds[1], y_bounds[0]],
                                   [x_bounds[1], y_bounds[1]], [x_bounds[0], y_bounds[1]]]})
            window = sensor_solid["coupling_face_mm"][0]
            yz.append({"name": sensor["name"], "kind": "sensor", "color": "#e1bee7",
                       "label": f"{sensor_label} window • {face}", "solid_component": sensor["name"],
                       "polygon": [[point[1], point[2]] for point in window], "window": True})

    # Orthographic overlays are explicit payload items, not an instruction to
    # a renderer to infer wrapping from a solid's silhouette.  Each polygon is
    # projected from the exact canonical face copied above.
    for overlay in overlay_components:
        for axes, target in (((0, 2), xz), ((0, 1), xy), ((1, 2), yz)):
            for face in overlay["faces_mm"]:
                target.append({
                    "name": overlay["name"], "kind": "mylar_overlay",
                    "label": overlay["label"], "color": MYLAR_DISPLAY_COLOR,
                    "edge_color": MYLAR_DISPLAY_EDGE, "hatch": MYLAR_DISPLAY_HATCH,
                    "stripe_pattern": "diagonal-gold", "reflectivity": overlay["reflectivity"],
                    "source_component": overlay["section_name"],
                    "source_face_index": overlay["source_face_index"],
                    "display_only": True,
                    "polygon": [[float(point[axes[0]]), float(point[axes[1]])] for point in face],
                })

    if not readout["valid"]:
        invalid_couplings.append({"face": readout.get("face"), "label": readout["message"], "reason": readout["code"]})
    states = {item["face"]: item["state"] for item in terminal_interfaces}
    excluded = [item["face"] for item in terminal_interfaces if item["state"] in {"mirror", "detector"}]
    legend_items = [
        ("open", "Open", "#455a64"), ("mirror", "Mirror", "#d62728"),
        ("detector", "Detector/Photosensor coupling", "#6a1b9a"),
    ]
    if model.mylar_enabled:
        legend_items.append(("mylar", f"Mylar overlay ({model.mylar_reflectivity * 100.0:g}%)", MYLAR_DISPLAY_COLOR))
    legend = [{"state": key, "label": value, "color": color} for key, value, color in legend_items]
    sensor_defaults_used = sensor_width_mm is None or sensor_height_mm is None
    status_message = readout["message"]
    if sensor_defaults_used and any(item.get("kind") == "sensor" for item in components):
        status_message += " Configured sensor defaults are being used."
    axis_payloads = {
        view_name: _orthographic_axis_payload(view_name, view_components)
        for view_name, view_components in (("xz", xz), ("xy", xy), ("yz", yz))
    }
    views = {
        "xz": {"title": "Longitudinal XZ", "x_label": "X (mm)", "y_label": "Z (mm)",
                "components": xz, "terminals": terminal_interfaces, "axes": axis_payloads["xz"]},
        "xy": {"title": "Transverse XY", "x_label": "X (mm)", "y_label": "Y (mm)",
                "components": xy, "terminals": terminal_interfaces, "axes": axis_payloads["xy"]},
        "yz": {"title": "Terminal YZ", "x_label": "Y (mm)", "y_label": "Z (mm)",
                "components": yz, "terminals": terminal_interfaces, "axes": axis_payloads["yz"]},
    }
    return {
        "model": "CompleteAssemblyScene", "assembly": model.name, "version": ASSEMBLY_GEOMETRY_VERSION,
        "impact": {"x_mm": float(model.impact_x_mm), "y_mm": float(model.impact_y_mm), "z": "derived"},
        "terminal_states": states, "terminal_interfaces": terminal_interfaces,
        "terminal_coupling_status": {face: model.terminal_coupling_status(face) for face in ("minus_x", "plus_x")},
        "invalid_couplings": invalid_couplings,
        "sensor_routing": {"mode": route, "derived_component": derived_component,
                           "automatic_position": bool(sensor_automatic_position),
                           "readout_terminal": readout.get("face"),
                           "guide_configuration": "valid" if model.has_light_guide else "absent",
                           "status": status_message,
                           "sensor_count": sum(1 for item in components if item.get("kind") == "sensor"),
                           "sensor_defaults": {"shape": sensor_solid["shape"] if any(item.get("kind") == "sensor" for item in components) else str(sensor_shape),
                                               "width_mm": 20.0 if sensor_width_mm is None else float(sensor_width_mm),
                                               "height_mm": 10.0 if sensor_height_mm is None else float(sensor_height_mm),
                                               "radius_mm": sensor_radius_mm, "semi_axis_a_mm": sensor_semi_axis_a_mm,
                                               "semi_axis_b_mm": sensor_semi_axis_b_mm},
                           "sensor_defaults_used": sensor_defaults_used,
                           "ambiguous": readout["code"] == "multiple_readout_terminals"},
        "optical_path": list(readout["optical_path"]),
        "components": components + overlay_components, "physical_components": components,
        "overlay_components": overlay_components,
        "mylar": {"enabled": bool(model.mylar_enabled), "reflectivity": float(model.mylar_reflectivity),
                   "color": MYLAR_DISPLAY_COLOR, "edge_color": MYLAR_DISPLAY_EDGE,
                   "hatch": MYLAR_DISPLAY_HATCH, "display_epsilon_mm": MYLAR_DISPLAY_EPSILON_MM,
                   "overlay_surface_count": len(overlay_components),
                   "excluded_terminal_faces": excluded, "external_only": True,
                   "physical_transport_unchanged": True},
        "views": views, "orthographic_views": views, "legend": legend,
        "summary": f"{model.name}: {len(sections)} geometry section(s); -X {TERMINAL_SURFACE_LABELS[states['minus_x']]}; +X {TERMINAL_SURFACE_LABELS[states['plus_x']]}; Mylar external-only",
        "guide_enabled": bool(model.guide_enabled), "guide_present": bool(model.has_light_guide),
        "guide_section_count": len(model.valid_guide_sections),
        "state_signature": (states["minus_x"], states["plus_x"], bool(model.has_light_guide), bool(model.mylar_enabled),
                             sum(1 for item in components if item.get("kind") == "sensor")),
    }


@dataclass(frozen=True)
class RigidAssemblyTransform:
    """One local-to-global rigid transform for an entire detector assembly.

    The convention is explicit and shared by every renderer: subtract the
    assembly rotation origin, apply the same right-handed rotation matrix to
    every point, then add the rotation origin and the signed translation.
    The default rotation origin is the assembly local/global origin.
    """

    translation_x_mm: float = 0.0
    translation_y_mm: float = 0.0
    translation_z_mm: float = 0.0
    rotation_axis: str = "None"
    rotation_angle_deg: float = 0.0
    rotation_origin_x_mm: float = 0.0
    rotation_origin_y_mm: float = 0.0
    rotation_origin_z_mm: float = 0.0

    def validate(self) -> None:
        if str(self.rotation_axis).strip().title() not in {"None", "X", "Y", "Z"}:
            raise ValueError("rotation_axis must be None, X, Y, or Z")
        for value in (
            self.translation_x_mm, self.translation_y_mm, self.translation_z_mm,
            self.rotation_angle_deg, self.rotation_origin_x_mm,
            self.rotation_origin_y_mm, self.rotation_origin_z_mm,
        ):
            if not math.isfinite(float(value)):
                raise ValueError("rigid assembly transform values must be finite")

    @property
    def normalized_axis(self) -> str:
        axis = str(self.rotation_axis).strip().title()
        return "None" if axis in {"", "None"} else axis

    @property
    def rotation_origin_mm(self) -> tuple[float, float, float]:
        return (
            float(self.rotation_origin_x_mm), float(self.rotation_origin_y_mm),
            float(self.rotation_origin_z_mm),
        )

    def rotation_matrix(self) -> np.ndarray:
        """Return the single right-handed matrix used for all assembly points."""

        self.validate()
        angle = math.radians(float(self.rotation_angle_deg))
        c, s = math.cos(angle), math.sin(angle)
        axis = self.normalized_axis
        if axis == "X":
            return np.asarray(((1.0, 0.0, 0.0), (0.0, c, -s), (0.0, s, c)), dtype=float)
        if axis == "Y":
            return np.asarray(((c, 0.0, s), (0.0, 1.0, 0.0), (-s, 0.0, c)), dtype=float)
        if axis == "Z":
            return np.asarray(((c, -s, 0.0), (s, c, 0.0), (0.0, 0.0, 1.0)), dtype=float)
        return np.eye(3, dtype=float)

    def apply_point_mm(self, point: Sequence[float]) -> tuple[float, float, float]:
        """Map one local millimetre point into global coordinates."""

        vector = np.asarray(tuple(point), dtype=float)
        if vector.shape != (3,) or not np.all(np.isfinite(vector)):
            raise ValueError("rigid transform points must be finite 3-vectors")
        origin = np.asarray(self.rotation_origin_mm, dtype=float)
        translation = np.asarray((self.translation_x_mm, self.translation_y_mm, self.translation_z_mm), dtype=float)
        result = origin + self.rotation_matrix().dot(vector - origin) + translation
        return tuple(float(value) for value in result)

    def apply_vector_mm(self, vector: Sequence[float]) -> tuple[float, float, float]:
        """Rotate a local direction without applying assembly translation."""

        values = np.asarray(tuple(vector), dtype=float)
        if values.shape != (3,) or not np.all(np.isfinite(values)):
            raise ValueError("rigid transform vectors must be finite 3-vectors")
        result = self.rotation_matrix().dot(values)
        return tuple(float(value) for value in result)


_RIGID_FACE_KEYS = (
    "faces_mm", "solid_faces_mm", "window_faces_mm", "coupling_face_mm",
    "interface_faces_mm", "readout_faces_mm", "mylar_faces_mm",
    "mylar_overlay_faces_mm", "mylar_source_faces_mm", "source_faces_mm",
)
_RIGID_POINT_KEYS = (
    "position_mm", "coupling_position_mm", "center_mm", "local_position_mm",
    "global_position_mm", "impact_marker_mm",
)


def _transform_face_list(value: Any, transform: RigidAssemblyTransform) -> Any:
    if not isinstance(value, (list, tuple)):
        return value
    # A face is a sequence of 3-vectors.  Leave malformed metadata untouched;
    # validation of physical faces remains the responsibility of the caller.
    if value and isinstance(value[0], (list, tuple)) and len(value[0]) == 3 and all(
        isinstance(number, (int, float, np.number)) for number in value[0]
    ):
        return [list(transform.apply_point_mm(point)) for point in value]
    return [_transform_face_list(item, transform) for item in value]


def transform_assembly_components(
    components: Sequence[Mapping[str, Any]], transform: RigidAssemblyTransform,
) -> list[dict[str, Any]]:
    """Apply one rigid transform to every physical component payload.

    This is the canonical global operation.  It intentionally handles the
    sensor's body/window/coupling faces, section faces, mirror/readout
    interfaces, and any future rendered Mylar face payload together.  No
    component receives a legacy or second transform.
    """

    transform.validate()
    result: list[dict[str, Any]] = []
    for source in components:
        item = deepcopy(dict(source))
        for key in _RIGID_FACE_KEYS:
            if key in item:
                item[key] = _transform_face_list(item[key], transform)
        for key in _RIGID_POINT_KEYS:
            if key in item and isinstance(item[key], (list, tuple)) and len(item[key]) == 3:
                item[key] = list(transform.apply_point_mm(item[key]))
        normal = item.get("normal_vector_mm")
        if isinstance(normal, (list, tuple)) and len(normal) == 3:
            item["normal_vector_mm"] = list(transform.apply_vector_mm(normal))
        nested = item.get("transform")
        if isinstance(nested, Mapping):
            nested_copy = deepcopy(dict(nested))
            for key in _RIGID_POINT_KEYS:
                if key in nested_copy and isinstance(nested_copy[key], (list, tuple)) and len(nested_copy[key]) == 3:
                    nested_copy[key] = list(transform.apply_point_mm(nested_copy[key]))
            if isinstance(nested_copy.get("normal_vector_mm"), (list, tuple)) and len(nested_copy["normal_vector_mm"]) == 3:
                nested_copy["normal_vector_mm"] = list(transform.apply_vector_mm(nested_copy["normal_vector_mm"]))
            item["transform"] = nested_copy
        # A rotated solid's old axis-aligned bounds are invalid.  Recompute
        # them from the already transformed body faces when available.
        points = [point for face in item.get("faces_mm", ()) for point in face]
        if points and all(len(point) == 3 for point in points):
            item["bounds_mm"] = [
                min(float(point[index]) for point in points) for index in range(3)
            ] + [max(float(point[index]) for point in points) for index in range(3)]
        result.append(item)
    return result


def transform_assembly_terminal_interfaces(
    terminals: Sequence[Mapping[str, Any]], transform: RigidAssemblyTransform,
) -> list[dict[str, Any]]:
    """Transform terminal metadata while retaining its local position."""

    transform.validate()
    result: list[dict[str, Any]] = []
    for source in terminals:
        item = deepcopy(dict(source))
        local = [
            float(item.get("x_mm", 0.0)), float(item.get("center_y_mm", 0.0)),
            float(item.get("center_z_mm", 0.0)),
        ]
        global_position = transform.apply_point_mm(local)
        item["local_position_mm"] = local
        item["global_position_mm"] = list(global_position)
        # Keep legacy fields useful to consumers that read terminal metadata
        # directly; they now describe the same global terminal point.
        item["x_mm"], item["center_y_mm"], item["center_z_mm"] = global_position
        item["global_x_mm"] = global_position[0]
        result.append(item)
    return result


def load_geometry(path: str | Path) -> GeometryModel:
    return GeometryModel.from_json(path)


__all__ = [
    "SECTION_SHAPES",
    "fused_silica_index",
    "MaterialSpec",
    "ActiveDetectorArea",
    "GeometrySection",
    "GeometryModel",
    "BarGeometry",
    "load_geometry",
    "ASSEMBLY_SECTION_SHAPES",
    "ASSEMBLY_AXES",
    "TERMINAL_SURFACE_LABELS",
    "TERMINAL_SURFACE_CHOICES",
    "ASSEMBLY_GEOMETRY_FORMAT",
    "ASSEMBLY_GEOMETRY_VERSION",
    "ASSEMBLY_GEOMETRY_TYPES",
    "MYLAR_DISPLAY_EPSILON_MM",
    "MYLAR_DISPLAY_COLOR",
    "MYLAR_DISPLAY_EDGE",
    "MYLAR_DISPLAY_HATCH",
    "normalize_terminal_surface",
    "terminal_surface_label",
    "terminal_surfaces_from_payload",
    "AssemblySection",
    "AssemblyGeometryModel",
    "assembly_geometry_to_dict",
    "assembly_geometry_from_dict",
    "sensor_solid_geometry_mm",
    "canonical_readout_frame_mm",
    "canonical_sensor_transform_mm",
    "canonical_section_faces_mm",
    "derive_mylar_overlay_surfaces",
    "build_complete_assembly_scene",
    "RigidAssemblyTransform",
    "transform_assembly_components",
    "transform_assembly_terminal_interfaces",
]
