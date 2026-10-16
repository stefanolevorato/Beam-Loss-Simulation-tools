#!/usr/bin/env python3
"""Monte Carlo Cherenkov transport and beam-monitor optical models.

The v31 Mylar-overlay update keeps this transport model local: the GUI owns
global assembly placement and this module continues to consume local physical
XYZ bounds and terminal-state contracts.  That prevents a global refresh from
creating a second sensor or changing local optical physics.

The backward-compatible default is a 0.10 m x 0.02 m x 0.01 m fused-silica
bar in x, y, z. ``SolidGeometry`` additionally transports tapered rectangular
prisms, finite cylinders, elliptical cylinders, and truncated cones.
``CompoundGeometry`` transports validated contiguous sections with local
Fresnel/Snell interfaces. ScintillatorLocalGeometry and
ScintillatorGlobalGeometry keep assembly-local solids separate from the
quartz-relative transform. Public model quantities use SI units.
"""

from __future__ import annotations

import argparse
import json
import math
import os
import tempfile
import time
from collections import Counter
from dataclasses import asdict, dataclass, field, replace
from pathlib import Path
from typing import Callable, Dict, List, Mapping, Optional, Sequence, Tuple, Union

import numpy as np


C_LIGHT = 299_792_458.0
ELEMENTARY_CHARGE = 1.602_176_634e-19
MEV = 1.0e6 * ELEMENTARY_CHARGE
NAN = float("nan")


class _SimulationStopped(Exception):
    """Internal cooperative-cancellation signal for a photon transport loop."""


@dataclass(frozen=True)
class BarGeometry:
    """Axis-aligned rectangular bar centered on the origin.

    x is the 100 mm bar length, y is the 20 mm width, and z is the 10 mm
    particle-crossing thickness. The particle entry face is z=+half_z and
    the detector face is x=+half_x, whose dimensions are 20 mm by 10 mm in
    y-z.
    """

    length_x_m: float = 0.100
    width_y_m: float = 0.020
    thickness_z_m: float = 0.010

    @property
    def half_x(self) -> float:
        return self.length_x_m / 2.0

    @property
    def half_y(self) -> float:
        return self.width_y_m / 2.0

    @property
    def half_z(self) -> float:
        return self.thickness_z_m / 2.0

    def contains(self, point: np.ndarray, tolerance: float = 1e-12) -> bool:
        x, y, z = point
        return (
            -self.half_x - tolerance <= x <= self.half_x + tolerance
            and -self.half_y - tolerance <= y <= self.half_y + tolerance
            and -self.half_z - tolerance <= z <= self.half_z + tolerance
        )


MATERIAL_ALIASES = {
    "fused silica": "fused_silica",
    "fused_silica": "fused_silica",
    "quartz": "fused_silica",
    "air": "air",
    "scintillator": "scintillator",
    "plastic scintillator": "scintillator",
    "light guide": "light_guide",
    "light_guide": "light_guide",
    "custom": "custom",
}
SUPPORTED_SOLID_SHAPES = {
    "box",
    "tapered_rectangular_prism",
    "cylinder",
    "elliptical_cylinder",
    "truncated_cone",
}


DETECTOR_RUN_MODES = ("cherenkov", "scintillator", "combined")

# Path coordinates are deliberately typed at the record boundary.  Physics
# transport can remain in an assembly-local frame while the GUI applies the
# one canonical local-to-global transform for rendering.
COORDINATE_FRAMES = {
    "geometry_m",
    "cherenkov_local_m",
    "scintillator_local_m",
    "global_m",
}


def normalize_coordinate_frame(value: object, default: str = "geometry_m") -> str:
    frame = str(value or default).strip().lower()
    aliases = {
        "local": "geometry_m",
        "assembly_local": "scintillator_local_m",
        "assembly_local_m": "scintillator_local_m",
        "global": "global_m",
        "global_mm": "global_m",
    }
    frame = aliases.get(frame, frame)
    if frame not in COORDINATE_FRAMES:
        raise ValueError(f"unsupported coordinate frame: {value}")
    return frame


def normalize_detector_run_mode(value: object, *, scintillator_enabled: bool = False) -> str:
    """Normalize the GUI/event-view detector participation mode.

    ``auto`` is retained as a safe compatibility alias for the historical
    default, but it is always normalized to Cherenkov-only.  Participation is
    selected by the normalized mode; the legacy checkbox argument is ignored
    so stale profile state cannot turn a run into a combined event.
    """

    mode = str(value or "auto").strip().lower().replace("-", "_")
    aliases = {
        "cherenkov_only": "cherenkov", "scintillator_only": "scintillator",
        "both": "combined", "two_assemblies": "combined",
        "cherenkov only": "cherenkov", "scintillator only": "scintillator",
        "all": "combined",
    }
    mode = aliases.get(mode, mode)
    if mode == "auto":
        return "cherenkov"
    if mode not in DETECTOR_RUN_MODES:
        raise ValueError(f"unsupported detector run mode: {value}")
    return mode


def _terminal_mode(value: object, default: str = "open") -> str:
    text = str(value if value is not None else default).strip().lower()
    aliases = {
        "open": "open", "mirror": "mirror", "detector": "detector",
        "detector/photosensor coupling": "detector",
        # v24 compatibility: the old guide terminal meant detector readout
        # through a configured guide; guide presence is now derived elsewhere.
        "guide": "detector", "optical/light-guide coupling": "detector", "light_guide": "detector",
    }
    mode = aliases.get(text, text)
    if mode not in {"open", "mirror", "detector"}:
        raise ValueError(f"unsupported terminal surface state: {value}")
    return mode


@dataclass(frozen=True)
class SolidGeometry:
    """A single convex editor solid used by the generalized optical transport.

    All lengths are SI metres.  The solid is centred at ``center_*_m``.
    ``tapered_rectangular_prism`` varies its Y/Z dimensions from the -X cap
    (start) to the +X cap (end).  Cylinder, elliptical-cylinder and
    truncated-cone length/radius parameters are measured along ``axis``.

    A ``SolidGeometry`` is one section of the transport model. Multiple
    sections are represented explicitly by ``CompoundGeometry`` below;
    callers should not silently flatten them to this class.
    """

    shape: str = "box"
    material: str = "fused_silica"
    center_x_m: float = 0.0
    center_y_m: float = 0.0
    center_z_m: float = 0.0
    length_x_m: float = 0.100
    width_y_m: float = 0.020
    thickness_z_m: float = 0.010
    width_y_start_m: float = 0.020
    width_y_end_m: float = 0.020
    thickness_z_start_m: float = 0.010
    thickness_z_end_m: float = 0.010
    axis: str = "Z"
    radius_m: float = 0.010
    length_m: float = 0.100
    semi_axis_a_m: float = 0.010
    semi_axis_b_m: float = 0.005
    radius_start_m: float = 0.010
    radius_end_m: float = 0.005
    # Explicit index is used for the first scintillator/light-guide model.
    # It is optional so all legacy quartz geometry constructors remain valid.
    refractive_index_value: Optional[float] = None
    # Passive canonical assembly metadata.  It is populated by
    # AssemblyGeometryModel and intentionally does not alter dimensions or
    # transport intersection math.
    terminal_interfaces: Tuple[Mapping[str, object], ...] = field(default_factory=tuple)
    # Passive copy of the canonical physical sensor body/window transform.
    # Transport still uses terminal state for photon interaction; this keeps
    # the optical adapter tied to the same readout geometry as the renderers.
    sensor_components: Tuple[Mapping[str, object], ...] = field(default_factory=tuple)

    def normalized_material(self) -> str:
        return MATERIAL_ALIASES.get(str(self.material).strip().lower(), "")

    def validate(self) -> None:
        if self.shape not in SUPPORTED_SOLID_SHAPES:
            raise ValueError(f"unsupported solid shape: {self.shape}")
        if self.normalized_material() not in {"fused_silica", "air", "scintillator", "light_guide", "custom"}:
            raise ValueError(f"unsupported solid material: {self.material}")
        if self.refractive_index_value is not None and (
            not math.isfinite(float(self.refractive_index_value)) or float(self.refractive_index_value) <= 0.0
        ):
            raise ValueError("refractive_index_value must be positive when supplied")
        if self.axis not in {"X", "Y", "Z"}:
            raise ValueError("solid axis must be X, Y, or Z")
        for key in ("center_x_m", "center_y_m", "center_z_m"):
            if not math.isfinite(float(getattr(self, key))):
                raise ValueError(f"{key} must be finite")
        fields = {
            "box": ("length_x_m", "width_y_m", "thickness_z_m"),
            "tapered_rectangular_prism": ("length_x_m", "width_y_start_m", "width_y_end_m", "thickness_z_start_m", "thickness_z_end_m"),
            "cylinder": ("radius_m", "length_m"),
            "elliptical_cylinder": ("semi_axis_a_m", "semi_axis_b_m", "length_m"),
            "truncated_cone": ("radius_start_m", "radius_end_m", "length_m"),
        }[self.shape]
        for key in fields:
            value = float(getattr(self, key))
            if not math.isfinite(value) or value <= 0.0:
                raise ValueError(f"{key} must be positive")

    @property
    def center_m(self) -> np.ndarray:
        return np.array([self.center_x_m, self.center_y_m, self.center_z_m], dtype=float)

    @property
    def bounds_min_m(self) -> np.ndarray:
        return self._bounds()[0]

    @property
    def bounds_max_m(self) -> np.ndarray:
        return self._bounds()[1]

    @property
    def half_x(self) -> float:
        return max(abs(float(self.bounds_min_m[0])), abs(float(self.bounds_max_m[0])))

    @property
    def half_y(self) -> float:
        return max(abs(float(self.bounds_min_m[1])), abs(float(self.bounds_max_m[1])))

    @property
    def half_z(self) -> float:
        return max(abs(float(self.bounds_min_m[2])), abs(float(self.bounds_max_m[2])))

    @property
    def length_x_extent_m(self) -> float:
        return float(self.bounds_max_m[0] - self.bounds_min_m[0])

    @property
    def width_y_extent_m(self) -> float:
        return float(self.bounds_max_m[1] - self.bounds_min_m[1])

    @property
    def thickness_z_extent_m(self) -> float:
        return float(self.bounds_max_m[2] - self.bounds_min_m[2])

    @property
    def characteristic_length_m(self) -> float:
        return max(self.length_x_extent_m, self.width_y_extent_m, self.thickness_z_extent_m)

    @property
    def boundary_tolerance_m(self) -> float:
        """Numerical tolerance used when classifying a solid boundary.

        This is deliberately separate from the transport nudge below.  A
        point may be accepted as being on a surface for validation purposes,
        but a photon must be moved farther than this tolerance after a hit or
        the next intersection query can see the same surface again.
        """

        return max(1e-10, self.characteristic_length_m * 1e-8)

    @property
    def transport_epsilon_m(self) -> float:
        """Distance for an oriented post-boundary transport probe."""

        return max(1e-12, min(1e-7, self.characteristic_length_m * 1e-7))

    @property
    def detector_portals_supported(self) -> bool:
        """Whether +/-X are planar terminal caps usable by end detectors."""

        return self.shape in {"box", "tapered_rectangular_prism"} or self.axis == "X"

    def _bounds(self) -> Tuple[np.ndarray, np.ndarray]:
        self.validate()
        cx, cy, cz = self.center_x_m, self.center_y_m, self.center_z_m
        if self.shape == "box":
            extents = (self.length_x_m / 2.0, self.width_y_m / 2.0, self.thickness_z_m / 2.0)
        elif self.shape == "tapered_rectangular_prism":
            extents = (self.length_x_m / 2.0, max(self.width_y_start_m, self.width_y_end_m) / 2.0, max(self.thickness_z_start_m, self.thickness_z_end_m) / 2.0)
        else:
            if self.shape == "cylinder":
                radius_a = radius_b = self.radius_m
            elif self.shape == "elliptical_cylinder":
                radius_a, radius_b = self.semi_axis_a_m, self.semi_axis_b_m
            else:
                radius_a = radius_b = max(self.radius_start_m, self.radius_end_m)
            half_length = self.length_m / 2.0
            extents = [radius_a, radius_b, radius_a]
            axis_index = {"X": 0, "Y": 1, "Z": 2}[self.axis]
            extents[axis_index] = half_length
            if self.axis == "X":
                extents = (half_length, radius_a, radius_b)
            elif self.axis == "Y":
                extents = (radius_a, half_length, radius_b)
            else:
                extents = (radius_a, radius_b, half_length)
        return (np.array([cx - extents[0], cy - extents[1], cz - extents[2]], dtype=float), np.array([cx + extents[0], cy + extents[1], cz + extents[2]], dtype=float))

    def refractive_index(self, wavelength_m: float | np.ndarray) -> float | np.ndarray:
        material = self.normalized_material()
        if self.refractive_index_value is not None:
            value = np.ones_like(np.asarray(wavelength_m, dtype=float)) * float(self.refractive_index_value)
            return float(value) if value.ndim == 0 else value
        if material == "air":
            value = np.ones_like(np.asarray(wavelength_m, dtype=float))
            return float(value) if value.ndim == 0 else value
        if material == "scintillator":
            # The canonical editor schema identifies plastic scintillator by
            # material but does not carry a separate dispersion table.  Keep
            # its documented default optical index explicit instead of
            # silently treating it as fused silica.
            value = np.ones_like(np.asarray(wavelength_m, dtype=float)) * 1.58
            return float(value) if value.ndim == 0 else value
        return fused_silica_index(wavelength_m)

    def contains(self, point: np.ndarray, tolerance: float = 1e-12) -> bool:
        self.validate()
        x, y, z = np.asarray(point, dtype=float) - self.center_m
        if self.shape == "box":
            return abs(x) <= self.length_x_m / 2.0 + tolerance and abs(y) <= self.width_y_m / 2.0 + tolerance and abs(z) <= self.thickness_z_m / 2.0 + tolerance
        if self.shape == "tapered_rectangular_prism":
            half_x = self.length_x_m / 2.0
            if abs(x) > half_x + tolerance:
                return False
            fraction = (x + half_x) / self.length_x_m
            hy = (self.width_y_start_m + fraction * (self.width_y_end_m - self.width_y_start_m)) / 2.0
            hz = (self.thickness_z_start_m + fraction * (self.thickness_z_end_m - self.thickness_z_start_m)) / 2.0
            return abs(y) <= hy + tolerance and abs(z) <= hz + tolerance
        q, u, v = self._axis_coordinates(x, y, z)
        if abs(q) > self.length_m / 2.0 + tolerance:
            return False
        if self.shape == "cylinder":
            return u * u + v * v <= self.radius_m * self.radius_m + tolerance
        if self.shape == "elliptical_cylinder":
            return (u / self.semi_axis_a_m) ** 2 + (v / self.semi_axis_b_m) ** 2 <= 1.0 + tolerance
        radius = self._cone_radius(q)
        return u * u + v * v <= radius * radius + tolerance

    def _axis_coordinates(self, x: float, y: float, z: float) -> Tuple[float, float, float]:
        if self.axis == "X":
            return x, y, z
        if self.axis == "Y":
            return y, x, z
        return z, x, y

    def _cone_radius(self, q: float) -> float:
        fraction = (q + self.length_m / 2.0) / self.length_m
        return self.radius_start_m + fraction * (self.radius_end_m - self.radius_start_m)

    def _candidate_parameters(self, point: np.ndarray, direction: np.ndarray) -> List[Tuple[float, str]]:
        p = np.asarray(point, dtype=float) - self.center_m
        d = np.asarray(direction, dtype=float)
        candidates: List[Tuple[float, str]] = []

        def add_plane(coordinate: float, component: float, boundary: float, name: str) -> None:
            if abs(component) > 1e-15:
                candidates.append(((boundary - coordinate) / component, name))

        if self.shape == "box":
            add_plane(p[0], d[0], -self.length_x_m / 2.0, "x-")
            add_plane(p[0], d[0], self.length_x_m / 2.0, "x+")
            add_plane(p[1], d[1], -self.width_y_m / 2.0, "y-")
            add_plane(p[1], d[1], self.width_y_m / 2.0, "y+")
            add_plane(p[2], d[2], -self.thickness_z_m / 2.0, "z-")
            add_plane(p[2], d[2], self.thickness_z_m / 2.0, "z+")
            return candidates
        if self.shape == "tapered_rectangular_prism":
            half_x = self.length_x_m / 2.0
            add_plane(p[0], d[0], -half_x, "x-")
            add_plane(p[0], d[0], half_x, "x+")
            sy = (self.width_y_end_m - self.width_y_start_m) / self.length_x_m / 2.0
            sz = (self.thickness_z_end_m - self.thickness_z_start_m) / self.length_x_m / 2.0
            hy0, hz0 = self.width_y_start_m / 2.0, self.thickness_z_start_m / 2.0
            for value, component, slope, offset, name in (
                (p[1], d[1], sy, hy0, "y+"), (-p[1], -d[1], sy, hy0, "y-"),
                (p[2], d[2], sz, hz0, "z+"), (-p[2], -d[2], sz, hz0, "z-"),
            ):
                coefficient = component - slope * d[0]
                constant = value - offset - slope * (p[0] + half_x)
                if abs(coefficient) > 1e-15:
                    candidates.append((-constant / coefficient, name))
            return candidates

        q, u, v = self._axis_coordinates(*p)
        dq, du, dv = self._axis_coordinates(*d)
        half_length = self.length_m / 2.0
        add_plane(q, dq, -half_length, "axis-")
        add_plane(q, dq, half_length, "axis+")
        if self.shape == "cylinder":
            a, b, slope = self.radius_m, self.radius_m, 0.0
        elif self.shape == "elliptical_cylinder":
            a, b, slope = self.semi_axis_a_m, self.semi_axis_b_m, 0.0
        else:
            a, b, slope = self.radius_start_m, self.radius_start_m, (self.radius_end_m - self.radius_start_m) / self.length_m
        if self.shape == "truncated_cone":
            radius_at_q = self._cone_radius(q)
            aa = du * du + dv * dv - (slope * dq) ** 2
            bb = 2.0 * (u * du + v * dv - radius_at_q * slope * dq)
            cc = u * u + v * v - radius_at_q * radius_at_q
        else:
            aa = (du / a) ** 2 + (dv / b) ** 2
            bb = 2.0 * (u * du / (a * a) + v * dv / (b * b))
            cc = (u / a) ** 2 + (v / b) ** 2 - 1.0
        if abs(aa) < 1e-18:
            if abs(bb) > 1e-15:
                candidates.append((-cc / bb, "side"))
        else:
            discriminant = bb * bb - 4.0 * aa * cc
            if discriminant >= -1e-14:
                root = math.sqrt(max(0.0, discriminant))
                candidates.extend([((-bb - root) / (2.0 * aa), "side"), ((-bb + root) / (2.0 * aa), "side")])
        return candidates

    def surface_normal(self, point: np.ndarray, surface: Optional[str] = None) -> np.ndarray:
        p = np.asarray(point, dtype=float) - self.center_m
        if surface == "x+": return np.array([1.0, 0.0, 0.0])
        if surface == "x-": return np.array([-1.0, 0.0, 0.0])
        if surface in {"axis+", "axis-"}:
            vector = np.zeros(3); vector[{"X": 0, "Y": 1, "Z": 2}[self.axis]] = 1.0 if surface == "axis+" else -1.0
            return vector
        if self.shape == "tapered_rectangular_prism":
            x, y, z = p; half_x = self.length_x_m / 2.0
            sy = (self.width_y_end_m - self.width_y_start_m) / self.length_x_m / 2.0
            sz = (self.thickness_z_end_m - self.thickness_z_start_m) / self.length_x_m / 2.0
            hy = (self.width_y_start_m + (x + half_x) * 2.0 * sy) / 2.0
            hz = (self.thickness_z_start_m + (x + half_x) * 2.0 * sz) / 2.0
            values = [(abs(y - hy), np.array([-sy, 1.0, 0.0])), (abs(y + hy), np.array([-sy, -1.0, 0.0])), (abs(z - hz), np.array([-sz, 0.0, 1.0])), (abs(z + hz), np.array([-sz, 0.0, -1.0]))]
            normal = min(values, key=lambda item: item[0])[1]
            return normal / np.linalg.norm(normal)
        if surface == "y+": return np.array([0.0, 1.0, 0.0])
        if surface == "y-": return np.array([0.0, -1.0, 0.0])
        if surface == "z+": return np.array([0.0, 0.0, 1.0])
        if surface == "z-": return np.array([0.0, 0.0, -1.0])
        q, u, v = self._axis_coordinates(*p)
        if self.shape == "cylinder":
            normal_uv = np.array([u, v])
        elif self.shape == "elliptical_cylinder":
            normal_uv = np.array([u / (self.semi_axis_a_m ** 2), v / (self.semi_axis_b_m ** 2)])
        else:
            radius = self._cone_radius(q)
            slope = (self.radius_end_m - self.radius_start_m) / self.length_m
            normal_uv = np.array([u, v])
        if self.shape == "truncated_cone":
            local = np.array([-slope * radius, u, v], dtype=float)
        else:
            local = np.array([0.0, normal_uv[0], normal_uv[1]], dtype=float)
        if self.axis == "X": normal = local
        elif self.axis == "Y": normal = np.array([local[1], local[0], local[2]])
        else: normal = np.array([local[1], local[2], local[0]])
        return normal / np.linalg.norm(normal)

    def ray_intersection(self, point: np.ndarray, direction: np.ndarray, min_distance: float = 1e-12) -> Tuple[float, str, np.ndarray]:
        """Return the first forward surface hit and its outward local normal."""

        direction = np.asarray(direction, dtype=float)
        direction = direction / np.linalg.norm(direction)
        candidates = []
        tolerance = self.boundary_tolerance_m
        # Do not use the relaxed boundary tolerance to decide whether the ray
        # starts inside.  After a surface reflection/escape the transport
        # point is intentionally nudged in the direction of travel; treating
        # a point just outside as inside is what caused the original
        # ``no forward boundary`` failure on tapered sides.
        start_inside = self.contains(np.asarray(point, dtype=float), tolerance=0.0)
        for distance, surface in self._candidate_parameters(point, direction):
            if distance <= min_distance or not math.isfinite(distance):
                continue
            hit = np.asarray(point, dtype=float) + distance * direction
            if self.contains(hit, tolerance=tolerance):
                if surface == "axis+" and self.axis == "X":
                    surface = "x+"
                elif surface == "axis-" and self.axis == "X":
                    surface = "x-"
                normal = self.surface_normal(hit, surface)
                dot = float(np.dot(direction, normal))
                if (start_inside and dot > -1e-8) or (not start_inside and dot < -1e-8):
                    candidates.append((distance, surface, normal))
        if not candidates:
            raise RuntimeError("solid has no forward boundary intersection")
        return min(candidates, key=lambda item: item[0])

    def ray_interval(self, origin: np.ndarray, direction: np.ndarray) -> Tuple[float, float]:
        """Return entry/exit parameters for a ray through this convex solid."""

        direction = np.asarray(direction, dtype=float)
        direction = direction / np.linalg.norm(direction)
        candidates = sorted({float(t) for t, _ in self._candidate_parameters(origin, direction) if math.isfinite(t)})
        epsilon = max(1e-10, self.characteristic_length_m * 1e-8)
        if self.contains(np.asarray(origin, dtype=float), tolerance=epsilon):
            future = [t for t in candidates if t > epsilon]
            if not future:
                raise ValueError("ray inside solid has no exit")
            return 0.0, future[0]
        entry = None
        exit_value = None
        for value in candidates:
            before = np.asarray(origin, dtype=float) + (value - epsilon) * direction
            after = np.asarray(origin, dtype=float) + (value + epsilon) * direction
            if entry is None and self.contains(after, tolerance=epsilon):
                entry = value
            elif entry is not None and self.contains(before, tolerance=epsilon):
                exit_value = value
                break
        if entry is None or exit_value is None or exit_value <= entry:
            raise ValueError("ray does not intersect the solid")
        return entry, exit_value

    def particle_track(self, particle: "Particle") -> Tuple[np.ndarray, float]:
        reference = np.array([particle.entry_x_m, particle.entry_y_m, self.bounds_max_m[2]], dtype=float)
        margin = max(4.0 * self.characteristic_length_m, 1e-6)
        origin = reference - particle.direction * margin
        entry_parameter, exit_parameter = self.ray_interval(origin, particle.direction)
        entry = origin + entry_parameter * particle.direction
        return entry, exit_parameter - entry_parameter

    def surface_rings(self, samples: int = 48) -> List[np.ndarray]:
        """Return compact surface rings for diagnostics and event plots."""

        cx, cy, cz = self.center_x_m, self.center_y_m, self.center_z_m
        if self.shape == "box":
            hx, hy, hz = self.length_x_m / 2.0, self.width_y_m / 2.0, self.thickness_z_m / 2.0
            return [np.asarray([(cx - hx, cy - hy, cz - hz), (cx + hx, cy - hy, cz - hz), (cx + hx, cy + hy, cz - hz), (cx - hx, cy + hy, cz - hz)]), np.asarray([(cx - hx, cy - hy, cz + hz), (cx + hx, cy - hy, cz + hz), (cx + hx, cy + hy, cz + hz), (cx - hx, cy + hy, cz + hz)])]
        if self.shape == "tapered_rectangular_prism":
            x0, x1 = cx - self.length_x_m / 2.0, cx + self.length_x_m / 2.0
            def ring(x, width, thickness):
                hy, hz = width / 2.0, thickness / 2.0
                return np.asarray([(x, cy - hy, cz - hz), (x, cy + hy, cz - hz), (x, cy + hy, cz + hz), (x, cy - hy, cz + hz)])
            return [ring(x0, self.width_y_start_m, self.thickness_z_start_m), ring(x1, self.width_y_end_m, self.thickness_z_end_m)]
        if self.shape == "cylinder":
            a = b = self.radius_m
        elif self.shape == "elliptical_cylinder":
            a, b = self.semi_axis_a_m, self.semi_axis_b_m
        else:
            a = b = self.radius_start_m
        half_length = self.length_m / 2.0
        def ring_at(q, first, second):
            points = []
            for index in range(samples):
                theta = 2.0 * math.pi * index / samples
                u, v = first * math.cos(theta), second * math.sin(theta)
                if self.axis == "X": points.append((cx + q, cy + u, cz + v))
                elif self.axis == "Y": points.append((cx + u, cy + q, cz + v))
                else: points.append((cx + u, cy + v, cz + q))
            return np.asarray(points)
        if self.shape == "truncated_cone":
            return [ring_at(-half_length, self.radius_start_m, self.radius_start_m), ring_at(half_length, self.radius_end_m, self.radius_end_m)]
        return [ring_at(-half_length, a, b), ring_at(half_length, a, b)]


@dataclass(frozen=True)
class CompoundGeometry:
    """A contiguous chain of editor solids along the global X axis.

    The first supported composite model is intentionally explicit: sections
    are ordered from -X to +X, their X intervals must touch exactly (within a
    numerical tolerance), and adjacent sections must have the same
    Y/Z cross-section at their shared plane.  This permits box-to-box,
    tapered-to-tapered, cylinder/cone/ellipse chains with a well-defined
    interface, while rejecting ledges, overlaps and gaps instead of silently
    treating them as a box.  Every section remains a real ``SolidGeometry``;
    photon transport uses the local section index and material.

    The particle track may cross any contiguous subset of the chain.  Optical
    transport supports internal Fresnel reflection/refraction and external
    +/-X detector portals on the terminal sections.  Optical interfaces on
    transverse ledges are deliberately outside this model and therefore are
    rejected by :meth:`validate`.
    """

    sections: Tuple[SolidGeometry, ...] = ()
    # Scintillator assemblies may deliberately couple different footprints.
    # Legacy CompoundGeometry keeps strict matching by default.
    allow_interface_mismatch: bool = False
    wrapped_external: bool = False
    mylar_reflectivity: float = 1.0
    terminal_interfaces: Tuple[Mapping[str, object], ...] = field(default_factory=tuple)
    sensor_components: Tuple[Mapping[str, object], ...] = field(default_factory=tuple)

    def __post_init__(self) -> None:
        object.__setattr__(self, "sections", tuple(self.sections))

    @staticmethod
    def _cross_section_signature(section: SolidGeometry, positive_x: bool) -> Tuple[str, float, float]:
        """Return the transverse footprint at an X cap as (kind, Y, Z)."""

        if section.shape == "box":
            return ("rectangle", section.width_y_m, section.thickness_z_m)
        if section.shape == "tapered_rectangular_prism":
            width = section.width_y_end_m if positive_x else section.width_y_start_m
            thickness = section.thickness_z_end_m if positive_x else section.thickness_z_start_m
            return ("rectangle", width, thickness)
        if section.shape == "cylinder":
            return ("ellipse", 2.0 * section.radius_m, 2.0 * section.radius_m)
        if section.shape == "elliptical_cylinder":
            return ("ellipse", 2.0 * section.semi_axis_a_m, 2.0 * section.semi_axis_b_m)
        if section.shape == "truncated_cone":
            radius = section.radius_end_m if positive_x else section.radius_start_m
            return ("ellipse", 2.0 * radius, 2.0 * radius)
        raise ValueError(f"unsupported compound section shape: {section.shape}")

    @staticmethod
    def _same_cross_section(left: Tuple[str, float, float], right: Tuple[str, float, float]) -> bool:
        return (
            left[0] == right[0]
            and math.isclose(left[1], right[1], rel_tol=1e-9, abs_tol=1e-12)
            and math.isclose(left[2], right[2], rel_tol=1e-9, abs_tol=1e-12)
        )

    def validate(self) -> None:
        if not self.sections:
            raise ValueError("compound geometry must contain at least one section")
        if not 0.0 <= float(self.mylar_reflectivity) <= 1.0:
            raise ValueError("mylar_reflectivity must be between 0 and 1")
        previous = None
        for index, section in enumerate(self.sections):
            if not isinstance(section, SolidGeometry):
                raise TypeError(f"compound section {index} is not SolidGeometry")
            section.validate()
            if section.shape not in {"box", "tapered_rectangular_prism"} and section.axis != "X":
                raise ValueError(
                    f"compound section {index} ({section.shape}) must use axis X; "
                    "only X-aligned sections have a supported +/-X interface"
                )
            if index and (
                not math.isclose(section.center_y_m, previous.center_y_m, rel_tol=0.0, abs_tol=1e-12)
                or not math.isclose(section.center_z_m, previous.center_z_m, rel_tol=0.0, abs_tol=1e-12)
            ):
                raise ValueError("compound sections must share the same Y/Z center")
            low, high = section.bounds_min_m[0], section.bounds_max_m[0]
            if previous is not None:
                scale = max(abs(low), abs(high), abs(previous.bounds_max_m[0]), 1e-6)
                tolerance = max(1e-12, scale * 1e-9)
                if low < previous.bounds_max_m[0] - tolerance:
                    raise ValueError(f"compound sections {index - 1} and {index} overlap along X")
                if abs(low - previous.bounds_max_m[0]) > tolerance:
                    raise ValueError(f"compound sections {index - 1} and {index} have a gap along X")
                left = self._cross_section_signature(previous, True)
                right = self._cross_section_signature(section, False)
                if not self.allow_interface_mismatch and not self._same_cross_section(left, right):
                    raise ValueError(
                        f"compound interface {index - 1}/{index} has incompatible Y/Z cross-sections"
                    )
            previous = section

    @property
    def bounds_min_m(self) -> np.ndarray:
        self.validate()
        return np.min(np.asarray([section.bounds_min_m for section in self.sections]), axis=0)

    @property
    def bounds_max_m(self) -> np.ndarray:
        self.validate()
        return np.max(np.asarray([section.bounds_max_m for section in self.sections]), axis=0)

    @property
    def half_x(self) -> float:
        return max(abs(float(self.bounds_min_m[0])), abs(float(self.bounds_max_m[0])))

    @property
    def half_y(self) -> float:
        return max(abs(float(self.bounds_min_m[1])), abs(float(self.bounds_max_m[1])))

    @property
    def half_z(self) -> float:
        return max(abs(float(self.bounds_min_m[2])), abs(float(self.bounds_max_m[2])))

    @property
    def length_x_extent_m(self) -> float:
        return float(self.bounds_max_m[0] - self.bounds_min_m[0])

    @property
    def width_y_extent_m(self) -> float:
        return float(self.bounds_max_m[1] - self.bounds_min_m[1])

    @property
    def thickness_z_extent_m(self) -> float:
        return float(self.bounds_max_m[2] - self.bounds_min_m[2])

    @property
    def characteristic_length_m(self) -> float:
        return max(self.length_x_extent_m, self.width_y_extent_m, self.thickness_z_extent_m)

    @property
    def boundary_tolerance_m(self) -> float:
        """Tolerance for shared-cap and boundary classification."""

        return max(1e-10, self.characteristic_length_m * 1e-8)

    @property
    def transport_epsilon_m(self) -> float:
        """Oriented nudge used after an interface or surface interaction."""

        return max(1e-12, min(1e-7, self.characteristic_length_m * 1e-7))

    @property
    def detector_portals_supported(self) -> bool:
        return bool(self.sections)

    def contains(self, point: np.ndarray, tolerance: float = 1e-12) -> bool:
        self.validate()
        return any(section.contains(point, tolerance=tolerance) for section in self.sections)

    def section_index_at(self, point: np.ndarray, direction: Optional[np.ndarray] = None) -> int:
        """Return the local section containing a point, resolving shared caps."""

        self.validate()
        point = np.asarray(point, dtype=float)
        candidates = [index for index, section in enumerate(self.sections) if section.contains(point, tolerance=self.boundary_tolerance_m)]
        if not candidates:
            raise ValueError("point is outside the compound geometry")
        if direction is not None and len(candidates) > 1:
            dx = float(np.asarray(direction, dtype=float)[0])
            if dx > 1e-12:
                return max(candidates)
            if dx < -1e-12:
                return min(candidates)
        return candidates[0]

    def resolve_transport_section(self, point: np.ndarray, direction: np.ndarray, current_index: int) -> Optional[int]:
        """Resolve the section after an oriented post-hit transport nudge.

        The shared X cap belongs to both sections within the numerical
        tolerance.  The forward probe gives the physically correct ownership
        after crossing it and also prevents a reflected/escaped photon from
        being sent back to the same surface indefinitely.  ``None`` means the
        photon has just left the compound through an external boundary.
        """

        self.validate()
        current_index = int(current_index)
        if not 0 <= current_index < len(self.sections):
            raise RuntimeError(
                f"compound transport state incoherent: invalid section={current_index}, "
                f"position={np.asarray(point, dtype=float).tolist()}, direction={np.asarray(direction, dtype=float).tolist()}"
            )
        point = np.asarray(point, dtype=float)
        direction = np.asarray(direction, dtype=float)
        norm = np.linalg.norm(direction)
        if norm <= 0.0 or not np.all(np.isfinite(direction)):
            raise RuntimeError(
                f"compound transport state incoherent: invalid direction={direction.tolist()}, "
                f"section={current_index}, position={point.tolist()}"
            )
        direction = direction / norm
        tolerance = self.boundary_tolerance_m
        probe = point + self.transport_epsilon_m * direction

        strict = [index for index, section in enumerate(self.sections) if section.contains(point, tolerance=0.0)]
        if current_index in strict:
            return current_index

        forward = [index for index, section in enumerate(self.sections) if section.contains(probe, tolerance=0.0)]
        if forward:
            if current_index in forward:
                return current_index
            dx = float(direction[0])
            if dx > 1e-12:
                return max(forward)
            if dx < -1e-12:
                return min(forward)
            return forward[0]

        tolerant = [index for index, section in enumerate(self.sections) if section.contains(point, tolerance=tolerance)]
        if tolerant:
            forward_tolerant = [index for index, section in enumerate(self.sections) if section.contains(probe, tolerance=tolerance)]
            if forward_tolerant:
                if current_index in forward_tolerant:
                    return current_index
                dx = float(direction[0])
                if dx > 1e-12:
                    return max(forward_tolerant)
                if dx < -1e-12:
                    return min(forward_tolerant)
                return forward_tolerant[0]

        # A boundary candidate at (or just behind) the point, with travel
        # directed out of the surface, is a real external escape.  Internal
        # caps are handled by the forward probe above.
        if self.nearby_exit_face(point, direction, current_index) is not None:
            return None

        raise RuntimeError(
            f"compound transport state incoherent: section={current_index}, "
            f"position={point.tolist()}, direction={direction.tolist()}, "
            f"strict_sections={strict}, tolerant_sections={tolerant}"
        )

    def nearby_exit_face(self, point: np.ndarray, direction: np.ndarray, section_index: int) -> Optional[str]:
        """Return a terminal face whose root is at/just behind ``point``."""

        section = self.sections[int(section_index)]
        direction = np.asarray(direction, dtype=float)
        direction = direction / np.linalg.norm(direction)
        limit = max(self.transport_epsilon_m * 2.0, self.boundary_tolerance_m * 2.0)
        candidates = []
        for distance, face in section._candidate_parameters(np.asarray(point, dtype=float), direction):
            if not math.isfinite(distance) or abs(distance) > limit:
                continue
            hit = np.asarray(point, dtype=float) + distance * direction
            normal = section.surface_normal(hit, face)
            if float(np.dot(direction, normal)) > 1e-8:
                adjacent = (
                    section_index + 1 if face == "x+" and section_index + 1 < len(self.sections)
                    else section_index - 1 if face == "x-" and section_index > 0
                    else None
                )
                if adjacent is None:
                    candidates.append((abs(distance), face))
        return min(candidates, key=lambda item: item[0])[1] if candidates else None

    def refractive_index_at(self, point: np.ndarray, wavelength_m: float, section_index: Optional[int] = None) -> float:
        index = self.section_index_at(point) if section_index is None else int(section_index)
        return float(self.sections[index].refractive_index(wavelength_m))

    def material_at(self, point: np.ndarray, section_index: Optional[int] = None) -> str:
        index = self.section_index_at(point) if section_index is None else int(section_index)
        return self.sections[index].normalized_material()

    def next_boundary(self, point: np.ndarray, direction: np.ndarray, section_index: int) -> Tuple[float, str, np.ndarray, bool, Optional[int]]:
        """Return the next boundary and whether it is an internal interface."""

        self.validate()
        section_index = int(section_index)
        section = self.sections[section_index]
        distance, face, normal = section.ray_intersection(point, direction)
        adjacent = None
        if face == "x+" and section_index + 1 < len(self.sections):
            adjacent = section_index + 1
        elif face == "x-" and section_index > 0:
            adjacent = section_index - 1
        return distance, face, normal, adjacent is not None, adjacent

    def ray_intersection(self, point: np.ndarray, direction: np.ndarray, min_distance: float = 1e-12) -> Tuple[float, str, np.ndarray]:
        """Return the first forward boundary of the compound union.

        For an inside point this is the boundary of its local section.  For
        an outside point the nearest section entry is returned.  Transport
        uses ``next_boundary`` when it needs the additional internal/external
        interface classification.
        """

        self.validate()
        if self.contains(np.asarray(point, dtype=float), tolerance=max(1e-10, self.characteristic_length_m * 1e-8)):
            index = self.section_index_at(point, direction)
            return self.next_boundary(point, direction, index)[:3]
        candidates = []
        for section in self.sections:
            try:
                distance, face, normal = section.ray_intersection(point, direction, min_distance=min_distance)
            except RuntimeError:
                continue
            candidates.append((distance, face, normal))
        if not candidates:
            raise RuntimeError("compound geometry has no forward boundary intersection")
        return min(candidates, key=lambda item: item[0])

    def ray_interval(self, origin: np.ndarray, direction: np.ndarray) -> Tuple[float, float]:
        intervals = self._particle_intervals(origin, direction)
        if not intervals:
            raise ValueError("ray does not intersect the compound geometry")
        return intervals[0][0], intervals[-1][1]

    def _particle_intervals(self, origin: np.ndarray, direction: np.ndarray) -> List[Tuple[float, float, int]]:
        self.validate()
        intervals = []
        for index, section in enumerate(self.sections):
            try:
                start, end = section.ray_interval(origin, direction)
            except ValueError:
                continue
            intervals.append((float(start), float(end), index))
        intervals.sort(key=lambda item: item[0])
        if intervals:
            tolerance = max(1e-10, self.characteristic_length_m * 1e-8)
            for left, right in zip(intervals, intervals[1:]):
                if right[0] - left[1] > tolerance:
                    raise ValueError("particle ray intersects compound sections with a gap")
        return intervals

    def particle_track_segments(self, particle: "Particle") -> List[Tuple[np.ndarray, float, int, float]]:
        reference = np.array([particle.entry_x_m, particle.entry_y_m, self.bounds_max_m[2]], dtype=float)
        margin = max(4.0 * self.characteristic_length_m, 1e-6)
        origin = reference - particle.direction * margin
        intervals = self._particle_intervals(origin, particle.direction)
        if not intervals:
            raise ValueError("particle track does not intersect the compound geometry")
        first_start = intervals[0][0]
        segments = []
        for start, end, index in intervals:
            entry = origin + start * particle.direction
            segments.append((entry, end - start, index, start - first_start))
        return segments

    def particle_track(self, particle: "Particle") -> Tuple[np.ndarray, float]:
        segments = self.particle_track_segments(particle)
        return segments[0][0], sum(segment[1] for segment in segments)

    def refractive_index(self, wavelength_m: float | np.ndarray) -> float | np.ndarray:
        materials = {section.normalized_material() for section in self.sections}
        if len(materials) != 1:
            raise ValueError("compound refractive_index requires a position because sections have different materials")
        return self.sections[0].refractive_index(wavelength_m)

    def surface_rings(self, samples: int = 48) -> List[np.ndarray]:
        rings = []
        for section in self.sections:
            rings.extend(section.surface_rings(samples=samples))
        return rings


@dataclass(frozen=True)
class ScintillatorAssembly:
    """Complete scintillator assembly in the explicit transport frame.

    ``transport_sections`` is the ordered physical compound, while the two
    index tuples retain ownership of active and guide sections.  The legacy
    ``scintillator``/``light_guide`` fields remain compatibility aliases for
    callers that construct a one-guide assembly directly.  A canonical GUI
    model must populate all sections through ``from_transport_geometry``;
    it is never reduced to a single displayed guide.
    """

    scintillator: SolidGeometry
    light_guide: Optional[SolidGeometry] = None
    gap_from_quartz_m: float = 0.0
    mylar_enabled: bool = True
    mylar_reflectivity: float = 0.95
    terminal_plus_surface: str = "detector"
    terminal_minus_surface: str = "open"
    coordinate_frame: str = "global_m"
    # The charged-particle impact is expressed in this assembly's transport
    # frame.  Keeping it with the physical assembly prevents the viewer from
    # inventing a marker at the quartz origin when the active volume is
    # translated or rotated for display.
    particle_impact_x_m: float = 0.0
    particle_impact_y_m: float = 0.0
    transport_sections: Tuple[SolidGeometry, ...] = field(default_factory=tuple)
    active_section_indices: Tuple[int, ...] = field(default_factory=tuple)
    guide_section_indices: Tuple[int, ...] = field(default_factory=tuple)
    section_names: Tuple[str, ...] = field(default_factory=tuple)

    @property
    def ordered_sections(self) -> Tuple[SolidGeometry, ...]:
        if self.transport_sections:
            return tuple(self.transport_sections)
        return (self.scintillator,) if self.light_guide is None else (self.scintillator, self.light_guide)

    @property
    def active_sections(self) -> Tuple[SolidGeometry, ...]:
        sections = self.ordered_sections
        indices = self.active_section_indices or (0,)
        return tuple(sections[index] for index in indices)

    @property
    def guide_sections(self) -> Tuple[SolidGeometry, ...]:
        sections = self.ordered_sections
        if self.guide_section_indices:
            return tuple(sections[index] for index in self.guide_section_indices)
        return () if self.light_guide is None else (self.light_guide,)

    @property
    def guide_names(self) -> Tuple[str, ...]:
        if len(self.section_names) == len(self.ordered_sections):
            return tuple(self.section_names[index] for index in self.guide_section_indices)
        return tuple(f"guide_{index + 1}" for index in range(len(self.guide_sections)))

    @classmethod
    def from_transport_geometry(
        cls, geometry: CompoundGeometry, *, active_section_indices: Sequence[int],
        guide_section_indices: Sequence[int], section_names: Sequence[str] = (),
        gap_from_quartz_m: float = 0.0, mylar_enabled: bool = True,
        mylar_reflectivity: float = 0.95, terminal_plus_surface: str = "detector",
        terminal_minus_surface: str = "open", coordinate_frame: str = "scintillator_local_m",
        particle_impact_x_m: float = 0.0, particle_impact_y_m: float = 0.0,
    ) -> "ScintillatorAssembly":
        """Adapt a validated complete compound without dropping sections."""

        if not isinstance(geometry, CompoundGeometry):
            raise TypeError("complete scintillator transport requires CompoundGeometry")
        geometry.validate()
        sections = tuple(geometry.sections)
        active_indices = tuple(int(index) for index in active_section_indices)
        guide_indices = tuple(int(index) for index in guide_section_indices)
        if not active_indices:
            raise ValueError("complete scintillator transport requires at least one active section")
        if not guide_indices and len(sections) != len(active_indices):
            raise ValueError("complete scintillator transport lost configured guide sections")
        assigned = active_indices + guide_indices
        if sorted(assigned) != list(range(len(sections))) or len(set(assigned)) != len(assigned):
            raise ValueError("active and guide section ownership does not cover the complete transport geometry")
        names = tuple(str(name) for name in section_names)
        if names and len(names) != len(sections):
            raise ValueError("section_names must match the complete transport geometry")
        return cls(
            sections[active_indices[0]], sections[guide_indices[0]] if guide_indices else None,
            gap_from_quartz_m=gap_from_quartz_m, mylar_enabled=mylar_enabled,
            mylar_reflectivity=mylar_reflectivity, terminal_plus_surface=terminal_plus_surface,
            terminal_minus_surface=terminal_minus_surface, coordinate_frame=coordinate_frame,
            particle_impact_x_m=particle_impact_x_m, particle_impact_y_m=particle_impact_y_m,
            transport_sections=sections, active_section_indices=active_indices,
            guide_section_indices=guide_indices, section_names=names,
        )

    def validate(self) -> None:
        if not isinstance(self.scintillator, SolidGeometry):
            raise TypeError("scintillator must be a SolidGeometry object")
        sections = self.ordered_sections
        if not sections:
            raise ValueError("scintillator assembly requires at least one transport section")
        for section in sections:
            if not isinstance(section, SolidGeometry):
                raise TypeError("scintillator transport sections must be SolidGeometry objects")
            section.validate()
            if section.shape not in {"box", "tapered_rectangular_prism"} and section.axis != "X":
                raise ValueError(
                    f"unsupported scintillator transport orientation: section shape={section.shape} axis={section.axis}; "
                    "only X-aligned sections have a supported sequential interface"
                )
        active_indices = self.active_section_indices or (0,)
        guide_indices = self.guide_section_indices or (() if self.light_guide is None else (1,))
        if sorted(active_indices + guide_indices) != list(range(len(sections))) or len(set(active_indices + guide_indices)) != len(sections):
            raise ValueError("scintillator transport section ownership is incomplete or duplicated")
        if any(index < 0 or index >= len(sections) for index in active_indices + guide_indices):
            raise ValueError("scintillator transport section index is out of range")
        if not math.isfinite(float(self.gap_from_quartz_m)) or self.gap_from_quartz_m < 0.0:
            raise ValueError("gap_from_quartz_m must be finite and non-negative")
        if not 0.80 <= float(self.mylar_reflectivity) <= 1.00:
            raise ValueError("mylar_reflectivity must be between 0.80 and 1.00")
        _terminal_mode(self.terminal_plus_surface, "detector")
        _terminal_mode(self.terminal_minus_surface, "open")
        normalize_coordinate_frame(self.coordinate_frame)
        for name, value in (("particle_impact_x_m", self.particle_impact_x_m), ("particle_impact_y_m", self.particle_impact_y_m)):
            if not math.isfinite(float(value)):
                raise ValueError(f"{name} must be finite")
        # CompoundGeometry performs the exact ordered-contiguity and shared
        # centre validation.  Do it here as well so direct assembly callers
        # cannot bypass the complete-geometry contract.
        CompoundGeometry(
            sections, allow_interface_mismatch=True,
            wrapped_external=bool(self.mylar_enabled),
            mylar_reflectivity=float(self.mylar_reflectivity),
        ).validate()
        detector_face = "plus_x" if _terminal_mode(self.terminal_plus_surface, "detector") == "detector" else "minus_x"
        if guide_indices:
            active_set, guide_set = set(active_indices), set(guide_indices)
            if detector_face == "plus_x" and not max(active_set) < min(guide_set):
                raise ValueError("configured guide sections must follow all active sections toward the +X readout")
            if detector_face == "minus_x" and not max(guide_set) < min(active_set):
                raise ValueError("configured guide sections must precede all active sections toward the -X readout")

    @property
    def transport_geometry(self) -> CompoundGeometry:
        self.validate()
        sections = self.ordered_sections
        return CompoundGeometry(
            sections,
            allow_interface_mismatch=True,
            wrapped_external=bool(self.mylar_enabled),
            mylar_reflectivity=float(self.mylar_reflectivity),
            terminal_interfaces=(
                {"face": "minus_x", "axis_face": "YZ", "state": _terminal_mode(self.terminal_minus_surface, "open")},
                {"face": "plus_x", "axis_face": "YZ", "state": _terminal_mode(self.terminal_plus_surface, "detector")},
            ),
        )

    @property
    def interface_faces(self) -> Tuple[str, ...]:
        if not self.guide_sections:
            return ()
        names = self.section_names if len(self.section_names) == len(self.ordered_sections) else tuple(f"section_{i + 1}" for i in range(len(self.ordered_sections)))
        active_indices = self.active_section_indices or (0,)
        guide_indices = self.guide_section_indices or ((1,) if self.light_guide is not None else ())
        return tuple(f"{names[left]}_{names[right]}" for left, right in zip(range(len(names) - 1), range(1, len(names))) if left in active_indices or right in guide_indices)

    @property
    def wrapped_faces(self) -> Tuple[str, ...]:
        if not self.mylar_enabled:
            return ()
        return ("scintillator_external",) if self.light_guide is None else ("scintillator_external", "light_guide_external")

    def relative_placement(self) -> Dict[str, object]:
        self.validate()
        return {
            "scintillator_x_mm": [self.scintillator.bounds_min_m[0] * 1e3, self.scintillator.bounds_max_m[0] * 1e3],
            "light_guide_x_mm": None if not self.guide_sections else [self.guide_sections[0].bounds_min_m[0] * 1e3, self.guide_sections[-1].bounds_max_m[0] * 1e3],
            "light_guide_sections": [[section.bounds_min_m[0] * 1e3, section.bounds_max_m[0] * 1e3] for section in self.guide_sections],
            "gap_from_quartz_mm": self.gap_from_quartz_m * 1e3,
            "coupling_face_wrapped": False,
            "readout_face_wrapped": False,
            "mylar_external_faces": list(self.wrapped_faces),
        }

    def preview(self) -> Dict[str, object]:
        """Return a deterministic, JSON-friendly assembly preview model."""

        self.validate()
        names = self.section_names if len(self.section_names) == len(self.ordered_sections) else tuple("scintillator" if index in (self.active_section_indices or (0,)) else f"light_guide_{index}" for index in range(len(self.ordered_sections)))
        guide_indices = self.guide_section_indices or ((1,) if self.light_guide is not None else ())
        readout_plus = _terminal_mode(self.terminal_plus_surface, "detector") == "detector"
        readout_name = (
            f"{names[guide_indices[-1] if readout_plus else guide_indices[0]]} {'+X' if readout_plus else '-X'}"
            if guide_indices else f"scintillator {'+X' if readout_plus else '-X'}"
        )
        return {
            "beam_axis": "+X assembly axis; particle beam remains the existing -Z track",
            "gap_from_quartz_mm": self.gap_from_quartz_m * 1e3,
            "mylar_enabled": bool(self.mylar_enabled),
            "mylar_reflectivity": float(self.mylar_reflectivity),
            "volumes": [
                {"name": names[index], "role": "guide" if index in guide_indices else "active", "shape": section.shape, "bounds_mm": (section.bounds_min_m * 1e3).tolist() + (section.bounds_max_m * 1e3).tolist()}
                for index, section in enumerate(self.ordered_sections)
            ],
            "interfaces": list(self.interface_faces),
            "readout_output_face": readout_name,
            "terminal_states": {"minus_x": _terminal_mode(self.terminal_minus_surface, "open"), "plus_x": _terminal_mode(self.terminal_plus_surface, "detector")},
        }


@dataclass(frozen=True)
class DetectorAssembly:
    """One detector instrument: scintillator plus two opposed light guides.

    The three solids are kept in the detector-local coordinate system and are
    sequential along X: left guide -> scintillator -> right guide.  A global
    placement is applied only by :class:`TwoDetectorBeamMonitor`, so changing
    an offset cannot change any local dimensions.
    """

    scintillator: SolidGeometry
    guide_left: Optional[SolidGeometry]
    guide_right: Optional[SolidGeometry]
    detector_id: int = 1
    mylar_enabled: bool = True
    mylar_reflectivity: float = 0.95
    # ``X`` preserves the original optical assembly. ``Z`` is the engineering
    # beam-monitor representation: the active volume is an XY plane with
    # thickness along Z and optional guides sit on its +/-Y (YZ) faces.
    active_axis: str = "X"

    def validate(self) -> None:
        solids = (self.guide_left, self.scintillator, self.guide_right)
        if not isinstance(self.scintillator, SolidGeometry):
            raise TypeError("detector assembly scintillator must be a SolidGeometry")
        axis = str(self.active_axis).upper()
        if axis not in {"X", "Z"}:
            raise ValueError("detector assembly active_axis must be X or Z")
        self.scintillator.validate()
        present = tuple(solid for solid in solids if solid is not None)
        if any(not isinstance(solid, SolidGeometry) for solid in present):
            raise TypeError("detector assembly guides must be SolidGeometry objects or None")
        for solid in present:
            solid.validate()
        if axis == "X":
            if any(solid.shape not in {"box", "tapered_rectangular_prism"} and solid.axis != "X" for solid in present + (self.scintillator,)):
                raise ValueError("two-detector X assembly volumes must be X-aligned")
        else:
            # In the Z-axis engineering model the active face is XY. Guides
            # are independent optional solids on the +/-Y faces; they are not
            # silently required and need not be part of optical CompoundGeometry.
            s_low, s_high = self.scintillator.bounds_min_m, self.scintillator.bounds_max_m
            for side, guide in (("left", self.guide_left), ("right", self.guide_right)):
                if guide is None:
                    continue
                g_low, g_high = guide.bounds_min_m, guide.bounds_max_m
                tol = max(1e-12, self.scintillator.characteristic_length_m * 1e-9)
                if g_high[2] < s_low[2] - tol or g_low[2] > s_high[2] + tol:
                    raise ValueError(f"{side} light guide must overlap the active volume in Z")
        if int(self.detector_id) < 1:
            raise ValueError("detector_id must be positive")
        if not 0.80 <= float(self.mylar_reflectivity) <= 1.00:
            raise ValueError("mylar_reflectivity must be between 0.80 and 1.00")
        if axis == "X":
            for left, right, name in (
                (self.guide_left, self.scintillator, "left guide/scintillator"),
                (self.scintillator, self.guide_right, "scintillator/right guide"),
            ):
                if left is None or right is None:
                    continue
                tolerance = max(1e-12, max(abs(left.bounds_max_m[0]), abs(right.bounds_min_m[0]), 1e-6) * 1e-9)
                if abs(left.bounds_max_m[0] - right.bounds_min_m[0]) > tolerance:
                    if left.bounds_max_m[0] > right.bounds_min_m[0] + tolerance:
                        raise ValueError(f"{name} overlap along X")
                    raise ValueError(f"{name} must be sequential with a zero interface gap")
            centers_yz = [(solid.center_y_m, solid.center_z_m) for solid in present + (self.scintillator,)]
            if any(not math.isclose(centers_yz[0][axis], center[axis], abs_tol=1e-12, rel_tol=0.0)
                   for center in centers_yz[1:] for axis in (0, 1)):
                raise ValueError("detector assembly X volumes must share Y/Z center")

    @property
    def bounds_min_m(self) -> np.ndarray:
        self.validate()
        return np.minimum.reduce([solid.bounds_min_m for solid in (self.guide_left, self.scintillator, self.guide_right) if solid is not None])

    @property
    def bounds_max_m(self) -> np.ndarray:
        self.validate()
        return np.maximum.reduce([solid.bounds_max_m for solid in (self.guide_left, self.scintillator, self.guide_right) if solid is not None])

    @property
    def transport_geometry(self) -> CompoundGeometry:
        self.validate()
        if self.active_axis == "Z":
            # The established photon transport is a sequential X model. Do
            # not pretend that the new two-plane engineering scene is already
            # a coupled multi-detector optical transport; retain the active
            # scintillator as a valid single-volume fallback.
            return CompoundGeometry((self.scintillator,), allow_interface_mismatch=True,
                                    wrapped_external=bool(self.mylar_enabled),
                                    mylar_reflectivity=float(self.mylar_reflectivity))
        return CompoundGeometry(
            tuple(solid for solid in (self.guide_left, self.scintillator, self.guide_right) if solid is not None),
            allow_interface_mismatch=True,
            wrapped_external=bool(self.mylar_enabled),
            mylar_reflectivity=float(self.mylar_reflectivity),
        )

    @property
    def interface_faces(self) -> Tuple[str, ...]:
        if self.active_axis == "Z":
            return tuple(name for name, solid in (
                ("guide_left_scintillator_Y-", self.guide_left),
                ("scintillator_guide_right_Y+", self.guide_right),
            ) if solid is not None)
        return tuple(name for name, solid in (
            ("guide_left_scintillator", self.guide_left),
            ("scintillator_guide_right", self.guide_right),
        ) if solid is not None)

    @property
    def readout_faces(self) -> Tuple[str, ...]:
        return tuple(name for name, solid in (
            ("guide_left_readout", self.guide_left),
            ("guide_right_readout", self.guide_right),
        ) if solid is not None)

    @property
    def wrapped_faces(self) -> Tuple[str, ...]:
        if not self.mylar_enabled:
            return ()
        return tuple(name for name, solid in (
            ("guide_left_external", self.guide_left),
            ("scintillator_external", self.scintillator),
            ("guide_right_external", self.guide_right),
        ) if solid is not None)

    def translated(self, dx_m: float = 0.0, dy_m: float = 0.0, dz_m: float = 0.0) -> "DetectorAssembly":
        """Return a placed copy; all local dimensions are preserved."""

        self.validate()
        return replace(
            self,
            scintillator=_translate_solid_geometry(self.scintillator, dx_m, dy_m, dz_m),
            guide_left=None if self.guide_left is None else _translate_solid_geometry(self.guide_left, dx_m, dy_m, dz_m),
            guide_right=None if self.guide_right is None else _translate_solid_geometry(self.guide_right, dx_m, dy_m, dz_m),
        )

    def preview(self) -> Dict[str, object]:
        self.validate()
        volumes = []
        for name, solid in (("guide_left", self.guide_left), ("scintillator", self.scintillator), ("guide_right", self.guide_right)):
            if solid is not None:
                volumes.append({"name": name, "shape": solid.shape,
                                "bounds_mm": (solid.bounds_min_m * 1e3).tolist() + (solid.bounds_max_m * 1e3).tolist()})
        return {
            "detector_id": int(self.detector_id),
            "active_axis": self.active_axis,
            "volumes": volumes,
            "interfaces": list(self.interface_faces),
            "readout_faces": list(self.readout_faces),
            "mylar_external_faces": list(self.wrapped_faces),
            "mylar_reflectivity": float(self.mylar_reflectivity),
        }


@dataclass(frozen=True)
class DetectorPlacement:
    """Signed global rigid transform for one detector assembly."""

    offset_x_m: float = 0.0
    offset_y_m: float = 0.0
    offset_z_m: float = 0.0
    rotation_axis: str = "None"
    rotation_angle_rad: float = 0.0

    def validate(self) -> None:
        for name in ("offset_x_m", "offset_y_m", "offset_z_m"):
            if not math.isfinite(float(getattr(self, name))):
                raise ValueError(f"{name} must be finite; signed values are allowed")
        if str(self.rotation_axis).upper() not in {"NONE", "X", "Y", "Z"}:
            raise ValueError("rotation_axis must be None, X, Y, or Z")
        if not math.isfinite(float(self.rotation_angle_rad)):
            raise ValueError("rotation_angle_rad must be finite")


@dataclass(frozen=True)
class TwoDetectorBeamMonitor:
    """Global beam monitor with two distinct detector assemblies.

    The Cherenkov radiator is followed by detector 1, an explicit physical
    separation, and detector 2.  Each detector has its own scintillator and
    two readout guides.  The existing single-assembly scintillation transport
    remains available; this class is the authoritative global placement,
    engineering-scene, and common-particle-path model for the two-detector
    instrument.
    """

    detector_1: DetectorAssembly
    detector_2: DetectorAssembly
    quartz_gap_m: float = 0.010
    detector_separation_m: float = 0.020
    placement_1: DetectorPlacement = field(default_factory=DetectorPlacement)
    placement_2: DetectorPlacement = field(default_factory=DetectorPlacement)
    beam_axis: str = "X"
    # Nominal active-plane centres.  They are explicit so Detector 2 is not
    # an implicit copy/translation of Detector 1.  Optional keeps old JSON
    # profiles compatible; legacy X profiles use the original sequential
    # placement path below.
    detector_1_z_m: Optional[float] = None
    detector_2_z_m: Optional[float] = None

    def validate(self, quartz_geometry=None) -> None:
        self.detector_1.validate()
        self.detector_2.validate()
        for name, value in (("quartz_gap_m", self.quartz_gap_m), ("detector_separation_m", self.detector_separation_m)):
            if not math.isfinite(float(value)) or float(value) < 0.0:
                raise ValueError(f"{name} must be finite and non-negative")
        self.placement_1.validate()
        self.placement_2.validate()
        axis = str(self.beam_axis).upper()
        if axis not in {"X", "Y", "Z"}:
            raise ValueError("global beam axis must be X, Y, or Z")
        if axis == "Z":
            if self.detector_1_z_m is None or self.detector_2_z_m is None:
                raise ValueError("Z-axis two-detector monitor requires detector_1_z_m and detector_2_z_m")
            if not all(math.isfinite(float(value)) for value in (self.detector_1_z_m, self.detector_2_z_m)):
                raise ValueError("detector Z positions must be finite")
            if float(self.detector_2_z_m) <= float(self.detector_1_z_m):
                raise ValueError("detector_2_z_m must be greater than detector_1_z_m")
        if quartz_geometry is not None:
            self.placed_detectors(quartz_geometry, _validate_only=True)

    def placed_detectors(self, quartz_geometry=None, *, _validate_only: bool = False) -> Tuple[DetectorAssembly, DetectorAssembly]:
        """Apply signed placements without mutating either local detector."""

        self.detector_1.validate()
        self.detector_2.validate()
        quartz = quartz_geometry or BarGeometry()
        quartz_low, quartz_high = _geometry_bounds(quartz)
        quartz_center = (quartz_low + quartz_high) / 2.0
        if str(self.beam_axis).upper() == "Z":
            # Explicit Z values are nominal active-plane centres.  A larger
            # physical gap/separation is authoritative and moves the placed
            # volume downstream without changing its local dimensions.
            d1_min_offset = float(self.detector_1.bounds_min_m[2] - self.detector_1.scintillator.center_z_m)
            d1_max_offset = float(self.detector_1.bounds_max_m[2] - self.detector_1.scintillator.center_z_m)
            d2_min_offset = float(self.detector_2.bounds_min_m[2] - self.detector_2.scintillator.center_z_m)
            minimum_d1_z = float(quartz_high[2]) + float(self.quartz_gap_m) - d1_min_offset
            effective_d1_z = max(float(self.detector_1_z_m), minimum_d1_z)
            d1_actual_center_z = effective_d1_z + self.placement_1.offset_z_m
            minimum_d2_z = d1_actual_center_z + d1_max_offset + float(self.detector_separation_m) - d2_min_offset
            effective_d2_z = max(float(self.detector_2_z_m), minimum_d2_z)
            d1_target = np.array([
                quartz_center[0] + self.placement_1.offset_x_m,
                quartz_center[1] + self.placement_1.offset_y_m,
                effective_d1_z + self.placement_1.offset_z_m,
            ])
            d2_target = np.array([
                quartz_center[0] + self.placement_2.offset_x_m,
                quartz_center[1] + self.placement_2.offset_y_m,
                effective_d2_z + self.placement_2.offset_z_m,
            ])
            d1 = self.detector_1.translated(*(d1_target - self.detector_1.scintillator.center_m))
            d2 = self.detector_2.translated(*(d2_target - self.detector_2.scintillator.center_m))
            tol = max(1e-12, self.detector_1.scintillator.characteristic_length_m * 1e-9)
            if d1.bounds_min_m[2] < quartz_high[2] + self.quartz_gap_m - tol:
                raise ValueError("detector 1 overlaps quartz or violates the quartz gap along Z")
            if d2.bounds_min_m[2] < d1.bounds_max_m[2] + self.detector_separation_m - tol:
                raise ValueError("detector 2 overlaps detector 1 or violates the detector separation along Z")
            return d1, d2
        detector_1_start = float(quartz_high[0]) + float(self.quartz_gap_m) + float(self.placement_1.offset_x_m)
        d1 = self.detector_1.translated(
            detector_1_start - self.detector_1.bounds_min_m[0],
            float(quartz_center[1]) + float(self.placement_1.offset_y_m) - self.detector_1.scintillator.center_y_m,
            float(quartz_center[2]) + float(self.placement_1.offset_z_m) - self.detector_1.scintillator.center_z_m,
        )
        detector_2_start = float(d1.bounds_max_m[0]) + float(self.detector_separation_m) + float(self.placement_2.offset_x_m)
        d2 = self.detector_2.translated(
            detector_2_start - self.detector_2.bounds_min_m[0],
            float(quartz_center[1]) + float(self.placement_2.offset_y_m) - self.detector_2.scintillator.center_y_m,
            float(quartz_center[2]) + float(self.placement_2.offset_z_m) - self.detector_2.scintillator.center_z_m,
        )
        tolerance = max(1e-12, abs(float(quartz_high[0])) * 1e-9)
        if d1.bounds_min_m[0] < quartz_high[0] - tolerance:
            raise ValueError("detector 1 signed offset overlaps the quartz radiator")
        if d2.bounds_min_m[0] < d1.bounds_max_m[0] - tolerance:
            raise ValueError("detector 2 signed offset overlaps detector 1")
        return d1, d2

    def particle_path(self, quartz_geometry=None) -> Dict[str, object]:
        """Return one common engineering particle path through both detectors."""

        quartz = quartz_geometry or BarGeometry()
        self.validate(quartz)
        d1, d2 = self.placed_detectors(quartz)
        q_low, q_high = _geometry_bounds(quartz)
        q_center = (q_low + q_high) / 2.0
        if str(self.beam_axis).upper() == "Z":
            points = [
                (q_center[0], q_center[1], q_low[2] - 0.020),
                (q_center[0], q_center[1], q_low[2]),
                (q_center[0], q_center[1], q_high[2]),
                (d1.scintillator.center_x_m, d1.scintillator.center_y_m, d1.scintillator.bounds_min_m[2]),
                (d1.scintillator.center_x_m, d1.scintillator.center_y_m, d1.scintillator.bounds_max_m[2]),
                (d2.scintillator.center_x_m, d2.scintillator.center_y_m, d2.scintillator.bounds_min_m[2]),
                (d2.scintillator.center_x_m, d2.scintillator.center_y_m, d2.scintillator.bounds_max_m[2]),
                (d2.scintillator.center_x_m, d2.scintillator.center_y_m, d2.bounds_max_m[2] + 0.020),
            ]
        else:
            points = [
                (q_low[0] - 0.020, q_center[1], q_center[2]),
                (q_low[0], q_center[1], q_center[2]),
                (q_high[0], q_center[1], q_center[2]),
                (d1.scintillator.bounds_min_m[0], d1.scintillator.center_y_m, d1.scintillator.center_z_m),
                (d1.scintillator.bounds_max_m[0], d1.scintillator.center_y_m, d1.scintillator.center_z_m),
                (d2.scintillator.bounds_min_m[0], d2.scintillator.center_y_m, d2.scintillator.center_z_m),
                (d2.scintillator.bounds_max_m[0], d2.scintillator.center_y_m, d2.scintillator.center_z_m),
                (d2.bounds_max_m[0] + 0.020, d2.scintillator.center_y_m, d2.scintillator.center_z_m),
            ]
        labels = (
            "particle entry", "quartz exit", "Detector 1 scintillator entry",
            "Detector 1 scintillator exit", "Detector 2 scintillator entry",
            "Detector 2 scintillator exit", "particle exit",
        )
        segments = [
            {"name": "entry_to_quartz", "start_mm": np.asarray(points[0]) * 1e3, "end_mm": np.asarray(points[1]) * 1e3},
            {"name": "quartz", "start_mm": np.asarray(points[1]) * 1e3, "end_mm": np.asarray(points[2]) * 1e3},
            {"name": "quartz_to_detector_1", "start_mm": np.asarray(points[2]) * 1e3, "end_mm": np.asarray(points[3]) * 1e3},
            {"name": "detector_1_scintillator", "start_mm": np.asarray(points[3]) * 1e3, "end_mm": np.asarray(points[4]) * 1e3},
            {"name": "detector_separation", "start_mm": np.asarray(points[4]) * 1e3, "end_mm": np.asarray(points[5]) * 1e3},
            {"name": "detector_2_scintillator", "start_mm": np.asarray(points[5]) * 1e3, "end_mm": np.asarray(points[6]) * 1e3},
            {"name": "exit", "start_mm": np.asarray(points[6]) * 1e3, "end_mm": np.asarray(points[7]) * 1e3},
        ]
        return {
            "points_mm": [list(np.asarray(point, dtype=float) * 1e3) for point in points],
            "labels": list(labels),
            "segments": [{key: (value.tolist() if isinstance(value, np.ndarray) else value) for key, value in segment.items()} for segment in segments],
            "order": "entry -> quartz -> Detector 1 XY active plane -> separation -> Detector 2 XY active plane -> exit" if str(self.beam_axis).upper() == "Z" else "entry -> quartz -> Detector 1 scintillator -> separation -> Detector 2 scintillator -> exit",
        }

    def relative_placement(self, quartz_geometry=None) -> Dict[str, object]:
        self.validate(quartz_geometry)
        d1, d2 = self.placed_detectors(quartz_geometry)
        return {
            "quartz_gap_mm": float(self.quartz_gap_m) * 1e3,
            "detector_separation_mm": float(self.detector_separation_m) * 1e3,
            "detector_1_signed_offset_mm": [self.placement_1.offset_x_m * 1e3, self.placement_1.offset_y_m * 1e3, self.placement_1.offset_z_m * 1e3],
            "detector_2_signed_offset_mm": [self.placement_2.offset_x_m * 1e3, self.placement_2.offset_y_m * 1e3, self.placement_2.offset_z_m * 1e3],
            "detector_1_rigid_rotation": {"axis": self.placement_1.rotation_axis, "angle_deg": math.degrees(self.placement_1.rotation_angle_rad)},
            "detector_2_rigid_rotation": {"axis": self.placement_2.rotation_axis, "angle_deg": math.degrees(self.placement_2.rotation_angle_rad)},
            "detector_1_z_mm": None if self.detector_1_z_m is None else float(self.detector_1_z_m) * 1e3,
            "detector_2_z_mm": None if self.detector_2_z_m is None else float(self.detector_2_z_m) * 1e3,
            "detector_1_bounds_mm": (d1.bounds_min_m * 1e3).tolist() + (d1.bounds_max_m * 1e3).tolist(),
            "detector_2_bounds_mm": (d2.bounds_min_m * 1e3).tolist() + (d2.bounds_max_m * 1e3).tolist(),
            "beam_axis": str(self.beam_axis).upper(),
        }

    def preview(self, quartz_geometry=None) -> Dict[str, object]:
        self.validate(quartz_geometry)
        d1, d2 = self.placed_detectors(quartz_geometry)
        return {
            "model": "TwoDetectorBeamMonitor",
            "beam_axis": str(self.beam_axis).upper(),
            "detector_1_z_mm": None if self.detector_1_z_m is None else float(self.detector_1_z_m) * 1e3,
            "detector_2_z_mm": None if self.detector_2_z_m is None else float(self.detector_2_z_m) * 1e3,
            "detectors": {"detector_1": d1.preview(), "detector_2": d2.preview()},
            "relative_placement": self.relative_placement(quartz_geometry),
            "particle_path": self.particle_path(quartz_geometry),
            "quartz_gap_mm": float(self.quartz_gap_m) * 1e3,
            "detector_separation_mm": float(self.detector_separation_m) * 1e3,
        }


@dataclass(frozen=True)
class DetectorSensor:
    """A readout sensor kept separate from active optical volumes.

    ``coupling_mode`` controls only the engineering connection: ``direct``
    couples to the selected active component, while ``light_guide`` requires
    one of the explicitly present guide components.  Response parameters are
    intentionally compact and are consumed by future multi-detector optical
    transport without changing the established Cherenkov detector API.
    """

    sensor_id: str = "sensor_1"
    component: str = "scintillator_1"
    coupling_mode: str = "direct"
    side: str = "Y+"
    position_x_m: float = 0.0
    position_y_m: float = 0.0
    position_z_m: float = 0.0
    shape: str = "rectangular"
    width_mm: float = 20.0
    height_mm: float = 10.0
    # v62 explicit readout dimensions.  width/height remain full extents for
    # v61 profile compatibility; circle uses radius and ellipse uses A/B.
    readout_radius_mm: Optional[float] = None
    readout_semiaxis_a_mm: Optional[float] = None
    readout_semiaxis_b_mm: Optional[float] = None
    quantum_efficiency: float = 1.0
    collection_efficiency: float = 1.0
    timing_sigma_s: float = 0.0

    @staticmethod
    def normalize_coupling_mode(value: object) -> str:
        """Accept canonical and legacy spellings used by stored profiles."""

        text = str(value or "direct").strip().lower().replace("-", "_").replace(" ", "_")
        aliases = {
            "direct": "direct", "active": "direct", "active_component": "direct",
            "light_guide": "light_guide", "guide": "light_guide", "through_guide": "light_guide",
            "complete_light_guide": "light_guide", "through_complete_light_guide": "light_guide",
        }
        return aliases.get(text, text)

    def validate(self, available_components: Optional[Sequence[str]] = None) -> None:
        if not str(self.sensor_id).strip():
            raise ValueError("sensor_id cannot be empty")
        if self.normalize_coupling_mode(self.coupling_mode) not in {"direct", "light_guide"}:
            raise ValueError("sensor coupling_mode must be direct or light_guide")
        if str(self.side).upper() not in {"Y-", "Y+", "Z-", "Z+", "X-", "X+"}:
            raise ValueError("sensor side must be one of X-/X+/Y-/Y+/Z-/Z+")
        if available_components is not None and self.component not in set(available_components):
            raise ValueError(f"sensor component is not available: {self.component}")
        if self.shape not in {"rectangular", "circular", "elliptical"}:
            raise ValueError("sensor shape must be rectangular, circular, or elliptical")
        for name in ("position_x_m", "position_y_m", "position_z_m", "timing_sigma_s"):
            if not math.isfinite(float(getattr(self, name))):
                raise ValueError(f"sensor {name} must be finite")
        for name in ("width_mm", "height_mm"):
            if not math.isfinite(float(getattr(self, name))) or float(getattr(self, name)) <= 0.0:
                raise ValueError(f"sensor {name} must be positive")
        if self.shape == "circular":
            value = min(self.width_mm, self.height_mm) / 2.0 if self.readout_radius_mm is None else float(self.readout_radius_mm)
            if not math.isfinite(value) or value <= 0.0:
                raise ValueError("sensor readout radius must be positive")
        if self.shape == "elliptical":
            for name, fallback in (("readout_semiaxis_a_mm", self.width_mm / 2.0), ("readout_semiaxis_b_mm", self.height_mm / 2.0)):
                value = fallback if getattr(self, name) is None else float(getattr(self, name))
                if not math.isfinite(value) or value <= 0.0:
                    raise ValueError(f"sensor {name} must be positive")
        for name in ("quantum_efficiency", "collection_efficiency"):
            value = float(getattr(self, name))
            if not 0.0 <= value <= 1.0:
                raise ValueError(f"sensor {name} must be between 0 and 1")
        if self.timing_sigma_s < 0.0:
            raise ValueError("sensor timing_sigma_s must be non-negative")

    def to_dict(self) -> Dict[str, object]:
        self.validate()
        payload = asdict(self)
        payload["position_x_mm"] = payload.pop("position_x_m") * 1e3
        payload["position_y_mm"] = payload.pop("position_y_m") * 1e3
        payload["position_z_mm"] = payload.pop("position_z_m") * 1e3
        payload["timing_sigma_ns"] = payload.pop("timing_sigma_s") * 1e9
        payload["readout_radius_mm"] = (float(self.readout_radius_mm) if self.readout_radius_mm is not None else min(self.width_mm, self.height_mm) / 2.0)
        payload["readout_semiaxes_mm"] = {
            "a": float(self.readout_semiaxis_a_mm) if self.readout_semiaxis_a_mm is not None else self.width_mm / 2.0,
            "b": float(self.readout_semiaxis_b_mm) if self.readout_semiaxis_b_mm is not None else self.height_mm / 2.0,
        }
        return payload

    @classmethod
    def from_dict(cls, payload: Mapping[str, object]) -> "DetectorSensor":
        def number(name, default=0.0):
            value = float(payload.get(name, default))
            if not math.isfinite(value):
                raise ValueError(f"sensor {name} must be finite")
            return value
        coupling_mode = cls.normalize_coupling_mode(payload.get("coupling_mode", "direct"))
        component = str(payload.get("component", "scintillator_1"))
        # Older profiles named individual guide components explicitly.  Those
        # names are now collapsed to the owning active component; actual
        # guide routing is derived from the assembly model.
        if component.startswith("guide_"):
            component = "scintillator_1"
        sensor = cls(
            sensor_id=str(payload.get("sensor_id", "sensor_1")),
            component=component,
            coupling_mode=coupling_mode,
            side=str(payload.get("side", "Y+")),
            position_x_m=number("position_x_mm") * 1e-3,
            position_y_m=number("position_y_mm") * 1e-3,
            position_z_m=number("position_z_mm") * 1e-3,
            shape=str(payload.get("shape", "rectangular")),
            width_mm=number("width_mm", 20.0),
            height_mm=number("height_mm", 10.0),
            readout_radius_mm=None if payload.get("readout_radius_mm") in (None, "") else number("readout_radius_mm"),
            readout_semiaxis_a_mm=(float(payload.get("readout_semiaxes_mm", {}).get("a"))
                                   if isinstance(payload.get("readout_semiaxes_mm"), Mapping) and payload.get("readout_semiaxes_mm", {}).get("a") not in (None, "") else None),
            readout_semiaxis_b_mm=(float(payload.get("readout_semiaxes_mm", {}).get("b"))
                                   if isinstance(payload.get("readout_semiaxes_mm"), Mapping) and payload.get("readout_semiaxes_mm", {}).get("b") not in (None, "") else None),
            quantum_efficiency=number("quantum_efficiency", 1.0),
            collection_efficiency=number("collection_efficiency", 1.0),
            timing_sigma_s=number("timing_sigma_ns", 0.0) * 1e-9,
        )
        sensor.validate()
        return sensor

@dataclass(frozen=True)
class ScintillatorLocalGeometry:
    """Local scintillator/light-guide solids without global placement.

    The two solids are expressed in one assembly-local coordinate system.
    Their dimensions, shapes, materials, and local sequential interface are
    independent of the quartz position and of the downstream gap.
    """

    scintillator: SolidGeometry
    light_guide: Optional[SolidGeometry] = None
    particle_impact_x_m: float = 0.0
    particle_impact_y_m: float = 0.0

    def validate(self) -> None:
        # A zero gap assembly validates the local coupling and deliberately
        # carries no quartz/global placement state.
        ScintillatorAssembly(self.scintillator, self.light_guide, gap_from_quartz_m=0.0).validate()

    def to_assembly(self, *, mylar_enabled: bool = True, mylar_reflectivity: float = 0.95) -> ScintillatorAssembly:
        self.validate()
        return ScintillatorAssembly(
            self.scintillator, self.light_guide, gap_from_quartz_m=0.0,
            mylar_enabled=mylar_enabled, mylar_reflectivity=mylar_reflectivity,
            coordinate_frame="scintillator_local_m",
            particle_impact_x_m=self.particle_impact_x_m,
            particle_impact_y_m=self.particle_impact_y_m,
        )

    def preview(self) -> Dict[str, object]:
        return self.to_assembly(mylar_enabled=False).preview()


@dataclass(frozen=True)
class ScintillatorGlobalGeometry:
    """Global placement of a local assembly relative to the quartz radiator."""

    gap_from_quartz_m: float = 0.010
    offset_x_m: float = 0.0
    offset_y_m: float = 0.0
    offset_z_m: float = 0.0
    beam_axis: str = "X"

    def validate(self) -> None:
        if not math.isfinite(float(self.gap_from_quartz_m)) or self.gap_from_quartz_m < 0.0:
            raise ValueError("global scintillator gap must be finite and non-negative")
        for name, value in (("offset_x_m", self.offset_x_m), ("offset_y_m", self.offset_y_m), ("offset_z_m", self.offset_z_m)):
            if not math.isfinite(float(value)):
                raise ValueError(f"global scintillator {name} must be finite")
        if str(self.beam_axis).upper() not in {"X", "Y", "Z"}:
            raise ValueError("global beam axis must be X, Y, or Z")


def _translate_solid_geometry(section: SolidGeometry, dx: float, dy: float, dz: float) -> SolidGeometry:
    """Return a dimension-identical solid with only its global center moved."""

    return replace(section, center_x_m=section.center_x_m + dx, center_y_m=section.center_y_m + dy, center_z_m=section.center_z_m + dz)


def compose_scintillator_assembly(
    local_geometry: ScintillatorLocalGeometry,
    global_geometry: ScintillatorGlobalGeometry,
    quartz_geometry: Optional[Union[BarGeometry, SolidGeometry, CompoundGeometry]] = None,
    *,
    mylar_enabled: bool = True,
    mylar_reflectivity: float = 0.95,
    terminal_plus_surface: str = "detector",
    terminal_minus_surface: str = "open",
) -> ScintillatorAssembly:
    """Apply only the global transform/gap to a validated local geometry."""

    local_geometry.validate()
    global_geometry.validate()
    quartz = quartz_geometry or BarGeometry()
    if isinstance(quartz, BarGeometry):
        quartz_max_x = quartz.half_x
        quartz_center_y = quartz_center_z = 0.0
    else:
        quartz_max_x = float(quartz.bounds_max_m[0])
        quartz_center_y = float((quartz.bounds_min_m[1] + quartz.bounds_max_m[1]) / 2.0)
        quartz_center_z = float((quartz.bounds_min_m[2] + quartz.bounds_max_m[2]) / 2.0)
    local_min_x = float(local_geometry.scintillator.bounds_min_m[0])
    dx = quartz_max_x + global_geometry.gap_from_quartz_m + global_geometry.offset_x_m - local_min_x
    dy = quartz_center_y + global_geometry.offset_y_m
    dz = quartz_center_z + global_geometry.offset_z_m
    scintillator = _translate_solid_geometry(local_geometry.scintillator, dx, dy, dz)
    guide = None if local_geometry.light_guide is None else _translate_solid_geometry(local_geometry.light_guide, dx, dy, dz)
    return ScintillatorAssembly(
        scintillator, guide, gap_from_quartz_m=global_geometry.gap_from_quartz_m,
        mylar_enabled=mylar_enabled, mylar_reflectivity=mylar_reflectivity,
        terminal_plus_surface=terminal_plus_surface, terminal_minus_surface=terminal_minus_surface,
        particle_impact_x_m=local_geometry.particle_impact_x_m,
        particle_impact_y_m=local_geometry.particle_impact_y_m,
    )


def assembly_preview(assembly: ScintillatorAssembly) -> Dict[str, object]:
    """Public headless preview helper used by tests and non-Tk clients."""

    return assembly.preview()


@dataclass(frozen=True)
class ScintillatorConfig:
    """Optical and generation controls for the first scintillator model.

    ``decay_time_s`` is the mean of the exponential emission-delay
    distribution.  A value of zero means that emission is instantaneous and
    no additional delay is added to the charged-particle time at the source.
    """

    enabled: bool = False
    sampling_mode: str = "fixed"
    photon_count: int = 100
    photon_yield_per_MeV: float = 10000.0
    deposited_energy_MeV: Optional[float] = None
    dEdx_MeV_per_m: Optional[float] = None
    decay_time_s: float = 2.0e-9
    emission_wavelength_min_m: float = 350e-9
    emission_wavelength_max_m: float = 550e-9
    readout_index: float = 1.0
    readout_acceptance_rad: float = math.pi / 2.0
    readout_quantum_efficiency: float = 1.0
    readout_collection_efficiency: float = 1.0
    absorption_length_m: Optional[float] = None
    max_bounces: int = 100

    def validate(self) -> None:
        if self.sampling_mode not in {"fixed", "physical", "poisson"}:
            raise ValueError("scintillator sampling_mode must be fixed or physical")
        if self.photon_count < 0 or self.max_bounces < 0:
            raise ValueError("scintillator photon_count and max_bounces must be non-negative")
        if self.photon_yield_per_MeV < 0.0 or not math.isfinite(self.photon_yield_per_MeV):
            raise ValueError("scintillator photon yield must be finite and non-negative")
        for name, value in (("deposited_energy_MeV", self.deposited_energy_MeV), ("dEdx_MeV_per_m", self.dEdx_MeV_per_m)):
            if value is not None and (not math.isfinite(float(value)) or float(value) < 0.0):
                raise ValueError(f"{name} must be finite and non-negative")
        if self.decay_time_s < 0.0 or not math.isfinite(self.decay_time_s):
            raise ValueError("scintillator decay_time_s must be finite and non-negative")
        if self.emission_wavelength_min_m <= 0.0 or self.emission_wavelength_max_m <= self.emission_wavelength_min_m:
            raise ValueError("scintillator emission wavelength range is invalid")
        if self.readout_index <= 0.0 or not 0.0 <= self.readout_acceptance_rad <= math.pi:
            raise ValueError("scintillator readout index/acceptance is invalid")
        for name, value in (("readout_quantum_efficiency", self.readout_quantum_efficiency), ("readout_collection_efficiency", self.readout_collection_efficiency)):
            if not 0.0 <= value <= 1.0:
                raise ValueError(f"{name} must be between 0 and 1")
        if self.absorption_length_m is not None and self.absorption_length_m <= 0.0:
            raise ValueError("scintillator absorption_length_m must be positive")


def sample_scintillation_delay_s(decay_time_s: float, rng: np.random.Generator) -> float:
    """Draw one scintillation emission delay in seconds.

    The parameter is the exponential mean (the GUI exposes the same value in
    nanoseconds).  The zero case is explicit so an instantaneous scintillator
    does not consume an RNG draw or add a numerical delay.
    """

    mean_s = float(decay_time_s)
    if not math.isfinite(mean_s) or mean_s < 0.0:
        raise ValueError("scintillator decay_time_s must be finite and non-negative")
    return 0.0 if mean_s == 0.0 else float(rng.exponential(mean_s))


@dataclass
class ScintillatorResult:
    assembly: ScintillatorAssembly
    config: ScintillatorConfig
    photon_records: List[PhotonRecord] = field(default_factory=list)
    deposited_energy_MeV: float = 0.0
    physical_expected_photons: float = 0.0
    track_length_m: float = 0.0
    stopped: bool = False
    source_row_count: int = 0
    readout_mismatch_count: int = 0
    readout_endpoint_tolerance_m: float = 0.0

    @property
    def requested_photons(self) -> int:
        return self.config.photon_count

    def summary(self) -> Dict[str, object]:
        records = self.photon_records
        detected = [record for record in records if record.detected]
        readout = [record for record in records if record.detector_transmitted]
        arrivals = [readout_arrival_time_ns(record) for record in records]
        arrivals = [value for value in arrivals if value is not None]
        guide_section_summary = []
        for position, (name, section) in enumerate(zip(self.assembly.guide_names, self.assembly.guide_sections)):
            guide_section_summary.append({
                "name": name,
                "shape": section.shape,
                "material": section.normalized_material(),
                "bounds_m": [float(value) for value in section.bounds_min_m] + [float(value) for value in section.bounds_max_m],
                "entries": sum(bool(position < len(record.guide_section_entry_counts) and record.guide_section_entry_counts[position]) for record in records),
                "crossings": sum(record.guide_section_entry_counts[position] if position < len(record.guide_section_entry_counts) else 0 for record in records),
                "path_segments_inside": sum(record.guide_section_path_segment_counts[position] if position < len(record.guide_section_path_segment_counts) else 0 for record in records),
            })
        return {
            "scintillator_enabled": self.config.enabled,
            "scintillator_sampling_mode": self.config.sampling_mode,
            "scintillator_effective_max_bounces": int(self.config.max_bounces),
            "scintillator_requested_photons": self.requested_photons,
            "scintillator_generated_photons": len(records),
            "scintillator_transported_photons": len(records),
            "scintillator_source_rows": self.source_row_count,
            "scintillator_accepted_records": len(records),
            "scintillator_excluded_records": self.readout_mismatch_count,
            "scintillator_deposited_energy_MeV": self.deposited_energy_MeV,
            "scintillator_physical_expected_photons": self.physical_expected_photons,
            "scintillator_photons_reaching_guide": sum(record.reached_guide for record in records),
            "scintillator_photons_reaching_readout": len(readout),
            "scintillator_readout_arrivals": len(arrivals),
            "scintillator_detected_photons": len(detected),
            "scintillator_mylar_reflections": sum(record.mylar_reflections for record in records),
            "scintillator_mylar_losses": sum(record.mylar_losses for record in records),
            "scintillator_absorbed": sum(record.outcome == "absorbed" for record in records),
            "scintillator_escaped": sum(record.outcome.startswith("escaped") for record in records),
            "scintillator_max_bounces": sum(record.outcome == "max_bounces" for record in records),
            "scintillator_readout_alignment_mismatches": self.readout_mismatch_count,
            "scintillator_readout_endpoint_tolerance_um": self.readout_endpoint_tolerance_m * 1e6,
            "scintillator_active_to_guide_crossings": sum(record.active_to_guide_crossings for record in records),
            "scintillator_guide_to_readout_crossings": sum(record.guide_to_readout_crossings for record in records),
            "scintillator_interface_transmissions": sum(record.interface_transmissions for record in records),
            "scintillator_interface_reflections": sum(record.interface_reflections for record in records),
            "scintillator_tir_reflections": sum(record.tir_reflections for record in records),
            "scintillator_guide_sections": guide_section_summary,
            "scintillator_configured_guide_section_count": len(guide_section_summary),
            "scintillator_mean_arrival_ns": float(np.mean(arrivals)) if arrivals else None,
            "scintillator_earliest_arrival_ns": min(arrivals) if arrivals else None,
            "scintillator_track_length_m": self.track_length_m,
            "scintillator_stopped": self.stopped,
        }


@dataclass(frozen=True)
class Particle:
    """Charged particle at the quartz entry face.

    Inclination is measured from the global/local +z beam axis. Azimuth is
    measured around +z, with azimuth zero in the +x direction and positive
    azimuth toward +y.  ``entry_x_m`` and ``entry_y_m`` select the location
    on the transverse entry face; the entry z face is selected from the
    direction (+half_z for a -z track, -half_z for a +z track).
    """

    mass_kg: float
    kinetic_energy_J: float
    entry_x_m: float = 0.0
    entry_y_m: float = 0.0
    inclination_rad: float = 0.0
    azimuth_rad: float = 0.0
    charge_e: float = 1.0
    dEdx_J_per_m: Optional[float] = None

    @property
    def entry_point_xy(self) -> np.ndarray:
        return np.array([self.entry_x_m, self.entry_y_m], dtype=float)

    @property
    def direction(self) -> np.ndarray:
        st = math.sin(self.inclination_rad)
        return np.array(
            [st * math.cos(self.azimuth_rad), st * math.sin(self.azimuth_rad), math.cos(self.inclination_rad)],
            dtype=float,
        )


@dataclass(frozen=True)
class OpticalConfig:
    wavelength_min_m: float = 300e-9
    wavelength_max_m: float = 650e-9
    photon_count: int = 1000
    absorption_length_m: Optional[float] = None
    rayleigh_scattering_length_m: Optional[float] = None
    surface_roughness_rad: float = 0.0
    max_bounces: int = 100
    detector_acceptance_rad: float = math.radians(30.0)
    detector_quantum_efficiency: float = 0.25
    detector_collection_efficiency: float = 0.80
    detector_index: float = 1.0
    # Canonical assembly terminal state; defaults preserve the historical
    # +X detector behavior.
    positive_end_mode: str = "detector"
    positive_mirror_reflectivity: float = 1.0
    opposite_end_mode: str = "open"
    opposite_mirror_reflectivity: float = 1.0
    opposite_detector_acceptance_rad: Optional[float] = None
    opposite_detector_quantum_efficiency: Optional[float] = None
    opposite_detector_collection_efficiency: Optional[float] = None
    opposite_detector_index: Optional[float] = None
    wavelength_grid_size: int = 256
    sampling_mode: str = "fixed"
    # Detector active areas are defined in the terminal Y/Z portal plane.
    # Width/height and semi-axes are SI metres; offsets are relative to the
    # centre of the corresponding terminal portal.
    detector_shape: str = "rectangular"
    detector_active_width_m: float = 0.020
    detector_active_height_m: float = 0.010
    detector_active_radius_m: float = 0.005
    detector_active_semi_axis_a_m: float = 0.010
    detector_active_semi_axis_b_m: float = 0.005
    detector_active_offset_y_m: float = 0.0
    detector_active_offset_z_m: float = 0.0
    detector_active_rotation_rad: float = 0.0
    opposite_detector_shape: str = "rectangular"
    opposite_detector_active_width_m: float = 0.020
    opposite_detector_active_height_m: float = 0.010
    opposite_detector_active_radius_m: float = 0.005
    opposite_detector_active_semi_axis_a_m: float = 0.010
    opposite_detector_active_semi_axis_b_m: float = 0.005
    opposite_detector_active_offset_y_m: float = 0.0
    opposite_detector_active_offset_z_m: float = 0.0
    opposite_detector_active_rotation_rad: float = 0.0


@dataclass
class PhotonRecord:
    """One generated and transported photon.

    ``arrival_time_s`` is populated only for photons transmitted through a
    detector portal. ``impact_position_m`` is populated only for detected
    photons and is NaN for other outcomes. ``detector_id`` is 1 for +x and 2
    for -x, while ``detector_face`` is ``plus_x`` or ``minus_x``.
    Detector diagnostic integer fields count detector attempts, TIR encounters,
    interface reflections, and post-transmission detector losses. TIR photons
    are reflected back into the bar and can later terminate for another reason.
    ``path_points_m`` is for event visualization and is not written as a ROOT
    scalar branch.
    """

    photon_id: int
    wavelength_m: float
    emission_position_m: Tuple[float, float, float]
    final_position_m: Tuple[float, float, float]
    impact_position_m: Tuple[float, float, float]
    path_length_m: float
    optical_time_s: float
    emission_time_s: float
    arrival_time_s: float
    beta_at_emission: float
    bounces: int
    scatter_count: int
    outcome: str
    # One final reason is stored on every generated record.  It is derived
    # from the actual terminal outcome, never from viewer geometry.
    outcome_reason: str = ""
    detector_id: int = 0
    detector_face: str = ""
    detected: bool = False
    detector_transmitted: bool = False
    detector_attempts: int = 0
    tir_plus_x: int = 0
    tir_minus_x: int = 0
    detector_interface_reflections: int = 0
    detector_angle_rejected: int = 0
    detector_active_area_rejected: int = 0
    detector_collection_loss: int = 0
    detector_qe_loss: int = 0
    detector_valid_sensor_hit: bool = False
    detector_incidence_angle_rad: float = NAN
    effective_detector_acceptance_rad: float = NAN
    effective_detector_quantum_efficiency: float = NAN
    effective_detector_collection_efficiency: float = NAN
    detector_window_shape: str = ""
    detector_window_width_m: float = NAN
    detector_window_height_m: float = NAN
    mirror_count: int = 0
    generated: bool = True
    path_points_m: Tuple[Tuple[float, float, float], ...] = field(default_factory=tuple)
    # Scintillator observables are additive fields; existing Cherenkov
    # records retain their prior defaults and meanings.
    photon_kind: str = "cherenkov"
    mylar_reflections: int = 0
    mylar_losses: int = 0
    reached_guide: bool = False
    reached_readout: bool = False
    path_coordinate_frame: str = "geometry_m"
    # Explicit section/interface accounting for Scintillator diagnostics.
    active_to_guide_crossings: int = 0
    guide_to_readout_crossings: int = 0
    interface_transmissions: int = 0
    interface_reflections: int = 0
    tir_reflections: int = 0
    # Tuples are aligned to ScintillatorAssembly.guide_sections.  They are
    # written from actual internal interface transitions and path midpoints,
    # never from viewer geometry or a synthetic guide counter.
    guide_section_entry_counts: Tuple[int, ...] = field(default_factory=tuple)
    guide_section_path_segment_counts: Tuple[int, ...] = field(default_factory=tuple)
    readout_alignment_valid: bool = True
    readout_endpoint_distance_m: float = NAN
    # Physical entry time into the configured finite photosensor/readout
    # window, before angular acceptance, collection, or QE response.  This
    # optional field is last to preserve positional construction compatibility.
    readout_arrival_time_s: float = NAN

    @property
    def time_s(self) -> float:
        """Backward-compatible alias for detector arrival time."""

        return self.arrival_time_s


def detector_diagnostic_counts(records: Sequence[PhotonRecord]) -> Dict[str, int]:
    """Count detector stages without conflating crossings and valid hits.

    ``detector_face_crossings`` is the number of transmitted detector portal
    crossings; ``detector_valid_sensor_hits`` is the subset inside the finite
    configured window.  Angle, collection, and QE losses are post-window
    losses.  Fresnel/TIR/interface reflections are reported separately.
    """

    records = list(records)
    counts = {
        "detector_attempts": sum(record.detector_attempts for record in records),
        "tir_plus_x": sum(record.tir_plus_x for record in records),
        "tir_minus_x": sum(record.tir_minus_x for record in records),
        "detector_interface_reflections": sum(record.detector_interface_reflections for record in records),
        "detector_angle_rejected": sum(record.detector_angle_rejected for record in records),
        "detector_active_area_rejected": sum(record.detector_active_area_rejected for record in records),
        "detector_collection_loss": sum(record.detector_collection_loss for record in records),
        "detector_qe_loss": sum(record.detector_qe_loss for record in records),
        "detector_valid_sensor_hits": sum(record.detector_valid_sensor_hit for record in records),
    }
    counts["detector_losses"] = (
        counts["detector_angle_rejected"]
        + counts["detector_active_area_rejected"]
        + counts["detector_collection_loss"]
        + counts["detector_qe_loss"]
    )
    return counts


def detector_outcome_reason_counts(records: Sequence[PhotonRecord]) -> Dict[str, int]:
    """Return mutually exclusive final outcome-reason counts."""

    counts = Counter(str(record.outcome_reason or "unspecified") for record in records)
    return dict(sorted(counts.items()))


def _finite_number(value: object) -> Optional[float]:
    """Return a finite float or ``None`` without inventing a measurement."""

    try:
        number = float(value)
    except (TypeError, ValueError):
        return None
    return number if math.isfinite(number) else None


def readout_arrival_time_ns(record: object) -> Optional[float]:
    """Return the configured sensor-window crossing time in nanoseconds.

    New records use ``readout_arrival_time_s``.  The fallback is deliberately
    limited to legacy records that explicitly carry ``detector_valid_sensor_hit``;
    it preserves older in-memory runs whose ``arrival_time_s`` was assigned at
    the same physical portal crossing without treating a bare detected flag as
    an arrival.
    """

    if not bool(getattr(record, "detector_valid_sensor_hit", False)) or not bool(getattr(record, "readout_alignment_valid", True)):
        return None
    value = _finite_number(getattr(record, "readout_arrival_time_s", NAN))
    if value is None:
        value = _finite_number(getattr(record, "arrival_time_s", NAN))
    return None if value is None else value * 1e9


def _owner_records(result: object, detector: str) -> List[object]:
    """Return only records owned by one detector, never a mixed stream."""

    detector = str(detector).strip().lower()
    if detector == "cherenkov":
        return list(getattr(result, "photon_records", ()) or ()) if detector_run_mode_for_result(result) in {"cherenkov", "combined"} else []
    if detector == "scintillator":
        scintillator = getattr(result, "scintillator", None)
        return list(getattr(scintillator, "photon_records", ()) or ()) if detector_run_mode_for_result(result) in {"scintillator", "combined"} else []
    raise ValueError(f"unknown detector owner: {detector}")


def detector_run_mode_for_result(result: object) -> str:
    """Normalize a result's detector participation mode without GUI imports."""

    return normalize_detector_run_mode(getattr(result, "detector_mode", "cherenkov"))


def _analysis_generated_count(result: object, detector: str, records: Sequence[object]) -> int:
    if detector == "cherenkov":
        return int(getattr(result, "generated_photons", len(records)) or 0)
    # ``source_row_count`` is the requested/source sample and can exceed the
    # records retained after cooperative stop; Results reports actual records
    # generated/transported for the completed result.
    return int(len(records))


def _analysis_detector(result: object, detector: str) -> Dict[str, object]:
    records = _owner_records(result, detector)
    times = [time for record in records if (time := readout_arrival_time_ns(record)) is not None]
    detected = [record for record in records if bool(getattr(record, "detected", False))]
    detected_times = [time for record in detected if (time := readout_arrival_time_ns(record)) is not None]
    return {
        "included": bool(records) or detector in ({"cherenkov", "combined"} if detector_run_mode_for_result(result) == "combined" else {detector_run_mode_for_result(result)}),
        "generated": _analysis_generated_count(result, detector, records),
        "arrivals": len(times),
        "detected": len(detected),
        "arrival_times_ns": times,
        "detected_times_ns": detected_times,
        "mean_arrival_ns": float(np.mean(times)) if times else None,
        "median_arrival_ns": float(np.median(times)) if times else None,
    }


def _analysis_histograms(detectors: Mapping[str, Mapping[str, object]], bins: int = 30) -> Dict[str, object]:
    all_times = [float(value) for item in detectors.values() for value in item.get("arrival_times_ns", ())]
    if not all_times:
        return {"bin_edges_ns": [], "bin_widths_ns": [], "histograms": {}}
    minimum, maximum = min(all_times), max(all_times)
    count = max(1, min(int(bins), max(1, len(all_times))))
    if math.isclose(minimum, maximum, rel_tol=0.0, abs_tol=1e-12):
        padding = max(0.5, abs(minimum) * 0.05, 1e-6)
        edges = np.linspace(minimum - padding, maximum + padding, count + 1)
    else:
        edges = np.linspace(minimum, maximum, count + 1)
    output = {}
    for detector, item in detectors.items():
        values = np.asarray(item.get("arrival_times_ns", ()), dtype=float)
        counts, _ = np.histogram(values, bins=edges)
        density = counts.astype(float)
        integral = float(np.sum(density * np.diff(edges)))
        if integral > 0.0:
            density /= integral
        output[detector] = {
            "counts": counts.astype(int).tolist(),
            "normalized_density_per_ns": density.tolist(),
            "normalization": "unit area: histogram integral equals 1 (probability density per ns)",
        }
    return {
        "bin_edges_ns": edges.tolist(),
        "bin_widths_ns": np.diff(edges).tolist(),
        "histograms": output,
    }


def _analysis_from_results(results: Sequence[SimulationResult], selected_event_id: Optional[int] = None, bins: int = 30) -> Dict[str, object]:
    results = list(results)
    if not results:
        return {"schema": "readout-arrival-analysis-v1", "units": "ns", "detectors": {}, "histograms": {"bin_edges_ns": [], "bin_widths_ns": [], "histograms": {}}}
    if selected_event_id is None:
        selected_event_id = int(results[0].event_id)
    selected = next((result for result in results if int(result.event_id) == int(selected_event_id)), None)
    if selected is None:
        raise KeyError(f"event_id {selected_event_id} is not present in the completed results")
    detectors = {detector: _analysis_detector(selected, detector) for detector in ("cherenkov", "scintillator")}
    aggregate = {detector: {"included": any(_analysis_detector(result, detector)["included"] for result in results), "generated": sum(int(_analysis_detector(result, detector)["generated"]) for result in results), "arrivals": 0, "detected": 0, "arrival_times_ns": [], "detected_times_ns": []} for detector in detectors}
    for result in results:
        for detector in aggregate:
            item = _analysis_detector(result, detector)
            aggregate[detector]["arrivals"] += int(item["arrivals"])
            aggregate[detector]["detected"] += int(item["detected"])
            aggregate[detector]["arrival_times_ns"].extend(item["arrival_times_ns"])
            aggregate[detector]["detected_times_ns"].extend(item["detected_times_ns"])
    for item in aggregate.values():
        values = item["arrival_times_ns"]
        item["mean_arrival_ns"] = float(np.mean(values)) if values else None
        item["median_arrival_ns"] = float(np.median(values)) if values else None
    return {
        "schema": "readout-arrival-analysis-v1",
        "units": "ns",
        "time_source": "readout_arrival_time_s: physical crossing of the finite configured sensor/readout window before angular acceptance, collection, and QE",
        "normalization": "unit area: each aggregate detector histogram is a probability density per ns with integral 1",
        "selected_event_id": int(selected.event_id),
        "selected": {"event_id": int(selected.event_id), "detectors": detectors, "histograms": _analysis_histograms(detectors, bins=bins)},
        "aggregate": {"detectors": aggregate, "histograms": _analysis_histograms(aggregate, bins=bins)},
    }


def selected_event_arrival_analysis(results: Sequence[SimulationResult], event_id: int, bins: int = 30) -> Dict[str, object]:
    """Return selected-event detector-owned arrival counts and common-bin data."""

    return _analysis_from_results(results, selected_event_id=event_id, bins=bins)["selected"]


def aggregate_arrival_analysis(results: Sequence[SimulationResult], bins: int = 30) -> Dict[str, object]:
    """Return detector-owned all-event arrivals and unit-area distributions."""

    return _analysis_from_results(results, bins=bins)["aggregate"]


def readout_arrival_analysis(results: Sequence[SimulationResult], event_id: Optional[int] = None, bins: int = 30) -> Dict[str, object]:
    """Return both selected-event and all-event readout arrival analyses."""

    return _analysis_from_results(results, selected_event_id=event_id, bins=bins)


def readout_arrival_analysis_from_rows(
    rows: Sequence[Mapping[str, object]],
    event_summaries: Sequence[Mapping[str, object]] = (),
    detector_mode: str = "cherenkov",
    event_id: Optional[int] = None,
    bins: int = 30,
) -> Dict[str, object]:
    """Rebuild Results-tab analysis from a compact persisted arrival file."""

    mode = normalize_detector_run_mode(detector_mode)
    summaries = list(event_summaries or ())
    event_ids = [int(item.get("event_id", index)) for index, item in enumerate(summaries) if isinstance(item, Mapping)]
    if event_id is None:
        event_id = event_ids[0] if event_ids else (int(rows[0].get("event_id", 0)) if rows else 0)
    rows = [item for item in rows if isinstance(item, Mapping)]
    def included(detector: str) -> bool:
        return mode == "combined" or mode == detector
    def detector_info(target_event_id: Optional[int]) -> Dict[str, Dict[str, object]]:
        selected_summary = next((item for item in summaries if int(item.get("event_id", -1)) == int(target_event_id)), {})
        output = {}
        for detector in ("cherenkov", "scintillator"):
            detector_rows = [item for item in rows if str(item.get("detector", "")) == detector and (target_event_id is None or int(item.get("event_id", -1)) == int(target_event_id))]
            times = [float(item["arrival_time_ns"]) for item in detector_rows if _finite_number(item.get("arrival_time_ns")) is not None]
            generated_key = "generated_photons" if detector == "cherenkov" else "scintillator_generated_photons"
            generated = int(selected_summary.get(generated_key, 0) or 0)
            detected_count = sum(bool(item.get("detected", False)) for item in detector_rows)
            output[detector] = {
                "included": included(detector), "generated": generated, "arrivals": len(times),
                "detected": detected_count, "arrival_times_ns": times,
                "detected_times_ns": [float(item["arrival_time_ns"]) for item in detector_rows if bool(item.get("detected", False)) and _finite_number(item.get("arrival_time_ns")) is not None],
                "mean_arrival_ns": float(np.mean(times)) if times else None,
                "median_arrival_ns": float(np.median(times)) if times else None,
            }
        return output
    selected = detector_info(event_id)
    aggregate = {}
    for detector in ("cherenkov", "scintillator"):
        detector_rows = [item for item in rows if str(item.get("detector", "")) == detector]
        times = [float(item["arrival_time_ns"]) for item in detector_rows if _finite_number(item.get("arrival_time_ns")) is not None]
        generated_key = "generated_photons" if detector == "cherenkov" else "scintillator_generated_photons"
        aggregate[detector] = {
            "included": included(detector),
            "generated": sum(int(item.get(generated_key, 0) or 0) for item in summaries),
            "arrivals": len(times), "detected": sum(bool(item.get("detected", False)) for item in detector_rows),
            "arrival_times_ns": times,
            "detected_times_ns": [float(item["arrival_time_ns"]) for item in detector_rows if bool(item.get("detected", False)) and _finite_number(item.get("arrival_time_ns")) is not None],
            "mean_arrival_ns": float(np.mean(times)) if times else None,
            "median_arrival_ns": float(np.median(times)) if times else None,
        }
    return {
        "schema": "readout-arrival-analysis-v1", "units": "ns",
        "time_source": "readout_arrival_time_s: physical crossing of the finite configured sensor/readout window before angular acceptance, collection, and QE",
        "normalization": "unit area: each aggregate detector histogram is a probability density per ns with integral 1",
        "selected_event_id": int(event_id),
        "selected": {"event_id": int(event_id), "detectors": selected, "histograms": _analysis_histograms(selected, bins=bins)},
        "aggregate": {"detectors": aggregate, "histograms": _analysis_histograms(aggregate, bins=bins)},
    }


def readout_arrival_rows(results: Sequence[SimulationResult]) -> List[Dict[str, object]]:
    """Serialize only finite in-window readout arrivals, with detector owner."""

    rows = []
    for result in results:
        for detector in ("cherenkov", "scintillator"):
            for record in _owner_records(result, detector):
                time_ns = readout_arrival_time_ns(record)
                if time_ns is None:
                    continue
                rows.append({
                    "event_id": int(result.event_id), "detector": detector,
                    "arrival_time_ns": float(time_ns), "detected": bool(getattr(record, "detected", False)),
                    "outcome": str(getattr(record, "outcome", "")),
                })
    return rows


def write_readout_arrivals(path: str | os.PathLike[str], results: Sequence[SimulationResult], detector_mode: Optional[str] = None) -> str:
    """Atomically write the compact persisted readout-arrival dataset."""

    output_path = Path(path).expanduser().resolve()
    output_path.parent.mkdir(parents=True, exist_ok=True)
    payload = {
        "schema": "readout-arrivals-v1",
        "detector_mode": normalize_detector_run_mode(detector_mode or (results[0].detector_mode if results else "cherenkov")),
        "units": "ns",
        "time_definition": "arrival_time_ns is the physical crossing of the configured finite active sensor/readout window before angular acceptance, collection, or QE response",
        "rows": readout_arrival_rows(results),
    }
    fd, temporary_path = tempfile.mkstemp(prefix=output_path.name + ".", suffix=".tmp", dir=str(output_path.parent))
    try:
        with os.fdopen(fd, "w", encoding="utf-8", newline="") as handle:
            json.dump(payload, handle, allow_nan=False, separators=(",", ":"))
            handle.write("\n")
            handle.flush()
            os.fsync(handle.fileno())
        os.replace(temporary_path, output_path)
    except Exception:
        try:
            os.unlink(temporary_path)
        except OSError:
            pass
        raise
    return str(output_path)


def scintillator_readout_endpoint_distance_m(record: PhotonRecord, assembly: ScintillatorAssembly) -> float:
    """Distance from a stored path endpoint to its physical readout face.

    The readout is the actual terminal portal of the configured active/guide
    compound (+X guide output or -X active terminal).  This deliberately does not use a separate guessed sensor
    coordinate: the editor sensor window is coupled to this same terminal
    face.  A valid endpoint is on the face plane and inside its real portal.
    """

    if not record.path_points_m:
        return math.inf
    positive = record.detector_face == "plus_x"
    terminal_candidates = assembly.guide_sections if assembly.guide_sections else assembly.active_sections
    terminal = terminal_candidates[-1] if positive else terminal_candidates[0]
    endpoint = np.asarray(record.path_points_m[-1], dtype=float)
    if endpoint.shape != (3,) or not np.all(np.isfinite(endpoint)):
        return math.inf
    face_x = float(terminal.bounds_max_m[0] if positive else terminal.bounds_min_m[0])
    plane_distance = abs(float(endpoint[0]) - face_x)
    # The transport's active-area test already uses the exact portal shape.
    # Recheck containment here so this function remains the invariant authority
    # even when called on a manually constructed PhotonRecord.
    portal = assembly.transport_geometry
    face = "plus_x" if positive else "minus_x"
    portal_kind, portal_width, portal_height, portal_a, portal_b = detector_portal_dimensions(portal, face)
    in_portal = active_area_contains(
        endpoint, portal,
        OpticalConfig(
            detector_active_width_m=float(portal_width),
            detector_active_height_m=float(portal_height),
            detector_shape="rectangular" if portal_kind == "rectangle" else "elliptical",
            detector_active_semi_axis_a_m=float(portal_width) / 2.0,
            detector_active_semi_axis_b_m=float(portal_height) / 2.0,
        ),
        face,
    )
    return plane_distance if in_portal else math.inf


def enforce_scintillator_readout_invariant(
    records: Sequence[PhotonRecord], assembly: ScintillatorAssembly,
) -> Tuple[List[PhotonRecord], int, float]:
    """Make detected accounting agree with the stored terminal trajectory.

    A mismatch is downgraded to an explicit non-detected outcome.  It is never
    hidden by drawing a synthetic marker at the sensor.
    """

    tolerance = max(1e-12, assembly.transport_geometry.characteristic_length_m * 1e-8)
    repaired: List[PhotonRecord] = []
    mismatches = 0
    for record in records:
        distance = scintillator_readout_endpoint_distance_m(record, assembly) if record.detector_transmitted else NAN
        expected_faces = {
            face for face, state in (
                ("plus_x", assembly.terminal_plus_surface),
                ("minus_x", assembly.terminal_minus_surface),
            ) if _terminal_mode(state, "open") == "detector"
        }
        endpoint_valid = bool(record.detector_transmitted and record.detector_face in expected_faces and math.isfinite(distance) and distance <= tolerance)
        detected_valid = bool(record.detected and endpoint_valid)
        if record.detected and not detected_valid:
            mismatches += 1
            repaired.append(replace(
                record,
                detected=False,
                outcome="readout_alignment_mismatch",
                outcome_reason="readout_alignment_mismatch",
                detector_valid_sensor_hit=False,
                readout_arrival_time_s=NAN,
                impact_position_m=_nan_position(),
                readout_alignment_valid=False,
                readout_endpoint_distance_m=float(distance),
            ))
        else:
            repaired.append(replace(
                record,
                readout_alignment_valid=(not record.detector_transmitted) or endpoint_valid,
                readout_endpoint_distance_m=float(distance),
            ))
    return repaired, mismatches, tolerance


@dataclass
class SimulationResult:
    geometry: Union[BarGeometry, SolidGeometry, CompoundGeometry]
    particle: Particle
    optical: OpticalConfig
    event_id: int = 0
    photon_records: List[PhotonRecord] = field(default_factory=list)
    track_length_m: float = 0.0
    beta_start: float = 0.0
    beta_end: float = 0.0
    threshold_start: bool = False
    threshold_end: bool = False
    generated_photons: int = 0
    physical_expected_photons: float = 0.0
    stopped: bool = False
    notes: List[str] = field(default_factory=list)
    scintillator: Optional[ScintillatorResult] = None
    detector_mode: str = "cherenkov"

    def summary(self) -> Dict[str, object]:
        outcomes = Counter(record.outcome for record in self.photon_records)
        detected = [record for record in self.photon_records if record.detected]
        crossings = [record for record in self.photon_records if record.detector_transmitted]
        plus_crossings = [record for record in crossings if record.detector_face == "plus_x"]
        minus_crossings = [record for record in crossings if record.detector_face == "minus_x"]
        diagnostic_counts = detector_diagnostic_counts(self.photon_records)
        total = len(self.photon_records)
        participates_cherenkov = self.detector_mode in {"cherenkov", "combined"}
        detector_response = {}
        for face, mode, acceptance, qe, collection, index in (
            ("plus_x", self.optical.positive_end_mode, self.optical.detector_acceptance_rad,
             self.optical.detector_quantum_efficiency, self.optical.detector_collection_efficiency,
             self.optical.detector_index),
            ("minus_x", self.optical.opposite_end_mode, *_opposite_detector_settings(self.optical)),
        ):
            if mode != "detector":
                continue
            shape, width, height, _radius, _a, _b, _oy, _oz, _rotation = _active_area_values(self.optical, face, self.geometry)
            detector_response[face] = {
                "acceptance_deg": math.degrees(float(acceptance)),
                "quantum_efficiency": float(qe), "collection_efficiency": float(collection),
                "coupling_index": float(index), "window_shape": shape,
                "window_width_mm": width * 1e3, "window_height_mm": height * 1e3,
            }
        summary = {
            "event_id": self.event_id,
            "detector_mode": self.detector_mode,
            "track_length_m": self.track_length_m,
            "track_length_mm": self.track_length_m * 1e3,
            "beta_start": self.beta_start,
            "beta_end": self.beta_end,
            "threshold_start": self.threshold_start,
            "threshold_end": self.threshold_end,
            "requested_photons": self.optical.photon_count if participates_cherenkov else 0,
            "requested_sample_photons": self.optical.photon_count if participates_cherenkov else 0,
            "requested_photons_ignored": participates_cherenkov and self.optical.sampling_mode == "poisson",
            "sampling_mode": self.optical.sampling_mode,
            "cherenkov_effective_max_bounces": int(self.optical.max_bounces) if participates_cherenkov else None,
            "sampling_description": (
                "fixed Monte Carlo test sample"
                if self.optical.sampling_mode == "fixed"
                else "physical Poisson sampling from the Frank-Tamm yield"
            ),
            "generated_photons": self.generated_photons,
            "stopped": self.stopped,
            "physical_expected_photons": self.physical_expected_photons if participates_cherenkov else 0.0,
            "poisson_mean": self.physical_expected_photons if participates_cherenkov and self.optical.sampling_mode == "poisson" else None,
            "transported_photons": total,
            "detected_photons": len(detected),
            "detected_plus_x": sum(record.detector_face == "plus_x" for record in detected),
            "detected_minus_x": sum(record.detector_face == "minus_x" for record in detected),
            "detector_face_crossings": len(crossings),
            "valid_sensor_hits": sum(record.detector_valid_sensor_hit for record in self.photon_records),
            "readout_arrivals": sum(readout_arrival_time_ns(record) is not None for record in self.photon_records),
            "detector_plus_x_crossings": len(plus_crossings),
            "detector_minus_x_crossings": len(minus_crossings),
            "detection_efficiency_among_crossings": len(detected) / len(crossings) if crossings else 0.0,
            "detection_efficiency_plus_x_crossings": sum(record.detected for record in plus_crossings) / len(plus_crossings) if plus_crossings else 0.0,
            "detection_efficiency_minus_x_crossings": sum(record.detected for record in minus_crossings) / len(minus_crossings) if minus_crossings else 0.0,
            "detector_hits": {
                "plus_x": sum(record.detector_face == "plus_x" for record in detected),
                "minus_x": sum(record.detector_face == "minus_x" for record in detected),
            },
            "detection_fraction_of_generated": (len(detected) / self.generated_photons if self.generated_photons else 0.0),
            "generation_fraction": (
                self.generated_photons / self.physical_expected_photons
                if participates_cherenkov and self.optical.sampling_mode == "poisson" and self.physical_expected_photons > 0.0
                else (self.generated_photons / self.optical.photon_count if self.optical.photon_count > 0 else 0.0)
            ) if participates_cherenkov else 0.0,
            "outcomes": dict(sorted(outcomes.items())),
            "outcome_reasons": detector_outcome_reason_counts(self.photon_records),
            "effective_detector_response": detector_response,
            "mean_detected_wavelength_nm": (float(np.mean([p.wavelength_m for p in detected])) * 1e9 if detected else None),
            "mean_detected_path_m": (float(np.mean([p.path_length_m for p in detected])) if detected else None),
            "mean_detected_path_mm": (float(np.mean([p.path_length_m for p in detected])) * 1e3 if detected else None),
            "mean_detected_time_ns": (float(np.mean([p.arrival_time_s for p in detected])) * 1e9 if detected else None),
            "mean_readout_arrival_ns": (float(np.mean([readout_arrival_time_ns(p) for p in self.photon_records if readout_arrival_time_ns(p) is not None])) if any(readout_arrival_time_ns(p) is not None for p in self.photon_records) else None),
            "mean_detected_bounces": (float(np.mean([p.bounces for p in detected])) if detected else None),
            "notes": list(self.notes),
        }
        summary.update(diagnostic_counts)
        if self.scintillator is not None:
            summary.update(self.scintillator.summary())
            scint_arrivals = [record.arrival_time_s for record in self.scintillator.photon_records if record.detected and math.isfinite(record.arrival_time_s)]
            summary["scintillator_cherenkov_coincidence_ns"] = (
                float(np.min(scint_arrivals)) * 1e9 - float(np.min([record.arrival_time_s for record in detected if math.isfinite(record.arrival_time_s)])) * 1e9
                if scint_arrivals and any(math.isfinite(record.arrival_time_s) for record in detected) else None
            )
        return summary


def fused_silica_index(wavelength_m: float | np.ndarray) -> float | np.ndarray:
    """Return fused-silica refractive index from the Malitson Sellmeier fit."""

    wavelength_um = np.asarray(wavelength_m, dtype=float) * 1e6
    if np.any(wavelength_um <= 0.0):
        raise ValueError("wavelength must be positive")
    b1, b2, b3 = 0.6961663, 0.4079426, 0.8974794
    c1, c2, c3 = 0.0684043**2, 0.1162414**2, 9.896161**2
    lam2 = wavelength_um * wavelength_um
    n2 = 1.0 + b1 * lam2 / (lam2 - c1) + b2 * lam2 / (lam2 - c2) + b3 * lam2 / (lam2 - c3)
    if np.any(n2 <= 0.0):
        raise ValueError("Sellmeier equation produced a non-positive n^2")
    value = np.sqrt(n2)
    return float(value) if value.ndim == 0 else value


def particle_beta(mass_kg: float, kinetic_energy_J: float) -> float:
    if mass_kg <= 0.0 or kinetic_energy_J < 0.0:
        raise ValueError("mass must be positive and kinetic energy must be non-negative")
    gamma = 1.0 + kinetic_energy_J / (mass_kg * C_LIGHT**2)
    return math.sqrt(max(0.0, 1.0 - 1.0 / (gamma * gamma)))


def direction_from_angles(inclination_rad: float, azimuth_rad: float) -> np.ndarray:
    st = math.sin(inclination_rad)
    value = np.array([st * math.cos(azimuth_rad), st * math.sin(azimuth_rad), math.cos(inclination_rad)], dtype=float)
    return value / np.linalg.norm(value)


def frank_tamm_spectrum(wavelength_m: np.ndarray, beta: float, charge_e: float = 1.0, refractive_index=None) -> np.ndarray:
    """Return dN/(dx d(lambda)) in photons per meter per meter."""

    wavelength_m = np.asarray(wavelength_m, dtype=float)
    index_function = fused_silica_index if refractive_index is None else refractive_index
    n = np.asarray(index_function(wavelength_m), dtype=float)
    allowed = beta * n > 1.0
    spectrum = np.zeros_like(wavelength_m)
    spectrum[allowed] = (
        2.0 * math.pi * (charge_e**2 * ELEMENTARY_CHARGE**2 / (4.0 * math.pi * 8.8541878128e-12))
        / (6.62607015e-34 * C_LIGHT)
        * (1.0 - 1.0 / (beta * beta * n[allowed] * n[allowed]))
        / (wavelength_m[allowed] ** 2)
    )
    return spectrum


def cherenkov_threshold(beta: float, wavelength_m: float = 450e-9, refractive_index=None) -> bool:
    index_function = fused_silica_index if refractive_index is None else refractive_index
    return bool(beta * float(index_function(wavelength_m)) > 1.0)


ACTIVE_AREA_SHAPES = {"rectangular", "circular", "elliptical"}


def normalize_active_area_shape(shape: str) -> str:
    """Normalize the public detector active-area shape names."""

    value = str(shape).strip().lower()
    aliases = {"rectangle": "rectangular", "box": "rectangular", "circle": "circular", "ellipse": "elliptical"}
    value = aliases.get(value, value)
    if value not in ACTIVE_AREA_SHAPES:
        raise ValueError("detector active-area shape must be rectangular, circular, or elliptical")
    return value


def _active_area_values(optical: OpticalConfig, detector_label: str, geometry=None) -> Tuple[str, float, float, float, float, float, float, float, float]:
    # A current complete assembly carries a separate photosensor/readout
    # component.  Its window, rather than the obsolete generic active-area
    # fields, is the transport aperture for that assembly.
    for item in getattr(geometry, "sensor_components", ()):
        if str(item.get("face", "")).strip().lower() != str(detector_label).strip().lower():
            continue
        dimensions = item.get("dimensions_mm", {})
        if not isinstance(dimensions, Mapping):
            continue
        if not all(float(dimensions.get(key, 0.0)) > 0.0 for key in ("width", "height")):
            # A malformed canonical component must not become a zero-sized
            # active window. Complete assemblies carry validated dimensions;
            # legacy callers retain the fallback below.
            continue
        width = float(dimensions.get("width", 0.0)) * 1e-3
        height = float(dimensions.get("height", 0.0)) * 1e-3
        shape = normalize_active_area_shape(str(item.get("shape", "rectangular")))
        radius_value = item.get("readout_radius_mm")
        axes_value = item.get("readout_semiaxes_mm")
        radius = (float(radius_value) * 1e-3 if radius_value not in (None, "") else min(width, height) / 2.0)
        if isinstance(axes_value, Mapping):
            semi_a = float(axes_value.get("a", width * 500.0)) * 1e-3
            semi_b = float(axes_value.get("b", height * 500.0)) * 1e-3
        else:
            semi_a, semi_b = width / 2.0, height / 2.0
        return (
            shape, width, height, radius,
            semi_a, semi_b, 0.0, 0.0, 0.0,
        )
    if detector_label == "plus_x":
        prefix = "detector_"
    elif detector_label == "minus_x":
        prefix = "opposite_detector_"
    else:
        raise ValueError("detector label must be plus_x or minus_x")
    return (
        normalize_active_area_shape(getattr(optical, prefix + "shape")),
        float(getattr(optical, prefix + "active_width_m")),
        float(getattr(optical, prefix + "active_height_m")),
        float(getattr(optical, prefix + "active_radius_m")),
        float(getattr(optical, prefix + "active_semi_axis_a_m")),
        float(getattr(optical, prefix + "active_semi_axis_b_m")),
        float(getattr(optical, prefix + "active_offset_y_m")),
        float(getattr(optical, prefix + "active_offset_z_m")),
        float(getattr(optical, prefix + "active_rotation_rad")),
    )


def _validate_active_area_values(optical: OpticalConfig, detector_label: str) -> None:
    shape, width, height, radius, semi_a, semi_b, offset_y, offset_z, rotation = _active_area_values(optical, detector_label)
    label = "+x" if detector_label == "plus_x" else "-x"
    for name, value in (("offset_y", offset_y), ("offset_z", offset_z), ("rotation", rotation)):
        if not math.isfinite(value):
            raise ValueError(f"{label} detector active-area {name} must be finite")
    dimensions = (("width", width), ("height", height), ("radius", radius), ("semi-axis a", semi_a), ("semi-axis b", semi_b))
    for name, value in dimensions:
        if not math.isfinite(value) or value <= 0.0:
            raise ValueError(f"{label} detector active-area {name} must be positive")


def _portal_cross_section(geometry: Union[BarGeometry, SolidGeometry, CompoundGeometry], positive_x: bool) -> Optional[Tuple[str, float, float, float, float]]:
    """Return portal kind, Y/Z dimensions and centre Y/Z for a terminal cap."""

    if isinstance(geometry, BarGeometry):
        return "rectangle", geometry.width_y_m, geometry.thickness_z_m, 0.0, 0.0
    if isinstance(geometry, CompoundGeometry):
        section = geometry.sections[-1] if positive_x else geometry.sections[0]
    elif isinstance(geometry, SolidGeometry):
        section = geometry
    else:
        return None
    if section.shape == "box":
        width, height, kind = section.width_y_m, section.thickness_z_m, "rectangle"
    elif section.shape == "tapered_rectangular_prism":
        width = section.width_y_end_m if positive_x else section.width_y_start_m
        height = section.thickness_z_end_m if positive_x else section.thickness_z_start_m
        kind = "rectangle"
    elif section.axis == "X" and section.shape == "cylinder":
        width = height = 2.0 * section.radius_m
        kind = "ellipse"
    elif section.axis == "X" and section.shape == "elliptical_cylinder":
        width, height, kind = 2.0 * section.semi_axis_a_m, 2.0 * section.semi_axis_b_m, "ellipse"
    elif section.axis == "X" and section.shape == "truncated_cone":
        radius = section.radius_end_m if positive_x else section.radius_start_m
        width = height = 2.0 * radius
        kind = "ellipse"
    else:
        return None
    return kind, float(width), float(height), float(section.center_y_m), float(section.center_z_m)


def detector_portal_dimensions(geometry: Union[BarGeometry, SolidGeometry, CompoundGeometry], face: str = "plus_x") -> Optional[Tuple[str, float, float, float, float]]:
    """Public portal helper: kind, width, height, centre-y, centre-z in metres."""

    if face not in {"plus_x", "minus_x"}:
        raise ValueError("face must be plus_x or minus_x")
    return _portal_cross_section(geometry, face == "plus_x")


def _active_boundary_points(shape: str, width: float, height: float, radius: float, semi_a: float, semi_b: float, rotation: float, samples: int = 720) -> np.ndarray:
    if shape == "rectangular":
        points = np.array([[-width / 2, -height / 2], [width / 2, -height / 2], [width / 2, height / 2], [-width / 2, height / 2]], dtype=float)
    else:
        angles = np.linspace(0.0, 2.0 * math.pi, max(32, samples), endpoint=False)
        if shape == "circular":
            points = np.column_stack((radius * np.cos(angles), radius * np.sin(angles)))
        else:
            points = np.column_stack((semi_a * np.cos(angles), semi_b * np.sin(angles)))
    c, s = math.cos(rotation), math.sin(rotation)
    rotation_matrix = np.array([[c, -s], [s, c]], dtype=float)
    return points @ rotation_matrix.T


def active_area_contains(position: Sequence[float], geometry: Union[BarGeometry, SolidGeometry, CompoundGeometry], optical: OpticalConfig, face: str = "plus_x") -> bool:
    """Return whether a transmitted portal impact lies inside its active area."""

    portal = _portal_cross_section(geometry, face == "plus_x")
    if portal is None:
        return False
    shape, width, height, radius, semi_a, semi_b, offset_y, offset_z, rotation = _active_area_values(optical, face, geometry)
    _, _, _, center_y, center_z = portal
    dy = float(position[1]) - center_y - offset_y
    dz = float(position[2]) - center_z - offset_z
    c, s = math.cos(rotation), math.sin(rotation)
    u, v = c * dy + s * dz, -s * dy + c * dz
    if shape == "rectangular":
        return abs(u) <= width / 2.0 + 1e-12 and abs(v) <= height / 2.0 + 1e-12
    if shape == "circular":
        return u * u + v * v <= radius * radius + 1e-24
    return (u / semi_a) ** 2 + (v / semi_b) ** 2 <= 1.0 + 1e-12


def _validate_active_area_against_portal(geometry: Union[BarGeometry, SolidGeometry, CompoundGeometry], optical: OpticalConfig, face: str) -> None:
    portal = _portal_cross_section(geometry, face == "plus_x")
    if portal is None:
        raise ValueError(f"{face} detector requires a planar X portal; the active solid has no compatible portal")
    _validate_active_area_values(optical, face)
    kind, portal_width, portal_height, center_y, center_z = portal
    shape, width, height, radius, semi_a, semi_b, offset_y, offset_z, rotation = _active_area_values(optical, face)
    boundary = _active_boundary_points(shape, width, height, radius, semi_a, semi_b, rotation)
    boundary[:, 0] += offset_y
    boundary[:, 1] += offset_z
    if kind == "rectangle":
        if np.any(np.abs(boundary[:, 0]) > portal_width / 2.0 + 1e-12) or np.any(np.abs(boundary[:, 1]) > portal_height / 2.0 + 1e-12):
            raise ValueError(f"{face} detector active area does not fit inside the rectangular portal ({portal_width * 1e3:g} x {portal_height * 1e3:g} mm)")
    else:
        normalized = (boundary[:, 0] / (portal_width / 2.0)) ** 2 + (boundary[:, 1] / (portal_height / 2.0)) ** 2
        if np.any(normalized > 1.0 + 1e-9):
            raise ValueError(f"{face} detector active area does not fit inside the elliptical portal ({portal_width * 1e3:g} x {portal_height * 1e3:g} mm)")
    # Offsets are relative to the actual terminal section centre.  Keep the
    # explicit variables above in the check so a non-zero portal centre is not
    # accidentally treated as a global origin.
    _ = center_y, center_z


def validate_detector_active_areas(geometry: Union[BarGeometry, SolidGeometry, CompoundGeometry], optical: OpticalConfig) -> None:
    """Validate configured active areas against all detector portals in use."""

    if optical.positive_end_mode == "detector":
        _validate_active_area_against_portal(geometry, optical, "plus_x")
    if optical.opposite_end_mode == "detector":
        _validate_active_area_against_portal(geometry, optical, "minus_x")


def _validate_inputs(geometry: Union[BarGeometry, SolidGeometry, CompoundGeometry], particle: Particle, optical: OpticalConfig) -> None:
    if isinstance(geometry, (SolidGeometry, CompoundGeometry)):
        geometry.validate()
        try:
            geometry.particle_track(particle)
        except ValueError as exc:
            raise ValueError(f"particle track does not intersect the solid: {exc}") from exc
    elif not isinstance(geometry, BarGeometry):
        raise TypeError("geometry must be BarGeometry, SolidGeometry, or CompoundGeometry")
    elif geometry.length_x_m <= 0 or geometry.width_y_m <= 0 or geometry.thickness_z_m <= 0:
        raise ValueError("bar dimensions must be positive")
    if not (-math.pi / 2.0 < particle.inclination_rad < math.pi / 2.0):
        raise ValueError("inclination_rad must be in (-pi/2, pi/2)")
    if isinstance(geometry, BarGeometry) and (abs(particle.entry_x_m) > geometry.half_x or abs(particle.entry_y_m) > geometry.half_y):
        raise ValueError("entry point is outside the z=+half_z face")
    if optical.wavelength_min_m <= 0 or optical.wavelength_max_m <= optical.wavelength_min_m:
        raise ValueError("wavelength range is invalid")
    if optical.photon_count < 0 or optical.max_bounces < 0 or optical.wavelength_grid_size < 16:
        raise ValueError("photon_count, max_bounces, and wavelength_grid_size are invalid")
    if optical.sampling_mode not in {"fixed", "poisson"}:
        raise ValueError("sampling_mode must be fixed or poisson")
    for name, value in (("detector_quantum_efficiency", optical.detector_quantum_efficiency), ("detector_collection_efficiency", optical.detector_collection_efficiency)):
        if not 0.0 <= value <= 1.0:
            raise ValueError(f"{name} must be between 0 and 1")
    for name, value in (("absorption_length_m", optical.absorption_length_m), ("rayleigh_scattering_length_m", optical.rayleigh_scattering_length_m)):
        if value is not None and value <= 0.0:
            raise ValueError(f"{name} must be positive when supplied")
    if optical.detector_index <= 0.0:
        raise ValueError("detector_index must be positive")
    if optical.positive_end_mode not in {"open", "mirror", "detector"}:
        raise ValueError("positive_end_mode must be open, mirror, or detector")
    if optical.opposite_end_mode not in {"open", "mirror", "detector"}:
        raise ValueError("opposite_end_mode must be open, mirror, or detector")
    if not 0.0 <= optical.positive_mirror_reflectivity <= 1.0:
        raise ValueError("positive_mirror_reflectivity must be between 0 and 1")
    if not 0.0 <= optical.opposite_mirror_reflectivity <= 1.0:
        raise ValueError("opposite_mirror_reflectivity must be between 0 and 1")
    if optical.opposite_detector_acceptance_rad is not None and not 0.0 <= optical.opposite_detector_acceptance_rad <= math.pi:
        raise ValueError("opposite_detector_acceptance_rad must be between 0 and pi")
    for name, value in (
        ("opposite_detector_quantum_efficiency", optical.opposite_detector_quantum_efficiency),
        ("opposite_detector_collection_efficiency", optical.opposite_detector_collection_efficiency),
    ):
        if value is not None and not 0.0 <= value <= 1.0:
            raise ValueError(f"{name} must be between 0 and 1")
    if optical.opposite_detector_index is not None and optical.opposite_detector_index <= 0.0:
        raise ValueError("opposite_detector_index must be positive")
    # Canonical complete assemblies validate their separate photosensor
    # against the actual terminal face before a run.  A second call to the
    # legacy generic portal validator would resurrect the obsolete 10 x 5 mm
    # error path.  Direct legacy geometry callers retain that validation.
    if not getattr(geometry, "sensor_components", ()):
        validate_detector_active_areas(geometry, optical)
    if particle.dEdx_J_per_m is not None and particle.dEdx_J_per_m < 0.0:
        raise ValueError("dEdx_J_per_m must be non-negative when supplied")


def _opposite_detector_settings(optical: OpticalConfig) -> Tuple[float, float, float, float]:
    """Return acceptance, QE, collection efficiency, and index for -x."""

    return (
        optical.opposite_detector_acceptance_rad if optical.opposite_detector_acceptance_rad is not None else optical.detector_acceptance_rad,
        optical.opposite_detector_quantum_efficiency if optical.opposite_detector_quantum_efficiency is not None else optical.detector_quantum_efficiency,
        optical.opposite_detector_collection_efficiency if optical.opposite_detector_collection_efficiency is not None else optical.detector_collection_efficiency,
        optical.opposite_detector_index if optical.opposite_detector_index is not None else optical.detector_index,
    )


def _entry_point(geometry: Union[BarGeometry, SolidGeometry, CompoundGeometry], particle: Particle) -> np.ndarray:
    if isinstance(geometry, (SolidGeometry, CompoundGeometry)):
        return geometry.particle_track(particle)[0]
    entry_z = -geometry.half_z if particle.direction[2] >= 0.0 else geometry.half_z
    return np.array([particle.entry_x_m, particle.entry_y_m, entry_z], dtype=float)


def _track_length(geometry: Union[BarGeometry, SolidGeometry, CompoundGeometry], particle: Particle) -> float:
    if isinstance(geometry, (SolidGeometry, CompoundGeometry)):
        return geometry.particle_track(particle)[1]
    point = _entry_point(geometry, particle)
    direction = particle.direction
    distances: List[float] = []
    for coordinate, direction_component, low, high in (
        (point[0], direction[0], -geometry.half_x, geometry.half_x),
        (point[1], direction[1], -geometry.half_y, geometry.half_y),
        (point[2], direction[2], -geometry.half_z, geometry.half_z),
    ):
        if direction_component > 1e-15:
            distances.append((high - coordinate) / direction_component)
        elif direction_component < -1e-15:
            distances.append((low - coordinate) / direction_component)
    positive = [value for value in distances if value > 1e-14]
    if not positive:
        raise ValueError("particle direction does not intersect the bar")
    return max(0.0, min(positive))


def _particle_beta_at_distance(particle: Particle, distance_m: float) -> float:
    energy = particle.kinetic_energy_J - (particle.dEdx_J_per_m or 0.0) * distance_m
    return particle_beta(particle.mass_kg, max(0.0, energy))


def _particle_time_to_distance(particle: Particle, distance_m: float) -> float:
    if distance_m <= 0.0:
        return 0.0
    steps = max(16, int(math.ceil(distance_m / 1e-4)))
    grid = np.linspace(0.0, distance_m, steps + 1)
    beta = np.array([_particle_beta_at_distance(particle, float(s)) for s in grid])
    inv_speed = 1.0 / np.maximum(beta, 1e-15) / C_LIGHT
    integral = np.trapezoid(inv_speed, grid) if hasattr(np, "trapezoid") else np.trapz(inv_speed, grid)
    return float(integral)


def _orthonormal_basis(direction: np.ndarray) -> Tuple[np.ndarray, np.ndarray]:
    reference = np.array([0.0, 0.0, 1.0]) if abs(direction[2]) < 0.9 else np.array([1.0, 0.0, 0.0])
    first = np.cross(direction, reference)
    first /= np.linalg.norm(first)
    second = np.cross(direction, first)
    return first, second


def _sample_unit_vector(rng: np.random.Generator) -> np.ndarray:
    vector = rng.normal(size=3)
    return vector / np.linalg.norm(vector)


def _fresnel_unpolarized(cos_i: float, n1: float, n2: float) -> Tuple[float, float]:
    cos_i = min(1.0, max(0.0, cos_i))
    sin_i2 = max(0.0, 1.0 - cos_i * cos_i)
    sin_t2 = (n1 / n2) ** 2 * sin_i2
    if sin_t2 >= 1.0:
        return 1.0, 0.0
    cos_t = math.sqrt(max(0.0, 1.0 - sin_t2))
    rs = ((n1 * cos_i - n2 * cos_t) / (n1 * cos_i + n2 * cos_t)) ** 2
    rp = ((n1 * cos_t - n2 * cos_i) / (n1 * cos_t + n2 * cos_i)) ** 2
    return 0.5 * (rs + rp), cos_t


def _refracted_direction(direction: np.ndarray, outward: np.ndarray, n1: float, n2: float, cos_t: float) -> np.ndarray:
    """Return the transmitted ray for a normal pointing from n1 to n2."""

    direction = np.asarray(direction, dtype=float)
    outward = np.asarray(outward, dtype=float)
    cos_i = max(0.0, float(np.dot(direction, outward)))
    eta = n1 / n2
    transmitted = eta * direction + (eta * cos_i - cos_t) * outward
    norm = np.linalg.norm(transmitted)
    if norm <= 0.0 or not np.all(np.isfinite(transmitted)):
        raise ValueError("invalid refracted ray direction")
    return transmitted / norm


def _rough_surface_normal(outward: np.ndarray, roughness_rad: float, rng: np.random.Generator) -> np.ndarray:
    if roughness_rad <= 0.0:
        return outward
    tangent_a, tangent_b = _orthonormal_basis(outward)
    perturbed = outward + roughness_rad * (rng.normal() * tangent_a + rng.normal() * tangent_b)
    perturbed /= np.linalg.norm(perturbed)
    return perturbed if np.dot(perturbed, outward) > 0.0 else outward


def _boundary_distance(point: np.ndarray, direction: np.ndarray, geometry: Union[BarGeometry, SolidGeometry, CompoundGeometry]) -> Tuple[float, str, np.ndarray]:
    if isinstance(geometry, (SolidGeometry, CompoundGeometry)):
        return geometry.ray_intersection(point, direction)
    candidates: List[Tuple[float, str, np.ndarray]] = []
    faces = (
        ("x+", geometry.half_x, 0, np.array([1.0, 0.0, 0.0])),
        ("x-", -geometry.half_x, 0, np.array([-1.0, 0.0, 0.0])),
        ("y+", geometry.half_y, 1, np.array([0.0, 1.0, 0.0])),
        ("y-", -geometry.half_y, 1, np.array([0.0, -1.0, 0.0])),
        ("z+", geometry.half_z, 2, np.array([0.0, 0.0, 1.0])),
        ("z-", -geometry.half_z, 2, np.array([0.0, 0.0, -1.0])),
    )
    for name, boundary, axis, normal in faces:
        component = direction[axis]
        distance = (boundary - point[axis]) / component if abs(component) > 1e-15 else math.inf
        if distance > 1e-12:
            candidates.append((distance, name, normal))
    if not candidates:
        raise RuntimeError("photon has no forward boundary intersection")
    return min(candidates, key=lambda item: item[0])


def _segment_event_distance(boundary_distance: float, absorption_length_m: Optional[float], scattering_length_m: Optional[float], rng: np.random.Generator) -> Tuple[float, Optional[str]]:
    distances = [(boundary_distance, None)]
    if absorption_length_m is not None:
        distances.append((float(rng.exponential(absorption_length_m)), "absorbed"))
    if scattering_length_m is not None:
        distances.append((float(rng.exponential(scattering_length_m)), "scattered"))
    return min(distances, key=lambda item: item[0])


def _nan_position() -> Tuple[float, float, float]:
    return (NAN, NAN, NAN)


def _transport_epsilon_for_geometry(geometry: Union[BarGeometry, SolidGeometry, CompoundGeometry]) -> float:
    """Return a post-hit nudge larger than the boundary classification tolerance."""

    if isinstance(geometry, (SolidGeometry, CompoundGeometry)):
        return geometry.transport_epsilon_m
    characteristic = max(geometry.length_x_m, geometry.width_y_m, geometry.thickness_z_m)
    return max(1e-12, min(1e-7, characteristic * 1e-7))


def _escaped_face_outcome(face: Optional[str]) -> str:
    """Map a terminal surface to the established escaped outcome names."""

    if not face:
        return "escaped"
    return f"escaped_{face}"


def _final_outcome_reason(outcome: str, detector_transmitted: bool, valid_sensor_hit: bool) -> str:
    """Map one terminal outcome to a stable, mutually exclusive reason."""

    text = str(outcome or "")
    if text in {"detected_plus_x", "detected_minus_x"}:
        return "detected"
    if text.endswith("_active_area_rejected"):
        return "sensor_window_miss"
    if text.endswith("_angle_rejected"):
        return "acceptance_reject"
    if text.endswith("_collection_loss"):
        return "collection_reject"
    if text.endswith("_qe_loss"):
        return "qe_reject"
    if text.startswith("tir_"):
        return "fresnel_tir"
    if detector_transmitted and not valid_sensor_hit:
        return "sensor_window_miss"
    if text == "absorbed":
        return "absorption"
    if text == "max_bounces":
        return "max_bounces"
    if text.startswith("mylar_loss"):
        return "mylar_absorption"
    if text.startswith("escaped"):
        return "boundary_escape"
    if text == "readout_alignment_mismatch":
        return "readout_alignment_mismatch"
    return text or "other"


def _transport_photon(photon_id: int, wavelength_m: float, emission_position: np.ndarray, initial_direction: np.ndarray, emission_time_s: float, beta_at_emission: float, geometry: Union[BarGeometry, SolidGeometry, CompoundGeometry], optical: OpticalConfig, rng: np.random.Generator, stop_event=None, photon_kind: str = "", status_callback: Optional[Callable[[str], None]] = None, coordinate_frame: str = "geometry_m", guide_section_indices: Sequence[int] = (), active_section_indices: Sequence[int] = ()) -> PhotonRecord:
    coordinate_frame = normalize_coordinate_frame(coordinate_frame)
    position = emission_position.astype(float).copy()
    characteristic_for_tolerance = (
        geometry.characteristic_length_m
        if isinstance(geometry, (SolidGeometry, CompoundGeometry))
        else max(geometry.length_x_m, geometry.width_y_m, geometry.thickness_z_m)
    )
    if not geometry.contains(position, tolerance=max(1e-10, characteristic_for_tolerance * 1e-8)):
        raise ValueError("photon emission position is outside the active solid")
    direction = initial_direction / np.linalg.norm(initial_direction)
    path = 0.0
    optical_time = 0.0
    bounces = 0
    scatters = 0
    outcome = "max_bounces"
    impact_position = _nan_position()
    detector_arrival_time = NAN
    detector_transmitted = False
    detector_id = 0
    detector_face = ""
    detector_attempts = 0
    tir_plus_x = 0
    tir_minus_x = 0
    detector_interface_reflections = 0
    detector_angle_rejected = 0
    detector_active_area_rejected = 0
    detector_collection_loss = 0
    detector_qe_loss = 0
    detector_valid_sensor_hit = False
    readout_arrival_time = NAN
    detector_incidence_angle_rad = NAN
    effective_detector_acceptance_rad = NAN
    effective_detector_quantum_efficiency = NAN
    effective_detector_collection_efficiency = NAN
    detector_window_shape = ""
    detector_window_width_m = NAN
    detector_window_height_m = NAN
    mirror_count = 0
    mylar_reflections = 0
    mylar_losses = 0
    reached_guide = False
    reached_readout = False
    active_to_guide_crossings = 0
    guide_to_readout_crossings = 0
    interface_transmissions = 0
    interface_reflections = 0
    tir_reflections = 0
    guide_section_entries: Counter = Counter()
    path_points: List[Tuple[float, float, float]] = [tuple(float(x) for x in position)]
    characteristic_length = geometry.characteristic_length_m if isinstance(geometry, (SolidGeometry, CompoundGeometry)) else max(geometry.length_x_m, geometry.width_y_m, geometry.thickness_z_m)
    # This epsilon is always applied as ``position + eps * direction``.  It
    # is intentionally larger than the boundary tolerance so that a photon
    # leaving a tapered side is not classified as still being on that side.
    eps = _transport_epsilon_for_geometry(geometry)
    section_index = geometry.section_index_at(position, direction) if isinstance(geometry, CompoundGeometry) else None
    n_inside = (
        float(geometry.refractive_index_at(position, wavelength_m, section_index))
        if isinstance(geometry, CompoundGeometry)
        else float(geometry.refractive_index(wavelength_m))
        if isinstance(geometry, SolidGeometry)
        else float(fused_silica_index(wavelength_m))
    )

    last_status = time.monotonic()
    for bounce_index in range(max(100, optical.max_bounces * 20 + 100)):
        if stop_event is not None and stop_event.is_set():
            raise _SimulationStopped
        if status_callback is not None and time.monotonic() - last_status >= 0.25:
            status_callback(f"optical transport heartbeat (photon {photon_id + 1}, bounce {bounce_index})")
            last_status = time.monotonic()
        if isinstance(geometry, CompoundGeometry):
            try:
                resolved_section = geometry.resolve_transport_section(position, direction, section_index)
            except RuntimeError as exc:
                raise RuntimeError(
                    f"compound photon transport state incoherent: section={section_index}, "
                    f"position={position.tolist()}, direction={direction.tolist()}: {exc}"
                ) from exc
            if resolved_section is None:
                outcome = _escaped_face_outcome(geometry.nearby_exit_face(position, direction, section_index))
                break
            section_index = resolved_section
            try:
                distance, face, outward, is_internal, adjacent_section = geometry.next_boundary(position, direction, section_index)
            except RuntimeError as exc:
                nearby_face = geometry.nearby_exit_face(position, direction, section_index)
                if nearby_face is not None:
                    outcome = _escaped_face_outcome(nearby_face)
                    break
                raise RuntimeError(
                    f"compound photon transport state incoherent: section={section_index}, "
                    f"position={position.tolist()}, direction={direction.tolist()}: {exc}"
                ) from exc
        else:
            distance, face, outward = _boundary_distance(position, direction, geometry)
            is_internal, adjacent_section = False, None
        event_distance, event = _segment_event_distance(distance, optical.absorption_length_m, optical.rayleigh_scattering_length_m, rng)
        position = position + event_distance * direction
        path += event_distance
        optical_time += event_distance * n_inside / C_LIGHT
        path_points.append(tuple(float(x) for x in position))
        if event == "absorbed":
            outcome = "absorbed"
            break
        if event == "scattered":
            scatters += 1
            direction = _sample_unit_vector(rng)
            position += eps * direction
            path_points.append(tuple(float(x) for x in position))
            continue

        cos_i = max(0.0, float(np.dot(direction, outward)))
        if is_internal:
            # This is an interface between two validated, adjacent sections.
            # Equal indices cross continuously.  Different indices use the
            # same Fresnel model as a detector portal, with Snell refraction
            # on transmission and specular reflection/TIR otherwise.
            n_out = float(geometry.refractive_index_at(position, wavelength_m, adjacent_section))
            equal_index = math.isclose(n_inside, n_out, rel_tol=1e-10, abs_tol=1e-12)
            reflectance, cos_t = _fresnel_unpolarized(cos_i, n_inside, n_out)
            transmit = equal_index or (cos_t > 0.0 and rng.random() >= reflectance)
            if transmit:
                interface_transmissions += 1
                if isinstance(geometry, CompoundGeometry) and adjacent_section in set(guide_section_indices):
                    reached_guide = True
                    guide_section_entries[int(adjacent_section)] += 1
                    if section_index in set(active_section_indices):
                        active_to_guide_crossings += 1
                direction = direction if equal_index else _refracted_direction(direction, outward, n_inside, n_out, cos_t)
                section_index = adjacent_section
                n_inside = n_out
                position += eps * direction
                path_points.append(tuple(float(x) for x in position))
                continue
            if bounces >= optical.max_bounces:
                outcome = "max_bounces"
                break
            interface_reflections += 1
            if cos_t <= 0.0:
                tir_reflections += 1
            reflected = direction - 2.0 * np.dot(direction, outward) * outward
            direction = reflected / np.linalg.norm(reflected)
            bounces += 1
            position += eps * direction
            path_points.append(tuple(float(x) for x in position))
            continue

        positive_mode = optical.positive_end_mode if face == "x+" else optical.opposite_end_mode
        n_out = optical.detector_index if face == "x+" and positive_mode == "detector" else 1.0
        detector_acceptance = optical.detector_acceptance_rad
        detector_qe = optical.detector_quantum_efficiency
        detector_collection = optical.detector_collection_efficiency
        detector_label = "plus_x"
        if face == "x-" and positive_mode == "detector":
            n_out = _opposite_detector_settings(optical)[3]
            detector_acceptance, detector_qe, detector_collection, _ = _opposite_detector_settings(optical)
            detector_label = "minus_x"
        is_detector_face = (face == "x+" and positive_mode == "detector") or (face == "x-" and positive_mode == "detector")
        reflectance, cos_t = _fresnel_unpolarized(cos_i, n_inside, n_out)
        if is_detector_face:
            detector_attempts += 1
            window_values = _active_area_values(optical, detector_label, geometry)
            detector_window_shape = str(window_values[0])
            detector_window_width_m = float(window_values[1])
            detector_window_height_m = float(window_values[2])
            effective_detector_acceptance_rad = float(detector_acceptance)
            effective_detector_quantum_efficiency = float(detector_qe)
            effective_detector_collection_efficiency = float(detector_collection)
            if cos_t <= 0.0:
                if detector_label == "plus_x":
                    tir_plus_x += 1
                else:
                    tir_minus_x += 1
                tir_reflections += 1
                outcome = f"tir_{detector_label}"
            elif rng.random() >= reflectance:
                detector_transmitted = True
                interface_transmissions += 1
                if isinstance(geometry, CompoundGeometry) and face in {"x+", "x-"}:
                    reached_readout = True
                    if section_index in set(guide_section_indices):
                        guide_to_readout_crossings += 1
                detector_id = 1 if detector_label == "plus_x" else 2
                detector_face = detector_label
                detector_arrival_time = emission_time_s + optical_time
                detector_position = tuple(float(x) for x in position)
                if not active_area_contains(detector_position, geometry, optical, detector_label):
                    detector_active_area_rejected += 1
                    outcome = f"detector_{detector_label}_active_area_rejected"
                else:
                    detector_valid_sensor_hit = True
                    # This is the physical active-window crossing.  The
                    # following angle/collection/QE branches may reject the
                    # photon, but they must not erase this arrival.
                    readout_arrival_time = detector_arrival_time
                    angle = math.acos(min(1.0, max(0.0, cos_i)))
                    detector_incidence_angle_rad = float(angle)
                if outcome.endswith("active_area_rejected"):
                    break
                if angle > detector_acceptance:
                    detector_angle_rejected += 1
                    outcome = f"detector_{detector_label}_angle_rejected"
                elif rng.random() >= detector_collection:
                    detector_collection_loss += 1
                    outcome = f"detector_{detector_label}_collection_loss"
                elif rng.random() >= detector_qe:
                    detector_qe_loss += 1
                    outcome = f"detector_{detector_label}_qe_loss"
                else:
                    outcome = f"detected_{detector_label}"
                if outcome.startswith("detected_"):
                    impact_position = detector_position
                break
            else:
                detector_interface_reflections += 1
                interface_reflections += 1

        # Mylar is an internal probabilistic reflector on every external
        # assembly surface except the coupling interface and +X readout.
        # Compound internal interfaces have already been handled above.
        terminal_state = None
        if isinstance(geometry, (SolidGeometry, CompoundGeometry)):
            for interface in getattr(geometry, "terminal_interfaces", ()):
                if interface.get("face") == ("plus_x" if face == "x+" else "minus_x" if face == "x-" else ""):
                    terminal_state = interface.get("state")
                    break
        mylar_face = (
            isinstance(geometry, CompoundGeometry) and geometry.wrapped_external and not is_detector_face
            and not (face in {"x+", "x-"} and terminal_state in {"mirror", "detector", "guide"})
        )
        if mylar_face and face != "x+":
            if rng.random() < geometry.mylar_reflectivity:
                mylar_reflections += 1
            else:
                mylar_losses += 1
                outcome = f"mylar_loss_{face}"
                break
        elif positive_mode == "mirror" and face in {"x+", "x-"}:
            mirror_reflectivity = optical.positive_mirror_reflectivity if face == "x+" else optical.opposite_mirror_reflectivity
            if rng.random() >= mirror_reflectivity:
                outcome = f"escaped_{face}"
                break
            mirror_count += 1
        elif positive_mode == "open" and face in {"x+", "x-"}:
            outcome = f"escaped_{face}"
            break
        elif not is_detector_face and reflectance < 1.0 and rng.random() >= reflectance:
            outcome = f"escaped_{face}"
            break

        if bounces >= optical.max_bounces:
            outcome = "max_bounces"
            break
        normal = _rough_surface_normal(outward, optical.surface_roughness_rad, rng)
        reflected = direction - 2.0 * np.dot(direction, normal) * normal
        if np.dot(reflected, outward) >= 0.0:
            reflected = direction - 2.0 * np.dot(direction, outward) * outward
        direction = reflected / np.linalg.norm(reflected)
        bounces += 1
        position += eps * direction
        path_points.append(tuple(float(x) for x in position))
    else:
        outcome = "max_bounces"

    guide_indices = tuple(int(index) for index in guide_section_indices)
    guide_path_segment_counts = []
    if isinstance(geometry, CompoundGeometry):
        for guide_index in guide_indices:
            count = 0
            section = geometry.sections[guide_index]
            for left, right in zip(path_points, path_points[1:]):
                midpoint = (np.asarray(left, dtype=float) + np.asarray(right, dtype=float)) / 2.0
                if section.contains(midpoint, tolerance=0.0):
                    count += 1
            guide_path_segment_counts.append(count)
    return PhotonRecord(
        photon_id=photon_id,
        wavelength_m=float(wavelength_m),
        emission_position_m=tuple(float(x) for x in emission_position),
        final_position_m=tuple(float(x) for x in position),
        impact_position_m=impact_position,
        path_length_m=float(path),
        optical_time_s=float(optical_time),
        emission_time_s=float(emission_time_s),
        arrival_time_s=float(detector_arrival_time),
        beta_at_emission=float(beta_at_emission),
        bounces=int(bounces),
        scatter_count=int(scatters),
        outcome=outcome,
        outcome_reason=_final_outcome_reason(outcome, detector_transmitted, detector_valid_sensor_hit),
        detector_id=detector_id,
        detector_face=detector_face,
        detected=outcome in {"detected_plus_x", "detected_minus_x"},
        detector_transmitted=detector_transmitted,
        detector_attempts=detector_attempts,
        tir_plus_x=tir_plus_x,
        tir_minus_x=tir_minus_x,
        detector_interface_reflections=detector_interface_reflections,
        detector_angle_rejected=detector_angle_rejected,
        detector_active_area_rejected=detector_active_area_rejected,
        detector_collection_loss=detector_collection_loss,
        detector_qe_loss=detector_qe_loss,
        detector_valid_sensor_hit=detector_valid_sensor_hit,
        detector_incidence_angle_rad=detector_incidence_angle_rad,
        effective_detector_acceptance_rad=effective_detector_acceptance_rad,
        effective_detector_quantum_efficiency=effective_detector_quantum_efficiency,
        effective_detector_collection_efficiency=effective_detector_collection_efficiency,
        detector_window_shape=detector_window_shape,
        detector_window_width_m=detector_window_width_m,
        detector_window_height_m=detector_window_height_m,
        mirror_count=mirror_count,
        path_points_m=tuple(path_points),
        photon_kind=photon_kind or ("scintillator" if isinstance(geometry, CompoundGeometry) and geometry.wrapped_external else "cherenkov"),
        mylar_reflections=mylar_reflections,
        mylar_losses=mylar_losses,
        reached_guide=reached_guide,
        reached_readout=reached_readout,
        path_coordinate_frame=coordinate_frame,
        active_to_guide_crossings=active_to_guide_crossings,
        guide_to_readout_crossings=guide_to_readout_crossings,
        interface_transmissions=interface_transmissions,
        interface_reflections=interface_reflections,
        tir_reflections=tir_reflections,
        guide_section_entry_counts=tuple(int(guide_section_entries.get(index, 0)) for index in guide_indices),
        guide_section_path_segment_counts=tuple(int(value) for value in guide_path_segment_counts),
        readout_arrival_time_s=float(readout_arrival_time),
    )


def _prepare_emission_sampler(particle: Particle, optical: OpticalConfig, track_length_m: float, geometry: Union[BarGeometry, SolidGeometry, CompoundGeometry] = None) -> Tuple[np.ndarray, np.ndarray, np.ndarray, np.ndarray, float]:
    sample_count = max(128, min(2048, optical.wavelength_grid_size))
    wavelengths = np.linspace(optical.wavelength_min_m, optical.wavelength_max_m, sample_count)
    position_count = 256
    positions = (np.arange(position_count, dtype=float) + 0.5) * track_length_m / position_count if track_length_m > 0.0 else np.array([0.0])
    energies = np.maximum(0.0, particle.kinetic_energy_J - (particle.dEdx_J_per_m or 0.0) * positions)
    betas = np.array([particle_beta(particle.mass_kg, energy) for energy in energies])
    if isinstance(geometry, CompoundGeometry):
        entry = geometry.particle_track(particle)[0]
        spectra = []
        for position, beta in zip(positions, betas):
            point = entry + float(position) * particle.direction
            section_index = geometry.section_index_at(point, particle.direction)
            local_section = geometry.sections[section_index]
            spectra.append(
                frank_tamm_spectrum(
                    wavelengths,
                    beta,
                    particle.charge_e,
                    refractive_index=local_section.refractive_index,
                )
            )
        spectra = np.asarray(spectra)
    else:
        refractive_index = geometry.refractive_index if isinstance(geometry, SolidGeometry) else fused_silica_index
        spectra = np.array([frank_tamm_spectrum(wavelengths, beta, particle.charge_e, refractive_index=refractive_index) for beta in betas])
    integrated_per_m = np.trapezoid(spectra, wavelengths, axis=1) if hasattr(np, "trapezoid") else np.trapz(spectra, wavelengths, axis=1)
    total_yield = float(np.trapezoid(integrated_per_m, positions) if hasattr(np, "trapezoid") else np.trapz(integrated_per_m, positions))
    return wavelengths, positions, betas, spectra, total_yield


def simulate(
    particle: Particle,
    geometry: Optional[Union[BarGeometry, SolidGeometry, CompoundGeometry]] = None,
    optical: Optional[OpticalConfig] = None,
    seed: int = 12345,
    event_id: int = 0,
    *,
    _progress_callback: Optional[Callable[[float], None]] = None,
    _status_callback: Optional[Callable[[str], None]] = None,
    stop_event=None,
    scintillator_assembly: Optional[ScintillatorAssembly] = None,
    scintillator_config: Optional[ScintillatorConfig] = None,
    detector_mode: str = "auto",
) -> SimulationResult:
    """Run one reproducible Monte Carlo event."""

    geometry = geometry or BarGeometry()
    optical = optical or OpticalConfig()
    mode = normalize_detector_run_mode(detector_mode)
    if mode in {"scintillator", "combined"} and (scintillator_assembly is None or scintillator_config is None):
        raise ValueError(f"{mode} mode requires a scintillator assembly and configuration")
    # Scintillator-only is a genuinely independent optical path.  Do not
    # validate or sample the Cherenkov geometry before entering that branch;
    # the shared particle/beam input is retained, while the detector-specific
    # track is derived from the scintillator volume itself.
    if mode == "scintillator":
        scintillator_assembly.validate()
        scintillator_config = replace(scintillator_config, enabled=True)
        scint_entry, scint_track_length = _scintillator_track(scintillator_assembly, particle)
        beta = particle_beta(particle.mass_kg, particle.kinetic_energy_J)
        result = SimulationResult(
            geometry=scintillator_assembly.transport_geometry, particle=particle, optical=optical or OpticalConfig(),
            event_id=event_id, track_length_m=scint_track_length, beta_start=beta, beta_end=beta,
            threshold_start=False, threshold_end=False, physical_expected_photons=0.0,
            detector_mode=mode,
            notes=[
                "Scintillator-only detector mode: Cherenkov generation and transport are disabled.",
                "The shared incident particle/beam configuration is retained for scintillator deposition.",
            ],
        )
        if _status_callback is not None:
            _status_callback(f"Generating scintillator photons for event {event_id + 1}")
        result.scintillator = simulate_scintillator(
            particle, scintillator_assembly, scintillator_config,
            seed=seed + 0x5EED, event_id=event_id, stop_event=stop_event,
            progress_callback=_progress_callback, status_callback=_status_callback,
        )
        result.stopped = result.scintillator.stopped
        if _progress_callback is not None:
            _progress_callback(1.0 if not result.stopped else 0.0)
        return result

    _validate_inputs(geometry, particle, optical)
    combined = mode == "combined"
    if combined:
        scintillator_config = replace(scintillator_config, enabled=True)
    detector_label = "combined detectors" if combined else "Cherenkov"
    if _status_callback is not None:
        _status_callback(f"Preparing {detector_label} event {event_id + 1}")
    rng = np.random.default_rng(seed)
    track_length = _track_length(geometry, particle)
    beta_start = particle_beta(particle.mass_kg, particle.kinetic_energy_J)
    final_energy = max(0.0, particle.kinetic_energy_J - (particle.dEdx_J_per_m or 0.0) * track_length)
    beta_end = particle_beta(particle.mass_kg, final_energy)
    wavelengths, positions, betas, spectra, total_yield = _prepare_emission_sampler(particle, optical, track_length, geometry=geometry)
    threshold_start = bool(np.any(spectra[0] > 0.0))
    threshold_end = bool(np.any(spectra[-1] > 0.0))
    result = SimulationResult(
        geometry=geometry, particle=particle, optical=optical, event_id=event_id,
        track_length_m=track_length, beta_start=beta_start, beta_end=beta_end,
        threshold_start=threshold_start, threshold_end=threshold_end,
        physical_expected_photons=total_yield,
        detector_mode=mode,
        notes=[
            (
                "Photon count is fixed by the requested Monte Carlo sample size; it is not a Poisson draw from total_yield."
                if optical.sampling_mode == "fixed"
                else "Generated photon count is drawn from Poisson(total_yield), with total_yield as the physical mean."
            ),
            "Rayleigh scattering is approximated as isotropic redistribution inside the active solid material.",
            "The detector +x face uses detector_index for Fresnel transmission when the active solid has a planar +x portal; other faces use the external air interface.",
            "Detector QE and collection efficiency apply only after interface transmission; QE=1 does not overcome total internal reflection.",
            f"The opposite -x end face mode is {optical.opposite_end_mode}.",
            "t=0 is the instant the particle crosses the active solid entry surface; +z is the forward beam direction.",
        ],
    )
    if isinstance(geometry, SolidGeometry):
        result.notes.append(f"Active solid: {geometry.shape}, material {geometry.normalized_material()}, axis {geometry.axis}.")
        if not geometry.detector_portals_supported:
            result.notes.append("The active solid has no planar +/-X terminal portal for the legacy end-detector logic; photons leaving curved/other-axis surfaces are not counted as detector crossings.")
    elif isinstance(geometry, CompoundGeometry):
        result.notes.append(
            "Active compound: " + "; ".join(
                f"section {index} {section.shape}/{section.normalized_material()}"
                for index, section in enumerate(geometry.sections)
            ) + ". Internal interfaces use local Fresnel/Snell transport."
        )
    if optical.sampling_mode == "poisson":
        sample_count = int(rng.poisson(total_yield)) if total_yield > 0.0 else 0
    else:
        sample_count = optical.photon_count
    if sample_count == 0 or total_yield <= 0.0:
        if sample_count > 0 or optical.sampling_mode == "fixed" and optical.photon_count > 0:
            result.notes.append("No wavelength in the requested range satisfies the Cherenkov condition along the track.")
        if _progress_callback is not None:
            _progress_callback(0.5 if combined else 1.0)
        if not combined:
            return result

    if _status_callback is not None:
        _status_callback(f"Transporting Cherenkov photons for event {event_id + 1}")
    if sample_count > 0 and total_yield > 0.0:
        position_probability = np.maximum(np.array([float(np.trapezoid(row, wavelengths) if hasattr(np, "trapezoid") else np.trapz(row, wavelengths)) for row in spectra]), 0.0)
        position_probability /= np.sum(position_probability)
        emission_indices = rng.choice(len(positions), size=sample_count, p=position_probability)
    else:
        emission_indices = ()
    particle_basis = _orthonormal_basis(particle.direction)
    last_status = time.monotonic()
    for sample_index, position_index in enumerate(emission_indices):
        if stop_event is not None and stop_event.is_set():
            result.stopped = True
            break
        local_beta = float(betas[position_index])
        spectrum = np.maximum(spectra[position_index], 0.0)
        spectrum_integral = float(np.trapezoid(spectrum, wavelengths) if hasattr(np, "trapezoid") else np.trapz(spectrum, wavelengths))
        if local_beta <= 0.0 or spectrum_integral <= 0.0:
            continue
        cdf = np.zeros_like(wavelengths)
        increments = 0.5 * (spectrum[:-1] + spectrum[1:]) * np.diff(wavelengths)
        cdf[1:] = np.cumsum(increments)
        cdf /= cdf[-1]
        wavelength = float(np.interp(rng.random(), cdf, wavelengths))
        # Use the active solid's optical model for the emission angle as well
        # as for the Frank--Tamm spectrum.  Keeping this tied to the geometry
        # is important for the explicitly supported ``air`` material and
        # prevents a material from being silently treated as fused silica.
        if isinstance(geometry, CompoundGeometry):
            emission_position = _entry_point(geometry, particle) + float(positions[position_index]) * particle.direction
            local_section_index = geometry.section_index_at(emission_position, particle.direction)
            n = float(geometry.refractive_index_at(emission_position, wavelength, local_section_index))
        else:
            n = float(geometry.refractive_index(wavelength)) if isinstance(geometry, SolidGeometry) else float(fused_silica_index(wavelength))
        cos_cherenkov = 1.0 / (local_beta * n)
        if cos_cherenkov >= 1.0:
            continue
        cone_angle = math.acos(cos_cherenkov)
        psi = rng.uniform(0.0, 2.0 * math.pi)
        photon_direction = cos_cherenkov * particle.direction + math.sin(cone_angle) * (math.cos(psi) * particle_basis[0] + math.sin(psi) * particle_basis[1])
        emission_distance = float(positions[position_index])
        emission_position = _entry_point(geometry, particle) + emission_distance * particle.direction
        emission_time_s = _particle_time_to_distance(particle, emission_distance)
        try:
            result.photon_records.append(_transport_photon(len(result.photon_records), wavelength, emission_position, photon_direction, emission_time_s, local_beta, geometry, optical, rng, stop_event=stop_event, status_callback=_status_callback, coordinate_frame="cherenkov_local_m"))
        except _SimulationStopped:
            result.stopped = True
            break
        if _progress_callback is not None and (sample_index == 0 or sample_index % 4 == 0 or sample_index + 1 == sample_count):
            fraction = (sample_index + 1) / sample_count
            _progress_callback(0.5 * fraction if combined else fraction)
        if _status_callback is not None and time.monotonic() - last_status >= 0.25:
            _status_callback(f"Transporting Cherenkov photons for event {event_id + 1}: photon {sample_index + 1}/{sample_count}")
            last_status = time.monotonic()
    result.generated_photons = len(result.photon_records)
    if mode == "combined":
        if _status_callback is not None:
            _status_callback(f"Generating scintillator photons for event {event_id + 1}")
        result.scintillator = simulate_scintillator(
            particle, scintillator_assembly, scintillator_config,
            seed=seed + 0x5EED, event_id=event_id, stop_event=stop_event,
            progress_callback=(lambda fraction: _progress_callback(0.5 + 0.5 * fraction)) if _progress_callback is not None else None,
            status_callback=_status_callback,
        )
        if result.scintillator.stopped:
            result.stopped = True
    if _progress_callback is not None:
        _progress_callback(1.0 if not result.stopped else min(1.0, len(result.photon_records) / max(sample_count, 1)))
    return result


def _scintillator_track(assembly: ScintillatorAssembly, particle: Particle) -> Tuple[np.ndarray, float]:
    """Return the first active entry and total active path length."""

    segments = _scintillator_track_segments(assembly, particle)
    if not segments:
        return assembly.active_sections[0].center_m.copy(), 0.0
    return segments[0][0].copy(), float(sum(item[1] for item in segments))


def _scintillator_track_segments(
    assembly: ScintillatorAssembly, particle: Particle,
) -> List[Tuple[np.ndarray, float, int, float]]:
    """Return actual charged-particle intervals through every active section.

    The reference point uses the shared Run particle x/y impact.  The
    assembly-local adapter sets that impact explicitly, so a global display
    offset or rotation cannot move the physical source independently of the
    beam configuration.
    """

    geometry = assembly.transport_geometry
    max_z = max(section.bounds_max_m[2] for section in assembly.active_sections)
    margin = max(4.0 * geometry.characteristic_length_m, 1e-6)
    reference = np.array([assembly.particle_impact_x_m, assembly.particle_impact_y_m, max_z], dtype=float)
    origin = reference - particle.direction * margin
    try:
        intervals = geometry._particle_intervals(origin, particle.direction)
    except ValueError:
        return []
    all_segments = [
        (origin + start * particle.direction, end - start, index, start - intervals[0][0])
        for start, end, index in intervals
    ]
    active_indices = set(assembly.active_section_indices or (0,))
    return [
        (entry.copy(), float(length), int(index), float(offset))
        for entry, length, index, offset in all_segments
        if int(index) in active_indices and float(length) > 0.0
    ]


def simulate_scintillator(
    particle: Particle,
    assembly: ScintillatorAssembly,
    config: ScintillatorConfig,
    seed: int = 12345,
    event_id: int = 0,
    stop_event=None,
    progress_callback: Optional[Callable[[float], None]] = None,
    status_callback: Optional[Callable[[str], None]] = None,
) -> ScintillatorResult:
    """Generate and transport scintillation photons through the assembly.

    Physical mode uses only the current simulator's constant dE/dx model:
    deposited energy is ``dEdx * scintillator track length`` unless an
    explicit deposited-energy override is supplied.  No detailed energy-loss
    or wavelength-spectrum model is implied.
    """

    assembly.validate()
    config.validate()
    track_segments = _scintillator_track_segments(assembly, particle)
    entry = track_segments[0][0].copy() if track_segments else assembly.active_sections[0].center_m.copy()
    track_length = float(sum(item[1] for item in track_segments))
    if config.deposited_energy_MeV is not None:
        deposited = float(config.deposited_energy_MeV)
    elif config.dEdx_MeV_per_m is not None:
        deposited = float(config.dEdx_MeV_per_m) * track_length
    elif particle.dEdx_J_per_m is not None:
        deposited = float(particle.dEdx_J_per_m) * track_length / MEV
    else:
        deposited = 0.0
    expected = deposited * config.photon_yield_per_MeV
    rng = np.random.default_rng(seed)
    count = int(rng.poisson(expected)) if config.sampling_mode in {"physical", "poisson"} and expected > 0.0 else int(config.photon_count)
    result = ScintillatorResult(
        assembly, config, deposited_energy_MeV=deposited,
        physical_expected_photons=expected, track_length_m=track_length,
        source_row_count=count,
    )
    if not config.enabled:
        if progress_callback is not None:
            progress_callback(1.0)
        return result
    transport_geometry = assembly.transport_geometry
    readout_face = "plus_x" if _terminal_mode(assembly.terminal_plus_surface, "detector") == "detector" else "minus_x"
    portal = detector_portal_dimensions(transport_geometry, readout_face)
    if portal is None:
        raise ValueError("scintillator readout requires a planar +X portal")
    _, portal_width, portal_height, _, _ = portal
    transport_optical = OpticalConfig(
        photon_count=count,
        wavelength_min_m=config.emission_wavelength_min_m,
        wavelength_max_m=config.emission_wavelength_max_m,
        absorption_length_m=config.absorption_length_m,
        max_bounces=config.max_bounces,
        detector_index=config.readout_index,
        detector_acceptance_rad=config.readout_acceptance_rad,
        detector_quantum_efficiency=config.readout_quantum_efficiency,
        detector_collection_efficiency=config.readout_collection_efficiency,
        detector_active_width_m=portal_width,
        detector_active_height_m=portal_height,
        detector_shape="rectangular" if portal[0] == "rectangle" else "elliptical",
        detector_active_semi_axis_a_m=portal_width / 2.0,
        detector_active_semi_axis_b_m=portal_height / 2.0,
        positive_end_mode=_terminal_mode(assembly.terminal_plus_surface, "detector"),
        opposite_end_mode=_terminal_mode(assembly.terminal_minus_surface, "open"),
        opposite_detector_index=config.readout_index,
        opposite_detector_acceptance_rad=config.readout_acceptance_rad,
        opposite_detector_quantum_efficiency=config.readout_quantum_efficiency,
        opposite_detector_collection_efficiency=config.readout_collection_efficiency,
        opposite_detector_shape="rectangular" if portal[0] == "rectangle" else "elliptical",
        opposite_detector_active_width_m=portal_width,
        opposite_detector_active_height_m=portal_height,
        opposite_detector_active_semi_axis_a_m=portal_width / 2.0,
        opposite_detector_active_semi_axis_b_m=portal_height / 2.0,
    )
    if count <= 0:
        if progress_callback is not None:
            progress_callback(1.0)
        return result
    if status_callback is not None:
        status_callback(f"Transporting scintillator photons for event {event_id + 1}")
    last_status = time.monotonic()
    segment_lengths = np.asarray([segment[1] for segment in track_segments], dtype=float)
    segment_probability = segment_lengths / np.sum(segment_lengths) if len(segment_lengths) and np.sum(segment_lengths) > 0.0 else None
    guide_indices = tuple(int(index) for index in (assembly.guide_section_indices or ((1,) if assembly.light_guide is not None else ())))
    active_indices = tuple(int(index) for index in assembly.active_section_indices or (0,))
    for index in range(count):
        if stop_event is not None and stop_event.is_set():
            result.stopped = True
            break
        if segment_probability is not None:
            segment_index = int(rng.choice(len(track_segments), p=segment_probability))
            segment_entry, segment_length, _, segment_offset = track_segments[segment_index]
            distance = float(rng.uniform(0.0, segment_length))
            position = segment_entry + distance * particle.direction
            particle_distance = float(segment_offset + distance)
        else:
            # Fixed diagnostic samples may intentionally use an explicit
            # deposited-energy override even when the shared beam misses the
            # active section.  Keep the source physically contained while
            # retaining zero charged-particle path length in that case.
            section = assembly.active_sections[0]
            position = section.center_m.copy()
            distance = 0.0
            particle_distance = 0.0
        direction = rng.normal(size=3)
        direction /= np.linalg.norm(direction)
        wavelength = rng.uniform(config.emission_wavelength_min_m, config.emission_wavelength_max_m)
        particle_time = _particle_time_to_distance(particle, particle_distance)
        decay = sample_scintillation_delay_s(config.decay_time_s, rng)
        try:
            record = _transport_photon(
                index, float(wavelength), position, direction, particle_time + decay,
                particle_beta(particle.mass_kg, particle.kinetic_energy_J), transport_geometry,
                transport_optical, rng, stop_event=stop_event,
                photon_kind="scintillator",
                status_callback=status_callback,
                coordinate_frame=assembly.coordinate_frame,
                guide_section_indices=guide_indices,
                active_section_indices=active_indices,
            )
        except _SimulationStopped:
            result.stopped = True
            break
        result.photon_records.append(record)
        if progress_callback is not None and (index == 0 or index % 4 == 0 or index + 1 == count):
            progress_callback((index + 1) / count)
        if status_callback is not None and time.monotonic() - last_status >= 0.25:
            status_callback(f"Transporting scintillator photons for event {event_id + 1}: photon {index + 1}/{count}")
            last_status = time.monotonic()
    # The endpoint invariant is enforced after transport, before summaries,
    # event tables, ROOT rows, or viewer scenes can observe the result.
    result.photon_records, result.readout_mismatch_count, result.readout_endpoint_tolerance_m = enforce_scintillator_readout_invariant(
        result.photon_records, assembly,
    )
    if progress_callback is not None:
        progress_callback(1.0 if not result.stopped else len(result.photon_records) / max(count, 1))
    return result


def simulate_events(
    particle: Particle,
    events: int,
    geometry: Optional[Union[BarGeometry, SolidGeometry, CompoundGeometry]] = None,
    optical: Optional[OpticalConfig] = None,
    seed: int = 12345,
    progress_callback: Optional[Callable[[int, int, float], None]] = None,
    status_callback: Optional[Callable[[str], None]] = None,
    stop_event=None,
    scintillator_assembly: Optional[ScintillatorAssembly] = None,
    scintillator_config: Optional[ScintillatorConfig] = None,
    detector_mode: str = "auto",
) -> List[SimulationResult]:
    """Simulate events with optional progress and cooperative cancellation.

    ``progress_callback(event_index, total_events, event_fraction)`` is called
    during photon transport and at the end of every completed event.  A
    supplied ``stop_event`` must provide ``is_set()``; cancellation returns
    the completed/partially transported results without raising and marks the
    interrupted result with ``stopped=True``.
    """

    if events < 1:
        raise ValueError("events must be at least 1")
    mode = normalize_detector_run_mode(
        detector_mode,
    )
    children = np.random.SeedSequence(seed).spawn(events)
    results: List[SimulationResult] = []
    for index, child in enumerate(children):
        if stop_event is not None and stop_event.is_set():
            if status_callback is not None:
                status_callback(f"Cancellation requested before event {index + 1}/{events}")
            break

        def event_status(message: str, event_index: int = index) -> None:
            if status_callback is not None:
                detector_text = {"cherenkov": "Cherenkov", "scintillator": "Scintillator", "combined": "combined detectors"}.get(mode, mode)
                status_callback(f"{detector_text} event {event_index + 1}/{events}: {message}")

        event_status("Starting")

        def event_progress(fraction: float, event_index: int = index) -> None:
            if progress_callback is not None:
                numeric_fraction = max(0.0, min(1.0, float(fraction)))
                progress_callback(event_index, events, numeric_fraction)

        result = simulate(
            particle,
            geometry,
            optical,
            seed=int(child.generate_state(1)[0]),
            event_id=index,
            _progress_callback=event_progress,
            _status_callback=event_status,
            stop_event=stop_event,
            scintillator_assembly=scintillator_assembly,
            scintillator_config=scintillator_config,
            detector_mode=mode,
        )
        results.append(result)
        if progress_callback is not None:
            progress_callback(index, events, 1.0)
        event_status("Completed" if not result.stopped else "Stopped")
        if result.stopped or (stop_event is not None and stop_event.is_set()):
            break
    return results


def aggregate_summary(results: Sequence[SimulationResult], detector_mode: Optional[str] = None) -> Dict[str, object]:
    """Return a JSON-friendly aggregate summary for one or more events."""

    results = list(results)
    outcomes = Counter(record.outcome for result in results for record in result.photon_records)
    detected = [record for result in results for record in result.photon_records if record.detected]
    crossings = [record for result in results for record in result.photon_records if record.detector_transmitted]
    plus_crossings = [record for record in crossings if record.detector_face == "plus_x"]
    minus_crossings = [record for record in crossings if record.detector_face == "minus_x"]
    mode = normalize_detector_run_mode(detector_mode or (results[0].detector_mode if results else "cherenkov"))
    participates_cherenkov = mode in {"cherenkov", "combined"}
    generated = sum(result.generated_photons for result in results) if participates_cherenkov else 0
    sampling_modes = sorted({result.optical.sampling_mode for result in results})
    sampling_mode = sampling_modes[0] if len(sampling_modes) == 1 and sampling_modes else "mixed"
    poisson_mean = (
        sum(result.physical_expected_photons for result in results if result.optical.sampling_mode == "poisson")
        if sampling_mode == "poisson" else None
    )
    diagnostic_counts = detector_diagnostic_counts(
        [record for result in results for record in result.photon_records]
    )
    scintillator_results = [result.scintillator for result in results if result.scintillator is not None]
    summary = {
        "events": len(results),
        "stopped": any(result.stopped for result in results),
        "detector_mode": mode,
        "sampling_mode": sampling_mode if participates_cherenkov else "n/a",
        "sampling_modes": sampling_modes,
        "sampling_description": (
            "fixed Monte Carlo test sample" if sampling_mode == "fixed"
            else "physical Poisson sampling from the Frank-Tamm yield" if sampling_mode == "poisson"
            else "mixed sampling modes"
        ),
        "requested_photons": sum(result.optical.photon_count for result in results) if participates_cherenkov else 0,
        "requested_sample_photons": sum(result.optical.photon_count for result in results) if participates_cherenkov else 0,
        "generated_photons": generated,
        "physical_expected_photons": sum(result.physical_expected_photons for result in results) if participates_cherenkov else 0.0,
        "poisson_mean": poisson_mean,
        "generation_fraction": (
            generated / poisson_mean if poisson_mean is not None and poisson_mean > 0.0
            else generated / sum(result.optical.photon_count for result in results)
            if sum(result.optical.photon_count for result in results) > 0 else 0.0
        ),
        "transported_photons": sum(len(result.photon_records) for result in results),
        "detected_photons": len(detected),
        "detected_plus_x": sum(record.detector_face == "plus_x" for record in detected),
        "detected_minus_x": sum(record.detector_face == "minus_x" for record in detected),
        "detector_face_crossings": len(crossings),
        "valid_sensor_hits": sum(record.detector_valid_sensor_hit for result in results for record in result.photon_records),
        "readout_arrivals": sum(readout_arrival_time_ns(record) is not None for result in results for record in result.photon_records),
        "detector_plus_x_crossings": len(plus_crossings),
        "detector_minus_x_crossings": len(minus_crossings),
        "detection_efficiency_among_crossings": len(detected) / len(crossings) if crossings else 0.0,
        "detection_efficiency_plus_x_crossings": sum(record.detected for record in plus_crossings) / len(plus_crossings) if plus_crossings else 0.0,
        "detection_efficiency_minus_x_crossings": sum(record.detected for record in minus_crossings) / len(minus_crossings) if minus_crossings else 0.0,
        "detector_hits": {
            "plus_x": sum(record.detector_face == "plus_x" for record in detected),
            "minus_x": sum(record.detector_face == "minus_x" for record in detected),
        },
        "detection_fraction_of_generated": len(detected) / generated if generated else 0.0,
        "outcomes": dict(sorted(outcomes.items())),
        "outcome_reasons": detector_outcome_reason_counts([record for result in results for record in result.photon_records]),
        "effective_max_bounces": {
            "cherenkov": int(results[0].optical.max_bounces) if participates_cherenkov and results else None,
            "scintillator": None,
        },
        "max_bounces_authority": "Run controls > Maximum photon bounces; one effective limit per participating detector",
        "mean_track_length_m": float(np.mean([result.track_length_m for result in results])) if results else None,
        "mean_track_length_mm": float(np.mean([result.track_length_m for result in results])) * 1e3 if results else None,
        "mean_detected_time_ns": float(np.mean([record.arrival_time_s for record in detected])) * 1e9 if detected else None,
        "mean_readout_arrival_ns": float(np.mean([readout_arrival_time_ns(record) for result in results for record in result.photon_records if readout_arrival_time_ns(record) is not None])) if any(readout_arrival_time_ns(record) is not None for result in results for record in result.photon_records) else None,
        "event_summaries": [result.summary() for result in results],
    }
    summary.update(diagnostic_counts)
    if scintillator_results and mode in {"scintillator", "combined"}:
        scint_summaries = [item.summary() for item in scintillator_results]
        summary["effective_max_bounces"]["scintillator"] = int(scintillator_results[0].config.max_bounces)
        for key in (
            "scintillator_requested_photons", "scintillator_generated_photons", "scintillator_transported_photons",
            "scintillator_source_rows", "scintillator_accepted_records", "scintillator_excluded_records",
            "scintillator_photons_reaching_guide", "scintillator_photons_reaching_readout", "scintillator_readout_arrivals", "scintillator_detected_photons",
            "scintillator_mylar_reflections", "scintillator_mylar_losses", "scintillator_absorbed",
            "scintillator_escaped", "scintillator_max_bounces", "scintillator_readout_alignment_mismatches",
            "scintillator_active_to_guide_crossings", "scintillator_guide_to_readout_crossings",
            "scintillator_interface_transmissions", "scintillator_interface_reflections", "scintillator_tir_reflections",
        ):
            summary[key] = sum(int(item.get(key, 0)) for item in scint_summaries)
        summary["scintillator_deposited_energy_MeV"] = sum(float(item.get("scintillator_deposited_energy_MeV", 0.0)) for item in scint_summaries)
        summary["scintillator_physical_expected_photons"] = sum(float(item.get("scintillator_physical_expected_photons", 0.0)) for item in scint_summaries)
        arrivals = [record.arrival_time_s for item in scintillator_results for record in item.photon_records if record.detected and math.isfinite(record.arrival_time_s)]
        summary["scintillator_mean_arrival_ns"] = float(np.mean(arrivals)) * 1e9 if arrivals else None
        summary["scintillator_earliest_arrival_ns"] = min(arrivals) * 1e9 if arrivals else None
        summary["scintillator_readout_endpoint_tolerance_um"] = max(
            float(item.get("scintillator_readout_endpoint_tolerance_um", 0.0)) for item in scint_summaries
        )
        guide_totals = {}
        for item in scint_summaries:
            for section in item.get("scintillator_guide_sections", ()):
                name = str(section.get("name", "guide"))
                target = guide_totals.setdefault(name, dict(section))
                target["entries"] = int(target.get("entries", 0)) + int(section.get("entries", 0))
                target["crossings"] = int(target.get("crossings", 0)) + int(section.get("crossings", 0))
                target["path_segments_inside"] = int(target.get("path_segments_inside", 0)) + int(section.get("path_segments_inside", 0))
        summary["scintillator_guide_sections"] = list(guide_totals.values())
        summary["scintillator_enabled"] = True
        summary["scintillator_event_summaries"] = scint_summaries
        cherenkov_arrivals = [record.arrival_time_s for result in results for record in result.photon_records if record.detected and math.isfinite(record.arrival_time_s)]
        summary["scintillator_cherenkov_coincidence_ns"] = (
            (min(arrivals) - min(cherenkov_arrivals)) * 1e9 if arrivals and cherenkov_arrivals else None
        )
    else:
        summary["scintillator_enabled"] = False
    return summary


def _as_results(result_or_results: SimulationResult | Sequence[SimulationResult]) -> List[SimulationResult]:
    if isinstance(result_or_results, SimulationResult):
        return [result_or_results]
    results = list(result_or_results)
    if not all(isinstance(result, SimulationResult) for result in results):
        raise TypeError("result_or_results must contain SimulationResult objects")
    return results


def _load_root_dependencies():
    try:
        import awkward as ak
        import uproot
    except ImportError as exc:
        raise ImportError(
            "ROOT output requires uproot and awkward. Install them with: "
            "python3 -m pip install uproot awkward"
        ) from exc
    return uproot, ak


def check_root_dependencies() -> str:
    """Return a preflight message or raise the exact ROOT dependency error."""

    uproot, awkward = _load_root_dependencies()
    return f"ROOT dependencies available: uproot {uproot.__version__}, awkward {awkward.__version__}"


def _root_event_entry_point(result: SimulationResult) -> np.ndarray:
    """Return a mode-owned entry point for optional ROOT event metadata."""

    mode = normalize_detector_run_mode(getattr(result, "detector_mode", "cherenkov"))
    if mode == "scintillator" and result.scintillator is not None:
        entry, _ = _scintillator_track(result.scintillator.assembly, result.particle)
        return np.asarray(entry, dtype=float)
    try:
        return _entry_point(result.geometry, result.particle)
    except ValueError as exc:
        raise RuntimeError(
            f"ROOT event metadata could not derive an entry point for detector_mode={mode}; "
            "simulation results remain valid, but this optional ROOT schema row cannot be written"
        ) from exc


def write_root(result_or_results: SimulationResult | Sequence[SimulationResult], filename: str) -> str:
    """Write verified Photons and Events TTrees and return the absolute path."""

    output_path = Path(filename).expanduser().resolve()
    temporary_path: Optional[str] = None
    try:
        uproot, ak = _load_root_dependencies()
        results = _as_results(result_or_results)
        def owned_records(result: SimulationResult) -> list[PhotonRecord]:
            mode = normalize_detector_run_mode(getattr(result, "detector_mode", "cherenkov"))
            records = list(result.photon_records) if mode in {"cherenkov", "combined"} else []
            if mode in {"scintillator", "combined"} and result.scintillator is not None:
                records.extend(result.scintillator.photon_records)
            return records
        records_by_event = {id(result): owned_records(result) for result in results}
        photon_rows = [(result.event_id, result, record) for result in results for record in records_by_event[id(result)]]
        photon_data = {
            "event_id": np.asarray([event_id for event_id, _, _ in photon_rows], dtype=np.int32),
            "photon_id": np.asarray([record.photon_id for _, _, record in photon_rows], dtype=np.int32),
            "generated": np.asarray([record.generated for _, _, record in photon_rows], dtype=np.bool_),
            "detected": np.asarray([record.detected for _, _, record in photon_rows], dtype=np.bool_),
            "detector_transmitted": np.asarray([record.detector_transmitted for _, _, record in photon_rows], dtype=np.bool_),
            "detector_valid_sensor_hit": np.asarray([record.detector_valid_sensor_hit for _, _, record in photon_rows], dtype=np.bool_),
            "readout_arrival": np.asarray([readout_arrival_time_ns(record) is not None for _, _, record in photon_rows], dtype=np.bool_),
            "detector_id": np.asarray([record.detector_id for _, _, record in photon_rows], dtype=np.int8),
            "detector_attempts": np.asarray([record.detector_attempts for _, _, record in photon_rows], dtype=np.int32),
            "tir_plus_x": np.asarray([record.tir_plus_x for _, _, record in photon_rows], dtype=np.int32),
            "tir_minus_x": np.asarray([record.tir_minus_x for _, _, record in photon_rows], dtype=np.int32),
            "detector_interface_reflections": np.asarray([record.detector_interface_reflections for _, _, record in photon_rows], dtype=np.int32),
            "detector_angle_rejected": np.asarray([record.detector_angle_rejected for _, _, record in photon_rows], dtype=np.int32),
            "detector_active_area_rejected": np.asarray([record.detector_active_area_rejected for _, _, record in photon_rows], dtype=np.int32),
            "detector_collection_loss": np.asarray([record.detector_collection_loss for _, _, record in photon_rows], dtype=np.int32),
            "detector_qe_loss": np.asarray([record.detector_qe_loss for _, _, record in photon_rows], dtype=np.int32),
            "detector_incidence_angle_rad": np.asarray([record.detector_incidence_angle_rad for _, _, record in photon_rows], dtype=np.float64),
            "effective_detector_acceptance_rad": np.asarray([record.effective_detector_acceptance_rad for _, _, record in photon_rows], dtype=np.float64),
            "effective_detector_quantum_efficiency": np.asarray([record.effective_detector_quantum_efficiency for _, _, record in photon_rows], dtype=np.float64),
            "effective_detector_collection_efficiency": np.asarray([record.effective_detector_collection_efficiency for _, _, record in photon_rows], dtype=np.float64),
            "detector_window_width_mm": np.asarray([record.detector_window_width_m * 1e3 for _, _, record in photon_rows], dtype=np.float64),
            "detector_window_height_mm": np.asarray([record.detector_window_height_m * 1e3 for _, _, record in photon_rows], dtype=np.float64),
            "detector_face": ak.Array([record.detector_face for _, _, record in photon_rows]),
            "outcome": ak.Array([record.outcome for _, _, record in photon_rows]),
            "photon_kind": ak.Array([record.photon_kind for _, _, record in photon_rows]),
            "path_coordinate_frame": ak.Array([record.path_coordinate_frame for _, _, record in photon_rows]),
            "wavelength_nm": np.asarray([record.wavelength_m * 1e9 for _, _, record in photon_rows], dtype=np.float64),
            "wavelength_m": np.asarray([record.wavelength_m for _, _, record in photon_rows], dtype=np.float64),
            "emission_x_m": np.asarray([record.emission_position_m[0] for _, _, record in photon_rows], dtype=np.float64),
            "emission_y_m": np.asarray([record.emission_position_m[1] for _, _, record in photon_rows], dtype=np.float64),
            "emission_z_m": np.asarray([record.emission_position_m[2] for _, _, record in photon_rows], dtype=np.float64),
            "emission_x_mm": np.asarray([record.emission_position_m[0] * 1e3 for _, _, record in photon_rows], dtype=np.float64),
            "emission_y_mm": np.asarray([record.emission_position_m[1] * 1e3 for _, _, record in photon_rows], dtype=np.float64),
            "emission_z_mm": np.asarray([record.emission_position_m[2] * 1e3 for _, _, record in photon_rows], dtype=np.float64),
            "impact_x_m": np.asarray([record.impact_position_m[0] for _, _, record in photon_rows], dtype=np.float64),
            "impact_y_m": np.asarray([record.impact_position_m[1] for _, _, record in photon_rows], dtype=np.float64),
            "impact_z_m": np.asarray([record.impact_position_m[2] for _, _, record in photon_rows], dtype=np.float64),
            "impact_x_mm": np.asarray([record.impact_position_m[0] * 1e3 for _, _, record in photon_rows], dtype=np.float64),
            "impact_y_mm": np.asarray([record.impact_position_m[1] * 1e3 for _, _, record in photon_rows], dtype=np.float64),
            "impact_z_mm": np.asarray([record.impact_position_m[2] * 1e3 for _, _, record in photon_rows], dtype=np.float64),
            "arrival_time_ns": np.asarray([record.arrival_time_s * 1e9 for _, _, record in photon_rows], dtype=np.float64),
            "readout_arrival_time_ns": np.asarray([readout_arrival_time_ns(record) if readout_arrival_time_ns(record) is not None else NAN for _, _, record in photon_rows], dtype=np.float64),
            "path_length_m": np.asarray([record.path_length_m for _, _, record in photon_rows], dtype=np.float64),
            "path_length_mm": np.asarray([record.path_length_m * 1e3 for _, _, record in photon_rows], dtype=np.float64),
            "bounces": np.asarray([record.bounces for _, _, record in photon_rows], dtype=np.int32),
            "mirror_count": np.asarray([record.mirror_count for _, _, record in photon_rows], dtype=np.int32),
            "scatter_count": np.asarray([record.scatter_count for _, _, record in photon_rows], dtype=np.int32),
            "particle_beta_start": np.asarray([result.beta_start for _, result, _ in photon_rows], dtype=np.float64),
            "particle_beta_end": np.asarray([result.beta_end for _, result, _ in photon_rows], dtype=np.float64),
            "particle_track_length_m": np.asarray([result.track_length_m for _, result, _ in photon_rows], dtype=np.float64),
            "particle_track_length_mm": np.asarray([result.track_length_m * 1e3 for _, result, _ in photon_rows], dtype=np.float64),
            "emission_time_ns": np.asarray([record.emission_time_s * 1e9 for _, _, record in photon_rows], dtype=np.float64),
            "beta_at_emission": np.asarray([record.beta_at_emission for _, _, record in photon_rows], dtype=np.float64),
            "scintillator_photon": np.asarray([record.photon_kind == "scintillator" for _, _, record in photon_rows], dtype=np.bool_),
            "scintillator_mylar_reflections": np.asarray([record.mylar_reflections for _, _, record in photon_rows], dtype=np.int32),
            "scintillator_mylar_losses": np.asarray([record.mylar_losses for _, _, record in photon_rows], dtype=np.int32),
            "scintillator_reached_guide": np.asarray([record.reached_guide for _, _, record in photon_rows], dtype=np.bool_),
            "scintillator_reached_readout": np.asarray([record.reached_readout for _, _, record in photon_rows], dtype=np.bool_),
            "scintillator_active_to_guide_crossings": np.asarray([record.active_to_guide_crossings for _, _, record in photon_rows], dtype=np.int32),
            "scintillator_guide_to_readout_crossings": np.asarray([record.guide_to_readout_crossings for _, _, record in photon_rows], dtype=np.int32),
            "scintillator_interface_transmissions": np.asarray([record.interface_transmissions for _, _, record in photon_rows], dtype=np.int32),
            "scintillator_interface_reflections": np.asarray([record.interface_reflections for _, _, record in photon_rows], dtype=np.int32),
            "scintillator_tir_reflections": np.asarray([record.tir_reflections for _, _, record in photon_rows], dtype=np.int32),
            "scintillator_readout_alignment_valid": np.asarray([record.readout_alignment_valid for _, _, record in photon_rows], dtype=np.bool_),
            "scintillator_readout_endpoint_distance_m": np.asarray([record.readout_endpoint_distance_m for _, _, record in photon_rows], dtype=np.float64),
        }
        photon_data["outcome_reason"] = ak.Array([record.outcome_reason for _, _, record in photon_rows])
        photon_data["detector_window_shape"] = ak.Array([record.detector_window_shape for _, _, record in photon_rows])
        photon_types = {name: ("string" if name in {"detector_face", "outcome", "outcome_reason", "detector_window_shape", "photon_kind", "path_coordinate_frame"} else values.dtype) for name, values in photon_data.items()}
        root_entries = [_root_event_entry_point(result) for result in results]
        event_data = {
            "event_id": np.asarray([result.event_id for result in results], dtype=np.int32),
            "sampling_mode": ak.Array([result.optical.sampling_mode for result in results]),
            "requested_photons": np.asarray([result.optical.photon_count for result in results], dtype=np.int32),
            "requested_sample_photons": np.asarray([result.optical.photon_count for result in results], dtype=np.int32),
            "generated_photons": np.asarray([result.generated_photons for result in results], dtype=np.int32),
            "transported_photons": np.asarray([len(records_by_event[id(result)]) for result in results], dtype=np.int32),
            "detected_photons": np.asarray([sum(record.detected for record in records_by_event[id(result)]) for result in results], dtype=np.int32),
            "detected_plus_x": np.asarray([sum(record.detected and record.detector_face == "plus_x" for record in records_by_event[id(result)]) for result in results], dtype=np.int32),
            "detected_minus_x": np.asarray([sum(record.detected and record.detector_face == "minus_x" for record in records_by_event[id(result)]) for result in results], dtype=np.int32),
            "detector_face_crossings": np.asarray([sum(record.detector_transmitted for record in records_by_event[id(result)]) for result in results], dtype=np.int32),
            "valid_sensor_hits": np.asarray([sum(record.detector_valid_sensor_hit for record in records_by_event[id(result)]) for result in results], dtype=np.int32),
            "readout_arrivals": np.asarray([sum(readout_arrival_time_ns(record) is not None for record in records_by_event[id(result)]) for result in results], dtype=np.int32),
            "detector_attempts": np.asarray([detector_diagnostic_counts(records_by_event[id(result)])["detector_attempts"] for result in results], dtype=np.int32),
            "tir_plus_x": np.asarray([detector_diagnostic_counts(records_by_event[id(result)])["tir_plus_x"] for result in results], dtype=np.int32),
            "tir_minus_x": np.asarray([detector_diagnostic_counts(records_by_event[id(result)])["tir_minus_x"] for result in results], dtype=np.int32),
            "detector_interface_reflections": np.asarray([detector_diagnostic_counts(records_by_event[id(result)])["detector_interface_reflections"] for result in results], dtype=np.int32),
            "detector_angle_rejected": np.asarray([detector_diagnostic_counts(records_by_event[id(result)])["detector_angle_rejected"] for result in results], dtype=np.int32),
            "detector_active_area_rejected": np.asarray([detector_diagnostic_counts(records_by_event[id(result)])["detector_active_area_rejected"] for result in results], dtype=np.int32),
            "detector_collection_loss": np.asarray([detector_diagnostic_counts(records_by_event[id(result)])["detector_collection_loss"] for result in results], dtype=np.int32),
            "detector_qe_loss": np.asarray([detector_diagnostic_counts(records_by_event[id(result)])["detector_qe_loss"] for result in results], dtype=np.int32),
            "detector_losses": np.asarray([detector_diagnostic_counts(records_by_event[id(result)])["detector_losses"] for result in results], dtype=np.int32),
            "physical_expected_photons": np.asarray([result.physical_expected_photons for result in results], dtype=np.float64),
            "poisson_mean": np.asarray([result.physical_expected_photons if result.optical.sampling_mode == "poisson" else NAN for result in results], dtype=np.float64),
            "generation_fraction": np.asarray([
                result.generated_photons / result.physical_expected_photons
                if result.optical.sampling_mode == "poisson" and result.physical_expected_photons > 0.0
                else result.generated_photons / result.optical.photon_count
                if result.optical.photon_count > 0 else 0.0
                for result in results
            ], dtype=np.float64),
            "entry_x_m": np.asarray([entry[0] for entry in root_entries], dtype=np.float64),
            "entry_y_m": np.asarray([entry[1] for entry in root_entries], dtype=np.float64),
            "entry_z_m": np.asarray([entry[2] for entry in root_entries], dtype=np.float64),
            "entry_x_mm": np.asarray([entry[0] * 1e3 for entry in root_entries], dtype=np.float64),
            "entry_y_mm": np.asarray([entry[1] * 1e3 for entry in root_entries], dtype=np.float64),
            "entry_z_mm": np.asarray([entry[2] * 1e3 for entry in root_entries], dtype=np.float64),
            "direction_x": np.asarray([result.particle.direction[0] for result in results], dtype=np.float64),
            "direction_y": np.asarray([result.particle.direction[1] for result in results], dtype=np.float64),
            "direction_z": np.asarray([result.particle.direction[2] for result in results], dtype=np.float64),
            "beta": np.asarray([result.beta_start for result in results], dtype=np.float64),
            "track_length_m": np.asarray([result.track_length_m for result in results], dtype=np.float64),
            "track_length_mm": np.asarray([result.track_length_m * 1e3 for result in results], dtype=np.float64),
            "scintillator_requested_photons": np.asarray([result.summary().get("scintillator_requested_photons", 0) for result in results], dtype=np.int32),
            "scintillator_generated_photons": np.asarray([result.summary().get("scintillator_generated_photons", 0) for result in results], dtype=np.int32),
            "scintillator_transported_photons": np.asarray([result.summary().get("scintillator_transported_photons", 0) for result in results], dtype=np.int32),
            "scintillator_photons_reaching_guide": np.asarray([result.summary().get("scintillator_photons_reaching_guide", 0) for result in results], dtype=np.int32),
            "scintillator_photons_reaching_readout": np.asarray([result.summary().get("scintillator_photons_reaching_readout", 0) for result in results], dtype=np.int32),
            "scintillator_mylar_reflections": np.asarray([result.summary().get("scintillator_mylar_reflections", 0) for result in results], dtype=np.int32),
            "scintillator_mylar_losses": np.asarray([result.summary().get("scintillator_mylar_losses", 0) for result in results], dtype=np.int32),
            "scintillator_detected_photons": np.asarray([result.summary().get("scintillator_detected_photons", 0) for result in results], dtype=np.int32),
            "scintillator_source_rows": np.asarray([result.summary().get("scintillator_source_rows", 0) for result in results], dtype=np.int32),
            "scintillator_accepted_records": np.asarray([result.summary().get("scintillator_accepted_records", 0) for result in results], dtype=np.int32),
            "scintillator_excluded_records": np.asarray([result.summary().get("scintillator_excluded_records", 0) for result in results], dtype=np.int32),
            "scintillator_readout_alignment_mismatches": np.asarray([result.summary().get("scintillator_readout_alignment_mismatches", 0) for result in results], dtype=np.int32),
            "scintillator_mean_arrival_ns": np.asarray([result.summary().get("scintillator_mean_arrival_ns") if result.summary().get("scintillator_mean_arrival_ns") is not None else NAN for result in results], dtype=np.float64),
            "scintillator_earliest_arrival_ns": np.asarray([result.summary().get("scintillator_earliest_arrival_ns") if result.summary().get("scintillator_earliest_arrival_ns") is not None else NAN for result in results], dtype=np.float64),
        }
        event_types = {name: ("string" if name == "sampling_mode" else values.dtype) for name, values in event_data.items()}
        metadata = {
            "coordinate_system": "x=100 mm length, y=20 mm width, z=10 mm thickness; entry z=+5 mm; end faces x=+/-50 mm",
            "time_reference": "t=0 is particle crossing of the quartz entry surface z=+half_z",
            "units": "coordinate branches have m and mm variants; wavelength nm; times ns; beta dimensionless",
            "tree_description": "Photons has one row per transported/generated photon; Events has one row per event",
            "detector_diagnostics": "TIR counters record reflected detector-face encounters; detector losses occur only after interface transmission",
            "readout_arrival_definition": "readout_arrival is a physical crossing of the configured finite active sensor/readout window before angular acceptance, collection, or QE response; readout_arrival_time_ns is the crossing time",
            "scintillator_path_invariant": "Scintillator PhotonRecord.path_points_m are in scintillator_local_m for local transport; detected endpoints are on the physical +X active/guide readout face and are transformed once by the event viewer",
            "sampling": "fixed uses photon_count as a controlled transport sample; poisson draws generated count from Poisson(physical_expected_photons)",
        }
        fd, temporary_path = tempfile.mkstemp(prefix=output_path.name + ".", suffix=".tmp", dir=str(output_path.parent))
        os.close(fd)
        with uproot.recreate(temporary_path) as root_file:
            photon_tree = root_file.mktree("Photons", photon_types)
            if photon_rows:
                photon_tree.extend(photon_data)
            event_tree = root_file.mktree("Events", event_types)
            if results:
                event_tree.extend(event_data)
            metadata_tree = root_file.mktree("Metadata", {"key": "string", "value": "string"})
            metadata_tree.extend({"key": ak.Array(list(metadata.keys())), "value": ak.Array(list(metadata.values()))})
        with uproot.open(temporary_path) as check_file:
            verified_keys = set(check_file.keys(cycle=False))
            if not {"Photons", "Events"}.issubset(verified_keys):
                raise RuntimeError("verification did not find Photons and Events TTrees")
            if check_file["Photons"].num_entries != len(photon_rows) or check_file["Events"].num_entries != len(results):
                raise RuntimeError("verified ROOT tree entry counts do not match the simulation results")
        os.replace(temporary_path, output_path)
        temporary_path = None
        absolute = str(output_path)
        # Both uproot file contexts above have exited before this diagnostic is
        # printed.  Flush it explicitly because the GUI is commonly launched
        # from a Windows console whose stdout may otherwise be buffered.
        print(f"ROOT output written: {absolute}", flush=True)
        return absolute
    except Exception as exc:
        if temporary_path:
            try:
                os.unlink(temporary_path)
            except OSError:
                pass
        platform_hint = (
            " Check that the destination is not open in ROOT/file preview software "
            "and that the Windows account can replace files in the run directory."
            if isinstance(exc, OSError) else ""
        )
        if isinstance(exc, ImportError):
            detail = (
                f"ROOT optional output unavailable for {output_path}: {exc}. "
                "Simulation results remain valid; PNG/JSON outputs are preserved. "
                "Install optional dependencies with: python3 -m pip install uproot awkward"
            )
        else:
            detail = (
                f"ROOT optional output serialization failed for {output_path}: {exc}.{platform_hint} "
                "Simulation results remain valid; PNG/JSON outputs are preserved."
            )
        raise RuntimeError(detail) from exc


def _geometry_bounds(geometry: Union[BarGeometry, SolidGeometry, CompoundGeometry]) -> Tuple[np.ndarray, np.ndarray]:
    if isinstance(geometry, (SolidGeometry, CompoundGeometry)):
        return geometry.bounds_min_m, geometry.bounds_max_m
    return np.array([-geometry.half_x, -geometry.half_y, -geometry.half_z], dtype=float), np.array([geometry.half_x, geometry.half_y, geometry.half_z], dtype=float)


def _set_axes_equal(ax, geometry: Union[BarGeometry, SolidGeometry, CompoundGeometry]) -> None:
    low, high = _geometry_bounds(geometry)
    ax.set_xlim(low[0] * 1e3, high[0] * 1e3)
    ax.set_ylim(low[1] * 1e3, high[1] * 1e3)
    ax.set_zlim(low[2] * 1e3, high[2] * 1e3)
    ax.set_box_aspect(tuple(np.maximum(high - low, 1e-12)))


def _draw_bar(ax, geometry: Union[BarGeometry, SolidGeometry, CompoundGeometry]) -> None:
    if isinstance(geometry, (SolidGeometry, CompoundGeometry)):
        rings = geometry.surface_rings()
        for ring in rings:
            closed = np.vstack([ring, ring[0]])
            ax.plot(closed[:, 0] * 1e3, closed[:, 1] * 1e3, closed[:, 2] * 1e3, color="#6b7280", linewidth=1.0)
        for ring_index in range(0, len(rings) - 1, 2):
            if len(rings[ring_index]) != len(rings[ring_index + 1]):
                continue
            stride = max(1, len(rings[ring_index]) // 12)
            for index in range(0, len(rings[ring_index]), stride):
                first, second = rings[ring_index][index], rings[ring_index + 1][index]
                ax.plot([first[0] * 1e3, second[0] * 1e3], [first[1] * 1e3, second[1] * 1e3], [first[2] * 1e3, second[2] * 1e3], color="#6b7280", linewidth=0.8)
        return
    from itertools import product

    corners = list(product((-geometry.half_x * 1e3, geometry.half_x * 1e3), (-geometry.half_y * 1e3, geometry.half_y * 1e3), (-geometry.half_z * 1e3, geometry.half_z * 1e3)))
    for i, first in enumerate(corners):
        for second in corners[i + 1:]:
            if sum(a != b for a, b in zip(first, second)) == 1:
                ax.plot([first[0], second[0]], [first[1], second[1]], [first[2], second[2]], color="#8ca6b8", linewidth=0.8, alpha=0.7)
    yy, zz = np.meshgrid([-geometry.half_y * 1e3, geometry.half_y * 1e3], [-geometry.half_z * 1e3, geometry.half_z * 1e3])
    for x, color in ((geometry.half_x * 1e3, "#f4b183"), (-geometry.half_x * 1e3, "#c5b0d5")):
        ax.plot_surface(np.full_like(yy, x), yy, zz, color=color, alpha=0.25, shade=False)


def _plot_canonical_event_scene(
    result: SimulationResult,
    output_path: str,
    event_scene: Mapping[str, object],
    display_photons: int = 40,
    show: bool = False,
    photon_message: Optional[str] = None,
) -> None:
    """Render the renderer-neutral global event scene.

    The GUI builds ``event_scene`` with the same canonical assembly builder as
    the interactive event viewer.  This function is intentionally only a
    consumer: it never reconstructs detector geometry, applies a second
    transform, or merges records from different detector owners.
    """

    import matplotlib.pyplot as plt
    from mpl_toolkits.mplot3d.art3d import Poly3DCollection

    mode = normalize_detector_run_mode(event_scene.get("mode", result.detector_mode))
    mode_labels = {
        "cherenkov": "Cherenkov only",
        "scintillator": "Scintillator only",
        "combined": "All",
    }
    fig = plt.figure(figsize=(10, 7))
    ax = fig.add_subplot(111, projection="3d")

    # Faces are already in global mm.  Keep body/guide geometry visible while
    # rendering sensor/readout components with a distinct, owner-specific
    # accent so a cloud of photons cannot be mistaken for a readout.
    component_styles = {
        "geometry": ("#5b9bd5", "-", 1.0, 0.10),
        "guide": ("#f0a35b", "-", 1.15, 0.12),
        "sensor": ("#2ca25f", "-", 1.5, 0.20),
        "readout_interface": ("#7b3294", "-", 1.5, 0.18),
        "mirror": ("#008c95", "--", 1.3, 0.14),
        "mylar_overlay": ("#d9a441", ":", 1.0, 0.08),
    }
    plotted_components = set()
    for component in event_scene.get("geometry_components", ()):
        kind = str(component.get("kind", "geometry"))
        color, linestyle, linewidth, alpha = component_styles.get(
            kind, (str(component.get("color", "#80b1d3")), "-", 1.0, 0.10)
        )
        assembly = str(component.get("assembly_type", ""))
        plotted_components.add(assembly)
        for face in component.get("faces_mm", ()):
            if len(face) < 2:
                continue
            polygon = np.asarray(face, dtype=float)
            closed = np.vstack([polygon, polygon[0]])
            ax.plot(closed[:, 0], closed[:, 1], closed[:, 2], color=color,
                    linestyle=linestyle, linewidth=linewidth, alpha=0.92)
            if kind in {"geometry", "guide", "sensor"} and len(face) >= 3:
                ax.add_collection3d(Poly3DCollection(
                    [polygon], facecolors=color, edgecolors="none", alpha=alpha
                ))

    # The scene builder has already applied detector participation filtering.
    # Photon lines retain their detector_type, so ownership is visible in the
    # static image and can never be inferred from outcome color alone.
    line_colors = {
        "cherenkov": {"detected_plus_x": "#d95f02", "detected_minus_x": "#e7298a", "escaped": "#1f78b4", "absorbed": "#b2182b", "max_bounces": "#7570b3", "losses": "#ff7f00", "unknown": "#6b7280"},
        "scintillator": {"detected_plus_x": "#1b9e77", "detected_minus_x": "#66a61e", "escaped": "#4daf4a", "absorbed": "#762a83", "max_bounces": "#386cb0", "losses": "#a6761d", "unknown": "#008c95"},
    }
    # ``build_interactive_event_scene`` applies the display budget independently
    # to each detector owner.  Keep the same per-owner contract defensively at
    # this consumer boundary, but never slice the combined list: a single
    # combined slice consumed the whole budget on Cherenkov in v48 and
    # silently removed every Scintillator trajectory.
    display_limit = max(0, int(display_photons))
    displayed = []
    displayed_by_owner = {}
    for line in event_scene.get("photon_lines", ()):
        owner = str(line.get("detector_type", line.get("assembly_type", "")))
        if displayed_by_owner.get(owner, 0) >= display_limit:
            continue
        displayed.append(line)
        displayed_by_owner[owner] = displayed_by_owner.get(owner, 0) + 1
    for line in displayed:
        points = np.asarray(line.get("points_mm", ()), dtype=float)
        if points.ndim != 2 or points.shape[0] < 2 or points.shape[1] < 3:
            continue
        detector_type = str(line.get("detector_type", line.get("assembly_type", "")))
        category = str(line.get("category", "unknown"))
        color = line_colors.get(detector_type, {}).get(category, "#6b7280")
        ax.plot(points[:, 0], points[:, 1], points[:, 2], color=color,
                alpha=0.68, linewidth=0.85)
        if bool(line.get("detected_marker")):
            ax.scatter([points[-1, 0]], [points[-1, 1]], [points[-1, 2]],
                       color=color, edgecolors="#111827", linewidths=0.55,
                       s=22, marker="D", depthshade=False)
    track = np.asarray(event_scene.get("track_mm", ()), dtype=float)
    if track.ndim == 2 and track.shape[0] >= 2:
        ax.plot(track[:, 0], track[:, 1], track[:, 2], color="black",
                linewidth=2.3, linestyle="--", label="shared particle track")
    for marker in event_scene.get("impact_markers", ()):
        point = marker.get("point_mm", ())
        if len(point) >= 3:
            ax.scatter([point[0]], [point[1]], [point[2]], color="#00a65a", s=26,
                       marker="o", depthshade=False)

    # Add only short ownership labels.  Full section names stay out of the
    # plot and remain available through the interactive viewer details.
    if "cherenkov" in plotted_components:
        ax.plot([], [], [], color="#d95f02", linewidth=2.0, label="Cherenkov paths")
    if "scintillator" in plotted_components:
        ax.plot([], [], [], color="#1b9e77", linewidth=2.0, label="Scintillator paths")
    if any(str(item.get("kind")) == "sensor" for item in event_scene.get("geometry_components", ())):
        ax.plot([], [], [], color="#2ca25f", linewidth=2.0, label="sensor/readout")
    if not displayed:
        ax.text2D(0.02, 0.96, photon_message or "No photon paths selected.",
                  transform=ax.transAxes, fontsize=9, color="#8a1c1c",
                  verticalalignment="top")

    bounds = event_scene.get("bounds_mm")
    if isinstance(bounds, (list, tuple)) and len(bounds) >= 6:
        lows = np.asarray(bounds[:3], dtype=float)
        highs = np.asarray(bounds[3:6], dtype=float)
    else:
        points = [point for component in event_scene.get("geometry_components", ())
                  for face in component.get("faces_mm", ()) for point in face]
        points.extend(point for line in event_scene.get("photon_lines", ())
                      for point in line.get("points_mm", ()))
        points.extend(event_scene.get("track_mm", ()))
        if not points:
            points = [[0.0, 0.0, 0.0]]
        values = np.asarray(points, dtype=float)
        lows, highs = np.min(values, axis=0), np.max(values, axis=0)
    span = np.maximum(highs - lows, 1.0)
    pad = np.maximum(span * 0.04, 0.5)
    lows, highs = lows - pad, highs + pad
    ax.set_xlim(lows[0], highs[0]); ax.set_ylim(lows[1], highs[1]); ax.set_zlim(lows[2], highs[2])
    ax.set_box_aspect(tuple(np.maximum(highs - lows, 1e-12)))
    ax.set_xlabel("x (mm)"); ax.set_ylabel("y (mm)"); ax.set_zlabel("z (mm)")
    per_detector = event_scene.get("per_detector", {})
    detector_labels = (("cherenkov", "C"), ("scintillator", "S"))
    displayed_parts = []
    available_parts = []
    for detector_type, short_label in detector_labels:
        if detector_type not in per_detector:
            continue
        details = per_detector[detector_type]
        displayed_parts.append(f"{short_label} {displayed_by_owner.get(detector_type, 0)}")
        available_parts.append(f"{short_label} {details.get('record_count', 0)}")
    displayed_text = " / ".join(displayed_parts) or "none"
    available_text = " / ".join(available_parts) or "none"
    ax.set_title(
        f"Event {result.event_id}: detector mode {mode_labels[mode]}, "
        f"displayed {displayed_text} / available {available_text}"
    )
    legend = ax.legend(loc="upper left", fontsize=8)
    if legend is not None:
        legend.set_title(mode_labels[mode])
    ax.text2D(0.02, 0.02, "Canonical global scene • paths and readouts keep detector ownership",
              transform=ax.transAxes, fontsize=8)
    fig.savefig(output_path, dpi=160, bbox_inches="tight")
    if show:
        plt.show()
    plt.close(fig)


def plot_event_display(
    result: SimulationResult,
    output_path: str,
    display_photons: int = 40,
    show: bool = False,
    particle_track_margin_m: float = 0.020,
    photon_message: Optional[str] = None,
    event_scene: Optional[Mapping[str, object]] = None,
) -> None:
    """Save a static 3D event display with representative photon paths."""

    if event_scene is not None:
        _plot_canonical_event_scene(
            result, output_path, event_scene, display_photons, show, photon_message
        )
        return

    import matplotlib.pyplot as plt

    colors = {
        "detected_plus_x": "#2ca02c",
        "detected_minus_x": "#d627d8",
        "absorbed": "#d62728",
        "max_bounces": "#9467bd",
    }
    fig = plt.figure(figsize=(10, 7))
    ax = fig.add_subplot(111, projection="3d")
    mode = normalize_detector_run_mode(result.detector_mode)
    draw_cherenkov = mode in {"cherenkov", "combined"}
    draw_scintillator = mode in {"scintillator", "combined"}
    if draw_cherenkov:
        _draw_bar(ax, result.geometry)
    if draw_scintillator and result.scintillator is not None:
        scintillator_geometry = result.scintillator.assembly.transport_geometry
        if scintillator_geometry is result.geometry:
            if not draw_cherenkov:
                _draw_bar(ax, scintillator_geometry)
        else:
            _draw_bar(ax, scintillator_geometry)
    try:
        entry = _entry_point(result.geometry, result.particle)
    except ValueError:
        if result.scintillator is None:
            raise
        entry, _ = _scintillator_track(result.scintillator.assembly, result.particle)
    exit_point = entry + result.track_length_m * result.particle.direction
    track_start = entry - particle_track_margin_m * result.particle.direction
    track_end = exit_point + particle_track_margin_m * result.particle.direction
    track_scale = 1e3
    ax.plot([track_start[0] * track_scale, track_end[0] * track_scale], [track_start[1] * track_scale, track_end[1] * track_scale], [track_start[2] * track_scale, track_end[2] * track_scale], color="black", linewidth=2.5, label="particle track")
    ax.scatter([entry[0] * track_scale], [entry[1] * track_scale], [entry[2] * track_scale], color="#00aa00", marker="o", s=32, label="entry")
    ax.scatter([exit_point[0] * track_scale], [exit_point[1] * track_scale], [exit_point[2] * track_scale], color="#aa0000", marker="s", s=32, label="exit")
    all_records = list(result.photon_records) if draw_cherenkov else []
    if draw_scintillator and result.scintillator is not None:
        all_records.extend(result.scintillator.photon_records)
    displayed_records = all_records[:max(0, display_photons)]
    for record in displayed_records:
        points = np.asarray(record.path_points_m)
        color = colors.get(record.outcome, "#1f77b4" if record.outcome.startswith("escaped") else "#ff7f0e")
        ax.plot(points[:, 0] * track_scale, points[:, 1] * track_scale, points[:, 2] * track_scale, color=color, alpha=0.65, linewidth=0.8)
        if record.mirror_count:
            low_x = _geometry_bounds(result.geometry)[0][0]
            mirror_points = points[np.isclose(points[:, 0], low_x, atol=2e-9)]
            if len(mirror_points):
                ax.scatter(mirror_points[:, 0] * track_scale, mirror_points[:, 1] * track_scale, mirror_points[:, 2] * track_scale, color="#111111", marker="x", s=18)
    visible_outcomes = [record.outcome for record in displayed_records]
    legend_specs = (
        ("detected +x", "#2ca02c", lambda outcome: outcome == "detected_plus_x"),
        ("detected -x", "#d627d8", lambda outcome: outcome == "detected_minus_x"),
        ("escaped", "#1f77b4", lambda outcome: outcome.startswith("escaped")),
        ("absorbed", "#d62728", lambda outcome: outcome == "absorbed"),
        ("max-bounces", "#9467bd", lambda outcome: outcome == "max_bounces"),
        ("detector/TIR loss", "#ff7f0e", lambda outcome: outcome.startswith("detector_") or outcome.startswith("tir_")),
    )
    for label, color, predicate in legend_specs:
        if any(predicate(outcome) for outcome in visible_outcomes):
            ax.plot([], [], [], color=color, label=label)
    if any(outcome not in {"detected_plus_x", "detected_minus_x", "absorbed", "max_bounces"} and not outcome.startswith(("escaped", "detector_", "tir_")) for outcome in visible_outcomes):
        ax.plot([], [], [], color="#ff7f0e", linestyle="--", label="other/unknown outcome")
    if any(record.mirror_count for record in displayed_records):
        ax.plot([], [], [], color="#111111", marker="x", linestyle="None", label="mirror interaction")
    ax.set_xlabel("x (mm)")
    ax.set_ylabel("y (mm)")
    ax.set_zlabel("z (mm)")
    detected_count = sum(record.detected for record in all_records)
    mode_label = {"cherenkov": "Cherenkov only", "scintillator": "Scintillator only", "combined": "All"}.get(mode, mode)
    ax.set_title(
        f"Event {result.event_id}: detector mode {mode_label}, generated {result.generated_photons} "
        f"(+ scintillator {len(result.scintillator.photon_records) if result.scintillator is not None else 0}), "
        f"displayed {len(displayed_records)}, transported {len(all_records)}, detected {detected_count}"
    )
    ax.text2D(
        0.02,
        0.02,
        f"Physical expected yield: {result.physical_expected_photons:.3g} photons/event",
        transform=ax.transAxes,
        fontsize=8,
    )
    if not displayed_records:
        ax.text2D(
            0.02,
            0.96,
            photon_message or "No photon paths selected.",
            transform=ax.transAxes,
            fontsize=9,
            color="#8a1c1c",
            verticalalignment="top",
        )
    geometries = [result.geometry] if draw_cherenkov else []
    if draw_scintillator and result.scintillator is not None:
        scintillator_geometry = result.scintillator.assembly.transport_geometry
        if scintillator_geometry is not result.geometry or not geometries:
            geometries.append(scintillator_geometry)
    lows = np.minimum.reduce([_geometry_bounds(item)[0] for item in geometries])
    highs = np.maximum.reduce([_geometry_bounds(item)[1] for item in geometries])
    ax.set_xlim(lows[0] * 1e3, highs[0] * 1e3)
    ax.set_ylim(lows[1] * 1e3, highs[1] * 1e3)
    ax.set_zlim(lows[2] * 1e3, highs[2] * 1e3)
    ax.set_box_aspect(tuple(np.maximum(highs - lows, 1e-12)))
    ax.legend(loc="upper left", fontsize=8)
    fig.savefig(output_path, dpi=160)
    if show:
        plt.show()
    plt.close(fig)


def plot_cumulative(results_or_result: SimulationResult | Sequence[SimulationResult], output_path: str, show: bool = False) -> None:
    """Save detector-face y-z impacts from all events."""

    import matplotlib.pyplot as plt

    results = _as_results(results_or_result)
    fig, ax = plt.subplots(figsize=(7, 6))
    geometry = results[0].geometry if results else BarGeometry()
    low, high = _geometry_bounds(geometry)
    scale = 1e3
    detector_portals = not isinstance(geometry, SolidGeometry) or geometry.detector_portals_supported
    if detector_portals:
        ax.add_patch(plt.Rectangle((low[1] * scale, low[2] * scale), (high[1] - low[1]) * scale, (high[2] - low[2]) * scale, fill=False, linewidth=2, edgecolor="#2ca02c", label="+x detector face"))
        ax.add_patch(plt.Rectangle((low[1] * scale, low[2] * scale), (high[1] - low[1]) * scale, (high[2] - low[2]) * scale, fill=False, linewidth=2, linestyle="--", edgecolor="#d627d8", label="-x detector face"))
    else:
        ax.text(0.5, 0.96, "No planar +/-X detector portal for this solid/axis", transform=ax.transAxes, ha="center", va="top", fontsize=8, color="#8a1c1c")
    points = [(record.impact_position_m[1] * scale, record.impact_position_m[2] * scale, record.arrival_time_s * 1e9, record.detector_face) for result in results for record in result.photon_records if record.detected and not math.isnan(record.impact_position_m[0])]
    if points:
        y = np.asarray([point[0] for point in points])
        z = np.asarray([point[1] for point in points])
        faces = [point[3] for point in points]
        plus = np.asarray([face == "plus_x" for face in faces])
        ax.scatter(y[plus], z[plus], color="#2ca02c", marker="o", s=24, alpha=0.8, label="hits on +x")
        ax.scatter(y[~plus], z[~plus], color="#d627d8", marker="s", s=24, alpha=0.8, label="hits on -x")
    else:
        ax.text(0.5, 0.5, "No detected photons on either end face\n0 hits", transform=ax.transAxes, ha="center", va="center")
    ax.set_xlim(low[1] * scale * 1.08, high[1] * scale * 1.08)
    ax.set_ylim(low[2] * scale * 1.08, high[2] * scale * 1.08)
    ax.set_aspect("equal")
    ax.set_xlabel("detector y (mm)")
    ax.set_ylabel("detector z (mm)")
    ax.set_title(f"Cumulative end-face detector impacts, {len(points)} hits, {len(results)} event(s)")
    ax.legend(loc="upper right", fontsize=8)
    fig.tight_layout()
    fig.savefig(output_path, dpi=160)
    if show:
        plt.show()
    plt.close(fig)


def plot_photons_per_event(results_or_result: SimulationResult | Sequence[SimulationResult], output_path: str, show: bool = False) -> None:
    """Plot generated, transported, and detected counts for either sampling mode."""

    import matplotlib.pyplot as plt

    results = _as_results(results_or_result)
    event_ids = np.asarray([result.event_id for result in results], dtype=int)
    generated = np.asarray([result.generated_photons for result in results], dtype=float)
    transported = np.asarray([len(result.photon_records) for result in results], dtype=float)
    detected = np.asarray([sum(record.detected for record in result.photon_records) for result in results], dtype=float)
    width = 0.25
    positions = np.arange(len(results), dtype=float)
    fig, ax = plt.subplots(figsize=(9, 5))
    ax.bar(positions - width, generated, width, label="generated", color="#4472c4")
    ax.bar(positions, transported, width, label="transported", color="#70ad47")
    ax.bar(positions + width, detected, width, label="detected", color="#d62728")
    mode = results[0].optical.sampling_mode if results else "fixed"
    requested = results[0].optical.photon_count if results else 0
    expected = float(np.mean([result.physical_expected_photons for result in results])) if results else 0.0
    if mode == "poisson":
        ax.axhline(expected, color="#444444", linestyle="--", linewidth=1, label=f"physical Poisson mean ({expected:.2f})")
    else:
        ax.axhline(requested, color="#444444", linestyle="--", linewidth=1, label=f"fixed test sample ({requested})")
    ax.set_xticks(positions, [str(event_id) for event_id in event_ids])
    ax.set_xlabel("event_id")
    ax.set_ylabel("photons per event")
    ax.set_title(f"Photon counts per event: {mode} sampling")
    explanation = (
        f"Poisson mode: generated count fluctuates around physical mean {expected:.3g}."
        if mode == "poisson" else
        f"Fixed mode: {requested} is a controlled transport sample, not physical yield. Mean Frank-Tamm yield: {expected:.3g}."
    )
    ax.text(0.01, 0.98, explanation, transform=ax.transAxes, va="top", fontsize=8)
    ax.legend(fontsize=8)
    fig.tight_layout()
    fig.savefig(output_path, dpi=160)
    if show:
        plt.show()
    plt.close(fig)


def plot_arrival_times(results_or_result: SimulationResult | Sequence[SimulationResult], output_path: str, show: bool = False) -> None:
    """Plot detector arrival times relative to particle entry."""

    import matplotlib.pyplot as plt

    results = _as_results(results_or_result)
    records = [record for result in results for record in result.photon_records if record.detector_transmitted and not math.isnan(record.arrival_time_s)]
    fig, ax = plt.subplots(figsize=(8, 5))
    if records:
        ax.hist([record.arrival_time_s * 1e9 for record in records], bins=min(30, max(8, len(records) // 3)), color="#4472c4")
        ax.set_xlabel("detector arrival time relative to entry (ns)")
        ax.set_ylabel("photons")
        ax.set_title(f"Detector arrival times, {len(records)} crossings")
    else:
        ax.text(0.5, 0.5, "No detector-face crossings", transform=ax.transAxes, ha="center", va="center")
        ax.set_title("Detector arrival times, 0 crossings")
        ax.set_axis_off()
    fig.tight_layout()
    fig.savefig(output_path, dpi=160)
    if show:
        plt.show()
    plt.close(fig)


def plot_wavelengths(results_or_result: SimulationResult | Sequence[SimulationResult], output_path: str, show: bool = False) -> None:
    """Plot wavelengths for detected photons."""

    import matplotlib.pyplot as plt

    results = _as_results(results_or_result)
    records = [record for result in results for record in result.photon_records if record.detected]
    fig, ax = plt.subplots(figsize=(8, 5))
    if records:
        ax.hist([record.wavelength_m * 1e9 for record in records], bins=20, color="#70ad47")
        ax.set_xlabel("detected wavelength (nm)")
        ax.set_ylabel("photons")
        ax.set_title(f"Detected wavelength distribution, {len(records)} photons")
    else:
        ax.text(0.5, 0.5, "No detected photons", transform=ax.transAxes, ha="center", va="center")
        ax.set_title("Detected wavelength distribution, 0 photons")
        ax.set_axis_off()
    fig.tight_layout()
    fig.savefig(output_path, dpi=160)
    if show:
        plt.show()
    plt.close(fig)


def plot_outcome_face_counts(results_or_result: SimulationResult | Sequence[SimulationResult], output_path: str, show: bool = False) -> None:
    """Plot terminal outcomes and detector-face crossing counts."""

    import matplotlib.pyplot as plt

    results = _as_results(results_or_result)
    outcome_counts = Counter(record.outcome for result in results for record in result.photon_records)
    labels = list(sorted(outcome_counts))
    values = [outcome_counts[label] for label in labels]
    fig, axes = plt.subplots(1, 2, figsize=(12, 5), constrained_layout=True)
    if labels:
        axes[0].barh(labels, values, color="#8064a2")
        axes[0].set_xlabel("photons")
        axes[0].set_title("Terminal outcome counts")
    else:
        axes[0].text(0.5, 0.5, "No transported photons", transform=axes[0].transAxes, ha="center", va="center")
        axes[0].set_axis_off()
    crossing_counts = Counter(record.detector_face for result in results for record in result.photon_records if record.detector_transmitted)
    face_labels = ["plus_x", "minus_x"]
    axes[1].bar(face_labels, [crossing_counts[label] for label in face_labels], color=["#2ca02c", "#d627d8"])
    axes[1].set_ylabel("crossings")
    axes[1].set_title("Detector-face crossings")
    fig.savefig(output_path, dpi=160)
    if show:
        plt.show()
    plt.close(fig)


def plot_result(result: SimulationResult, output_path: str) -> None:
    """Save the compact outcome and detected-spectrum diagnostic."""

    import matplotlib.pyplot as plt

    records = result.photon_records
    fig, axes = plt.subplots(1, 2, figsize=(11, 4.2), constrained_layout=True)
    counts = Counter(record.outcome for record in records)
    axes[0].bar(list(counts.keys()), list(counts.values()), color="#4472c4")
    axes[0].set_ylabel("Photons")
    axes[0].set_title("Transport outcomes")
    axes[0].tick_params(axis="x", rotation=45)
    detected = [record for record in records if record.detected]
    if detected:
        axes[1].hist([record.wavelength_m * 1e9 for record in detected], bins=16, color="#70ad47")
        axes[1].set_xlabel("Detected wavelength (nm)")
        axes[1].set_ylabel("Photons")
    else:
        axes[1].text(0.5, 0.5, "No detected photons", ha="center", va="center")
        axes[1].set_axis_off()
    axes[1].set_title("Detected spectrum")
    fig.savefig(output_path, dpi=150)
    plt.close(fig)


def _cli_parser() -> argparse.ArgumentParser:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--self-test", action="store_true", help="run internal consistency tests and exit")
    parser.add_argument("--seed", type=int, default=12345)
    parser.add_argument("--events", type=int, default=1, help="number of independent events")
    parser.add_argument("--mass-mev-c2", type=float, default=0.51099895)
    parser.add_argument("--kinetic-mev", type=float, default=2.0)
    parser.add_argument("--entry-x-mm", type=float, default=0.0)
    parser.add_argument("--entry-y-mm", type=float, default=0.0, help="entry y coordinate on the z=+5 mm face")
    parser.add_argument("--inclination-deg", type=float, default=0.0, help="angle from +z; positive tilts toward +x when azimuth is zero")
    parser.add_argument("--azimuth-deg", type=float, default=0.0, help="azimuth around +z, zero toward +x, positive toward +y")
    parser.add_argument("--charge-e", type=float, default=1.0)
    parser.add_argument("--wavelength-min-nm", type=float, default=300.0)
    parser.add_argument("--wavelength-max-nm", type=float, default=650.0)
    parser.add_argument("--photon-count", type=int, default=1000)
    parser.add_argument("--sampling-mode", choices=("fixed", "poisson"), default="fixed", help="fixed Monte Carlo test sample or physical Poisson sampling from the Frank-Tamm yield")
    parser.add_argument("--dE-dx-MeV-per-m", type=float, default=None)
    parser.add_argument("--absorption-length-m", type=float, default=None)
    parser.add_argument("--rayleigh-length-m", type=float, default=None)
    parser.add_argument("--surface-roughness-deg", type=float, default=0.0)
    parser.add_argument("--max-bounces", type=int, default=100)
    parser.add_argument("--detector-acceptance-deg", type=float, default=30.0)
    parser.add_argument("--detector-qe", type=float, default=0.25)
    parser.add_argument("--collection-efficiency", type=float, default=0.80)
    parser.add_argument("--detector-index", type=float, default=1.0)
    parser.add_argument("--opposite-end-mode", choices=("open", "mirror", "detector"), default="open", help="behavior at the -x end face")
    parser.add_argument("--opposite-mirror-reflectivity", type=float, default=1.0, help="-x mirror reflectivity in mirror mode")
    parser.add_argument("--opposite-detector-acceptance-deg", type=float, default=None, help="-x detector acceptance, default inherited from +x")
    parser.add_argument("--opposite-detector-qe", type=float, default=None, help="-x detector QE, default inherited from +x")
    parser.add_argument("--opposite-collection-efficiency", type=float, default=None, help="-x detector collection efficiency, default inherited from +x")
    parser.add_argument("--opposite-detector-index", type=float, default=None, help="-x detector outside index, default inherited from +x")
    parser.add_argument("--plot", type=str, default=None, help="optional compact diagnostic PNG")
    parser.add_argument("--plot-event", type=str, default=None, help="optional static 3D event display PNG")
    parser.add_argument("--plot-cumulative", type=str, default=None, help="optional cumulative detector-impact PNG")
    parser.add_argument("--plot-photons-per-event", type=str, default=None, help="optional generated/transported/detected per-event PNG")
    parser.add_argument("--event-display-event", type=int, default=0, help="zero-based event index for --plot-event")
    parser.add_argument("--display-photons", type=int, default=40, help="maximum photon paths in the 3D display")
    parser.add_argument("--particle-track-margin-mm", type=float, default=20.0, help="track extension before entry and after exit")
    parser.add_argument("--show", action="store_true", help="show requested matplotlib plots interactively")
    parser.add_argument("--root-output", type=str, default=None, help="optional ROOT output filename")
    parser.add_argument("--check-root", action="store_true", help="check uproot and awkward availability and exit")
    parser.add_argument("--json", action="store_true", help="print summary as JSON")
    return parser


def _self_test() -> None:
    geometry = BarGeometry()
    assert math.isclose(geometry.length_x_m, 0.10)
    assert math.isclose(geometry.width_y_m, 0.020)
    assert math.isclose(geometry.thickness_z_m, 0.010)
    assert math.isclose(geometry.half_y, 0.010)
    assert math.isclose(geometry.half_z, 0.005)
    assert math.isclose(float(fused_silica_index(450e-9)), 1.46557, rel_tol=2e-4)
    assert not cherenkov_threshold(0.5)
    assert cherenkov_threshold(0.9)
    mass = 0.51099895 * MEV / C_LIGHT**2
    particle = Particle(mass_kg=mass, kinetic_energy_J=2.0 * MEV)
    assert np.allclose(particle.direction, [0.0, 0.0, 1.0])
    assert np.allclose(_entry_point(geometry, particle), [0.0, 0.0, -geometry.half_z])
    assert math.isclose(_track_length(geometry, particle), geometry.thickness_z_m)
    optical = OpticalConfig(photon_count=24, max_bounces=20, absorption_length_m=0.5)
    first = simulate(particle, geometry, optical, seed=7)
    second = simulate(particle, geometry, optical, seed=7)
    assert first.summary() == second.summary()
    assert first.generated_photons == 24
    assert all(record.path_length_m >= 0.0 for record in first.photon_records)
    assert all(record.emission_time_s >= 0.0 for record in first.photon_records)
    assert any(record.emission_time_s > 0.0 for record in first.photon_records)
    for record in first.photon_records:
        if record.detected:
            assert record.arrival_time_s > record.emission_time_s
            assert math.isclose(record.impact_position_m[0], geometry.half_x, abs_tol=1e-9)
        else:
            assert all(math.isnan(value) for value in record.impact_position_m)
        if record.detector_transmitted:
            assert record.arrival_time_s > record.emission_time_s
        else:
            assert math.isnan(record.arrival_time_s)

    # Regression geometry from the GUI screenshot: a 100 x 20 x 10 mm box
    # followed by a 20 mm tapered section whose Y width changes 20 -> 10 mm.
    # The common cap is a true internal interface, not two external surfaces.
    compound_box = SolidGeometry(
        shape="box", material="fused_silica", center_x_m=0.0,
        length_x_m=0.100, width_y_m=0.020, thickness_z_m=0.010,
    )
    compound_taper = SolidGeometry(
        shape="tapered_rectangular_prism", material="fused_silica", center_x_m=0.060,
        length_x_m=0.020, width_y_start_m=0.020, width_y_end_m=0.010,
        thickness_z_start_m=0.010, thickness_z_end_m=0.010,
    )
    compound = CompoundGeometry((compound_box, compound_taper))
    compound.validate()
    compound_optical = OpticalConfig(
        photon_count=100, max_bounces=20, detector_index=1.4,
        detector_acceptance_rad=math.pi / 2.0,
        detector_quantum_efficiency=1.0, detector_collection_efficiency=1.0,
        detector_active_width_m=0.010, detector_active_height_m=0.010,
    )
    compound_result = simulate(particle, compound, compound_optical, seed=67)
    assert compound_result.generated_photons == 100
    assert all(not record.outcome.startswith("transport_error") for record in compound_result.photon_records)

    # Explicit +X and oblique rays both cross box -> tapered without a
    # spurious Fresnel bounce when both sections are fused silica.
    interface_optical = OpticalConfig(
        max_bounces=40, detector_index=1.4, detector_acceptance_rad=math.pi / 2.0,
        detector_quantum_efficiency=1.0, detector_collection_efficiency=1.0,
    )
    for ray_direction in (np.array([1.0, 0.0, 0.0]), np.array([0.99, 0.05, 0.05])):
        ray_record = _transport_photon(
            0, 450e-9, np.array([0.0, 0.0, 0.0]), ray_direction,
            0.0, 0.99, compound, interface_optical, np.random.default_rng(71),
        )
        assert ray_record.detector_transmitted
        assert ray_record.detector_face == "plus_x"
        assert ray_record.bounces == 0

    # A point nudged just outside a tapered terminal side has no future
    # intersection.  It is an external escape, not a transport exception.
    air_box = SolidGeometry(
        shape="box", material="air", center_x_m=0.0,
        length_x_m=0.100, width_y_m=0.020, thickness_z_m=0.010,
    )
    air_taper = SolidGeometry(
        shape="tapered_rectangular_prism", material="air", center_x_m=0.060,
        length_x_m=0.020, width_y_start_m=0.020, width_y_end_m=0.010,
        thickness_z_start_m=0.010, thickness_z_end_m=0.010,
    )
    air_compound = CompoundGeometry((air_box, air_taper))
    boundary_point = np.array([0.0528, -0.0093, 0.0])
    local_normal = air_taper.surface_normal(boundary_point, "y-")
    escaped_start = boundary_point + 0.5e-9 * local_normal
    escaped_record = _transport_photon(
        0, 450e-9, escaped_start, local_normal, 0.0, 0.99,
        air_compound, OpticalConfig(max_bounces=10), np.random.default_rng(72),
    )
    assert escaped_record.outcome == "escaped_y-"

    # Different materials retain Fresnel/Snell and TIR at the same internal
    # cap.  The helper checks cover transmitted Snell refraction and critical
    # angle handling; the direct ray exercises an actual internal TIR bounce.
    reflectance, cos_t = _fresnel_unpolarized(0.8, 1.0, 1.5)
    assert 0.0 <= reflectance < 1.0 and cos_t > 0.0
    refracted = _refracted_direction(np.array([0.8, 0.0, 0.6]), np.array([1.0, 0.0, 0.0]), 1.0, 1.5, cos_t)
    assert math.isclose(np.linalg.norm(refracted), 1.0, rel_tol=1e-12)
    tir_reflectance, tir_cos_t = _fresnel_unpolarized(0.5, 1.5, 1.0)
    assert tir_reflectance == 1.0 and tir_cos_t == 0.0
    different_material = CompoundGeometry((
        compound_box,
        SolidGeometry(
            shape="tapered_rectangular_prism", material="air", center_x_m=0.060,
            length_x_m=0.020, width_y_start_m=0.020, width_y_end_m=0.010,
            thickness_z_start_m=0.010, thickness_z_end_m=0.010,
        ),
    ))
    tir_record = _transport_photon(
        0, 450e-9, np.array([0.049, 0.0, 0.0]), np.array([0.5, 0.0, math.sqrt(0.75)]),
        0.0, 0.99, different_material, interface_optical, np.random.default_rng(73),
    )
    assert tir_record.bounces >= 1

    # Real geometry validation remains strict for gaps, overlaps, and unknown
    # materials; the boundary fix must not mask editor configuration errors.
    with_gap = SolidGeometry(
        shape="tapered_rectangular_prism", center_x_m=0.061, length_x_m=0.020,
        width_y_start_m=0.020, width_y_end_m=0.010,
        thickness_z_start_m=0.010, thickness_z_end_m=0.010,
    )
    with_overlap = SolidGeometry(
        shape="tapered_rectangular_prism", center_x_m=0.059, length_x_m=0.020,
        width_y_start_m=0.020, width_y_end_m=0.010,
        thickness_z_start_m=0.010, thickness_z_end_m=0.010,
    )
    for invalid_geometry, expected in (
        (CompoundGeometry((compound_box, with_gap)), "gap"),
        (CompoundGeometry((compound_box, with_overlap)), "overlap"),
        (CompoundGeometry((compound_box, SolidGeometry(center_x_m=0.060, material="unknown"))), "unsupported solid material"),
    ):
        try:
            invalid_geometry.validate()
        except (TypeError, ValueError) as exc:
            assert expected in str(exc)
        else:
            raise AssertionError(f"invalid compound was accepted: {expected}")

    air_tir = simulate(
        particle,
        optical=OpticalConfig(
            photon_count=64,
            detector_index=1.0,
            detector_acceptance_rad=math.pi / 2.0,
            detector_quantum_efficiency=1.0,
            detector_collection_efficiency=1.0,
            max_bounces=60,
        ),
        seed=41,
    )
    assert air_tir.summary()["tir_plus_x"] > 0
    assert any(record.tir_plus_x > 0 for record in air_tir.photon_records)
    assert air_tir.summary()["detector_interface_reflections"] == 0
    opposite_air_tir = simulate(
        particle,
        optical=OpticalConfig(
            photon_count=64,
            opposite_end_mode="detector",
            opposite_detector_index=1.0,
            opposite_detector_acceptance_rad=math.pi / 2.0,
            opposite_detector_quantum_efficiency=1.0,
            opposite_detector_collection_efficiency=1.0,
            max_bounces=60,
        ),
        seed=45,
    )
    assert opposite_air_tir.summary()["tir_minus_x"] > 0
    coupled = simulate(
        particle,
        optical=OpticalConfig(
            photon_count=64,
            detector_index=1.4,
            detector_acceptance_rad=math.pi / 2.0,
            detector_quantum_efficiency=1.0,
            detector_collection_efficiency=1.0,
        ),
        seed=43,
    )
    assert coupled.summary()["detector_face_crossings"] > 0
    assert coupled.summary()["detected_plus_x"] > 0
    below = simulate(Particle(mass_kg=mass, kinetic_energy_J=0.05 * MEV), optical=OpticalConfig(photon_count=10), seed=1)
    assert below.generated_photons == 0
    poisson_optical = OpticalConfig(photon_count=24, sampling_mode="poisson")
    poisson_first = simulate(particle, optical=poisson_optical, seed=52)
    poisson_second = simulate(particle, optical=poisson_optical, seed=52)
    assert poisson_first.summary() == poisson_second.summary()
    assert poisson_first.generated_photons == len(poisson_first.photon_records)
    assert poisson_first.summary()["sampling_mode"] == "poisson"
    assert math.isclose(poisson_first.summary()["poisson_mean"], poisson_first.physical_expected_photons)
    poisson_zero = simulate(
        Particle(mass_kg=mass, kinetic_energy_J=0.05 * MEV),
        optical=OpticalConfig(photon_count=24, sampling_mode="poisson"),
        seed=53,
    )
    assert poisson_zero.generated_photons == 0
    assert poisson_zero.photon_records == []
    assert poisson_zero.summary()["physical_expected_photons"] == 0.0
    multi = simulate_events(particle, 2, optical=OpticalConfig(photon_count=3), seed=19)
    assert [result.event_id for result in multi] == [0, 1]
    assert multi[0].photon_records[0].photon_id == 0
    mirror_result = simulate(
        particle,
        optical=OpticalConfig(photon_count=32, max_bounces=40, opposite_end_mode="mirror", opposite_mirror_reflectivity=1.0),
        seed=31,
    )
    assert any(record.mirror_count > 0 for record in mirror_result.photon_records)
    assert any(any(point[0] <= geometry.half_x * -0.999999 for point in record.path_points_m) for record in mirror_result.photon_records if record.mirror_count > 0)
    opposite_detector_result = simulate(
        particle,
        optical=OpticalConfig(
            photon_count=64,
            max_bounces=200,
            opposite_end_mode="detector",
            opposite_detector_index=1.4,
            opposite_detector_acceptance_rad=math.pi / 2.0,
            opposite_detector_quantum_efficiency=1.0,
            opposite_detector_collection_efficiency=1.0,
        ),
        seed=37,
    )
    minus_hits = [record for record in opposite_detector_result.photon_records if record.detector_face == "minus_x" and record.detected]
    assert minus_hits
    assert all(math.isclose(record.impact_position_m[0], -geometry.half_x, abs_tol=1e-9) for record in minus_hits)
    assert opposite_detector_result.summary()["detected_minus_x"] == len(minus_hits)

    # Beam-monitor extension checks: separate volumes, shape combinations,
    # strict gap/Mylar validation, fixed and explicit constant-dE/dx modes.
    scintillator = SolidGeometry(
        shape="box", material="scintillator", center_x_m=0.010,
        length_x_m=0.020, width_y_m=0.020, thickness_z_m=0.010,
        refractive_index_value=1.58,
    )
    light_guide = SolidGeometry(
        shape="cylinder", material="light_guide", axis="X", center_x_m=0.040,
        length_m=0.040, radius_m=0.010, refractive_index_value=1.49,
    )
    assembly = ScintillatorAssembly(scintillator, light_guide, gap_from_quartz_m=0.002, mylar_reflectivity=0.80)
    assembly.validate()
    local_model = ScintillatorLocalGeometry(
        replace(scintillator, center_x_m=0.010), replace(light_guide, center_x_m=0.040)
    )
    global_model = ScintillatorGlobalGeometry(gap_from_quartz_m=0.002)
    composed_model = compose_scintillator_assembly(local_model, global_model, geometry)
    assert math.isclose(local_model.scintillator.length_x_m, scintillator.length_x_m)
    assert math.isclose(composed_model.gap_from_quartz_m, global_model.gap_from_quartz_m)
    assert not math.isclose(composed_model.scintillator.center_x_m, local_model.scintillator.center_x_m)
    assert assembly.relative_placement()["coupling_face_wrapped"] is False
    assert "scintillator_light_guide" not in assembly.wrapped_faces
    try:
        ScintillatorAssembly(scintillator, light_guide, gap_from_quartz_m=-1.0).validate()
    except ValueError:
        pass
    else:
        raise AssertionError("negative scintillator gap was accepted")
    try:
        ScintillatorAssembly(scintillator, light_guide, mylar_reflectivity=0.79).validate()
    except ValueError:
        pass
    else:
        raise AssertionError("invalid Mylar reflectivity was accepted")
    fixed_scintillator = simulate_scintillator(
        particle, assembly, ScintillatorConfig(enabled=True, photon_count=8, sampling_mode="fixed", max_bounces=8), seed=77,
    )
    assert len(fixed_scintillator.photon_records) == 8
    physical_scintillator = simulate_scintillator(
        particle, assembly,
        ScintillatorConfig(enabled=True, sampling_mode="physical", photon_yield_per_MeV=3.0, deposited_energy_MeV=2.0, max_bounces=8),
        seed=78,
    )
    assert physical_scintillator.deposited_energy_MeV == 2.0
    assert math.isclose(physical_scintillator.physical_expected_photons, 6.0)
    assert set(physical_scintillator.summary()) >= {"scintillator_mylar_reflections", "scintillator_mylar_losses"}

    # Two-detector engineering model: each detector owns one scintillator and
    # two distinct guide/readout solids, while placements remain signed.
    left_guide = replace(light_guide, center_x_m=-0.020)
    right_guide = replace(light_guide, center_x_m=0.040)
    detector_1 = DetectorAssembly(scintillator, left_guide, right_guide, detector_id=1, mylar_reflectivity=0.80)
    detector_2 = DetectorAssembly(
        replace(scintillator, material="scintillator", refractive_index_value=1.57),
        replace(left_guide, shape="cylinder"), replace(right_guide, shape="elliptical_cylinder"),
        detector_id=2, mylar_reflectivity=1.00,
    )
    monitor = TwoDetectorBeamMonitor(
        detector_1, detector_2, quartz_gap_m=0.002, detector_separation_m=0.010,
        placement_1=DetectorPlacement(offset_x_m=-0.001, offset_y_m=-0.002, offset_z_m=0.001),
        placement_2=DetectorPlacement(offset_x_m=-0.003, offset_y_m=0.002, offset_z_m=-0.001),
    )
    monitor.validate(geometry)
    placed_1, placed_2 = monitor.placed_detectors(geometry)
    assert len(placed_1.preview()["volumes"]) == 3 and len(placed_2.preview()["volumes"]) == 3
    assert math.isclose(placed_1.scintillator.length_x_m, detector_1.scintillator.length_x_m)
    assert math.isclose(placed_2.scintillator.length_x_m, detector_2.scintillator.length_x_m)
    path = monitor.particle_path(geometry)
    assert "Detector 1 scintillator" in path["order"] and "Detector 2 scintillator" in path["order"]
    assert len(path["points_mm"]) >= 8
    assert monitor.relative_placement(geometry)["detector_1_signed_offset_mm"][0] < 0.0

    # Definitive beam-monitor representation: two independent XY active
    # planes at distinct Z values, with optional Y-face guides and signed
    # lateral placement. This is an engineering scene model; established
    # single-assembly optical transport remains the physics path above.
    z_active = SolidGeometry(shape="box", material="scintillator", length_x_m=0.040, width_y_m=0.030, thickness_z_m=0.010, refractive_index_value=1.58)
    z_guide_a = SolidGeometry(shape="box", material="light_guide", center_y_m=-0.025, length_x_m=0.040, width_y_m=0.020, thickness_z_m=0.010, refractive_index_value=1.49)
    z_guide_b = SolidGeometry(shape="box", material="light_guide", center_y_m=0.025, length_x_m=0.040, width_y_m=0.020, thickness_z_m=0.010, refractive_index_value=1.49)
    z_detector_1 = DetectorAssembly(z_active, z_guide_a, z_guide_b, detector_id=1, active_axis="Z")
    z_detector_2 = DetectorAssembly(replace(z_active, refractive_index_value=1.57), None, z_guide_b, detector_id=2, active_axis="Z")
    z_monitor = TwoDetectorBeamMonitor(
        z_detector_1, z_detector_2, quartz_gap_m=0.010, detector_separation_m=0.030,
        placement_1=DetectorPlacement(offset_x_m=-0.002, offset_y_m=0.001),
        placement_2=DetectorPlacement(offset_x_m=0.003, offset_y_m=-0.001),
        beam_axis="Z", detector_1_z_m=0.020, detector_2_z_m=0.070,
    )
    z_monitor.validate(geometry)
    z1, z2 = z_monitor.placed_detectors(geometry)
    assert z1.active_axis == z2.active_axis == "Z" and z1.guide_left is not None and z2.guide_left is None
    assert math.isclose(z1.scintillator.width_y_m, z_detector_1.scintillator.width_y_m)
    assert z_monitor.relative_placement(geometry)["detector_1_z_mm"] != z_monitor.relative_placement(geometry)["detector_2_z_mm"]
    z_path = z_monitor.particle_path(geometry)
    assert "XY active plane" in z_path["order"] and len(z_path["points_mm"]) >= 8
    DetectorSensor(component="scintillator_1", coupling_mode="direct").validate({"scintillator_1"})
    DetectorSensor(component="scintillator_1", coupling_mode="light_guide").validate({"scintillator_1"})
    try:
        with tempfile.TemporaryDirectory(prefix="quartz_cherenkov_self_test_") as temporary_directory:
            root_output = Path(temporary_directory) / "quartz_cherenkov_self_test.root"
            root_path = write_root(multi, root_output)
            import uproot
            with uproot.open(root_path) as root_file:
                assert "Photons" in root_file and "Events" in root_file
                assert root_file["Events"].num_entries == len(multi)
        print("self-test: ROOT writing and verification available")
    except RuntimeError as exc:
        assert "uproot and awkward" in str(exc)
        print("self-test: ROOT writing skipped (optional uproot/awkward unavailable)")
    print("self-test: PASS")


def main(argv: Optional[Sequence[str]] = None) -> int:
    args = _cli_parser().parse_args(argv)
    if args.check_root:
        try:
            print(check_root_dependencies())
        except RuntimeError as exc:
            raise SystemExit(str(exc)) from exc
        return 0
    if args.self_test:
        _self_test()
        return 0
    if args.events < 1:
        raise SystemExit("--events must be at least 1")
    if args.particle_track_margin_mm < 0.0:
        raise SystemExit("--particle-track-margin-mm must be non-negative")
    particle = Particle(
        mass_kg=args.mass_mev_c2 * MEV / C_LIGHT**2,
        kinetic_energy_J=args.kinetic_mev * MEV,
        entry_x_m=args.entry_x_mm * 1e-3,
        entry_y_m=args.entry_y_mm * 1e-3,
        inclination_rad=math.radians(args.inclination_deg),
        azimuth_rad=math.radians(args.azimuth_deg),
        charge_e=args.charge_e,
        dEdx_J_per_m=None if args.dE_dx_MeV_per_m is None else args.dE_dx_MeV_per_m * MEV,
    )
    optical = OpticalConfig(
        wavelength_min_m=args.wavelength_min_nm * 1e-9,
        wavelength_max_m=args.wavelength_max_nm * 1e-9,
        photon_count=args.photon_count,
        sampling_mode=args.sampling_mode,
        absorption_length_m=args.absorption_length_m,
        rayleigh_scattering_length_m=args.rayleigh_length_m,
        surface_roughness_rad=math.radians(args.surface_roughness_deg),
        max_bounces=args.max_bounces,
        detector_acceptance_rad=math.radians(args.detector_acceptance_deg),
        detector_quantum_efficiency=args.detector_qe,
        detector_collection_efficiency=args.collection_efficiency,
        detector_index=args.detector_index,
        opposite_end_mode=args.opposite_end_mode,
        opposite_mirror_reflectivity=args.opposite_mirror_reflectivity,
        opposite_detector_acceptance_rad=None if args.opposite_detector_acceptance_deg is None else math.radians(args.opposite_detector_acceptance_deg),
        opposite_detector_quantum_efficiency=args.opposite_detector_qe,
        opposite_detector_collection_efficiency=args.opposite_collection_efficiency,
        opposite_detector_index=args.opposite_detector_index,
    )
    results = simulate_events(particle, args.events, optical=optical, seed=args.seed)
    if args.plot:
        plot_result(results[0], args.plot)
    if args.plot_event:
        if not 0 <= args.event_display_event < len(results):
            raise SystemExit(f"--event-display-event must be in [0, {len(results) - 1}]")
        plot_event_display(results[args.event_display_event], args.plot_event, args.display_photons, args.show, args.particle_track_margin_mm * 1e-3)
    if args.plot_cumulative:
        plot_cumulative(results, args.plot_cumulative, args.show)
    if args.plot_photons_per_event:
        plot_photons_per_event(results, args.plot_photons_per_event, args.show)
    if args.root_output:
        write_root(results, args.root_output)
    print(json.dumps(aggregate_summary(results), indent=2, sort_keys=True))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
