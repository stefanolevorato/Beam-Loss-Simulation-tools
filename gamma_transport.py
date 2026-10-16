#!/usr/bin/env python3
"""First usable gamma-transport stage for the frozen v62 beam monitor.

This module is deliberately separate from the established optical transport.
It is a small, reproducible photon Monte Carlo for the 12.5 keV--2.5 MeV
profile used by the project.  It includes Compton scattering, photoelectric
absorption, pair production above threshold, straight-line CSDA transport of
the resulting electrons/positrons, and local Cherenkov/scintillation light
production from those charged tracks.

The interaction coefficients are coarse effective mass-attenuation tables for
fused silica and plastic/PMMA-like materials.  They are a first-stage model,
not a Geant4 replacement: atomic relaxation, fluorescence/Auger cascades,
multiple-scattering electron transport, optical photon ray transport, and
the detector/electronics response are intentionally outside this module.
Every result carries that limitation in its ``model_status`` field.
"""

from __future__ import annotations

import argparse
import json
import math
from collections import Counter, defaultdict
from dataclasses import asdict, dataclass, field
from pathlib import Path
from typing import Any, Iterable, Mapping, Optional, Sequence

import numpy as np

try:  # Package and script execution are both useful in the sandbox.
    from . import quartz_cherenkov_sim as sim
except ImportError:  # pragma: no cover - exercised by the CLI script path
    import quartz_cherenkov_sim as sim


ME_C2_KEV = 510.99895
PAIR_THRESHOLD_KEV = 2.0 * ME_C2_KEV
MODEL_STATUS = "first_stage_effective_cross_sections_no_geant4"

# Energies are MeV; values are cm^2/g.  These are intentionally compact
# effective tables (log-log interpolation) rather than silently pretending to
# be element-by-element evaluated databases.
_ENERGY_MEV = np.array([0.0125, 0.025, 0.05, 0.10, 0.20, 0.30, 0.50, 0.80, 1.0, 1.5, 2.0, 2.5])
_TABLES = {
    "fused_silica": {
        "photoelectric": np.array([7.2, 1.05, 0.18, 0.026, 0.0056, 0.0022, 0.00062, 0.00023, 0.00015, 0.000075, 0.000048, 0.000034]),
        "compton": np.array([0.115, 0.105, 0.094, 0.082, 0.069, 0.061, 0.050, 0.041, 0.037, 0.031, 0.027, 0.024]),
        "pair": np.array([0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.00005, 0.00030, 0.00055, 0.00135, 0.00205, 0.00270]),
    },
    "plastic": {
        "photoelectric": np.array([2.4, 0.34, 0.058, 0.0085, 0.0018, 0.00070, 0.00020, 0.000075, 0.000050, 0.000025, 0.000016, 0.000011]),
        "compton": np.array([0.105, 0.097, 0.087, 0.076, 0.064, 0.057, 0.047, 0.039, 0.035, 0.030, 0.026, 0.023]),
        "pair": np.array([0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.00002, 0.00012, 0.00022, 0.00055, 0.00082, 0.00105]),
    },
}
_DENSITY_G_CM3 = {"fused_silica": 2.20, "plastic": 1.032, "pmma": 1.18, "air": 0.0012, "custom": 1.0}


def _json_number(value: Any) -> float:
    number = float(value)
    if not math.isfinite(number):
        raise ValueError(f"non-finite numeric value: {value!r}")
    return number


def _unit(v: Sequence[float]) -> np.ndarray:
    arr = np.asarray(v, dtype=float)
    norm = float(np.linalg.norm(arr))
    if norm <= 0.0 or not np.all(np.isfinite(arr)):
        raise ValueError("direction must be finite and non-zero")
    return arr / norm


def _basis(direction: np.ndarray) -> tuple[np.ndarray, np.ndarray]:
    direction = _unit(direction)
    helper = np.array([0.0, 0.0, 1.0]) if abs(float(direction[2])) < 0.9 else np.array([0.0, 1.0, 0.0])
    first = np.cross(direction, helper)
    first /= np.linalg.norm(first)
    second = np.cross(direction, first)
    return first, second


def _random_direction(rng: np.random.Generator) -> np.ndarray:
    return _unit(rng.normal(size=3))


def _material_key(material: str) -> str:
    text = str(material).strip().lower().replace(" ", "_")
    if text in {"fused_silica", "quartz"}:
        return "fused_silica"
    if text in {"scintillator", "plastic_scintillator", "light_guide", "pmma", "plastic"}:
        return "plastic"
    if text == "air":
        return "air"
    return "custom"


def mass_attenuation_coefficients(energy_keV: float, material: str) -> dict[str, float]:
    """Return effective mass coefficients in cm^2/g for one material."""

    energy_mev = max(float(energy_keV) / 1000.0, float(_ENERGY_MEV[0]))
    key = _material_key(material)
    if key == "air":
        # Air is only a transport gap here; its interactions are negligible
        # over the dimensions of the supplied assemblies.
        return {"photoelectric": 0.0, "compton": 0.0, "pair": 0.0}
    table = _TABLES["fused_silica" if key == "fused_silica" else "plastic"]
    out: dict[str, float] = {}
    for name, values in table.items():
        positive = np.maximum(values, 1e-12)
        value = math.exp(float(np.interp(math.log(energy_mev), np.log(_ENERGY_MEV), np.log(positive))))
        out[name] = 0.0 if name == "pair" and energy_keV < PAIR_THRESHOLD_KEV else value
    return out


def electron_range_g_cm2(kinetic_energy_keV: float) -> float:
    """Approximate electron CSDA range using a compact empirical relation."""

    energy_mev = max(float(kinetic_energy_keV) / 1000.0, 1e-6)
    # Legacy expression retained pending an authoritative range reference.
    # See gamma_range_model_diagnostic.md; do not treat the attribution or
    # output as reference-validated. The lower-energy clamp avoids the
    # logarithmic singularity in the current expression.
    value = 0.412 * energy_mev ** 1.265 - 0.0954 * math.log(max(energy_mev, 0.001))
    return max(1e-5, value)


def _beta_from_kinetic_keV(kinetic_keV: float) -> float:
    gamma = 1.0 + max(float(kinetic_keV), 0.0) / ME_C2_KEV
    return math.sqrt(max(0.0, 1.0 - 1.0 / (gamma * gamma)))


def _sample_compton_cosine(energy_keV: float, rng: np.random.Generator) -> float:
    """Rejection sample the Klein-Nishina angular density in cos(theta)."""

    epsilon = max(float(energy_keV), 1e-9) / ME_C2_KEV
    for _ in range(10000):
        cosine = float(rng.uniform(-1.0, 1.0))
        ratio = 1.0 / (1.0 + epsilon * (1.0 - cosine))
        sine2 = max(0.0, 1.0 - cosine * cosine)
        weight = ratio * ratio * (ratio + 1.0 / ratio - sine2)
        if float(rng.random()) < weight / 2.0:
            return cosine
    return 1.0


def _scatter_direction(direction: np.ndarray, cosine: float, rng: np.random.Generator) -> np.ndarray:
    first, second = _basis(direction)
    sine = math.sqrt(max(0.0, 1.0 - cosine * cosine))
    phi = float(rng.uniform(0.0, 2.0 * math.pi))
    return _unit(cosine * direction + sine * (math.cos(phi) * first + math.sin(phi) * second))


@dataclass(frozen=True)
class SourceLine:
    energy_keV: float
    relative_weight: float
    label: str = "user line"


@dataclass(frozen=True)
class GammaSource:
    bins: tuple[tuple[float, float], ...]
    lines: tuple[SourceLine, ...] = ()
    continuum_weight_scale: float = 1.0
    source_profile_path: str = ""
    source_direction: tuple[float, float, float] = (0.0, 0.0, 1.0)

    @classmethod
    def from_profile(cls, path: str | Path, line_weights: Optional[Mapping[float, float]] = None) -> "GammaSource":
        payload = json.loads(Path(path).read_text(encoding="utf-8"))
        if payload.get("schema") != "gamma-source-profile-v1":
            raise ValueError("unsupported gamma source profile schema")
        spectrum = payload.get("spectrum", {})
        bins = []
        for row in spectrum.get("bins", []):
            energy = _json_number(row["energy_center_keV"])
            weight = _json_number(row["relative_weight"])
            if energy > 0.0 and weight > 0.0:
                bins.append((energy, weight))
        if not bins:
            raise ValueError("gamma source profile has no positive continuum bins")
        provided = {float(k): _json_number(v) for k, v in (line_weights or {}).items()}
        lines = []
        for row in payload.get("line_components", []):
            energy = _json_number(row["energy_keV"])
            if energy in provided:
                weight = provided[energy]
                if weight < 0.0:
                    raise ValueError(f"line weight for {energy:g} keV must be non-negative")
                if weight > 0.0:
                    lines.append(SourceLine(energy, weight, str(row.get("label", "user line"))))
        geometry = payload.get("source_geometry", {})
        direction = tuple(float(value) for value in geometry.get("direction_global", (0.0, 0.0, 1.0)))
        _unit(direction)
        return cls(
            tuple(bins), tuple(lines), float(payload.get("continuum_component", {}).get("weight_scale", 1.0)),
            str(path), tuple(direction),
        )

    def sample_energy_keV(self, rng: np.random.Generator) -> tuple[float, str]:
        continuum_weight = max(0.0, self.continuum_weight_scale) * sum(weight for _, weight in self.bins)
        line_weight = sum(line.relative_weight for line in self.lines)
        if line_weight > 0.0 and float(rng.random()) < line_weight / (continuum_weight + line_weight):
            weights = np.asarray([line.relative_weight for line in self.lines], dtype=float)
            index = int(rng.choice(len(self.lines), p=weights / weights.sum()))
            return self.lines[index].energy_keV, "line"
        weights = np.asarray([weight for _, weight in self.bins], dtype=float)
        index = int(rng.choice(len(self.bins), p=weights / weights.sum()))
        return self.bins[index][0], "continuum"


@dataclass(frozen=True)
class TransportVolume:
    name: str
    assembly: str
    component: str
    solid: Any
    material: str
    scintillation_yield_photons_per_MeV: float = 10000.0

    @property
    def density_g_cm3(self) -> float:
        return _DENSITY_G_CM3[_material_key(self.material)]


@dataclass(frozen=True)
class RigidTransformedSolid:
    """Ray-query adapter for a solid placed by one rigid global transform.

    The frozen v62 optical geometry stores solids in an assembly-local frame.
    Gamma transport needs the same solids in the global +Z frame.  Queries are
    transformed into the local frame, so the underlying v62 shape algorithms
    remain the geometry authority and no approximate replacement shape is
    introduced.  The bounds are conservative AABBs of the transformed local
    bounds; exact intersections still come from the wrapped solid.
    """

    base: Any
    rotation_matrix: tuple[tuple[float, float, float], ...] = ((1.0, 0.0, 0.0), (0.0, 1.0, 0.0), (0.0, 0.0, 1.0))
    translation_m: tuple[float, float, float] = (0.0, 0.0, 0.0)
    rotation_origin_m: tuple[float, float, float] = (0.0, 0.0, 0.0)

    def _matrix(self) -> np.ndarray:
        return np.asarray(self.rotation_matrix, dtype=float)

    def _translation(self) -> np.ndarray:
        return np.asarray(self.translation_m, dtype=float)

    def _origin(self) -> np.ndarray:
        return np.asarray(self.rotation_origin_m, dtype=float)

    def _to_local_point(self, point: Sequence[float]) -> np.ndarray:
        matrix = self._matrix()
        return matrix.T.dot(np.asarray(point, dtype=float) - self._origin() - self._translation()) + self._origin()

    def _to_local_direction(self, direction: Sequence[float]) -> np.ndarray:
        return self._matrix().T.dot(np.asarray(direction, dtype=float))

    @property
    def bounds_min_m(self) -> np.ndarray:
        low = np.asarray(self.base.bounds_min_m, dtype=float)
        high = np.asarray(self.base.bounds_max_m, dtype=float)
        corners = np.asarray([
            (x, y, z)
            for x in (low[0], high[0])
            for y in (low[1], high[1])
            for z in (low[2], high[2])
        ])
        transformed = (self._matrix().dot((corners - self._origin()).T)).T + self._origin() + self._translation()
        return np.min(transformed, axis=0)

    @property
    def bounds_max_m(self) -> np.ndarray:
        low = np.asarray(self.base.bounds_min_m, dtype=float)
        high = np.asarray(self.base.bounds_max_m, dtype=float)
        corners = np.asarray([
            (x, y, z)
            for x in (low[0], high[0])
            for y in (low[1], high[1])
            for z in (low[2], high[2])
        ])
        transformed = (self._matrix().dot((corners - self._origin()).T)).T + self._origin() + self._translation()
        return np.max(transformed, axis=0)

    @property
    def characteristic_length_m(self) -> float:
        return float(self.base.characteristic_length_m)

    def validate(self) -> None:
        matrix = self._matrix()
        if matrix.shape != (3, 3) or not np.all(np.isfinite(matrix)):
            raise ValueError("rigid gamma solid rotation matrix must be finite 3x3")
        if not np.allclose(matrix.T.dot(matrix), np.eye(3), atol=1e-9, rtol=0.0):
            raise ValueError("rigid gamma solid rotation matrix must be orthonormal")
        if not np.all(np.isfinite(self._translation())) or not np.all(np.isfinite(self._origin())):
            raise ValueError("rigid gamma solid transform must be finite")
        self.base.validate()

    def ray_interval(self, origin: Sequence[float], direction: Sequence[float]) -> tuple[float, float]:
        return self.base.ray_interval(self._to_local_point(origin), self._to_local_direction(direction))

    def ray_intersection(self, origin: Sequence[float], direction: Sequence[float]):
        distance, face, normal = self.base.ray_intersection(self._to_local_point(origin), self._to_local_direction(direction))
        global_normal = self._matrix().dot(np.asarray(normal, dtype=float))
        return distance, face, global_normal


def rigid_transform_solid(
    solid: Any,
    rotation_matrix: Sequence[Sequence[float]],
    translation_m: Sequence[float] = (0.0, 0.0, 0.0),
    rotation_origin_m: Sequence[float] = (0.0, 0.0, 0.0),
) -> RigidTransformedSolid:
    """Wrap one v62 solid in a validated local-to-global rigid transform."""

    wrapped = RigidTransformedSolid(
        solid,
        tuple(tuple(float(value) for value in row) for row in rotation_matrix),
        tuple(float(value) for value in translation_m),
        tuple(float(value) for value in rotation_origin_m),
    )
    wrapped.validate()
    return wrapped


@dataclass
class GammaTrack:
    track_id: int
    kind: str
    parent_track_id: Optional[int]
    energy_keV: float
    position_m: tuple[float, float, float]
    direction: tuple[float, float, float]
    volume: str
    range_mm: float
    path_mm: float
    deposited_energy_keV: float
    escaped_energy_keV: float
    beta: float
    segments: list["ChargedTrackSegment"] = field(default_factory=list)


@dataclass
class ChargedTrackSegment:
    """One in-material part of a straight charged-particle path."""

    start_position_m: tuple[float, float, float]
    end_position_m: tuple[float, float, float]
    volume: str
    assembly: str
    component: str
    material: str
    path_mm: float
    kinetic_energy_start_keV: float
    kinetic_energy_end_keV: float
    deposited_energy_keV: float


@dataclass
class GammaRaySegment:
    """One straight gamma flight segment retained for event visualization."""

    track_id: int
    parent_track_id: Optional[int]
    energy_keV: float
    start_position_m: tuple[float, float, float]
    end_position_m: tuple[float, float, float]
    volume: str
    interaction: str


@dataclass
class EnergyDeposit:
    event_id: int
    volume: str
    assembly: str
    component: str
    material: str
    interaction: str
    position_m: tuple[float, float, float]
    energy_keV: float
    incident_gamma_keV: float
    track_id: int
    secondary_track_ids: tuple[int, ...] = ()
    charged_track_id: Optional[int] = None
    interaction_site: bool = True


@dataclass
class LightProduction:
    event_id: int
    volume: str
    assembly: str
    mechanism: str
    parent_track_id: int
    deposited_energy_keV: float
    path_mm: float
    expected_photons: float
    sampled_photons: int
    wavelength_range_nm: tuple[float, float]
    transport_status: str = "produced_local_only"


@dataclass
class GammaEvent:
    event_id: int
    source_energy_keV: float
    source_component: str
    deposits: list[EnergyDeposit] = field(default_factory=list)
    charged_tracks: list[GammaTrack] = field(default_factory=list)
    gamma_tracks: list[GammaRaySegment] = field(default_factory=list)
    light: list[LightProduction] = field(default_factory=list)
    interaction_counts: Counter = field(default_factory=Counter)
    escaped_gamma_energy_keV: float = 0.0
    # Gamma energy discarded below the transport threshold without creating a
    # charged track. Absorbed photoelectric energy is carried by that track.
    terminated_gamma_energy_keV: float = 0.0

    def energy_accounting(self) -> dict[str, float]:
        deposited = sum(item.energy_keV for item in self.deposits)
        escaped_charged = sum(item.escaped_energy_keV for item in self.charged_tracks)
        accounted = deposited + escaped_charged + self.escaped_gamma_energy_keV + self.terminated_gamma_energy_keV
        return {
            "source_energy_keV": self.source_energy_keV,
            "deposited_energy_keV": deposited,
            "escaped_charged_energy_keV": escaped_charged,
            "escaped_gamma_energy_keV": self.escaped_gamma_energy_keV,
            "terminated_gamma_energy_keV": self.terminated_gamma_energy_keV,
            "balance_residual_keV": self.source_energy_keV - accounted,
        }


@dataclass
class GammaTransportResult:
    source: GammaSource
    events: list[GammaEvent]
    model_status: str = MODEL_STATUS
    stopped: bool = False
    volumes: tuple[TransportVolume, ...] = ()
    requested_events: int = 0
    notes: list[str] = field(default_factory=lambda: [
        "Effective coarse mass-attenuation tables are used for fused silica and plastic/PMMA-like materials.",
        "Electrons and positrons use straight-line CSDA range/deposition; multiple scattering and detailed relaxation are not modeled.",
        "Charged stopping uses the legacy unvalidated range as one initial areal-range budget with uniform loss per g/cm^2 across configured material segments.",
        "Cherenkov/scintillation photons are counted locally from charged secondaries; no optical-photon or fiber transport is performed.",
        "Line components remain disabled unless explicit user relative weights are supplied.",
    ])

    def summary(self) -> dict[str, Any]:
        deposits = [item for event in self.events for item in event.deposits]
        light = [item for event in self.events for item in event.light]
        tracks = [item for event in self.events for item in event.charged_tracks]
        by_volume: dict[str, dict[str, Any]] = {}
        for volume in self.volumes:
            by_volume[volume.name] = {
                "assembly": volume.assembly, "component": volume.component,
                "material": volume.material, "interactions": {}, "deposited_energy_keV": 0.0,
                "deposited_energy_MeV": 0.0,
            }
        for item in deposits:
            row = by_volume.setdefault(item.volume, {"assembly": item.assembly, "component": item.component, "material": item.material, "interactions": {}, "deposited_energy_keV": 0.0})
            row["deposited_energy_keV"] += item.energy_keV
            if item.interaction_site:
                row["interactions"][item.interaction] = row["interactions"].get(item.interaction, 0) + 1
        for row in by_volume.values():
            row["deposited_energy_MeV"] = row["deposited_energy_keV"] / 1000.0
        by_assembly: dict[str, dict[str, Any]] = {}
        for name, row in by_volume.items():
            target = by_assembly.setdefault(row["assembly"], {"deposited_energy_keV": 0.0, "interactions": {}, "expected_cherenkov_photons": 0.0, "expected_scintillation_photons": 0.0})
            target["deposited_energy_keV"] += row["deposited_energy_keV"]
            for mechanism in row["interactions"]:
                target["interactions"][mechanism] = target["interactions"].get(mechanism, 0) + row["interactions"][mechanism]
        for item in light:
            by_assembly[item.assembly][f"expected_{item.mechanism}_photons"] += item.expected_photons
        interaction_totals: Counter = Counter()
        for event in self.events:
            interaction_totals.update(event.interaction_counts)
        total_deposited_keV = sum(item.energy_keV for item in deposits)
        escaped_charged_keV = sum(track.escaped_energy_keV for track in tracks)
        escaped_gamma_keV = sum(event.escaped_gamma_energy_keV for event in self.events)
        terminated_gamma_keV = sum(event.terminated_gamma_energy_keV for event in self.events)
        source_energy_keV = sum(event.source_energy_keV for event in self.events)
        accounted_energy_keV = total_deposited_keV + escaped_charged_keV + escaped_gamma_keV + terminated_gamma_keV
        return {
            "model_status": self.model_status,
            "stopped": self.stopped,
            "events": len(self.events),
            "event_accounting": {
                "requested": int(self.requested_events or len(self.events)),
                "completed": len(self.events),
                "failed": 0,
                "excluded": 0,
                "stopped": bool(self.stopped),
                "reasons": {
                    "failed": "none",
                    "excluded": "not applicable: gamma stage has no sensor acceptance/detection filter",
                },
            },
            "source_profile": self.source.source_profile_path,
            "source_direction_global": list(self.source.source_direction),
            "active_line_energies_keV": [line.energy_keV for line in self.source.lines],
            "active_line_weights": [line.relative_weight for line in self.source.lines],
            "total_deposited_energy_keV": total_deposited_keV,
            "escaped_charged_energy_keV": escaped_charged_keV,
            "energy_accounting": {
                "source_energy_keV": source_energy_keV,
                "deposited_energy_keV": total_deposited_keV,
                "escaped_charged_energy_keV": escaped_charged_keV,
                "escaped_gamma_energy_keV": escaped_gamma_keV,
                "terminated_gamma_energy_keV": terminated_gamma_keV,
                "balance_residual_keV": source_energy_keV - accounted_energy_keV,
            },
            "total_charged_secondary_tracks": len(tracks),
            "total_expected_cherenkov_photons": sum(item.expected_photons for item in light if item.mechanism == "cherenkov"),
            "total_expected_scintillation_photons": sum(item.expected_photons for item in light if item.mechanism == "scintillation"),
            "interaction_counts": dict(sorted(interaction_totals.items())),
            "by_assembly": by_assembly,
            "by_volume": by_volume,
            "escaped_gamma_energy_keV": escaped_gamma_keV,
            "notes": list(self.notes),
        }

    def to_dict(self) -> dict[str, Any]:
        payload = {"summary": self.summary(), "model_status": self.model_status, "notes": self.notes, "events": []}
        for event in self.events:
            payload["events"].append({
                "event_id": event.event_id,
                "source_energy_keV": event.source_energy_keV,
                "source_component": event.source_component,
                "deposits": [asdict(item) for item in event.deposits],
                "charged_tracks": [asdict(item) for item in event.charged_tracks],
                "gamma_tracks": [asdict(item) for item in event.gamma_tracks],
                "light": [asdict(item) for item in event.light],
                "interaction_counts": dict(event.interaction_counts),
                "escaped_gamma_energy_keV": event.escaped_gamma_energy_keV,
                "terminated_gamma_energy_keV": event.terminated_gamma_energy_keV,
                "energy_accounting": event.energy_accounting(),
            })
        return payload


def volumes_from_v62_geometry(quartz_geometry: Any, monitor: Any, quartz_name: str = "Cherenkov quartz") -> list[TransportVolume]:
    """Adapt the frozen v62 geometry and two-detector monitor to gamma volumes."""

    volumes: list[TransportVolume] = []
    if isinstance(quartz_geometry, sim.BarGeometry):
        # The optical simulator keeps its historical BarGeometry as a small
        # compatibility type; gamma transport needs the richer ray interval
        # contract, so adapt it without changing the optical model.
        quartz_geometry = sim.SolidGeometry(
            shape="box", material="fused_silica", length_x_m=quartz_geometry.length_x_m,
            width_y_m=quartz_geometry.width_y_m, thickness_z_m=quartz_geometry.thickness_z_m,
        )
    sections = tuple(getattr(quartz_geometry, "sections", (quartz_geometry,)))
    for index, solid in enumerate(sections):
        volumes.append(TransportVolume(quartz_name if len(sections) == 1 else f"{quartz_name} section {index + 1}", "cherenkov", "quartz", solid, "fused_silica", 0.0))
    if monitor is None:
        return volumes
    placed = monitor.placed_detectors(quartz_geometry)
    for detector_index, detector in enumerate(placed, 1):
        volumes.append(TransportVolume(f"Detector {detector_index} scintillator", f"detector_{detector_index}", "scintillator", detector.scintillator, detector.scintillator.material))
        for side, guide in (("left", detector.guide_left), ("right", detector.guide_right)):
            if guide is not None:
                volumes.append(TransportVolume(f"Detector {detector_index} {side} guide", f"detector_{detector_index}", f"guide_{side}", guide, guide.material, 0.0))
    return volumes


def _ray_interval(solid: Any, origin: np.ndarray, direction: np.ndarray) -> Optional[tuple[float, float]]:
    try:
        start, end = solid.ray_interval(origin, direction)
    except (AttributeError, RuntimeError, ValueError):
        return None
    if not (math.isfinite(float(start)) and math.isfinite(float(end))):
        return None
    end = float(end)
    start = float(start)
    if end <= max(start, 0.0):
        return None
    return max(start, 0.0), end


def _next_volume(volumes: Sequence[TransportVolume], position: np.ndarray, direction: np.ndarray) -> Optional[tuple[TransportVolume, float, float]]:
    candidates = []
    for volume in volumes:
        interval = _ray_interval(volume.solid, position, direction)
        if interval is not None:
            candidates.append((interval[0], volume, interval[1]))
    if not candidates:
        return None
    entry, volume, exit_distance = min(candidates, key=lambda item: item[0])
    return volume, entry, exit_distance


def _charged_track(
    event: GammaEvent, volume: TransportVolume, position: np.ndarray, direction: np.ndarray,
    kinetic_keV: float, kind: str, parent_track_id: Optional[int], track_id: int,
    incident_gamma_keV: float, rng: np.random.Generator, max_light_samples: int,
    volumes: Optional[Sequence[TransportVolume]] = None, interaction: str = "charged_deposition",
) -> tuple[GammaTrack, int]:
    """Transport a straight CSDA track through each configured solid it enters.

    The existing range expression supplies one initial areal-range budget.
    Energy loss is uniform per g/cm^2 over that budget, so a density change
    changes the distance traveled in each material without resetting the
    particle's range at a boundary. Gaps are treated as vacuum and do not
    consume the budget. Once no configured solid remains ahead, the residual
    kinetic energy is recorded as escaped.
    """

    kinetic_keV = max(0.0, float(kinetic_keV))
    direction = _unit(direction)
    start_position = np.asarray(position, dtype=float)
    range_g_cm2 = electron_range_g_cm2(kinetic_keV) if kinetic_keV > 0.0 else 0.0
    nominal_range_mm = range_g_cm2 / max(volume.density_g_cm3, 1e-12) * 10.0
    configured_volumes = tuple(volumes) if volumes is not None else (volume,)
    remaining_areal_range = range_g_cm2
    remaining_energy = kinetic_keV
    path_m = 0.0
    deposited = 0.0
    segments: list[ChargedTrackSegment] = []
    beta = _beta_from_kinetic_keV(kinetic_keV)
    cursor = start_position.copy()
    max_steps = max(32, len(configured_volumes) * 12)
    for _ in range(max_steps):
        if remaining_energy <= 1e-12 or remaining_areal_range <= 1e-15:
            remaining_energy = 0.0
            break
        query_origin = cursor
        hit = _next_volume(configured_volumes, query_origin, direction)
        if hit is None:
            break
        segment_volume, entry, exit_distance = hit
        interval_length = max(0.0, float(exit_distance) - float(entry))
        segment_start = query_origin + direction * max(0.0, float(entry))
        density = max(segment_volume.density_g_cm3, 0.0)
        if interval_length <= 1e-12 or density <= 0.0:
            cursor = query_origin + direction * (max(float(exit_distance), float(entry)) + 1e-9)
            continue
        areal_rate_g_cm2_per_m = density * 100.0
        stopping_distance = remaining_areal_range / areal_rate_g_cm2_per_m
        segment_path_m = min(interval_length, stopping_distance)
        if segment_path_m <= 0.0:
            remaining_energy = 0.0
            break
        areal_loss = min(remaining_areal_range, segment_path_m * areal_rate_g_cm2_per_m)
        reaches_stop = areal_loss >= remaining_areal_range - 1e-14
        energy_loss = remaining_energy if reaches_stop else min(remaining_energy, kinetic_keV * areal_loss / max(range_g_cm2, 1e-30))
        segment_energy_start = remaining_energy
        remaining_energy = max(0.0, remaining_energy - energy_loss)
        remaining_areal_range = max(0.0, remaining_areal_range - areal_loss)
        segment_end = segment_start + direction * segment_path_m
        segment = ChargedTrackSegment(
            tuple(float(value) for value in segment_start),
            tuple(float(value) for value in segment_end),
            segment_volume.name, segment_volume.assembly, segment_volume.component,
            segment_volume.material, segment_path_m * 1000.0, segment_energy_start,
            remaining_energy, energy_loss,
        )
        segments.append(segment)
        path_m += segment_path_m
        deposited += energy_loss
        if energy_loss > 0.0:
            event.deposits.append(EnergyDeposit(
                event.event_id, segment_volume.name, segment_volume.assembly,
                segment_volume.component, segment_volume.material, interaction,
                tuple(float(value) for value in segment_end), energy_loss,
                incident_gamma_keV, parent_track_id if parent_track_id is not None else -1,
                (), track_id, False,
            ))
            material = _material_key(segment_volume.material)
            segment_beta = _beta_from_kinetic_keV(segment_energy_start)
            if material == "fused_silica" and segment_beta > 1.0 / 1.47:
                wavelengths = np.linspace(300e-9, 650e-9, 256)
                spectrum = sim.frank_tamm_spectrum(wavelengths, segment_beta, charge_e=1.0)
                integral = float(np.trapezoid(spectrum, wavelengths) if hasattr(np, "trapezoid") else np.trapz(spectrum, wavelengths))
                expected = max(0.0, integral * segment_path_m)
                sampled = int(rng.poisson(expected)) if expected > 0.0 else 0
                event.light.append(LightProduction(event.event_id, segment_volume.name, segment_volume.assembly, "cherenkov", track_id, energy_loss, segment_path_m * 1000.0, expected, min(sampled, max_light_samples), (300.0, 650.0)))
            if material == "plastic":
                expected = energy_loss / 1000.0 * segment_volume.scintillation_yield_photons_per_MeV
                sampled = int(rng.poisson(expected)) if expected > 0.0 else 0
                event.light.append(LightProduction(event.event_id, segment_volume.name, segment_volume.assembly, "scintillation", track_id, energy_loss, segment_path_m * 1000.0, expected, min(sampled, max_light_samples), (350.0, 550.0)))
        if reaches_stop:
            remaining_energy = 0.0
            break
        cursor = segment_end + direction * 1e-9

    track = GammaTrack(
        track_id, kind, parent_track_id, kinetic_keV,
        tuple(float(value) for value in start_position), tuple(float(value) for value in direction),
        volume.name, nominal_range_mm, path_m * 1000.0, deposited,
        max(0.0, remaining_energy), beta, segments,
    )
    event.charged_tracks.append(track)
    next_track_id = track_id + 1
    return track, next_track_id


def _interaction_for_gamma(volume: TransportVolume, volumes: Sequence[TransportVolume], energy_keV: float, direction: np.ndarray, position: np.ndarray, gamma_track_id: int, event: GammaEvent, rng: np.random.Generator, pending: list[dict[str, Any]], next_track_id: int, max_light_samples: int) -> int:
    coefficients = mass_attenuation_coefficients(energy_keV, volume.material)
    density = volume.density_g_cm3
    linear = {key: value * density * 100.0 for key, value in coefficients.items()}  # cm^-1 -> m^-1
    total = sum(linear.values())
    if total <= 0.0:
        event.escaped_gamma_energy_keV += energy_keV
        return next_track_id
    draw = float(rng.random()) * total
    interaction = "photoelectric"
    cumulative = 0.0
    for candidate in ("compton", "photoelectric", "pair"):
        cumulative += linear[candidate]
        if draw <= cumulative:
            interaction = candidate
            break
    event.interaction_counts[interaction] += 1
    secondary_ids: list[int] = []
    if interaction == "photoelectric":
        electron_direction = _random_direction(rng)
        track, next_track_id = _charged_track(event, volume, position, electron_direction, energy_keV, "electron", gamma_track_id, next_track_id, energy_keV, rng, max_light_samples, volumes, interaction)
        secondary_ids.append(track.track_id)
        event.deposits.append(EnergyDeposit(event.event_id, volume.name, volume.assembly, volume.component, volume.material, interaction, tuple(position), 0.0, energy_keV, gamma_track_id, tuple(secondary_ids), None, True))
    elif interaction == "compton":
        cosine = _sample_compton_cosine(energy_keV, rng)
        scattered_energy = energy_keV / (1.0 + (energy_keV / ME_C2_KEV) * (1.0 - cosine))
        scattered_direction = _scatter_direction(direction, cosine, rng)
        electron_energy = max(0.0, energy_keV - scattered_energy)
        first, second = _basis(direction)
        momentum = energy_keV * direction - scattered_energy * scattered_direction
        electron_direction = _unit(momentum) if np.linalg.norm(momentum) > 1e-12 else _random_direction(rng)
        track, next_track_id = _charged_track(event, volume, position, electron_direction, electron_energy, "electron", gamma_track_id, next_track_id, energy_keV, rng, max_light_samples, volumes, interaction)
        secondary_ids.append(track.track_id)
        event.deposits.append(EnergyDeposit(event.event_id, volume.name, volume.assembly, volume.component, volume.material, interaction, tuple(position), 0.0, energy_keV, gamma_track_id, tuple(secondary_ids), None, True))
        if scattered_energy > 1.0:
            pending.append({"energy_keV": scattered_energy, "position": position + scattered_direction * 1e-10, "direction": scattered_direction})
        else:
            event.terminated_gamma_energy_keV += scattered_energy
    else:  # pair production; positron annihilation gammas are queued locally.
        kinetic_total = max(0.0, energy_keV - PAIR_THRESHOLD_KEV)
        electron_energy = kinetic_total * float(rng.random())
        positron_energy = kinetic_total - electron_energy
        for kind, kinetic in (("electron", electron_energy), ("positron", positron_energy)):
            track, next_track_id = _charged_track(event, volume, position, _random_direction(rng), kinetic, kind, gamma_track_id, next_track_id, energy_keV, rng, max_light_samples, volumes, interaction)
            secondary_ids.append(track.track_id)
        event.deposits.append(EnergyDeposit(event.event_id, volume.name, volume.assembly, volume.component, volume.material, interaction, tuple(position), 0.0, energy_keV, gamma_track_id, tuple(secondary_ids), None, True))
        annihilation_direction = _random_direction(rng)
        pending.append({"energy_keV": ME_C2_KEV, "position": position + annihilation_direction * 1e-10, "direction": annihilation_direction})
        pending.append({"energy_keV": ME_C2_KEV, "position": position - annihilation_direction * 1e-10, "direction": -annihilation_direction})
    return next_track_id


def simulate_event(
    source: GammaSource, volumes: Sequence[TransportVolume], event_id: int, rng: np.random.Generator,
    max_interactions: int = 1000, max_light_samples: int = 10000,
) -> GammaEvent:
    energy, component = source.sample_energy_keV(rng)
    event = GammaEvent(event_id, energy, component)
    pending = [{
        "energy_keV": energy,
        "position": np.array((0.0, 0.0, min(float(v.solid.bounds_min_m[2]) for v in volumes) - 1e-6)),
        "direction": _unit(source.source_direction),
        "gamma_track_id": 0,
    }] if volumes else []
    next_track_id = 1
    next_gamma_track_id = 1
    interactions = 0
    while pending and interactions < max_interactions:
        gamma = pending.pop(0)
        gamma_energy = float(gamma["energy_keV"])
        position = np.asarray(gamma["position"], dtype=float)
        direction = _unit(gamma["direction"])
        gamma_track_id = int(gamma.get("gamma_track_id", 0))
        parent_gamma_track_id = gamma.get("parent_gamma_track_id")
        while gamma_energy > 1.0:
            hit = _next_volume(volumes, position, direction)
            if hit is None:
                event.escaped_gamma_energy_keV += gamma_energy
                break
            volume, entry, exit_distance = hit
            material = _material_key(volume.material)
            coeff = mass_attenuation_coefficients(gamma_energy, material)
            total = sum(coeff.values()) * volume.density_g_cm3 * 100.0
            if total <= 0.0:
                event.escaped_gamma_energy_keV += gamma_energy
                break
            interaction_distance = float(rng.exponential(1.0 / total))
            if interaction_distance >= exit_distance - entry:
                segment_start = position.copy()
                segment_end = position + direction * max(exit_distance, entry + 1e-10)
                event.gamma_tracks.append(GammaRaySegment(
                    gamma_track_id, parent_gamma_track_id, gamma_energy,
                    tuple(float(value) for value in segment_start),
                    tuple(float(value) for value in segment_end), volume.name, "transmitted",
                ))
                position = segment_end
                # Continue from the volume exit, with a tiny forward nudge to
                # avoid selecting the same boundary repeatedly.
                position = position + direction * 1e-10
                continue
            segment_start = position.copy()
            position = position + direction * (entry + interaction_distance)
            interactions += 1
            pending_count = len(pending)
            next_track_id = _interaction_for_gamma(volume, volumes, gamma_energy, direction, position, gamma_track_id, event, rng, pending, next_track_id, max_light_samples)
            interaction = event.deposits[-1].interaction if event.deposits else "interaction"
            event.gamma_tracks.append(GammaRaySegment(
                gamma_track_id, parent_gamma_track_id, gamma_energy,
                tuple(float(value) for value in segment_start),
                tuple(float(value) for value in position), volume.name, interaction,
            ))
            for queued_gamma in pending[pending_count:]:
                queued_gamma["gamma_track_id"] = next_gamma_track_id
                queued_gamma["parent_gamma_track_id"] = gamma_track_id
                next_gamma_track_id += 1
            break
    if interactions >= max_interactions:
        event.escaped_gamma_energy_keV += sum(float(item["energy_keV"]) for item in pending)
    return event


def simulate(
    source: GammaSource, volumes: Sequence[TransportVolume], events: int = 1, seed: int = 12345,
    max_interactions: int = 1000, max_light_samples: int = 10000,
    progress_callback=None, status_callback=None, stop_event=None,
) -> GammaTransportResult:
    if not volumes:
        raise ValueError("gamma transport requires at least one geometry volume")
    if events < 1:
        raise ValueError("events must be at least one")
    for volume in volumes:
        volume.solid.validate()
    children = np.random.SeedSequence(seed).spawn(events)
    results = []
    stopped = False
    for index, child in enumerate(children):
        if stop_event is not None and stop_event.is_set():
            stopped = True
            break
        if status_callback is not None:
            status_callback(f"Gamma transport: event {index + 1}/{events}")
        results.append(simulate_event(source, volumes, index, np.random.default_rng(child), max_interactions, max_light_samples))
        if progress_callback is not None:
            progress_callback(index, events, 1.0)
        if stop_event is not None and stop_event.is_set():
            stopped = True
            break
    return GammaTransportResult(source, results, stopped=stopped, volumes=tuple(volumes), requested_events=events)


def _default_geometry() -> tuple[Any, Any]:
    quartz = sim.SolidGeometry(shape="box", material="fused_silica", length_x_m=0.02, width_y_m=0.02, thickness_z_m=0.01)
    def detector(z: float, material: str = "scintillator") -> sim.DetectorAssembly:
        active = sim.SolidGeometry(shape="box", material=material, center_z_m=z, length_x_m=0.02, width_y_m=0.02, thickness_z_m=0.01)
        return sim.DetectorAssembly(active, None, None, active_axis="Z")
    monitor = sim.TwoDetectorBeamMonitor(detector(0.020), detector(0.055), beam_axis="Z", detector_1_z_m=0.020, detector_2_z_m=0.055)
    return quartz, monitor


def _parse_line_weight(value: str) -> tuple[float, float]:
    try:
        energy, weight = value.split("=", 1)
        return float(energy), float(weight)
    except ValueError as exc:
        raise argparse.ArgumentTypeError("line weight must be ENERGY_KEV=RELATIVE_WEIGHT") from exc


def _write_plot(result: GammaTransportResult, path: Path) -> None:
    import matplotlib.pyplot as plt
    summary = result.summary()
    labels = list(summary["interaction_counts"])
    values = [summary["interaction_counts"][label] for label in labels]
    assemblies = list(summary["by_assembly"])
    deposited = [summary["by_assembly"][label]["deposited_energy_keV"] / 1000.0 for label in assemblies]
    fig, axes = plt.subplots(1, 2, figsize=(10, 4.2), constrained_layout=True)
    axes[0].bar(labels or ["none"], values or [0], color="#4472c4")
    axes[0].set_title("Gamma interaction counts")
    axes[0].set_ylabel("interactions")
    axes[0].tick_params(axis="x", rotation=25)
    axes[1].bar(assemblies or ["none"], deposited or [0], color="#70ad47")
    axes[1].set_title("Deposited energy by assembly")
    axes[1].set_ylabel("MeV")
    axes[1].tick_params(axis="x", rotation=25)
    fig.savefig(path, dpi=160)
    plt.close(fig)


def main(argv: Optional[Sequence[str]] = None) -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--source-profile", type=Path, required=False, help="approved gamma_source_profile_collimated_v1.json")
    parser.add_argument("--line-weight", action="append", default=[], type=_parse_line_weight, metavar="KEV=WEIGHT", help="enable a source line with an explicit relative weight; repeatable")
    parser.add_argument("--events", type=int, default=100)
    parser.add_argument("--seed", type=int, default=12345)
    parser.add_argument("--json-output", type=Path, default=None)
    parser.add_argument("--plot", type=Path, default=None)
    parser.add_argument("--self-test", action="store_true")
    args = parser.parse_args(argv)
    if args.self_test:
        assert mass_attenuation_coefficients(511.0, "fused_silica")["compton"] > 0.0
        assert mass_attenuation_coefficients(511.0, "fused_silica")["pair"] == 0.0
        assert mass_attenuation_coefficients(1461.0, "fused_silica")["pair"] > 0.0
        print("gamma_transport self-test passed")
        return 0
    if args.source_profile is None:
        parser.error("--source-profile is required unless --self-test is used")
    weights = dict(args.line_weight)
    source = GammaSource.from_profile(args.source_profile, weights)
    quartz, monitor = _default_geometry()
    result = simulate(source, volumes_from_v62_geometry(quartz, monitor), events=args.events, seed=args.seed)
    payload = result.to_dict()
    if args.json_output:
        args.json_output.parent.mkdir(parents=True, exist_ok=True)
        args.json_output.write_text(json.dumps(payload, indent=2, sort_keys=True) + "\n", encoding="utf-8")
    if args.plot:
        args.plot.parent.mkdir(parents=True, exist_ok=True)
        _write_plot(result, args.plot)
    print(json.dumps(result.summary(), indent=2, sort_keys=True))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
