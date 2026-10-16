#!/usr/bin/env python3
"""Tkinter workflow for the quartz Cherenkov and beam-monitor simulator.

Run with ``python3 quartz_cherenkov_gui.py``. The GUI uses the corrected
100 mm x 20 mm x 10 mm geometry and delegates CompoundGeometry transport,
including oriented boundary hand-off and local tapered-surface normals, to
quartz_cherenkov_sim.py. It writes one timestamped directory containing a
ROOT file when optional ROOT dependencies are available, PNG plots, and a JSON
summary. The main notebook keeps local scintillator geometry, scintillator
detector response, and global beam-monitor placement as independent root-level
tabs. The latter is the assembly view of the Cherenkov radiator plus
downstream scintillator/light-guide instrument; it is not a synonym for
scintillator.
"""

from __future__ import annotations

import json
import hashlib
import math
import os
import queue
import subprocess
import sys
import threading
import time
from copy import copy, deepcopy
from dataclasses import asdict, dataclass, is_dataclass, replace
from datetime import datetime
from pathlib import Path
from types import MappingProxyType
from typing import Callable, Dict, Mapping, Optional, Sequence, Tuple

import matplotlib
import numpy as np

matplotlib.use("Agg")

try:
    import tkinter as tk
    from tkinter import filedialog, messagebox, ttk, font as tkfont
    TK_AVAILABLE = True
    TK_IMPORT_ERROR: Optional[Exception] = None
except ImportError as exc:
    TK_AVAILABLE = False
    TK_IMPORT_ERROR = exc

    class _TkUnavailable:
        Tk = object

    tk = _TkUnavailable()
    filedialog = None
    messagebox = None
    ttk = None
    tkfont = None


def _has_instance_attr(instance, name: str) -> bool:
    """Check GUI instance state without invoking Tkinter's fallback lookup.

    ``tkinter.Misc.__getattr__`` forwards unknown names to ``self.tk``. A
    headless probe made with ``Tk.__new__`` has no interpreter, so ``hasattr``
    on such an object can recurse while trying to find ``tk``. These checks
    concern instance state only; inspect the instance dictionary directly.
    """

    try:
        namespace = object.__getattribute__(instance, "__dict__")
    except (AttributeError, TypeError):
        return False
    return name in namespace


_MISSING_INSTANCE_ATTR = object()


def _get_instance_attr(instance, name: str, default=_MISSING_INSTANCE_ATTR):
    """Read optional GUI state without falling through Tkinter.__getattr__."""

    try:
        namespace = object.__getattribute__(instance, "__dict__")
    except (AttributeError, TypeError):
        if default is _MISSING_INSTANCE_ATTR:
            raise AttributeError(name) from None
        return default
    if name in namespace:
        return namespace[name]
    if default is not _MISSING_INSTANCE_ATTR:
        return default
    raise AttributeError(name)

TK_ARROW_LAST = getattr(tk, "LAST", "last")

import quartz_cherenkov_sim as simulator
import gamma_transport
from quartz_geometry import (
    AssemblyGeometryModel, AssemblySection, RigidAssemblyTransform,
    MYLAR_DISPLAY_COLOR, MYLAR_DISPLAY_EDGE, MYLAR_DISPLAY_HATCH,
    TERMINAL_SURFACE_CHOICES, TERMINAL_SURFACE_LABELS,
    assembly_geometry_from_dict, assembly_geometry_to_dict,
    complete_detector_assembly_geometry_from_dict,
    complete_detector_assembly_geometry_to_dict,
    build_complete_assembly_scene, normalize_terminal_surface,
    resolve_sensor_readout_dimensions_mm,
    canonical_section_faces_mm, derive_mylar_overlay_surfaces,
    sensor_solid_geometry_mm,
    terminal_surface_label,
    transform_assembly_components, transform_assembly_terminal_interfaces,
)


APP_VERSION = "0.46.0"
APP_BUILD_ID = "build-20260929-v78-gamma-charged-boundaries"
MAX_DISPLAYED_EVENT_ROWS = 2000

# The embedded global preview uses one unambiguous particle-track color.  It
# is intentionally different from the green local-coordinate axis used by
# older local previews and from the global XYZ triad below.
GLOBAL_PARTICLE_PATH_COLOR = "#0b8f55"

# The embedded preview and detached Global Assembly viewer deliberately share
# one camera contract.  A triad is allowed to be translated to a readable
# viewport corner, but its directions must always be the projected global
# basis under this exact yaw/pitch transform.
GLOBAL_PREVIEW_DEFAULT_YAW = 0.65
GLOBAL_PREVIEW_DEFAULT_PITCH = 0.55
GLOBAL_PREVIEW_AXIS_LENGTH_PX = 30.0


def global_assembly_camera_basis(yaw: float, pitch: float) -> Dict[str, tuple]:
    """Return the camera-space projections of canonical global unit axes.

    ``screen_y`` is positive upward, matching the renderer's y inversion.
    The returned ``depth`` is only for painter ordering; it is intentionally
    not a physical transform and never changes global scene coordinates.
    """

    cy, sy = math.cos(float(yaw)), math.sin(float(yaw))
    cp, sp = math.cos(float(pitch)), math.sin(float(pitch))

    def transform(vector: Sequence[float]) -> tuple[float, float, float]:
        x, y, z = (float(value) for value in vector)
        rotated_x = cy * x - sy * y
        rotated_y = sy * x + cy * y
        return rotated_x, cp * rotated_y - sp * z, sp * rotated_y + cp * z

    return {
        axis: transform(direction)
        for axis, direction in (
            ("X", (1.0, 0.0, 0.0)),
            ("Y", (0.0, 1.0, 0.0)),
            ("Z", (0.0, 0.0, 1.0)),
        )
    }


def global_assembly_camera_project(
    point: Sequence[float], center: Sequence[float], scale: float,
    width: float, height: float, yaw: float, pitch: float,
) -> tuple[float, float, float]:
    """Project one global/display point with the shared orbit camera."""

    basis = global_assembly_camera_basis(yaw, pitch)
    relative = tuple(float(point[index]) - float(center[index]) for index in range(3))
    # The basis columns are the transformed world axes; using the same dot
    # products for points and vectors makes path/triad agreement exact.
    screen_x = sum(relative[index] * basis[axis][0] for index, axis in enumerate(("X", "Y", "Z")))
    screen_y = sum(relative[index] * basis[axis][1] for index, axis in enumerate(("X", "Y", "Z")))
    depth = sum(relative[index] * basis[axis][2] for index, axis in enumerate(("X", "Y", "Z")))
    return float(width) / 2.0 + screen_x * float(scale), float(height) / 2.0 - screen_y * float(scale), depth


def global_assembly_orientation_triad(yaw: float, pitch: float) -> Dict[str, tuple]:
    """Return 2-D camera-projected global XYZ vectors for an orientation gizmo."""

    return {
        axis: (float(vector[0]), float(vector[1]))
        for axis, vector in global_assembly_camera_basis(yaw, pitch).items()
    }


# A GUI run has one worker thread, but keeping the writer synchronized also
# makes direct/headless callers and future status producers safe.  Each line
# is opened, flushed, and closed immediately: an antivirus scanner or a
# Windows file previewer can therefore never hold the run log open for the
# whole simulation.
_RUN_LOG_LOCK = threading.RLock()


def _ttk_background(widget, style_name: str = "TLabelframe", fallback: str = "#f0f0f0") -> str:
    """Return a valid Tk background without querying an unsupported ttk option."""

    safe_fallback = str(fallback or "white").strip() or "white"
    if ttk is None:
        return safe_fallback
    try:
        color = ttk.Style(widget).lookup(style_name, "background")
    except Exception:
        color = ""
    color = str(color or "").strip()
    return color or safe_fallback


def _append_persistent_run_log(path: Path, message: object) -> None:
    """Append one flushed line to a run-local log without recursive logging."""

    try:
        text = str(message).rstrip()
        with _RUN_LOG_LOCK:
            path.parent.mkdir(parents=True, exist_ok=True)
            with path.open("a", encoding="utf-8", newline="") as handle:
                handle.write(f"{datetime.now().astimezone().isoformat(timespec='milliseconds')} {text}\n")
                handle.flush()
                # This is intentionally best-effort.  fsync is useful on
                # Windows when a process is interrupted, but some networked
                # or removable filesystems reject it even after a successful
                # flush; never turn logging failure into a simulation failure.
                try:
                    os.fsync(handle.fileno())
                except OSError:
                    pass
    except Exception:
        # Logging must not recurse into the worker's error route.
        pass


def _json_safe(value):
    """Convert ordinary/numpy/dataclass values to strict JSON-safe values."""

    if is_dataclass(value):
        return _json_safe(asdict(value))
    if isinstance(value, Mapping):
        return {str(key): _json_safe(item) for key, item in value.items()}
    if isinstance(value, (tuple, list)):
        return [_json_safe(item) for item in value]
    if isinstance(value, np.ndarray):
        return _json_safe(value.tolist())
    if isinstance(value, np.generic):
        return _json_safe(value.item())
    if isinstance(value, Path):
        return str(value)
    if isinstance(value, float) and not math.isfinite(value):
        return None
    return value


def _safe_file_fingerprint(path_value) -> Dict[str, str]:
    """Return display-safe path, mtime, and SHA-256 values."""

    unavailable = "unavailable"
    try:
        path = Path(path_value).expanduser().resolve()
    except (OSError, TypeError, ValueError):
        return {"path": unavailable, "mtime": unavailable, "sha256": unavailable}
    result = {"path": str(path), "mtime": unavailable, "sha256": unavailable}
    try:
        result["mtime"] = datetime.fromtimestamp(
            path.stat().st_mtime
        ).astimezone().isoformat(timespec="seconds")
    except (OSError, OverflowError, ValueError):
        pass
    try:
        digest = hashlib.sha256()
        with path.open("rb") as stream:
            for chunk in iter(lambda: stream.read(1024 * 1024), b""):
                digest.update(chunk)
        result["sha256"] = digest.hexdigest()
    except (OSError, TypeError, ValueError):
        pass
    return result


def source_diagnostic_text() -> str:
    """Build the visible provenance diagnostic for the current process."""

    gui = _safe_file_fingerprint(__file__)
    simulator_path = getattr(simulator, "__file__", None)
    sim = _safe_file_fingerprint(simulator_path)
    gui_hash = gui["sha256"][:16] if gui["sha256"] != "unavailable" else "unavailable"
    sim_hash = sim["sha256"][:16] if sim["sha256"] != "unavailable" else "unavailable"
    return (
        f"Loaded GUI file: {gui['path']}\n"
        f"GUI version/build: {APP_VERSION} / {APP_BUILD_ID}\n"
        f"GUI source mtime: {gui['mtime']}\n"
        f"GUI SHA-256: {gui_hash}\n"
        f"Loaded simulator file: {sim['path']}\n"
        f"Simulator SHA-256: {sim_hash}\n"
        f"Simulator mtime: {sim['mtime']}"
    )


PARTICLE_MASSES_MEV_C2 = {
    "electron": 0.51099895,
    "muon": 105.6583755,
    "charged pion": 139.57039,
}

SAMPLING_MODE_LABELS = {
    "fixed": "Fixed Monte Carlo test sample",
    "poisson": "Physical Poisson from Frank-Tamm",
}
SAMPLING_MODE_BY_LABEL = {label: mode for mode, label in SAMPLING_MODE_LABELS.items()}
SCINTILLATION_SAMPLING_MODE_LABELS = {
    "fixed": "Fixed Monte Carlo test sample",
    "physical": "Physical yield from deposited energy",
}
SCINTILLATION_SAMPLING_MODE_BY_LABEL = {label: mode for mode, label in SCINTILLATION_SAMPLING_MODE_LABELS.items()}
SCINTILLATOR_DECAY_TIME_DEFAULT_NS = 2.0
DETECTOR_RUN_MODE_LABELS = {
    "cherenkov": "Cherenkov only",
    "scintillator": "Scintillator only",
    "combined": "All",
}
DETECTOR_RUN_MODE_BY_LABEL = {label: mode for mode, label in DETECTOR_RUN_MODE_LABELS.items()}
GAMMA_RUN_MODE = "gamma"
GAMMA_RUN_MODE_LABEL = "Gamma transport"
# Keep the v62 detector-mode mapping stable for optical callers and old
# regression checks; Gamma is an explicit additional Run choice at the GUI
# layer and never enters simulator.normalize_detector_run_mode().
RUN_MODE_BY_LABEL = {**DETECTOR_RUN_MODE_BY_LABEL, GAMMA_RUN_MODE_LABEL: GAMMA_RUN_MODE}
RUN_MODE_COMBO_VALUES = tuple(DETECTOR_RUN_MODE_LABELS.values()) + (GAMMA_RUN_MODE_LABEL,)


def detector_run_mode(value: object) -> str:
    """Return the one normalized detector-participation mode used by runs."""

    text = str(value or "").strip()
    if RUN_MODE_BY_LABEL.get(text, text).lower() == GAMMA_RUN_MODE:
        return GAMMA_RUN_MODE
    return simulator.normalize_detector_run_mode(DETECTOR_RUN_MODE_BY_LABEL.get(text, text or "cherenkov"))


def detector_run_mode_label(value: object) -> str:
    mode = detector_run_mode(value)
    return GAMMA_RUN_MODE_LABEL if mode == GAMMA_RUN_MODE else DETECTOR_RUN_MODE_LABELS[mode]


def detector_run_mode_from_values(values: Mapping[str, object]) -> str:
    """Migrate an old value store once, then use only the canonical mode.

    New stores contain the visible Run choice.  The only compatibility use of
    ``scintillator_enabled`` is for v42/older stores that have no explicit
    choice (or still say ``auto``); that value is immediately normalized and
    is never retained as active run state.
    """

    raw = _value_text(values, "detector_run_mode", "")
    if raw and raw.strip().lower() not in {"auto", "automatic"}:
        return detector_run_mode(RUN_MODE_BY_LABEL.get(raw, raw))
    if raw.strip().lower() in {"auto", "automatic"} and "scintillator_enabled" in values:
        return "combined" if _bool_text(values, "scintillator_enabled") else "cherenkov"
    return "cherenkov"
ROOT_NOTEBOOK_TAB_LABELS = (
    "Run",
    "Results",
    "Gamma Results",
    "Cherenkov",
    "Scintillator",
    "Global Assembly 3D",
)

ASSEMBLY_SUBTAB_LABELS = {
    "Cherenkov": ("Cherenkov Geometry Editor", "Cherenkov Sensor"),
    "Scintillator": ("Scintillator Geometry Editor", "Scintillator Sensor"),
}


def responsive_three_view_panel_bounds(width: float, height: float, *, margin: float = 2.0,
                                       gutter: float = 4.0) -> Tuple[Tuple[float, float, float, float], ...]:
    """Return three equal-width view rectangles for the live canvas size.

    The returned rectangles intentionally leave only the requested outer
    margin and the two requested inter-panel gutters.  In particular, the
    final rectangle ends at ``width - margin`` rather than at a historical
    fixed scene width.  This small pure helper is also used by headless
    layout checks so the resize contract does not depend on Tk.
    """

    canvas_width = max(float(width), 0.0)
    canvas_height = max(float(height), 0.0)
    outer_margin = max(float(margin), 0.0)
    outer_margin = min(outer_margin, canvas_width / 2.0)
    remaining_width = max(canvas_width - 2.0 * outer_margin, 0.0)
    inter_panel_gutter = min(max(float(gutter), 0.0), remaining_width / 2.0)
    panel_width = max((remaining_width - 2.0 * inter_panel_gutter) / 3.0, 0.0)
    top = min(outer_margin, canvas_height)
    bottom = max(top, canvas_height - outer_margin)
    return tuple(
        (
            outer_margin + index * (panel_width + inter_panel_gutter),
            top,
            outer_margin + index * (panel_width + inter_panel_gutter) + panel_width,
            bottom,
        )
        for index in range(3)
    )


def run_status_panel_layout_contract(panel_width: float = 430.0, *, text: str = "") -> Dict[str, object]:
    """Describe the fixed Run-status cells used by the Tk layout.

    This pure contract keeps status strings from becoming geometry requests:
    both labels are assigned the same finite width and wrap within the status
    panel.  It is intentionally testable without constructing a Tk display.
    """

    width = max(float(panel_width), 240.0)
    label_width = max(width - 14.0, 160.0)
    # Tk's wraplength is pixels; this estimate is only a deterministic
    # diagnostic for long-string regression tests, not a renderer substitute.
    chars_per_line = max(16, int(label_width / 7.0))
    value = str(text)
    all_lines = tuple(value[index:index + chars_per_line] for index in range(0, len(value), chars_per_line)) or ("",)
    lines = all_lines[:2]
    return {
        "panel_width_px": width,
        "status_cell_width_px": label_width,
        "progress_cell_width_px": label_width,
        "wrapped_lines": lines,
        "status_row_height_px": 34,
        "progress_row": 2,
        "row_heights_px": (32, 34, 22, 34, 30, 28),
        "status_text_clipped": len(all_lines) > len(lines),
        "dynamic_wrap_height": False,
        "fixed_width_cells": True,
        "overlap_with_simulation_parameters": False,
    }


# The global and event viewers deliberately share this small, renderer-neutral
# label policy.  Keeping it independent of Tk makes the placement contract
# testable at several canvas sizes and camera projections without a display.
_CALLOUT_KIND_NAMES = {
    "geometry": {"cherenkov": "quartz", "scintillator": "scintillator"},
    "guide": "guide",
    "sensor": "sensor",
    "readout_interface": "sensor",
    "mirror": "mirror",
    "mylar_overlay": "Mylar",
    "particle_path": "track",
}


def detector_callout_label(component: Mapping[str, object]) -> str:
    """Return a short detector-prefixed label for a 3D scene item.

    Full section names are diagnostic detail, not canvas annotations.  The
    prefix is intentionally stable (C/S) so combined scenes cannot make a
    viewer guess which detector owns a track, sensor, or mirror.
    """

    assembly = str(component.get("assembly_type", component.get("detector_type", ""))).lower()
    prefix = "C" if assembly.startswith("cher") or assembly in {"1", "detector_1"} else "S"
    kind = str(component.get("kind", "geometry")).lower()
    name = _CALLOUT_KIND_NAMES.get(kind)
    if isinstance(name, Mapping):
        name = name.get("cherenkov" if prefix == "C" else "scintillator", "body")
    if not name:
        name = "body"
    return f"{prefix} {name}"


def scene_legend_rect(width: float, height: float) -> Tuple[float, float, float, float]:
    """Return the compact legend card below the reserved title band.

    The old v52 card started at y=8 while the title also occupied the top
    band, so the two text layers could intersect.  Keep the card widths
    responsive so the legend and measurement cards also remain side by side
    on a narrow headless or popup canvas.
    """

    canvas_width = max(float(width), 1.0)
    legend_width = min(292.0, max(178.0, (canvas_width - 24.0) * 0.46))
    legend_height = 120.0
    return 8.0, 60.0, 8.0 + legend_width, 60.0 + legend_height


def scene_measurement_rect(width: float, height: float) -> Tuple[float, float, float, float]:
    """Reserve the measurement card beside the legend below the title band."""

    canvas_width = max(float(width), 1.0)
    legend_width = min(292.0, max(178.0, (canvas_width - 24.0) * 0.46))
    card_width = max(178.0, canvas_width - legend_width - 24.0)
    card_height = 120.0
    return canvas_width - card_width - 8.0, 60.0, canvas_width - 8.0, 60.0 + card_height


def scene_reserved_rectangles(width: float, height: float, *, sidebar_width: float = 0.0) -> Tuple[Tuple[float, float, float, float], ...]:
    """Return stable title, card, sidebar, and bottom reservations.

    These are screen-space rectangles, not point placeholders: the same
    reservations are passed to every collision-aware callout allocator.
    """

    canvas_width, canvas_height = max(float(width), 1.0), max(float(height), 1.0)
    sidebar = max(float(sidebar_width), 0.0)
    return (
        scene_legend_rect(canvas_width, canvas_height),
        scene_measurement_rect(canvas_width - sidebar, canvas_height),
        (0.0, 0.0, canvas_width, 52.0),
        (0.0, max(0.0, canvas_height - 58.0), canvas_width, canvas_height),
        (max(0.0, canvas_width - sidebar), 0.0, canvas_width, canvas_height) if sidebar else (0.0, 0.0, 0.0, 0.0),
    )


def _event_viewer_default_measure(text: str, font_size: float) -> float:
    """Conservative headless text-width estimate for the event header."""

    return max(1.0, 0.58 * float(font_size) * len(str(text)))


def _event_viewer_ellipsis(text: str, max_width: float, measure: Callable[[str], float]) -> Tuple[str, bool]:
    """Fit one line to a measured width, retaining an explicit ellipsis."""

    value = str(text).replace("\n", " ").strip()
    if measure(value) <= max_width:
        return value, False
    suffix = "…"
    if measure(suffix) > max_width:
        return "", True
    end = value
    while end and measure(end + suffix) > max_width:
        end = end[:-1]
    return (end.rstrip() + suffix) if end.rstrip() else suffix, True


def _event_viewer_wrap_status(
    text: str, max_width: float, measure: Callable[[str], float], max_lines: int = 2,
) -> Tuple[Tuple[str, ...], bool]:
    """Word-wrap a status message into a bounded number of measured lines."""

    value = " ".join(str(text).split())
    if not value:
        return ("",), False
    words = value.split(" ")
    lines = []
    current = ""
    for word in words:
        candidate = word if not current else f"{current} {word}"
        if measure(candidate) <= max_width:
            current = candidate
            continue
        if current:
            lines.append(current)
            current = ""
        # A single very long token still gets wrapped instead of extending
        # the status band beyond the fixed layout contract.
        while measure(word) > max_width and len(word) > 1:
            cut = len(word)
            while cut > 1 and measure(word[:cut]) > max_width:
                cut -= 1
            lines.append(word[:cut])
            word = word[cut:]
        current = word
    if current:
        lines.append(current)
    if not lines:
        lines = [""]
    clipped = len(lines) > max_lines
    if clipped:
        kept = list(lines[:max_lines])
        tail = " ".join(lines[max_lines - 1:])
        kept[-1], _ = _event_viewer_ellipsis(tail, max_width, measure)
        lines = kept
    return tuple(lines[:max_lines]), clipped


def interactive_event_view_layout(
    width: float, height: float, title: str, status: str, *,
    title_measure: Optional[Callable[[str], float]] = None,
    status_measure: Optional[Callable[[str], float]] = None,
    title_line_height: float = 18.0,
    status_line_height: float = 13.0,
    status_max_lines: int = 2,
) -> Dict[str, object]:
    """Return non-overlapping event-viewer header, legend, and scene bands.

    The callbacks are Tk font metrics in the live renderer and deterministic
    estimates in headless checks.  Status text is deliberately capped at two
    lines; ``full_status`` remains available to the scrollable details pane.
    Geometry is fitted to ``scene_rect`` below the legend, so a resize cannot
    move solids, tracks, or photon paths under the header cards.
    """

    canvas_width = max(float(width), 1.0)
    canvas_height = max(float(height), 1.0)
    title_measure = title_measure or (lambda value: _event_viewer_default_measure(value, 12.0))
    status_measure = status_measure or (lambda value: _event_viewer_default_measure(value, 9.0))
    text_width = max(canvas_width - 16.0, 24.0)
    title_text, title_clipped = _event_viewer_ellipsis(str(title), text_width, title_measure)
    title_height = max(16.0, float(title_line_height))
    status_height = max(11.0, float(status_line_height))
    title_rect = (
        max(0.0, (canvas_width - min(title_measure(title_text), text_width)) / 2.0),
        5.0,
        min(canvas_width, (canvas_width + min(title_measure(title_text), text_width)) / 2.0),
        5.0 + title_height,
    )
    status_top = title_rect[3] + 4.0
    status_lines, status_clipped = _event_viewer_wrap_status(
        status, text_width, status_measure, max(1, int(status_max_lines)),
    )
    status_rect = (8.0, status_top, max(8.0, canvas_width - 8.0), status_top + status_height * len(status_lines))
    header_bottom = status_rect[3] + 6.0

    legend_width = min(292.0, max(178.0, (canvas_width - 24.0) * 0.46))
    legend_height = min(120.0, max(64.0, canvas_height * 0.34))
    legend_rect = (8.0, header_bottom + 8.0, min(canvas_width - 8.0, 8.0 + legend_width), header_bottom + 8.0 + legend_height)
    bottom_reserved = min(28.0, max(0.0, canvas_height * 0.12))
    scene_rect = (
        8.0,
        legend_rect[3] + 8.0,
        max(8.0, canvas_width - 8.0),
        max(legend_rect[3] + 8.0, canvas_height - bottom_reserved),
    )
    header_rect = (0.0, 0.0, canvas_width, header_bottom)
    return {
        "canvas_size": (canvas_width, canvas_height),
        "title": title_text,
        "full_title": str(title),
        "title_clipped": title_clipped,
        "title_rect": title_rect,
        "status": "\n".join(status_lines),
        "status_lines": status_lines,
        "full_status": str(status),
        "status_clipped": status_clipped,
        "status_rect": status_rect,
        "header_rect": header_rect,
        "legend_rect": legend_rect,
        "scene_rect": scene_rect,
        "bottom_reserved": bottom_reserved,
        "occupied_rectangles": (header_rect, legend_rect, (0.0, canvas_height - bottom_reserved, canvas_width, canvas_height)),
    }


def _rectangles_overlap(first, second, gap: float = 2.0) -> bool:
    """Return true when two screen-space rectangles overlap with a gap."""

    left_a, top_a, right_a, bottom_a = (float(value) for value in first)
    left_b, top_b, right_b, bottom_b = (float(value) for value in second)
    return not (
        right_a + gap <= left_b or right_b + gap <= left_a or
        bottom_a + gap <= top_b or bottom_b + gap <= top_a
    )


def collision_aware_callout_layout(items: Sequence[Mapping[str, object]], width: float, height: float,
                                   *, occupied_rectangles: Sequence[Sequence[float]] = (),
                                   show_labels: bool = True, text_measure: Optional[Callable] = None) -> Dict[str, object]:
    """Place short callouts without screen-space rectangle collisions.

    ``items`` contain ``target`` (the projected geometry point), ``label``
    and optional ``detail``/``component_id`` fields.  Candidate offsets and
    iteration order are fixed, so the result is deterministic for a camera
    state.  Every accepted callout has a leader back to its target; items that
    cannot fit are returned in ``omitted`` for the fixed legend/sidebar.
    """

    canvas_width, canvas_height = max(float(width), 1.0), max(float(height), 1.0)
    occupied = [tuple(float(value) for value in rect) for rect in occupied_rectangles]
    accepted = []
    omitted = []
    # Text anchors are top-left boxes.  The order favours the readable outer
    # side of a component, then walks a deterministic clockwise ring.
    offsets = (
        (14.0, -22.0), (14.0, 8.0), (-14.0, -22.0), (-14.0, 8.0),
        (24.0, -42.0), (24.0, 28.0), (-24.0, -42.0), (-24.0, 28.0),
        (42.0, -22.0), (42.0, 8.0), (-42.0, -22.0), (-42.0, 8.0),
        (0.0, -42.0), (0.0, 28.0), (58.0, -52.0), (58.0, 38.0),
        (-58.0, -52.0), (-58.0, 38.0),
    )

    def measure(text: str) -> Tuple[float, float]:
        if text_measure is not None:
            measured = text_measure(text)
            if len(measured) == 2:
                return max(float(measured[0]), 12.0), max(float(measured[1]), 10.0)
        return max(24.0, 7.0 * len(str(text)) + 8.0), 16.0

    # Mylar and repeated section labels are useful geometry detail, but
    # dozens of identical canvas labels make the two required sensors
    # unreadable.  Collapse only these low-priority families into one
    # deterministic summary item; every source component ID remains attached
    # to that summary and omitted details remain available in the legend.
    prepared = []
    grouped = {}
    for source in items:
        kind = str(source.get("kind", "geometry"))
        assembly = str(source.get("assembly_type", ""))
        label = str(source.get("label") or detector_callout_label(source))
        if kind in {"mylar_overlay", "geometry", "guide"}:
            grouped.setdefault((assembly, kind, label), []).append(dict(source))
        else:
            prepared.append(dict(source))
    for (assembly, kind, label), group in grouped.items():
        if len(group) == 1:
            prepared.append(group[0])
            continue
        representative = dict(group[0])
        count = len(group)
        representative["label"] = f"{label} sections ({count})"
        representative["detail"] = f"{label}: {count} sections aggregated; individual IDs retained"
        representative["component_ids"] = tuple(str(item.get("component_id", item.get("name", ""))) for item in group)
        representative["aggregated_count"] = count
        prepared.append(representative)

    # Keep sensor/readout callouts first so a finite canvas always spends
    # space on the two detector identities before low-priority geometry.
    prepared.sort(key=lambda item: (0 if str(item.get("kind", "")).lower() in {"sensor", "readout_interface"} else 1, str(item.get("component_id", ""))))
    for sequence, source in enumerate(prepared):
        item = dict(source)
        target = item.get("target", item.get("anchor", (0.0, 0.0)))
        target_x, target_y = float(target[0]), float(target[1])
        label = str(item.get("label") or detector_callout_label(item))
        detail = str(item.get("detail", item.get("source_label", label)))
        if not show_labels:
            omitted.append({**item, "label": label, "detail": detail, "reason": "labels hidden"})
            continue
        text_width, text_height = measure(label)
        placed = None
        for offset_x, offset_y in offsets:
            left = target_x + offset_x if offset_x >= 0 else target_x + offset_x - text_width
            top = target_y + offset_y
            right, bottom = left + text_width, top + text_height
            if left < 4.0 or top < 4.0 or right > canvas_width - 4.0 or bottom > canvas_height - 4.0:
                continue
            rect = (left, top, right, bottom)
            if any(_rectangles_overlap(rect, other) for other in occupied):
                continue
            placed = {
                **item,
                "sequence": sequence,
                "label": label,
                "detail": detail,
                "target": (target_x, target_y),
                "position": (left, top),
                "anchor": "nw",
                "rect": rect,
                "leader_start": (left if offset_x >= 0 else right, top + text_height / 2.0),
            }
            break
        if placed is None:
            omitted.append({**item, "label": label, "detail": detail, "reason": "no free callout space"})
            continue
        occupied.append(placed["rect"])
        accepted.append(placed)
    return {
        "callouts": tuple(accepted),
        "omitted": tuple(omitted),
        "occupied_rectangles": tuple(occupied),
        "canvas_size": (canvas_width, canvas_height),
    }


# Explicit descriptive alias for callers/tests that prefer the longer name.
collision_aware_3d_label_layout = collision_aware_callout_layout


def compact_scene_legend_lines(items: Sequence[Mapping[str, object]], omitted: Sequence[Mapping[str, object]] = (), detector_mode: str = "combined") -> Tuple[str, ...]:
    """Build compact, count-aware legend/sidebar lines.

    Repeated Mylar/solid faces are grouped by detector and component family;
    sensor, mirror, and readout entries remain individually named.  The
    interactive details pane is scrollable, so grouping reduces collisions
    without silently dropping the component family or its count.
    """

    mode = detector_run_mode(detector_mode)
    if mode == "cherenkov":
        lines = [
            "MODE: CHERENKOV ONLY",
            "CHERENKOV: C quartz / C light guide",
            "READOUT: C photosensor / terminal",
            "PARTICLE PATH: black dashed, shared",
            "MIRRORS: teal terminal faces",
            "MYLAR: gold dashed overlays",
            "PHOTON PATH COLORS: Cherenkov outcomes",
        ]
    elif mode == "scintillator":
        lines = [
            "MODE: SCINTILLATOR ONLY",
            "SCINTILLATOR: S active / S light guide",
            "READOUT: S sensor / guide output",
            "PARTICLE PATH: black dashed, shared",
            "MIRRORS: teal terminal faces",
            "MYLAR: gold dashed overlays",
            "PHOTON PATH COLORS: scintillator outcomes",
        ]
    else:
        lines = [
            "MODE: ALL",
            "CHERENKOV: C quartz / C light guide",
            "SCINTILLATOR: S active / S light guide",
            "READOUT: C/S sensors and terminals",
            "PARTICLE PATH: black dashed, shared",
            "MIRRORS: teal terminal faces",
            "MYLAR: gold dashed overlays",
            "PHOTON PATH COLORS: detector-specific outcomes",
        ]
    def detail_lines(source_items: Sequence[Mapping[str, object]]) -> list[str]:
        grouped: Dict[tuple[str, str, str], list[Mapping[str, object]]] = {}
        for item in source_items:
            label = detector_callout_label(item)
            kind = str(item.get("kind", "geometry")).lower()
            detail = str(item.get("detail", item.get("label", "component")))
            # Keep readouts/sensors/mirrors separate; repeated geometry and
            # Mylar are the noisy families that should be counted.
            family = kind if kind in {"mylar_overlay", "geometry", "guide"} else f"{kind}:{item.get('component_id', label)}:{detail}"
            grouped.setdefault((label, family, kind), []).append(item)
        lines_out = []
        for (label, _family, _kind), group in grouped.items():
            first = group[0]
            details = []
            for item in group:
                value = str(item.get("detail", item.get("label", "component")))
                if value not in details:
                    details.append(value)
            detail = "; ".join(details)
            suffix = f" (×{len(group)})" if len(group) > 1 else ""
            lines_out.append(f"  {label} — {detail}{suffix}")
        return lines_out

    if items:
        lines.append("COMPONENT DETAILS:")
        lines.extend(detail_lines(items))
    if omitted:
        lines.append("CALLOUTS IN LEGEND (space):")
        lines.extend(detail_lines(omitted))
    return tuple(lines)

GEOMETRY_SHAPE_LABELS = {
    "box": "box / rectangular",
    "tapered_rectangular_prism": "tapered rectangular prism",
    "cylinder": "cylinder",
    "elliptical_cylinder": "elliptical cylinder",
    "truncated_cone": "truncated cone",
}
GEOMETRY_SHAPE_BY_LABEL = {label: key for key, label in GEOMETRY_SHAPE_LABELS.items()}
GEOMETRY_SHAPE_ALIASES = {
    "rectangular": "box",
    "rectangular_prism": "box",
    "box / rectangular": "box",
    "tapered box": "tapered_rectangular_prism",
    "tapered rectangular prism": "tapered_rectangular_prism",
}
GEOMETRY_AXES = ("X", "Y", "Z")
DETECTOR_ACTIVE_AREA_SHAPES = ("rectangular", "circular", "elliptical")
DETECTOR_ACTIVE_AREA_SHAPE_LABELS = {
    "rectangular": "rectangular",
    "circular": "circular",
    "elliptical": "elliptical",
}

# Event viewer filters.  ``photon_paths`` is deliberately an explicit
# "show everything" category; the other categories are a union, so a user
# can turn it off and inspect only selected terminal outcomes.  The strings
# are independent of the simulator's exact outcome spelling and are mapped
# by ``classify_photon_record`` below.
EVENT_VIEW_CATEGORY_LABELS = {
    "photon_paths": "Photon paths (all)",
    "detected": "Detected",
    "detected_plus_x": "Detected +x",
    "detected_minus_x": "Detected -x",
    "losses": "Losses (detector/TIR)",
    "escaped": "Escaped",
    "absorbed": "Absorbed",
    "max_bounces": "Max bounces",
    "unknown": "Other/unknown",
}


def _record_value(record, name: str, default=None):
    """Read a photon record field from dataclass- or mapping-like records."""

    if isinstance(record, Mapping):
        return record.get(name, default)
    return getattr(record, name, default)


def classify_photon_record(record) -> str:
    """Return one stable viewer category for a simulator photon record.

    The simulator currently emits strings such as ``detected_plus_x``,
    ``detector_plus_x_angle_rejected``, ``tir_plus_x``, ``escaped_x-``,
    ``absorbed`` and ``max_bounces``.  This helper intentionally uses the
    record's existing outcome and detector-face fields and has an explicit
    ``unknown`` fallback; it never infers a new physical outcome.
    """

    outcome = str(_record_value(record, "outcome", "") or "").strip().lower()
    face = str(_record_value(record, "detector_face", "") or "").strip().lower()
    detected = bool(_record_value(record, "detected", False))
    if outcome == "detected_plus_x" or face in {"plus_x", "+x", "x+"} and detected:
        return "detected_plus_x"
    if outcome == "detected_minus_x" or face in {"minus_x", "-x", "x-"} and detected:
        return "detected_minus_x"
    if outcome == "absorbed" or outcome.startswith("absorbed"):
        return "absorbed"
    if outcome == "max_bounces" or outcome.startswith("max_bounce"):
        return "max_bounces"
    if outcome.startswith("escaped") or outcome in {"escape", "escaped"}:
        return "escaped"
    if outcome.startswith(("detector_", "tir_")) or "_loss" in outcome or outcome.endswith("_rejected"):
        return "losses"
    # A detected flag without a recognized face is still a detected record,
    # but cannot be truthfully assigned to +x or -x.
    if detected:
        return "unknown"
    return "unknown"


def photon_record_categories(record) -> set[str]:
    """Return all checkbox categories matching one photon record."""

    category = classify_photon_record(record)
    categories = {category}
    if category in {"detected_plus_x", "detected_minus_x"}:
        categories.add("detected")
    return categories


def filter_photon_records(records, selected_categories) -> list:
    """Filter records for event viewers without changing simulation results.

    ``photon_paths`` selects every record.  With that category disabled, the
    remaining checkboxes form a union.  An empty selection intentionally
    returns no photon records while leaving geometry and the particle track
    available to the viewers.
    """

    selected = {str(value) for value in (selected_categories or ())}
    records = list(records or ())
    if "photon_paths" in selected:
        return records
    if not selected:
        return []
    return [record for record in records if photon_record_categories(record) & selected]


@dataclass(frozen=True)
class RunSettings:
    output_directory: Path
    events: int
    photons_per_event: int
    sampling_mode: str
    kinetic_mev: float
    particle_type: str
    mass_mev_c2: float
    entry_x_mm: float
    entry_y_mm: float
    inclination_deg: float
    azimuth_deg: float
    seed: int
    detector_index: float
    detector_acceptance_deg: float
    detector_qe: float
    collection_efficiency: float
    opposite_end_mode: str
    mirror_reflectivity: float
    opposite_detector_index: Optional[float]
    opposite_detector_acceptance_deg: Optional[float]
    opposite_detector_qe: Optional[float]
    opposite_collection_efficiency: Optional[float]
    detector_shape: str
    detector_active_width_mm: float
    detector_active_height_mm: float
    detector_active_radius_mm: float
    detector_active_semi_axis_a_mm: float
    detector_active_semi_axis_b_mm: float
    detector_active_offset_y_mm: float
    detector_active_offset_z_mm: float
    detector_active_rotation_deg: float
    opposite_detector_shape: str
    opposite_detector_active_width_mm: float
    opposite_detector_active_height_mm: float
    opposite_detector_active_radius_mm: float
    opposite_detector_active_semi_axis_a_mm: float
    opposite_detector_active_semi_axis_b_mm: float
    opposite_detector_active_offset_y_mm: float
    opposite_detector_active_offset_z_mm: float
    opposite_detector_active_rotation_deg: float
    track_margin_mm: float
    display_photons: int
    scintillator_shape: str
    scintillator_length_x_mm: float
    scintillator_width_y_mm: float
    scintillator_thickness_z_mm: float
    scintillator_width_y_start_mm: float
    scintillator_width_y_end_mm: float
    scintillator_thickness_z_start_mm: float
    scintillator_thickness_z_end_mm: float
    scintillator_axis: str
    scintillator_radius_mm: float
    scintillator_length_mm: float
    scintillator_semi_axis_a_mm: float
    scintillator_semi_axis_b_mm: float
    scintillator_radius_start_mm: float
    scintillator_radius_end_mm: float
    scintillator_center_x_mm: float
    scintillator_center_y_mm: float
    scintillator_center_z_mm: float
    scintillator_material: str
    scintillator_refractive_index: float
    guide_shape: str
    guide_length_x_mm: float
    guide_width_y_mm: float
    guide_thickness_z_mm: float
    guide_width_y_start_mm: float
    guide_width_y_end_mm: float
    guide_thickness_z_start_mm: float
    guide_thickness_z_end_mm: float
    guide_axis: str
    guide_radius_mm: float
    guide_length_mm: float
    guide_semi_axis_a_mm: float
    guide_semi_axis_b_mm: float
    guide_radius_start_mm: float
    guide_radius_end_mm: float
    guide_center_x_mm: float
    guide_center_y_mm: float
    guide_center_z_mm: float
    guide_material: str
    guide_refractive_index: float
    scintillator_gap_mm: float
    scintillator_sampling_mode: str
    scintillator_photon_count: int
    scintillator_yield_per_MeV: float
    scintillator_deposited_energy_MeV: Optional[float]
    scintillator_dEdx_MeV_per_m: Optional[float]
    scintillator_decay_time_ns: float
    scintillator_wavelength_min_nm: float
    scintillator_wavelength_max_nm: float
    scintillator_readout_index: float
    scintillator_readout_acceptance_deg: float
    scintillator_readout_qe: float
    scintillator_readout_collection_efficiency: float
    scintillator_absorption_length_mm: Optional[float]
    scintillator_max_bounces: int
    maximum_photon_bounces: int
    mylar_enabled: bool
    mylar_reflectivity_percent: float
    global_offset_x_mm: float
    global_offset_y_mm: float
    global_offset_z_mm: float
    global_beam_axis: str
    detector_separation_mm: float
    detector2_offset_x_mm: float
    detector2_offset_y_mm: float
    detector2_offset_z_mm: float
    detector1_guide_left_shape: str
    detector1_guide_right_shape: str
    detector2_scintillator_shape: str
    detector2_scintillator_material: str
    detector2_scintillator_refractive_index: float
    detector2_guide_left_shape: str
    detector2_guide_right_shape: str
    detector_1_z_mm: float
    detector_2_z_mm: float
    detector1_guide_left_enabled: bool
    detector1_guide_right_enabled: bool
    detector2_guide_left_enabled: bool
    detector2_guide_right_enabled: bool
    sensor_profile_values: Dict[str, str]
    detector_run_mode: str = "auto"
    assembly_1_rotation_axis: str = "None"
    assembly_1_rotation_angle_deg: float = 0.0
    assembly_2_rotation_axis: str = "None"
    assembly_2_rotation_angle_deg: float = 0.0
    gamma_source_profile_path: str = ""
    gamma_line_weight_511_keV: float = 0.0
    gamma_line_weight_835_keV: float = 0.0
    gamma_line_weight_1173_keV: float = 0.0
    gamma_line_weight_1275_keV: float = 0.0
    gamma_line_weight_1332_keV: float = 0.0
    gamma_line_weight_1461_keV: float = 0.0
    gamma_max_interactions: int = 1000
    gamma_max_light_samples: int = 10000


@dataclass(frozen=True)
class RunSnapshot:
    """Detached, read-only run inputs handed from Tk to the worker.

    The GUI owns the live Tk variables and mutable editor models.  A worker
    receives only this plain-data snapshot and reconstructs its private
    simulation inputs, so no worker operation can accidentally call Tk or
    observe a model while the editor is changing.
    """

    settings: Mapping[str, object]
    assembly_models: Mapping[str, Mapping[str, object]]
    transport_geometry: Mapping[str, object]


def _freeze_snapshot_data(value):
    """Recursively make serialized run data immutable for the worker."""

    if isinstance(value, Mapping):
        return MappingProxyType({str(key): _freeze_snapshot_data(item) for key, item in value.items()})
    if isinstance(value, (list, tuple)):
        return tuple(_freeze_snapshot_data(item) for item in value)
    return value


def _thaw_snapshot_data(value):
    """Turn immutable snapshot containers back into ordinary constructor data."""

    if isinstance(value, Mapping):
        return {str(key): _thaw_snapshot_data(item) for key, item in value.items()}
    if isinstance(value, tuple):
        return [_thaw_snapshot_data(item) for item in value]
    return value


def _transport_geometry_snapshot(geometry) -> Mapping[str, object]:
    """Serialize the detached transport geometry without retaining GUI state."""

    if isinstance(geometry, simulator.BarGeometry):
        payload = asdict(geometry)
        payload["geometry_type"] = "bar"
    elif isinstance(geometry, simulator.CompoundGeometry):
        payload = asdict(geometry)
        payload["geometry_type"] = "compound"
    elif isinstance(geometry, simulator.SolidGeometry):
        payload = asdict(geometry)
        payload["geometry_type"] = "solid"
    else:
        raise TypeError(f"unsupported transport geometry type: {type(geometry).__name__}")
    return _freeze_snapshot_data(_json_safe(payload))


def _transport_geometry_from_snapshot(snapshot: Mapping[str, object]):
    """Rebuild private simulator geometry from the worker-only snapshot."""

    payload = _thaw_snapshot_data(snapshot)
    geometry_type = payload.pop("geometry_type", None)
    if geometry_type == "bar":
        return simulator.BarGeometry(**payload)
    if geometry_type == "solid":
        payload["terminal_interfaces"] = tuple(dict(item) for item in payload.get("terminal_interfaces", ()))
        payload["sensor_components"] = tuple(dict(item) for item in payload.get("sensor_components", ()))
        return simulator.SolidGeometry(**payload)
    if geometry_type == "compound":
        payload["sections"] = tuple(
            simulator.SolidGeometry(
                **{
                    **dict(section),
                    "terminal_interfaces": tuple(dict(item) for item in section.get("terminal_interfaces", ())),
                    "sensor_components": tuple(dict(item) for item in section.get("sensor_components", ())),
                }
            )
            for section in payload.get("sections", ())
        )
        payload["terminal_interfaces"] = tuple(dict(item) for item in payload.get("terminal_interfaces", ()))
        payload["sensor_components"] = tuple(dict(item) for item in payload.get("sensor_components", ()))
        return simulator.CompoundGeometry(**payload)
    raise TypeError(f"unsupported transport geometry snapshot type: {geometry_type}")


def _make_run_snapshot(settings: RunSettings, geometry, assembly_models: Mapping[str, AssemblyGeometryModel]) -> RunSnapshot:
    """Capture all mutable GUI-owned run inputs on the Tk thread."""

    settings_data = _freeze_snapshot_data(_json_safe(asdict(settings)))
    mode = detector_run_mode(settings.detector_run_mode)
    # The v41 phase contract named this boundary "Building assemblies";
    # retain that phrase in source/provenance while adding the selected mode.
    participating = {"cherenkov", "scintillator"} if mode in {"combined", GAMMA_RUN_MODE} else {mode}
    assembly_data = {
        key: _json_safe(assembly_geometry_to_dict(deepcopy(model), key))
        for key, model in assembly_models.items()
        if key in participating
    }
    return RunSnapshot(
        settings=settings_data,
        assembly_models=_freeze_snapshot_data(assembly_data),
        transport_geometry=_transport_geometry_snapshot(deepcopy(geometry)),
    )


def _settings_from_snapshot(snapshot: RunSnapshot) -> RunSettings:
    payload = _thaw_snapshot_data(snapshot.settings)
    payload["output_directory"] = Path(payload["output_directory"])
    payload["sensor_profile_values"] = dict(payload.get("sensor_profile_values", {}))
    return RunSettings(**payload)


def _assembly_models_from_snapshot(snapshot: RunSnapshot) -> Dict[str, AssemblyGeometryModel]:
    return {
        key: assembly_geometry_from_dict(_thaw_snapshot_data(payload), key)
        for key, payload in snapshot.assembly_models.items()
    }


GeometrySection = AssemblySection


SCINTILLATOR_PROFILE_FORMAT = "quartz-cherenkov-scintillator-assembly"
SCINTILLATOR_PROFILE_VERSION = 1
# Populated immediately after default_values is defined; this early binding
# keeps the profile helpers close to the assembly model without forward calls.
SCINTILLATOR_PROFILE_KEYS = ()


def _bool_text(values: Mapping[str, object], key: str) -> bool:
    value = _value_text(values, key, "0").lower()
    if value in {"1", "true", "yes", "on", "enabled"}:
        return True
    if value in {"0", "false", "no", "off", "disabled", ""}:
        return False
    raise ValueError(f"{key} must be a boolean")


def _optional_float(values: Mapping[str, object], key: str) -> Optional[float]:
    value = _value_text(values, key, "")
    if not value:
        return None
    try:
        number = float(value)
    except (TypeError, ValueError) as exc:
        raise ValueError(f"{key} must be numeric when supplied") from exc
    if not math.isfinite(number):
        raise ValueError(f"{key} must be finite")
    return number


def _assembly_section_from_values(values: Mapping[str, object], prefix: str) -> simulator.SolidGeometry:
    def number(suffix: str, default: float = 0.0) -> float:
        text = _value_text(values, f"{prefix}_{suffix}", str(default))
        try:
            value = float(text)
        except (TypeError, ValueError) as exc:
            raise ValueError(f"{prefix}_{suffix} must be numeric") from exc
        if not math.isfinite(value):
            raise ValueError(f"{prefix}_{suffix} must be finite")
        return value
    shape = _value_text(values, f"{prefix}_shape", "box").lower()
    aliases = {"rectangular": "box", "rectangular_prism": "box", "tapered box": "tapered_rectangular_prism"}
    shape = aliases.get(shape, shape)
    material = _value_text(values, f"{prefix}_material", "scintillator")
    material_key = {"plastic scintillator": "scintillator", "light guide": "light_guide"}.get(material.lower(), material.lower())
    geometry = simulator.SolidGeometry(
        shape=shape, material=material_key,
        center_x_m=number("center_x_mm") * 1e-3, center_y_m=number("center_y_mm") * 1e-3, center_z_m=number("center_z_mm") * 1e-3,
        length_x_m=number("length_x_mm", 20.0) * 1e-3, width_y_m=number("width_y_mm", 20.0) * 1e-3, thickness_z_m=number("thickness_z_mm", 10.0) * 1e-3,
        width_y_start_m=number("width_y_start_mm", 20.0) * 1e-3, width_y_end_m=number("width_y_end_mm", 20.0) * 1e-3,
        thickness_z_start_m=number("thickness_z_start_mm", 10.0) * 1e-3, thickness_z_end_m=number("thickness_z_end_mm", 10.0) * 1e-3,
        axis=_value_text(values, f"{prefix}_axis", "X").upper(), radius_m=number("radius_mm", 10.0) * 1e-3, length_m=number("length_mm", 20.0) * 1e-3,
        semi_axis_a_m=number("semi_axis_a_mm", 10.0) * 1e-3, semi_axis_b_m=number("semi_axis_b_mm", 5.0) * 1e-3,
        radius_start_m=number("radius_start_mm", 10.0) * 1e-3, radius_end_m=number("radius_end_mm", 8.0) * 1e-3,
        refractive_index_value=number("refractive_index", 1.5),
    )
    geometry.validate()
    return geometry


def scintillator_assembly_from_values(values: Mapping[str, object]) -> simulator.ScintillatorAssembly:
    """Build a composed assembly for compatibility with existing callers.

    New callers should use the three explicit builders below.  This wrapper
    keeps the previous API and uses the legacy default quartz reference.
    """

    local = scintillator_local_geometry_from_values(values)
    global_geometry = scintillator_global_geometry_from_values(values)
    reflectivity = float(_value_text(values, "mylar_reflectivity_percent", "95")) / 100.0
    plus_surface = _value_text(values, "scintillator_terminal_plus_surface", "Detector/Photosensor coupling")
    minus_surface = _value_text(values, "scintillator_terminal_minus_surface", "Open")
    if normalize_terminal_surface(plus_surface) == "guide" and not _scintillator_guide_configured(values):
        plus_surface = "Open"
    if normalize_terminal_surface(minus_surface) == "guide" and not _scintillator_guide_configured(values):
        minus_surface = "Open"
    return simulator.compose_scintillator_assembly(
        local, global_geometry, simulator.BarGeometry(),
        mylar_enabled=_bool_text(values, "mylar_enabled"), mylar_reflectivity=reflectivity,
        terminal_plus_surface=plus_surface, terminal_minus_surface=minus_surface,
    )


def scintillator_assembly_from_model(
    model: AssemblyGeometryModel, *, values: Optional[Mapping[str, object]] = None,
    particle: Optional[simulator.Particle] = None,
) -> simulator.ScintillatorAssembly:
    """Convert the canonical complete assembly to the transport config.

    Every valid active and guide row is retained in physical X order.  This
    adapter is the run authority; the legacy value-field adapter remains only
    for compatibility and one-guide profiles.
    """

    if not isinstance(model, AssemblyGeometryModel):
        raise TypeError("scintillator transport requires the canonical AssemblyGeometryModel")
    model.validate()
    if model.guide_enabled and not model.guide_sections:
        raise ValueError("Scintillator guide is enabled but no guide section is configured")
    if model.guide_enabled and len(model.valid_guide_sections) != len(model.guide_sections):
        invalid = next(index for index, section in enumerate(model.guide_sections) if not model._section_is_valid(section))
        raise ValueError(f"Scintillator guide section {invalid + 1} has an unsupported shape, axis, or material")
    active = [(section.to_simulation_geometry(), "active", index, section.name) for index, section in enumerate(model.sections)]
    guides = [(section.to_simulation_geometry(), "guide", index, section.name) for index, section in enumerate(model.valid_guide_sections)] if model.guide_enabled else []
    ordered = sorted(active + guides, key=lambda item: float(item[0].bounds_min_m[0]))
    if not ordered:
        raise ValueError("Scintillator canonical assembly has no transport sections")
    solids = tuple(item[0] for item in ordered)
    active_indices = tuple(position for position, item in enumerate(ordered) if item[1] == "active")
    guide_indices = tuple(position for position, item in enumerate(ordered) if item[1] == "guide")
    transport_geometry = simulator.CompoundGeometry(
        solids, allow_interface_mismatch=True, wrapped_external=bool(model.mylar_enabled),
        mylar_reflectivity=float(model.mylar_reflectivity),
        terminal_interfaces=(
            {"face": "minus_x", "axis_face": "YZ", "state": model.terminal_state("minus_x")},
            {"face": "plus_x", "axis_face": "YZ", "state": model.terminal_state("plus_x")},
        ),
    )
    # The editable assembly-local impact is the transport source of truth.
    # The shared global beam contract validates that this point agrees with
    # the Run particle trajectory when both detectors participate; falling
    # back to particle.entry_x/y here was the source of a local/global split.
    impact_x = float(model.impact_x_mm) * 1e-3
    impact_y = float(model.impact_y_mm) * 1e-3
    return simulator.ScintillatorAssembly.from_transport_geometry(
        transport_geometry, active_section_indices=active_indices,
        guide_section_indices=guide_indices,
        section_names=tuple(item[3] for item in ordered),
        mylar_enabled=bool(model.mylar_enabled), mylar_reflectivity=float(model.mylar_reflectivity),
        terminal_plus_surface=model.terminal_plus_surface,
        terminal_minus_surface=model.terminal_minus_surface,
        coordinate_frame="scintillator_local_m", particle_impact_x_m=impact_x,
        particle_impact_y_m=impact_y,
    )


def _localize_section(section: simulator.SolidGeometry, start_x_m: float) -> simulator.SolidGeometry:
    """Move a solid so its X minimum is ``start_x_m`` without changing size."""

    centered = replace(section, center_x_m=0.0, center_y_m=0.0, center_z_m=0.0)
    return replace(centered, center_x_m=start_x_m - float(centered.bounds_min_m[0]))


def scintillator_local_geometry_from_values(values: Mapping[str, object]) -> simulator.ScintillatorLocalGeometry:
    """Build only local solids; no quartz gap or global placement is read."""

    # The local beam-monitor convention is centered on the active volume in
    # X/Y.  Guides remain sequential toward +X, but impact_x=0 is the active
    # center rather than its -X face.  v44 localized the active low face to
    # zero, making the default beam impact a boundary source and leaving the
    # canonical editor marker in a different local frame.
    scintillator_template = _assembly_section_from_values(values, "scintillator")
    scintillator = _localize_section(
        scintillator_template,
        -0.5 * scintillator_template.length_x_extent_m,
    )
    guide = None
    if _scintillator_guide_configured(values):
        guide_template = _localize_section(_assembly_section_from_values(values, "guide"), 0.0)
        guide_start = float(scintillator.bounds_max_m[0])
        guide = _localize_section(guide_template, guide_start)
    local = simulator.ScintillatorLocalGeometry(
        scintillator, guide,
        particle_impact_x_m=float(_value_text(values, "scintillator_impact_x_mm", "0")) * 1e-3,
        particle_impact_y_m=float(_value_text(values, "scintillator_impact_y_mm", "0")) * 1e-3,
    )
    local.validate()
    return local


def _scintillator_guide_configured(values: Mapping[str, object]) -> bool:
    """Apply the explicit enable-and-configured rule to legacy value stores."""

    if not _bool_text(values, "scintillator_guide_enabled"):
        return False
    shape = _value_text(values, "guide_shape", "").strip().lower()
    if shape in {"", "none", "disabled", "off"}:
        return False
    try:
        guide = _assembly_section_from_values(values, "guide")
        guide.validate()
    except (TypeError, ValueError):
        return False
    return True


def scintillator_global_geometry_from_values(values: Mapping[str, object]) -> simulator.ScintillatorGlobalGeometry:
    """Build only global placement/axis values from the single global source."""

    global_geometry = simulator.ScintillatorGlobalGeometry(
        gap_from_quartz_m=float(_value_text(values, "scintillator_gap_mm", "0")) * 1e-3,
        offset_x_m=float(_value_text(values, "global_offset_x_mm", "0")) * 1e-3,
        offset_y_m=float(_value_text(values, "global_offset_y_mm", "0")) * 1e-3,
        offset_z_m=float(_value_text(values, "global_offset_z_mm", "0")) * 1e-3,
        beam_axis=_value_text(values, "global_beam_axis", "Z").upper(),
    )
    global_geometry.validate()
    return global_geometry


def scintillator_global_assembly_from_values(values: Mapping[str, object], quartz_geometry=None) -> simulator.ScintillatorAssembly:
    local = scintillator_local_geometry_from_values(values)
    global_geometry = scintillator_global_geometry_from_values(values)
    plus_surface = _value_text(values, "scintillator_terminal_plus_surface", "Detector/Photosensor coupling")
    minus_surface = _value_text(values, "scintillator_terminal_minus_surface", "Open")
    if normalize_terminal_surface(plus_surface) == "guide" and not _scintillator_guide_configured(values):
        plus_surface = "Open"
    if normalize_terminal_surface(minus_surface) == "guide" and not _scintillator_guide_configured(values):
        minus_surface = "Open"
    return simulator.compose_scintillator_assembly(
        local, global_geometry, quartz_geometry or simulator.BarGeometry(),
        mylar_enabled=_bool_text(values, "mylar_enabled"),
        mylar_reflectivity=float(_value_text(values, "mylar_reflectivity_percent", "95")) / 100.0,
        terminal_plus_surface=plus_surface, terminal_minus_surface=minus_surface,
    )


def _solid_shape_variant(solid: simulator.SolidGeometry, values: Mapping[str, object], key: str, *, material_key: Optional[str] = None, refractive_key: Optional[str] = None) -> simulator.SolidGeometry:
    """Return a local solid variant without changing its dimensional fields."""

    shape = _value_text(values, key, solid.shape).lower()
    shape = GEOMETRY_SHAPE_ALIASES.get(shape, shape)
    material = _value_text(values, material_key, solid.material) if material_key else solid.material
    refractive = solid.refractive_index_value
    if refractive_key:
        try:
            refractive = float(_value_text(values, refractive_key, str(refractive if refractive is not None else 1.5)))
        except (TypeError, ValueError) as exc:
            raise ValueError(f"{refractive_key} must be numeric") from exc
    return replace(solid, shape=shape, material=material, refractive_index_value=refractive)


def _detector_assembly_from_local(
    local: simulator.ScintillatorLocalGeometry,
    values: Mapping[str, object],
    detector_id: int,
    *,
    left_shape_key: str,
    right_shape_key: str,
    scintillator_shape_key: Optional[str] = None,
    scintillator_material_key: Optional[str] = None,
    scintillator_refractive_key: Optional[str] = None,
    axis: str = "Z",
) -> simulator.DetectorAssembly:
    """Build one independent local detector.

    The legacy X branch remains sequential for compatibility.  The new Z
    branch keeps the active solid in the local XY plane and places optional
    guides on its Y-/Y+ faces; placement is applied later by the global model.
    """

    scintillator = local.scintillator
    if scintillator_shape_key:
        scintillator = _solid_shape_variant(
            scintillator, values, scintillator_shape_key,
            material_key=scintillator_material_key,
            refractive_key=scintillator_refractive_key,
        )
        scintillator = _localize_section(scintillator, 0.0)
    axis = str(axis).upper()
    if axis == "Z":
        if scintillator.shape in {"cylinder", "elliptical_cylinder", "truncated_cone"}:
            scintillator = replace(scintillator, axis="Z")
        scintillator = replace(scintillator, center_x_m=0.0, center_y_m=0.0, center_z_m=0.0)
        def z_guide(shape_key: str, enabled_key: str, side: float) -> Optional[simulator.SolidGeometry]:
            if not _bool_text(values, enabled_key):
                return None
            guide = _solid_shape_variant(local.light_guide, values, shape_key)
            # Box/tapered solids retain their XY footprint and active Z span.
            # Rotational shapes are represented with their axis along Y so the
            # guide is a side volume and still has a finite Z footprint.
            if guide.shape in {"cylinder", "elliptical_cylinder", "truncated_cone"}:
                guide = replace(guide, axis="Y")
            guide_half_y = max(guide.width_y_extent_m, 1e-6) / 2.0
            guide_y = scintillator.bounds_max_m[1] + guide_half_y if side > 0 else scintillator.bounds_min_m[1] - guide_half_y
            guide = replace(guide, center_x_m=scintillator.center_x_m,
                            center_y_m=guide_y,
                            center_z_m=scintillator.center_z_m)
            return guide
        left = z_guide(left_shape_key, f"detector{detector_id}_guide_left_enabled", -1.0)
        right = z_guide(right_shape_key, f"detector{detector_id}_guide_right_enabled", 1.0)
    else:
        right = left = None
        if local.light_guide is not None and _bool_text(values, f"detector{detector_id}_guide_right_enabled"):
            right = _solid_shape_variant(local.light_guide, values, right_shape_key)
            right = _localize_section(right, scintillator.bounds_max_m[0])
        if local.light_guide is not None and _bool_text(values, f"detector{detector_id}_guide_left_enabled"):
            left = _solid_shape_variant(local.light_guide, values, left_shape_key)
            left = _localize_section(left, scintillator.bounds_min_m[0] - left.length_x_extent_m)
    detector = simulator.DetectorAssembly(
        scintillator=scintillator,
        guide_left=left,
        guide_right=right,
        detector_id=detector_id,
        mylar_enabled=_bool_text(values, "mylar_enabled"),
        mylar_reflectivity=float(_value_text(values, "mylar_reflectivity_percent", "95")) / 100.0,
        active_axis=axis,
    )
    detector.validate()
    return detector


def two_detector_beam_monitor_from_values(values: Mapping[str, object], quartz_geometry=None) -> simulator.TwoDetectorBeamMonitor:
    """Build the explicit two-detector global model from GUI values.

    Existing local scintillator/guide fields remain the shared local template;
    each detector receives separate solid instances.  Detector 2 has explicit
    scintillator material/index and left/right guide-shape controls.
    """

    quartz = quartz_geometry or simulator.BarGeometry()
    local = scintillator_local_geometry_from_values(values)
    # The side-guide controls are subordinate to the shared Scintillator
    # option and its configured local guide section.  This prevents stale
    # legacy side flags from resurrecting guide volumes after a disable.
    effective_values = dict(values)
    if not _scintillator_guide_configured(values):
        for key in ("detector1_guide_left_enabled", "detector1_guide_right_enabled", "detector2_guide_left_enabled", "detector2_guide_right_enabled"):
            effective_values[key] = "0"
    detector_1 = _detector_assembly_from_local(
        local, effective_values, 1,
        left_shape_key="detector1_guide_left_shape",
        right_shape_key="detector1_guide_right_shape",
        axis=_value_text(effective_values, "global_beam_axis", "Z").upper(),
    )
    detector_2 = _detector_assembly_from_local(
        local, effective_values, 2,
        left_shape_key="detector2_guide_left_shape",
        right_shape_key="detector2_guide_right_shape",
        scintillator_shape_key="detector2_scintillator_shape",
        scintillator_material_key="detector2_scintillator_material",
        scintillator_refractive_key="detector2_scintillator_refractive_index",
        axis=_value_text(effective_values, "global_beam_axis", "Z").upper(),
    )
    def number(key: str, fallback: str = "0") -> float:
        try:
            value = float(_value_text(values, key, fallback))
        except (TypeError, ValueError) as exc:
            raise ValueError(f"{key} must be numeric") from exc
        if not math.isfinite(value):
            raise ValueError(f"{key} must be finite")
        return value
    beam_axis = _value_text(values, "global_beam_axis", "Z").upper()
    monitor = simulator.TwoDetectorBeamMonitor(
        detector_1=detector_1,
        detector_2=detector_2,
        quartz_gap_m=number("scintillator_gap_mm", "0") * 1e-3,
        detector_separation_m=number("detector_separation_mm", "25") * 1e-3,
        placement_1=simulator.DetectorPlacement(
            offset_x_m=number("global_offset_x_mm") * 1e-3,
            offset_y_m=number("global_offset_y_mm") * 1e-3,
            offset_z_m=number("global_offset_z_mm") * 1e-3,
            rotation_axis=_value_text(values, "assembly_1_rotation_axis", "None"),
            rotation_angle_rad=math.radians(number("assembly_1_rotation_angle_deg", "0")),
        ),
        placement_2=simulator.DetectorPlacement(
            offset_x_m=number("detector2_offset_x_mm") * 1e-3,
            offset_y_m=number("detector2_offset_y_mm") * 1e-3,
            offset_z_m=number("detector2_offset_z_mm") * 1e-3,
            rotation_axis=_value_text(values, "assembly_2_rotation_axis", "None"),
            rotation_angle_rad=math.radians(number("assembly_2_rotation_angle_deg", "0")),
        ),
        beam_axis=beam_axis,
        detector_1_z_m=number("detector_1_z_mm", "20") * 1e-3 if beam_axis == "Z" else None,
        detector_2_z_m=number("detector_2_z_mm", "55") * 1e-3 if beam_axis == "Z" else None,
    )
    monitor.validate(quartz)
    return monitor


def scintillator_config_from_values(values: Mapping[str, object], *, max_bounces_override: Optional[int] = None) -> simulator.ScintillatorConfig:
    """Build the independent scintillation generation/readout configuration."""

    mode_value = _value_text(values, "scintillator_sampling_mode", SCINTILLATION_SAMPLING_MODE_LABELS["fixed"])
    mode = SCINTILLATION_SAMPLING_MODE_BY_LABEL.get(mode_value, mode_value.lower())
    def number(key: str, default: float = 0.0) -> float:
        text = _value_text(values, key, str(default))
        try:
            value = float(text)
        except (TypeError, ValueError) as exc:
            raise ValueError(f"{key} must be numeric") from exc
        if not math.isfinite(value):
            raise ValueError(f"{key} must be finite")
        return value
    # Detector participation is owned by RunSettings.detector_run_mode.  The
    # config remains enabled as a physical-parameter object so old profiles
    # cannot silently disable a selected scintillator run.
    if max_bounces_override is not None:
        try:
            numeric_override = float(max_bounces_override)
            if not math.isfinite(numeric_override) or numeric_override != int(numeric_override):
                raise ValueError("not an integer")
            max_bounces = int(numeric_override)
        except (TypeError, ValueError) as exc:
            raise ValueError("maximum photon bounces override must be an integer") from exc
        if max_bounces <= 0:
            raise ValueError("maximum photon bounces override must be positive")
    else:
        max_bounces = int(number("scintillator_max_bounces", 100))
    config = simulator.ScintillatorConfig(
        enabled=True, sampling_mode=mode,
        photon_count=int(number("scintillator_photon_count", 100)), photon_yield_per_MeV=number("scintillator_yield_per_MeV", 10000.0),
        deposited_energy_MeV=_optional_float(values, "scintillator_deposited_energy_MeV"), dEdx_MeV_per_m=_optional_float(values, "scintillator_dEdx_MeV_per_m"),
        decay_time_s=scintillator_decay_time_s_from_ns(number("scintillator_decay_time_ns", SCINTILLATOR_DECAY_TIME_DEFAULT_NS)),
        emission_wavelength_min_m=number("scintillator_wavelength_min_nm", 350.0) * 1e-9, emission_wavelength_max_m=number("scintillator_wavelength_max_nm", 550.0) * 1e-9,
        readout_index=number("scintillator_readout_index", 1.0), readout_acceptance_rad=math.radians(number("scintillator_readout_acceptance_deg", 90.0)),
        readout_quantum_efficiency=number("scintillator_readout_qe", 1.0), readout_collection_efficiency=number("scintillator_readout_collection_efficiency", 1.0),
        absorption_length_m=None if _optional_float(values, "scintillator_absorption_length_mm") is None else _optional_float(values, "scintillator_absorption_length_mm") * 1e-3,
        max_bounces=max_bounces,
    )
    config.validate()
    return config


def scintillator_decay_time_s_from_ns(decay_time_ns: object) -> float:
    """Validate the GUI value and convert nanoseconds to simulator seconds."""

    try:
        value_ns = float(decay_time_ns)
    except (TypeError, ValueError) as exc:
        raise ValueError("Scintillation decay time must be numeric") from exc
    if not math.isfinite(value_ns) or value_ns < 0.0:
        raise ValueError("Scintillation decay time must be finite and non-negative")
    return value_ns * 1e-9


def validate_detector_characteristics_values(values: Mapping[str, object]) -> simulator.ScintillatorConfig:
    """Validate response plus the detector-owned Mylar optical setting."""

    config = scintillator_config_from_values(values)
    reflectivity = float(_value_text(values, "mylar_reflectivity_percent", "95"))
    if not math.isfinite(reflectivity) or not 80.0 <= reflectivity <= 100.0:
        raise ValueError("Mylar reflectivity must be between 80 and 100 percent")
    return config


def scintillator_profile_from_values(values: Mapping[str, object]) -> Dict[str, object]:
    """Return a versioned assembly profile; validation precedes serialization."""

    assembly = scintillator_assembly_from_values(values)
    config = scintillator_config_from_values(values)
    payload = {key: _value_text(values, key, "") for key in SCINTILLATOR_PROFILE_KEYS}
    payload.update({
        "format": SCINTILLATOR_PROFILE_FORMAT, "version": SCINTILLATOR_PROFILE_VERSION,
        "assembly_preview": assembly.preview(), "config": asdict(config),
        # New grouped sections make the legacy all-in-one profile explicit
        # while keeping its original flat keys readable by old versions.
        "local_geometry": local_geometry_profile_from_values(values),
        "detector_characteristics": detector_characteristics_profile_from_values(values),
        "global_geometry": global_geometry_profile_from_values(values),
    })
    return payload


def apply_scintillator_profile(values: Mapping[str, object], profile: Mapping[str, object]) -> Dict[str, object]:
    """Validate a profile and apply all fields atomically."""

    if not isinstance(profile, Mapping) or profile.get("format") != SCINTILLATOR_PROFILE_FORMAT or profile.get("version") != SCINTILLATOR_PROFILE_VERSION:
        raise ValueError("Invalid or unsupported scintillator assembly profile")
    missing = [key for key in SCINTILLATOR_PROFILE_KEYS if key not in profile]
    if missing:
        raise ValueError(f"Scintillator profile is missing keys: {', '.join(missing)}")
    incoming = dict(profile)
    candidate = dict(values)
    for key in SCINTILLATOR_PROFILE_KEYS:
        if key in incoming and key in candidate:
            candidate[key] = incoming[key]
    # Full validation happens against the candidate, before any mutation.
    scintillator_assembly_from_values(candidate)
    scintillator_config_from_values(candidate)
    old = {key: _value_text(values, key, "") for key in SCINTILLATOR_PROFILE_KEYS if key in values}
    try:
        for key, value in candidate.items():
            if key not in old:
                continue
            target = values[key]
            display = str(value)
            target.set(display) if hasattr(target, "set") and callable(target.set) else values.__setitem__(key, display)  # type: ignore[index]
    except Exception:
        for key, value in old.items():
            target = values[key]
            target.set(value) if hasattr(target, "set") and callable(target.set) else values.__setitem__(key, value)  # type: ignore[index]
        raise
    return scintillator_profile_from_values(values)


def save_scintillator_profile(path: Path | str, values: Mapping[str, object]) -> Dict[str, object]:
    profile = scintillator_profile_from_values(values)
    Path(path).expanduser().resolve().write_text(json.dumps(profile, indent=2, sort_keys=True) + "\n", encoding="utf-8")
    return profile


def load_scintillator_profile(path: Path | str, values: Mapping[str, object]) -> Dict[str, object]:
    profile = json.loads(Path(path).expanduser().resolve().read_text(encoding="utf-8"))
    return apply_scintillator_profile(values, profile)


def scintillator_assembly_preview(values: Mapping[str, object]) -> Dict[str, object]:
    return scintillator_assembly_from_values(values).preview()


def assembly_preview_from_values(values: Mapping[str, object]) -> Dict[str, object]:
    """Compatibility-friendly name for headless preview callers."""

    return scintillator_assembly_preview(values)


SCINTILLATOR_LOCAL_PROFILE_FORMAT = "quartz-cherenkov-scintillator-geometry"
SCINTILLATOR_DETECTOR_PROFILE_FORMAT = "quartz-cherenkov-scintillator-detector"
SCINTILLATOR_GLOBAL_PROFILE_FORMAT = "quartz-cherenkov-global-geometry"
SCINTILLATOR_COMPONENT_PROFILE_VERSION = 1
SCINTILLATOR_LOCAL_PROFILE_KEYS = (
    "scintillator_shape", "scintillator_length_x_mm", "scintillator_width_y_mm", "scintillator_thickness_z_mm",
    "scintillator_width_y_start_mm", "scintillator_width_y_end_mm", "scintillator_thickness_z_start_mm", "scintillator_thickness_z_end_mm",
    "scintillator_axis", "scintillator_radius_mm", "scintillator_length_mm", "scintillator_semi_axis_a_mm", "scintillator_semi_axis_b_mm",
    "scintillator_radius_start_mm", "scintillator_radius_end_mm", "scintillator_material", "scintillator_refractive_index",
    "guide_shape", "guide_length_x_mm", "guide_width_y_mm", "guide_thickness_z_mm", "guide_width_y_start_mm", "guide_width_y_end_mm",
    "guide_thickness_z_start_mm", "guide_thickness_z_end_mm", "guide_axis", "guide_radius_mm", "guide_length_mm", "guide_semi_axis_a_mm",
    "guide_semi_axis_b_mm", "guide_radius_start_mm", "guide_radius_end_mm", "guide_material", "guide_refractive_index",
)
SCINTILLATOR_DETECTOR_PROFILE_KEYS = (
    "scintillator_sampling_mode", "scintillator_photon_count", "scintillator_yield_per_MeV",
    "scintillator_deposited_energy_MeV", "scintillator_dEdx_MeV_per_m", "scintillator_decay_time_ns", "scintillator_wavelength_min_nm",
    "scintillator_wavelength_max_nm", "scintillator_readout_index", "scintillator_readout_acceptance_deg", "scintillator_readout_qe",
    "scintillator_readout_collection_efficiency", "scintillator_absorption_length_mm", "scintillator_max_bounces", "mylar_enabled",
    "mylar_reflectivity_percent",
)
SCINTILLATOR_GLOBAL_PROFILE_KEYS = ("scintillator_gap_mm", "global_offset_x_mm", "global_offset_y_mm", "global_offset_z_mm", "global_beam_axis")
TWO_DETECTOR_GLOBAL_PROFILE_KEYS = (
    "detector_separation_mm", "detector_1_z_mm", "detector_2_z_mm",
    "detector2_offset_x_mm", "detector2_offset_y_mm", "detector2_offset_z_mm",
    "detector1_guide_left_enabled", "detector1_guide_right_enabled",
    "detector2_guide_left_enabled", "detector2_guide_right_enabled",
    "detector1_guide_left_shape", "detector1_guide_right_shape",
    "detector2_scintillator_shape", "detector2_scintillator_material", "detector2_scintillator_refractive_index",
    "detector2_guide_left_shape", "detector2_guide_right_shape",
)


def _component_profile(values: Mapping[str, object], fmt: str, keys: Sequence[str], validator: Callable[[Mapping[str, object]], object]) -> Dict[str, object]:
    payload = {"format": fmt, "version": SCINTILLATOR_COMPONENT_PROFILE_VERSION}
    payload["values"] = {key: _value_text(values, key, "") for key in keys}
    validator({**values, **payload["values"]})
    return payload


def local_geometry_profile_from_values(values: Mapping[str, object]) -> Dict[str, object]:
    return _component_profile(values, SCINTILLATOR_LOCAL_PROFILE_FORMAT, SCINTILLATOR_LOCAL_PROFILE_KEYS, scintillator_local_geometry_from_values)


def validate_local_geometry_profile(profile: Mapping[str, object]) -> Dict[str, object]:
    return _validate_component_profile(profile, SCINTILLATOR_LOCAL_PROFILE_FORMAT, SCINTILLATOR_LOCAL_PROFILE_KEYS)


def detector_characteristics_profile_from_values(values: Mapping[str, object]) -> Dict[str, object]:
    return _component_profile(values, SCINTILLATOR_DETECTOR_PROFILE_FORMAT, SCINTILLATOR_DETECTOR_PROFILE_KEYS, validate_detector_characteristics_values)


def validate_detector_characteristics_profile(profile: Mapping[str, object]) -> Dict[str, object]:
    return _validate_component_profile(profile, SCINTILLATOR_DETECTOR_PROFILE_FORMAT, SCINTILLATOR_DETECTOR_PROFILE_KEYS)


def global_geometry_profile_from_values(values: Mapping[str, object]) -> Dict[str, object]:
    payload = _component_profile(values, SCINTILLATOR_GLOBAL_PROFILE_FORMAT, SCINTILLATOR_GLOBAL_PROFILE_KEYS, scintillator_global_geometry_from_values)
    candidate = dict(values)
    candidate.update({key: _value_text(values, key, default_values().get(key, "")) for key in TWO_DETECTOR_GLOBAL_PROFILE_KEYS})
    monitor = two_detector_beam_monitor_from_values(candidate)
    payload["values"].update({key: _value_text(candidate, key, "") for key in TWO_DETECTOR_GLOBAL_PROFILE_KEYS})
    payload["two_detector"] = {
        "model": "TwoDetectorBeamMonitor",
        "detector_1": {
            "signed_offset_mm": [_value_text(candidate, "global_offset_x_mm", "0"), _value_text(candidate, "global_offset_y_mm", "0"), _value_text(candidate, "global_offset_z_mm", "0")],
            "guide_left_shape": _value_text(candidate, "detector1_guide_left_shape", "box"),
            "guide_right_shape": _value_text(candidate, "detector1_guide_right_shape", "box"),
        },
        "detector_2": {
            "signed_offset_mm": [_value_text(candidate, "detector2_offset_x_mm", "0"), _value_text(candidate, "detector2_offset_y_mm", "0"), _value_text(candidate, "detector2_offset_z_mm", "0")],
            "scintillator_shape": _value_text(candidate, "detector2_scintillator_shape", "box"),
            "scintillator_material": _value_text(candidate, "detector2_scintillator_material", "scintillator"),
            "scintillator_refractive_index": _value_text(candidate, "detector2_scintillator_refractive_index", "1.58"),
            "guide_left_shape": _value_text(candidate, "detector2_guide_left_shape", "box"),
            "guide_right_shape": _value_text(candidate, "detector2_guide_right_shape", "box"),
        },
        "detector_separation_mm": _value_text(candidate, "detector_separation_mm", "25"),
        "beam_axis": _value_text(candidate, "global_beam_axis", "Z"),
        "detector_1_z_mm": _value_text(candidate, "detector_1_z_mm", "20"),
        "detector_2_z_mm": _value_text(candidate, "detector_2_z_mm", "55"),
        "relative_placement": monitor.relative_placement(),
    }
    payload["active_geometry"] = {
        "model": "two independent XY active planes",
        "detector_1": payload["two_detector"]["detector_1"],
        "detector_2": payload["two_detector"]["detector_2"],
    }
    payload["global_placement"] = {
        "beam_axis": _value_text(candidate, "global_beam_axis", "Z"),
        "detector_1_z_mm": _value_text(candidate, "detector_1_z_mm", "20"),
        "detector_2_z_mm": _value_text(candidate, "detector_2_z_mm", "55"),
        "signed_offsets_mm": {
            "detector_1": [_value_text(candidate, "global_offset_x_mm", "0"), _value_text(candidate, "global_offset_y_mm", "0"), _value_text(candidate, "global_offset_z_mm", "0")],
            "detector_2": [_value_text(candidate, "detector2_offset_x_mm", "0"), _value_text(candidate, "detector2_offset_y_mm", "0"), _value_text(candidate, "detector2_offset_z_mm", "0")],
        },
    }
    payload["sensors"] = sensor_profile_from_values(candidate)
    return payload


SENSOR_PROFILE_KEYS = tuple(
    key for index in (1, 2) for key in (
        f"sensor_{index}_component", f"sensor_{index}_coupling", f"sensor_{index}_side",
        f"sensor_{index}_position_x_mm", f"sensor_{index}_position_y_mm", f"sensor_{index}_position_z_mm",
        f"sensor_{index}_automatic_position",
        f"sensor_{index}_shape", f"sensor_{index}_width_mm", f"sensor_{index}_height_mm",
        f"sensor_{index}_radius_mm", f"sensor_{index}_semi_axis_a_mm", f"sensor_{index}_semi_axis_b_mm",
        f"sensor_{index}_index", f"sensor_{index}_acceptance_deg", f"sensor_{index}_qe",
        f"sensor_{index}_collection", f"sensor_{index}_timing_sigma_ns",
        f"sensor_{index}_ideal_active",
    )
)


def normalize_sensor_coupling(value: object) -> str:
    """Normalize current and legacy direct/guide routing values."""

    return simulator.DetectorSensor.normalize_coupling_mode(value)


def derived_sensor_component(assembly: str, coupling: object) -> str:
    """Return the read-only component status shown in each Sensor tab."""

    route = normalize_sensor_coupling(coupling)
    if route == "light_guide":
        return "Complete light guide"
    if assembly == "cherenkov":
        return "Active quartz radiator"
    if assembly == "scintillator":
        return "Active scintillator"
    raise ValueError("assembly must be cherenkov or scintillator")


def _set_value(values: Mapping[str, object], key: str, value: object) -> None:
    """Set a Tk variable or plain mapping entry for headless callers."""

    target = values[key]
    if hasattr(target, "set") and callable(target.set):
        target.set(str(value))
    else:
        values[key] = str(value)  # type: ignore[index]


def ideal_sensor_reference(assembly: str) -> Dict[str, float]:
    """Return the established ideal optical reference for one detector."""

    if assembly == "cherenkov":
        return {"coupling_index": 1.46, "acceptance_deg": 90.0,
                "quantum_efficiency": 1.0, "collection_efficiency": 1.0,
                "timing_sigma_ns": 0.0}
    if assembly == "scintillator":
        return {"coupling_index": 1.49, "acceptance_deg": 90.0,
                "quantum_efficiency": 1.0, "collection_efficiency": 1.0,
                "timing_sigma_ns": 0.0}
    raise ValueError("assembly must be cherenkov or scintillator")


def apply_ideal_sensor_settings(values: Mapping[str, object], assembly: str) -> Dict[str, object]:
    """Apply ideal response values without touching geometry or terminal state.

    The scintillator readout fields predate the detector-specific sensor
    profile, so both names are updated as a compatibility bridge.  The
    resulting profile remains detector-specific and the ideal-active marker
    is independent for Cherenkov and Scintillator.
    """

    index = "1" if assembly == "cherenkov" else "2" if assembly == "scintillator" else None
    if index is None:
        raise ValueError("assembly must be cherenkov or scintillator")
    ideal = ideal_sensor_reference(assembly)
    updates = {
        f"sensor_{index}_index": ideal["coupling_index"],
        f"sensor_{index}_acceptance_deg": ideal["acceptance_deg"],
        f"sensor_{index}_qe": ideal["quantum_efficiency"],
        f"sensor_{index}_collection": ideal["collection_efficiency"],
        f"sensor_{index}_timing_sigma_ns": ideal["timing_sigma_ns"],
        f"sensor_{index}_ideal_active": "1",
    }
    if assembly == "scintillator":
        updates.update({
            "scintillator_readout_index": ideal["coupling_index"],
            "scintillator_readout_acceptance_deg": ideal["acceptance_deg"],
            "scintillator_readout_qe": ideal["quantum_efficiency"],
            "scintillator_readout_collection_efficiency": ideal["collection_efficiency"],
        })
    for key, value in updates.items():
        if key in values:
            _set_value(values, key, value)
    return {key: _value_text(values, key, str(value)) for key, value in updates.items()}


def sensor_profile_from_values(values: Mapping[str, object], assembly_models: Optional[Mapping[str, AssemblyGeometryModel]] = None) -> Dict[str, object]:
    """Serialize the detector sensors separately from active geometry."""
    sensors = []
    components = {"scintillator_1", "scintillator_2", "quartz"}
    for index in (1, 2):
        assembly = "cherenkov" if index == 1 else "scintillator"
        coupling = normalize_sensor_coupling(_value_text(values, f"sensor_{index}_coupling", "direct"))
        canonical_model = assembly_models.get(assembly) if assembly_models is not None else None
        if canonical_model is not None:
            canonical_scene = complete_local_assembly_scene(canonical_model, values, assembly=assembly)
            coupling = str(canonical_scene.get("sensor_routing", {}).get("mode", "none"))
            if coupling == "none":
                coupling = "direct"
        # Keep old component fields readable, but make the actual component a
        # derived implementation detail rather than a user routing choice.
        component = "quartz" if index == 1 else "scintillator_2"
        def number(suffix, default="0"):
            return float(_value_text(values, f"sensor_{index}_{suffix}", default))
        dimensions = sensor_readout_dimensions_from_values(values, str(index))
        sensor = simulator.DetectorSensor(
            sensor_id=f"sensor_{index}", component=component,
            coupling_mode=coupling, side=_value_text(values, f"sensor_{index}_side", "Y+"),
            position_x_m=number("position_x_mm") * 1e-3, position_y_m=number("position_y_mm") * 1e-3, position_z_m=number("position_z_mm", "20" if index == 1 else "55") * 1e-3,
            shape=_value_text(values, f"sensor_{index}_shape", "rectangular"), width_mm=float(dimensions["width_mm"]), height_mm=float(dimensions["height_mm"]),
            readout_radius_mm=float(dimensions["radius_mm"]), readout_semiaxis_a_mm=float(dimensions["semi_axis_a_mm"]), readout_semiaxis_b_mm=float(dimensions["semi_axis_b_mm"]),
            quantum_efficiency=number("qe", "1"), collection_efficiency=number("collection", "1"), timing_sigma_s=number("timing_sigma_ns") * 1e-9,
        )
        sensor.validate(components)
        payload = sensor.to_dict()
        payload["derived_component"] = derived_sensor_component(assembly, coupling)
        payload["ideal_active"] = _bool_text(values, f"sensor_{index}_ideal_active")
        payload["coupling_index"] = number("index", "1.0")
        payload["acceptance_deg"] = number("acceptance_deg", "90.0")
        # Automatic placement is canonical terminal placement, not a free
        # point.  Keep the serialized readout metadata tied to the same
        # physical transform consumed by all local/global renderers.
        if assembly_models is not None:
            assembly_name = "cherenkov" if index == 1 else "scintillator"
            model = assembly_models.get(assembly_name)
            if model is not None:
                scene = complete_local_assembly_scene(model, values, assembly=assembly_name)
                matching = [item for item in scene["components"] if item.get("kind") == "sensor"]
                payload["canonical_sensor_present"] = bool(matching)
                payload["canonical_invalid_couplings"] = list(scene.get("invalid_couplings", ()))
                if matching:
                    canonical = matching[0]
                    payload["position_mm"] = list(canonical["position_mm"])
                    payload["coupling_position_mm"] = list(canonical["coupling_position_mm"])
                    payload["coupled_terminal"] = canonical.get("coupled_terminal")
                    payload["coupled_target"] = canonical.get("coupled_target")
                    payload["transform"] = dict(canonical.get("transform", {}))
                else:
                    payload["transform"] = None
                    payload["coupling_position_mm"] = None
        payload["automatic_position"] = _bool_text(values, f"sensor_{index}_automatic_position")
        sensors.append(payload)
    return {"format": "quartz-cherenkov-sensors", "version": 1, "sensors": sensors}


def apply_sensor_profile(values: Mapping[str, object], profile: Mapping[str, object]) -> Dict[str, object]:
    """Atomically validate and apply a two-sensor profile."""
    if not isinstance(profile, Mapping) or profile.get("format") != "quartz-cherenkov-sensors" or profile.get("version") != 1:
        raise ValueError("Invalid sensor profile")
    incoming = profile.get("sensors")
    if not isinstance(incoming, Sequence) or len(incoming) != 2:
        raise ValueError("Sensor profile must contain exactly two sensors")
    candidate = dict(values)
    for index, payload in enumerate(incoming, 1):
        sensor = simulator.DetectorSensor.from_dict(payload)
        coupling = normalize_sensor_coupling(payload.get("coupling_mode", "direct"))
        for target, source, default in (
            (f"sensor_{index}_component", "component", "scintillator_1" if index == 1 else "scintillator_2"),
            (f"sensor_{index}_coupling", "coupling_mode", "direct"), (f"sensor_{index}_side", "side", "Y+"),
            (f"sensor_{index}_shape", "shape", "rectangular"),
        ):
            if target.endswith("_coupling"):
                candidate[target] = coupling
            elif target.endswith("_component"):
                candidate[target] = "quartz" if index == 1 else "scintillator_2"
            elif source in payload:
                candidate[target] = payload[source]
        for target, source in ((f"sensor_{index}_position_x_mm", "position_x_mm"), (f"sensor_{index}_position_y_mm", "position_y_mm"), (f"sensor_{index}_position_z_mm", "position_z_mm"), (f"sensor_{index}_width_mm", "width_mm"), (f"sensor_{index}_height_mm", "height_mm"), (f"sensor_{index}_qe", "quantum_efficiency"), (f"sensor_{index}_collection", "collection_efficiency"), (f"sensor_{index}_timing_sigma_ns", "timing_sigma_ns")):
            candidate[target] = payload.get(source, 0.0)
        candidate[f"sensor_{index}_radius_mm"] = payload.get("readout_radius_mm", min(float(payload.get("width_mm", 20.0)), float(payload.get("height_mm", 10.0))) / 2.0)
        axes = payload.get("readout_semiaxes_mm", {})
        candidate[f"sensor_{index}_semi_axis_a_mm"] = axes.get("a", float(payload.get("width_mm", 20.0)) / 2.0) if isinstance(axes, Mapping) else float(payload.get("width_mm", 20.0)) / 2.0
        candidate[f"sensor_{index}_semi_axis_b_mm"] = axes.get("b", float(payload.get("height_mm", 10.0)) / 2.0) if isinstance(axes, Mapping) else float(payload.get("height_mm", 10.0)) / 2.0
        candidate[f"sensor_{index}_automatic_position"] = "1" if bool(payload.get("automatic_position", _bool_text(values, f"sensor_{index}_automatic_position"))) else "0"
        candidate[f"sensor_{index}_ideal_active"] = "1" if bool(payload.get("ideal_active", _bool_text(values, f"sensor_{index}_ideal_active"))) else "0"
    sensor_profile_from_values(candidate)
    old = {key: _value_text(values, key, "") for key in SENSOR_PROFILE_KEYS if key in values}
    try:
        for key in SENSOR_PROFILE_KEYS:
            if key in values:
                value = str(candidate[key])
                target = values[key]
                target.set(value) if hasattr(target, "set") and callable(target.set) else values.__setitem__(key, value)  # type: ignore[index]
    except Exception:
        for key, value in old.items():
            target = values[key]
            target.set(value) if hasattr(target, "set") and callable(target.set) else values.__setitem__(key, value)  # type: ignore[index]
        raise
    return sensor_profile_from_values(values)


def _global_profile_incoming_values(profile: Mapping[str, object]) -> Dict[str, object]:
    """Flatten new two-detector sections while accepting old global profiles."""

    raw = profile.get("values")
    if not isinstance(raw, Mapping):
        raise ValueError(f"{SCINTILLATOR_GLOBAL_PROFILE_FORMAT} profile is missing values")
    incoming = dict(raw)
    placement = profile.get("global_placement", {})
    if isinstance(placement, Mapping):
        for key in ("beam_axis", "detector_1_z_mm", "detector_2_z_mm"):
            if key in placement and key not in incoming:
                incoming[key if key != "beam_axis" else "global_beam_axis"] = placement[key]
        offsets = placement.get("signed_offsets_mm", {})
        if isinstance(offsets, Mapping):
            for detector_name, prefix in (("detector_1", "global"), ("detector_2", "detector2")):
                values = offsets.get(detector_name)
                if isinstance(values, Sequence) and len(values) == 3:
                    for axis, value in zip(("x", "y", "z"), values):
                        incoming.setdefault(f"{prefix}_offset_{axis}_mm", value)
    for key in ("detector_1_z_mm", "detector_2_z_mm"):
        if key not in incoming and key in profile:
            incoming[key] = profile[key]
    two_detector = profile.get("two_detector", {})
    if isinstance(two_detector, Mapping):
        if "detector_separation_mm" in two_detector and "detector_separation_mm" not in incoming:
            incoming["detector_separation_mm"] = two_detector["detector_separation_mm"]
        for section_name, prefix in (("detector_1", "detector1"), ("detector_2", "detector2")):
            section = two_detector.get(section_name, {})
            if not isinstance(section, Mapping):
                continue
            offsets = section.get("signed_offset_mm")
            if isinstance(offsets, (list, tuple)) and len(offsets) == 3:
                target_prefix = "global" if section_name == "detector_1" else prefix
                for axis, value in zip(("x", "y", "z"), offsets):
                    target_key = f"{target_prefix}_offset_{axis}_mm"
                    if target_key not in incoming:
                        incoming[target_key] = value
            mapping = {
                "guide_left_shape": f"{prefix}_guide_left_shape",
                "guide_right_shape": f"{prefix}_guide_right_shape",
                "scintillator_shape": "detector2_scintillator_shape",
                "scintillator_material": "detector2_scintillator_material",
                "scintillator_refractive_index": "detector2_scintillator_refractive_index",
            }
            for source, target in mapping.items():
                if source in section and target not in incoming:
                    incoming[target] = section[source]
    return incoming


def validate_global_geometry_profile(profile: Mapping[str, object]) -> Dict[str, object]:
    base = _validate_component_profile(profile, SCINTILLATOR_GLOBAL_PROFILE_FORMAT, SCINTILLATOR_GLOBAL_PROFILE_KEYS)
    candidate = default_values()
    candidate.update(_global_profile_incoming_values(profile))
    two_detector_beam_monitor_from_values(candidate)
    return base


def _validate_component_profile(profile: Mapping[str, object], fmt: str, keys: Sequence[str]) -> Dict[str, object]:
    if not isinstance(profile, Mapping) or profile.get("format") != fmt or profile.get("version") != SCINTILLATOR_COMPONENT_PROFILE_VERSION:
        raise ValueError(f"Invalid {fmt} profile")
    incoming = profile.get("values")
    if not isinstance(incoming, Mapping) or any(key not in incoming for key in keys):
        raise ValueError(f"{fmt} profile is missing values")
    candidate = default_values()
    candidate.update({key: incoming[key] for key in keys})
    if fmt == SCINTILLATOR_LOCAL_PROFILE_FORMAT:
        scintillator_local_geometry_from_values(candidate)
    elif fmt == SCINTILLATOR_DETECTOR_PROFILE_FORMAT:
        validate_detector_characteristics_values(candidate)
    else:
        scintillator_global_geometry_from_values(candidate)
    return {"format": fmt, "version": SCINTILLATOR_COMPONENT_PROFILE_VERSION, "values": dict(incoming)}


def _apply_component_profile(values: Mapping[str, object], profile: Mapping[str, object], fmt: str, keys: Sequence[str], validator: Callable[[Mapping[str, object]], object]) -> Dict[str, object]:
    if not isinstance(profile, Mapping) or profile.get("format") != fmt or profile.get("version") != SCINTILLATOR_COMPONENT_PROFILE_VERSION:
        raise ValueError(f"Invalid {fmt} profile")
    incoming = profile.get("values")
    if not isinstance(incoming, Mapping) or any(key not in incoming for key in keys):
        raise ValueError(f"{fmt} profile is missing values")
    candidate = dict(values)
    candidate.update({key: incoming[key] for key in keys})
    validator(candidate)
    old = {key: _value_text(values, key, "") for key in keys if key in values}
    try:
        for key in keys:
            if key not in values:
                continue
            target = values[key]
            value = str(incoming[key])
            if hasattr(target, "set") and callable(target.set):
                target.set(value)
            else:
                values[key] = value  # type: ignore[index]
    except Exception:
        for key, value in old.items():
            target = values[key]
            if hasattr(target, "set") and callable(target.set):
                target.set(value)
            else:
                values[key] = value  # type: ignore[index]
        raise
    return {"format": fmt, "version": SCINTILLATOR_COMPONENT_PROFILE_VERSION, "values": {key: _value_text(values, key, "") for key in keys}}


def apply_local_geometry_profile(values: Mapping[str, object], profile: Mapping[str, object]) -> Dict[str, object]:
    return _apply_component_profile(values, profile, SCINTILLATOR_LOCAL_PROFILE_FORMAT, SCINTILLATOR_LOCAL_PROFILE_KEYS, scintillator_local_geometry_from_values)


def apply_detector_characteristics_profile(values: Mapping[str, object], profile: Mapping[str, object]) -> Dict[str, object]:
    return _apply_component_profile(values, profile, SCINTILLATOR_DETECTOR_PROFILE_FORMAT, SCINTILLATOR_DETECTOR_PROFILE_KEYS, validate_detector_characteristics_values)


def apply_global_geometry_profile(values: Mapping[str, object], profile: Mapping[str, object]) -> Dict[str, object]:
    if not isinstance(profile, Mapping) or profile.get("format") != SCINTILLATOR_GLOBAL_PROFILE_FORMAT or profile.get("version") != SCINTILLATOR_COMPONENT_PROFILE_VERSION:
        raise ValueError(f"Invalid {SCINTILLATOR_GLOBAL_PROFILE_FORMAT} profile")
    incoming = _global_profile_incoming_values(profile)
    if any(key not in incoming for key in SCINTILLATOR_GLOBAL_PROFILE_KEYS):
        raise ValueError(f"{SCINTILLATOR_GLOBAL_PROFILE_FORMAT} profile is missing values")
    candidate = dict(values)
    candidate.update(incoming)
    scintillator_global_geometry_from_values(candidate)
    two_detector_beam_monitor_from_values(candidate)
    apply_keys = tuple(dict.fromkeys(SCINTILLATOR_GLOBAL_PROFILE_KEYS + TWO_DETECTOR_GLOBAL_PROFILE_KEYS))
    old = {key: _value_text(values, key, "") for key in apply_keys if key in values}
    try:
        for key in apply_keys:
            if key not in values or key not in incoming:
                continue
            target = values[key]
            value = str(incoming[key])
            target.set(value) if hasattr(target, "set") and callable(target.set) else values.__setitem__(key, value)  # type: ignore[index]
    except Exception:
        for key, value in old.items():
            target = values[key]
            target.set(value) if hasattr(target, "set") and callable(target.set) else values.__setitem__(key, value)  # type: ignore[index]
        raise
    return global_geometry_profile_from_values(values)


def save_component_profile(path: Path | str, profile: Mapping[str, object]) -> None:
    Path(path).expanduser().resolve().write_text(json.dumps(profile, indent=2, sort_keys=True) + "\n", encoding="utf-8")


def load_component_profile(path: Path | str, values: Mapping[str, object], apply_function: Callable[[Mapping[str, object], Mapping[str, object]], Dict[str, object]]) -> Dict[str, object]:
    profile = json.loads(Path(path).expanduser().resolve().read_text(encoding="utf-8"))
    return apply_function(values, profile)


def global_assembly_preview(values: Mapping[str, object], quartz_geometry=None) -> Dict[str, object]:
    quartz = quartz_geometry or simulator.BarGeometry()
    monitor = two_detector_beam_monitor_from_values(values, quartz)
    result = monitor.preview(quartz)
    result["volumes"] = [volume for detector in result["detectors"].values() for volume in detector["volumes"]]
    # Preserve the single-assembly preview keys for existing consumers while
    # adding the explicit two-detector sections above.
    result["gap_from_quartz_mm"] = result["quartz_gap_mm"]
    result["mylar_enabled"] = _bool_text(values, "mylar_enabled")
    result["mylar_reflectivity"] = float(_value_text(values, "mylar_reflectivity_percent", "95")) / 100.0
    guide_present = any(volume.get("name", "").startswith("guide") for volume in result["volumes"])
    result["guide_present"] = guide_present
    result["interfaces"] = (["detector_1 guide/scintillator", "detector_1 scintillator/guide", "detector_2 guide/scintillator", "detector_2 scintillator/guide"] if guide_present else [])
    result["readout_output_face"] = "detector 1/2 active +X faces" if not guide_present else "detector 1/2 left and right guide readout faces"
    result.update({
        "quartz_bounds_mm": (np.asarray(simulator._geometry_bounds(quartz)) * 1e3).tolist() if hasattr(simulator, "_geometry_bounds") else None,
        "global_geometry": asdict(scintillator_global_geometry_from_values(values)),
        "beam_axis": _value_text(values, "global_beam_axis", "Z").upper(),
        "quartz_gap_mm": float(_value_text(values, "scintillator_gap_mm", "0")),
    })
    return result


def local_geometry_scene(values: Mapping[str, object]) -> Dict[str, object]:
    """Return headless-testable local solids, rings, and complete faces."""

    local = scintillator_local_geometry_from_values(values)
    solids = (local.scintillator,) if local.light_guide is None else (local.scintillator, local.light_guide)
    volume_rings = [
        _surface_rings_mm(solid)
        for solid in solids
    ]
    volume_faces = [
        solid_surface_faces_mm(solid)
        for solid in solids
    ]
    points = [point for faces in volume_faces for face in faces for point in face]
    lows = tuple(min(point[index] for point in points) for index in range(3))
    highs = tuple(max(point[index] for point in points) for index in range(3))
    def x_plane_indices(faces, x_value):
        tolerance = max(1e-9, abs(float(x_value)) * 1e-9)
        return {
            index for index, face in enumerate(faces)
            if max(point[0] for point in face) - min(point[0] for point in face) <= tolerance
            and abs(sum(point[0] for point in face) / len(face) - float(x_value)) <= tolerance
        }
    scintillator_high_x = max(point[0] for point in volume_faces[0][1])
    scintillator_coupling = x_plane_indices(volume_faces[0], scintillator_high_x)
    guide_coupling = set()
    guide_readout = set()
    if local.light_guide is not None:
        guide_low_x = min(point[0] for point in volume_faces[1][0])
        guide_high_x = max(point[0] for face in volume_faces[1] for point in face)
        guide_coupling = x_plane_indices(volume_faces[1], guide_low_x)
        guide_readout = x_plane_indices(volume_faces[1], guide_high_x)
    has_guide = local.light_guide is not None
    # Adapt the legacy value-store entry point to the same canonical assembly
    # model used by the v30 editor.  The solids above remain the physical
    # payload; this adapter exists only so old local-preview callers receive
    # the new explicit overlay metadata as well.
    def assembly_section(solid, name: str, material: str) -> AssemblySection:
        mm = lambda value: float(value) * 1e3
        return AssemblySection(
            name=name, shape=solid.shape, material=material, axis=solid.axis,
            length_x_mm=mm(solid.length_x_m), width_y_mm=mm(solid.width_y_m),
            thickness_z_mm=mm(solid.thickness_z_m),
            width_y_start_mm=mm(solid.width_y_start_m), width_y_end_mm=mm(solid.width_y_end_m),
            thickness_z_start_mm=mm(solid.thickness_z_start_m), thickness_z_end_mm=mm(solid.thickness_z_end_m),
            radius_mm=mm(solid.radius_m), length_mm=mm(solid.length_m),
            semi_axis_a_mm=mm(solid.semi_axis_a_m), semi_axis_b_mm=mm(solid.semi_axis_b_m),
            radius_start_mm=mm(solid.radius_start_m), radius_end_mm=mm(solid.radius_end_m),
            center_x_mm=mm(solid.center_x_m), center_y_mm=mm(solid.center_y_m), center_z_mm=mm(solid.center_z_m),
        )
    active_section = assembly_section(local.scintillator, "Scintillator active", "scintillator")
    guide_section = assembly_section(local.light_guide, "Light guide", "light guide") if has_guide else None
    overlay_model = AssemblyGeometryModel(
        "Scintillator", [active_section], [guide_section] if guide_section is not None else [], has_guide,
        _value_text(values, "scintillator_terminal_plus_surface", "Detector/Photosensor coupling"),
        _value_text(values, "scintillator_terminal_minus_surface", "Open"),
        _bool_text(values, "mylar_enabled"),
        float(_value_text(values, "mylar_reflectivity_percent", "95")) / 100.0,
    )
    mylar_overlays = derive_mylar_overlay_surfaces(overlay_model)
    overlay_by_section = {}
    for item in mylar_overlays:
        overlay_by_section.setdefault(int(item["section_index"]), []).append(item)
    return {
        "scintillator_rings_mm": volume_rings[0],
        "guide_rings_mm": volume_rings[1] if has_guide else [],
        "scintillator_faces_mm": volume_faces[0],
        "guide_faces_mm": volume_faces[1] if has_guide else [],
        "bounds_mm": [*lows, *highs],
        "views": ("xy", "xz", "yz", "isometric"),
        "coupling_interface_x_mm": float(local.scintillator.bounds_max_m[0]) * 1e3,
        "readout_face": "light_guide +X" if has_guide else "scintillator +X",
        "mylar_external_surfaces": bool(_bool_text(values, "mylar_enabled")),
        "mylar_enabled": bool(_bool_text(values, "mylar_enabled")),
        "mylar_reflectivity": float(_value_text(values, "mylar_reflectivity_percent", "95")) / 100.0,
        "mylar_overlay_surfaces_mm": [face for item in mylar_overlays for face in item["faces_mm"]],
        "mylar_overlay_components": mylar_overlays,
        "mylar_overlay_style": {"color": MYLAR_DISPLAY_COLOR, "edge_color": MYLAR_DISPLAY_EDGE, "hatch": MYLAR_DISPLAY_HATCH},
        "scintillator_mylar_face_indices": sorted({int(item["source_face_index"]) for item in overlay_by_section.get(0, ())}),
        "guide_mylar_face_indices": sorted({int(item["source_face_index"]) for item in overlay_by_section.get(1, ())}) if has_guide else [],
        "scintillator_coupling_face_indices": sorted(scintillator_coupling),
        "guide_coupling_face_indices": sorted(guide_coupling),
        "guide_readout_face_indices": sorted(guide_readout),
        "interfaces_unwrapped": ("scintillator_light_guide", "light_guide_readout") if has_guide else (),
        "guide_enabled": bool(_bool_text(values, "scintillator_guide_enabled")),
        "guide_present": bool(has_guide),
    }


def local_geometry_preview(values: Mapping[str, object]) -> Dict[str, object]:
    """Compatibility alias for clients that call previews by this name."""

    return local_geometry_scene(values)


def _surface_rings_mm(geometry) -> list:
    if isinstance(geometry, (simulator.SolidGeometry, simulator.CompoundGeometry)):
        return [[tuple(float(value) * 1e3 for value in point) for point in ring] for ring in geometry.surface_rings()]
    hx, hy, hz = geometry.half_x * 1e3, geometry.half_y * 1e3, geometry.half_z * 1e3
    return [[(-hx, -hy, -hz), (hx, -hy, -hz), (hx, hy, -hz), (-hx, hy, -hz)], [(-hx, -hy, hz), (hx, -hy, hz), (hx, hy, hz), (-hx, hy, hz)]]


def solid_surface_faces_mm(geometry, samples: int = 32) -> list:
    """Build complete polygonal faces from a solid's two terminal rings.

    SolidGeometry.surface_rings remains the transport/diagnostic API and is
    intentionally unchanged. This renderer-only helper adds two cap polygons
    plus a closed strip of lateral faces between corresponding ring samples.
    """

    if not isinstance(geometry, simulator.SolidGeometry):
        raise TypeError("solid_surface_faces_mm expects a SolidGeometry")
    if samples < 4:
        raise ValueError("surface face samples must be at least 4")
    rings = [
        [tuple(float(value) * 1e3 for value in point) for point in ring]
        for ring in geometry.surface_rings(samples=int(samples))
    ]
    if len(rings) != 2 or any(len(ring) < 3 for ring in rings):
        raise ValueError("solid surface rings must contain two polygonal caps")
    first, second = rings
    faces = [first, list(reversed(second))]
    count = min(len(first), len(second))
    for index in range(count):
        next_index = (index + 1) % count
        faces.append([first[index], first[next_index], second[next_index], second[index]])
    return faces


def geometry_surface_faces_mm(geometry, samples: int = 32) -> list:
    """Return complete renderer faces for a bar, solid, or compound geometry."""

    if isinstance(geometry, simulator.SolidGeometry):
        return solid_surface_faces_mm(geometry, samples=samples)
    if isinstance(geometry, simulator.CompoundGeometry):
        faces = []
        for section in geometry.sections:
            faces.extend(solid_surface_faces_mm(section, samples=samples))
        return faces
    if isinstance(geometry, simulator.BarGeometry):
        hx, hy, hz = geometry.half_x * 1e3, geometry.half_y * 1e3, geometry.half_z * 1e3
        return [
            [(-hx, -hy, -hz), (-hx, -hy, hz), (-hx, hy, hz), (-hx, hy, -hz)],
            [(hx, -hy, -hz), (hx, hy, -hz), (hx, hy, hz), (hx, -hy, hz)],
            [(-hx, -hy, -hz), (hx, -hy, -hz), (hx, -hy, hz), (-hx, -hy, hz)],
            [(-hx, hy, -hz), (-hx, hy, hz), (hx, hy, hz), (hx, hy, -hz)],
            [(-hx, -hy, -hz), (-hx, hy, -hz), (hx, hy, -hz), (hx, -hy, -hz)],
            [(-hx, -hy, hz), (hx, -hy, hz), (hx, hy, hz), (-hx, hy, hz)],
        ]
    raise TypeError("unsupported geometry for renderer faces")


def _x_plane_face_indices(faces, x_value) -> set[int]:
    tolerance = max(1e-9, abs(float(x_value)) * 1e-9)
    return {
        index for index, face in enumerate(faces)
        if max(point[0] for point in face) - min(point[0] for point in face) <= tolerance
        and abs(sum(point[0] for point in face) / len(face) - float(x_value)) <= tolerance
    }


def _axis_plane_face_indices(faces, axis_index: int, value: float) -> set[int]:
    tolerance = max(1e-9, abs(float(value)) * 1e-9)
    return {
        index for index, face in enumerate(faces)
        if max(point[axis_index] for point in face) - min(point[axis_index] for point in face) <= tolerance
        and abs(sum(point[axis_index] for point in face) / len(face) - float(value)) <= tolerance
    }


def _terminal_face_mm_axis(solid: simulator.SolidGeometry, axis_index: int, high: bool, faces=None) -> list:
    faces = list(faces if faces is not None else solid_surface_faces_mm(solid))
    value = (solid.bounds_max_m[axis_index] if high else solid.bounds_min_m[axis_index]) * 1e3
    indices = _axis_plane_face_indices(faces, axis_index, value)
    if indices:
        return faces[min(indices)]
    # A curved/discretized side may not contain an exact plane; select the
    # polygon nearest the requested terminal coordinate.
    return min(faces, key=lambda face: abs(sum(point[axis_index] for point in face) / len(face) - value))


def _global_assembly_render_point(point, scene: Mapping[str, object]) -> tuple:
    """Map actual coordinates to the readable global viewer coordinates."""

    x, y, z = (float(value) for value in point)
    if scene.get("display_mode") == "schematic/not-to-scale":
        q_high = float(scene["quartz_high_x_mm"])
        d1_low = float(scene.get("detector1_low_x_mm", scene["assembly_low_x_mm"]))
        d1_high = float(scene.get("detector1_high_x_mm", scene.get("assembly_high_x_mm", d1_low)))
        d2_low = float(scene.get("detector2_low_x_mm", d1_high))
        gap = max(float(scene.get("gap_mm", 0.0)), 0.0)
        separation = max(float(scene.get("detector_separation_mm", scene.get("assembly_separation_mm", 0.0))), 0.0)
        display_gap = float(scene["display_gap_mm"])
        display_separation = float(scene.get("display_detector_separation_mm", display_gap))
        gap_delta = max(gap - display_gap, 0.0)
        separation_delta = max(separation - display_separation, 0.0)
        if x > q_high and gap > display_gap and x < d1_low:
            x = q_high + (x - q_high) * display_gap / gap
        elif x >= d1_low and x < d1_high:
            x -= gap_delta
        elif x >= d1_high and x < d2_low:
            x -= gap_delta
            x = d1_high - gap_delta + (x - (d1_high - gap_delta)) * display_separation / max(separation, 1e-12)
        elif x >= d2_low:
            x -= gap_delta + separation_delta
    return x, y, z


def global_assembly_fit(scene: Mapping[str, object], width: float, height: float) -> Dict[str, object]:
    """Return a camera fit covering solids, overlays, readouts, and beam ends."""

    faces = [face for component in scene.get("render_components", ()) for face in component["faces_mm"]]
    if not faces:
        faces = scene["quartz_faces_mm"] + scene["scintillator_faces_mm"] + scene["guide_faces_mm"]
    points = [
        _global_assembly_render_point(point, scene)
        for face in faces for point in face
    ]
    for component in scene.get("render_components", ()):
        for face in tuple(component.get("interface_faces_mm", ())) + tuple(component.get("readout_faces_mm", ())):
            points.extend(_global_assembly_render_point(point, scene) for point in face)
    points.extend(
        _global_assembly_render_point(point, scene)
        for point in scene.get("display_particle_path_mm", scene.get("particle_path_mm", ()))
    )
    lows = tuple(min(point[index] for point in points) for index in range(3))
    highs = tuple(max(point[index] for point in points) for index in range(3))
    center = tuple((lows[index] + highs[index]) / 2.0 for index in range(3))
    extent = max(max(highs[index] - lows[index] for index in range(3)), 1.0)
    margin = 72.0
    scale = min(max(float(width) - margin, 120.0), max(float(height) - margin, 120.0)) * 0.43 / extent
    return {
        "center_mm": center,
        "bounds_mm": [*lows, *highs],
        "extent_mm": extent,
        "scale_px_per_mm": scale,
        "display_mode": scene.get("display_mode", "true-scale"),
    }


def _global_preview_raw_point(
    point, yaw: float = GLOBAL_PREVIEW_DEFAULT_YAW,
    pitch: float = GLOBAL_PREVIEW_DEFAULT_PITCH,
) -> tuple:
    """Project a display point with the embedded viewer's camera contract."""

    projected = global_assembly_camera_project(
        point, (0.0, 0.0, 0.0), 1.0, 0.0, 0.0, yaw, pitch,
    )
    return projected[0], projected[1]


def global_preview_axis_vectors(
    yaw: float = GLOBAL_PREVIEW_DEFAULT_YAW,
    pitch: float = GLOBAL_PREVIEW_DEFAULT_PITCH,
) -> Dict[str, tuple]:
    """Return projected global XYZ unit vectors for the embedded triad."""

    return global_assembly_orientation_triad(yaw, pitch)


def global_assembly_preview_layout(
    scene: Mapping[str, object], width: float, height: float,
    *, yaw: float = GLOBAL_PREVIEW_DEFAULT_YAW,
    pitch: float = GLOBAL_PREVIEW_DEFAULT_PITCH,
) -> Dict[str, object]:
    """Return a headless-testable, canvas-sized layout for the global preview.

    Rendered component faces, readout/interface faces, and both physical beam
    endpoints participate in the fit.  Labels and the corner gizmo are
    reserved outside the plot bounds, so they cannot cover or shrink the
    assembly.  ``yaw``/``pitch`` are the same camera state used by the popup.
    """

    canvas_width = max(float(width), 420.0)
    canvas_height = max(float(height), 260.0)
    component_specs = tuple(
        (component["label"], component["faces_mm"])
        for component in scene.get("render_components", ())
    ) or (
        ("QUARTZ", scene["quartz_faces_mm"]),
        ("SCINTILLATOR", scene["scintillator_faces_mm"]),
        ("LIGHT GUIDE", scene["guide_faces_mm"]),
    )
    raw_components = []
    all_points = []
    for label, faces in component_specs:
        projected_faces = []
        for face in faces:
            projected = [
                _global_preview_raw_point(
                    _global_assembly_render_point(point, scene), yaw, pitch,
                )
                for point in face
            ]
            projected_faces.append(projected)
            all_points.extend(projected)
        raw_components.append((label, projected_faces))
    # Readout/interface overlays are part of the visible scene even when a
    # component's solid faces do not include their full terminal polygon.
    for component in scene.get("render_components", ()):
        for face in tuple(component.get("interface_faces_mm", ())) + tuple(component.get("readout_faces_mm", ())):
            all_points.extend(
                _global_preview_raw_point(_global_assembly_render_point(point, scene), yaw, pitch)
                for point in face
            )
    display_path = scene.get("display_particle_path_mm", scene.get("particle_path_mm", ()))
    all_points.extend(
        _global_preview_raw_point(_global_assembly_render_point(point, scene), yaw, pitch)
        for point in display_path
    )
    if not all_points:
        raise ValueError("global assembly has no renderable faces")
    min_x = min(point[0] for point in all_points)
    max_x = max(point[0] for point in all_points)
    min_y = min(point[1] for point in all_points)
    max_y = max(point[1] for point in all_points)
    span_x = max(max_x - min_x, 1.0)
    span_y = max(max_y - min_y, 1.0)
    # The title and axis strip are outside the physical drawing area.  The
    # former v55 legend/measurement reservations are deliberately gone: the
    # main preview has a dedicated details pane instead of canvas overlays.
    plot_left, plot_right = 28.0, canvas_width - 28.0
    plot_top, plot_bottom = 46.0, canvas_height - 42.0
    available_width = max(plot_right - plot_left, 240.0)
    available_height = max(plot_bottom - plot_top, 140.0)
    # 0.80 targets a large but bounded 2D drawing with room for annotations.
    scale = min(available_width / span_x, available_height / span_y) * 0.80
    center_x = (min_x + max_x) / 2.0
    center_y = (min_y + max_y) / 2.0
    draw_center_x = (plot_left + plot_right) / 2.0
    draw_center_y = (plot_top + plot_bottom) / 2.0 + 4.0

    def place(point):
        return (
            draw_center_x + (point[0] - center_x) * scale,
            draw_center_y - (point[1] - center_y) * scale,
        )

    components = []
    for label, faces in raw_components:
        placed_faces = [[place(point) for point in face] for face in faces]
        points = [point for face in placed_faces for point in face]
        components.append({
            "label": label,
            "faces": placed_faces,
            "bounds": (
                min(point[0] for point in points), min(point[1] for point in points),
                max(point[0] for point in points), max(point[1] for point in points),
            ),
        })
    geometry_bounds = (
        min(item["bounds"][0] for item in components),
        min(item["bounds"][1] for item in components),
        max(item["bounds"][2] for item in components),
        max(item["bounds"][3] for item in components),
    )
    placed_fit_points = [place(point) for point in all_points]
    content_bounds = (
        min(point[0] for point in placed_fit_points), min(point[1] for point in placed_fit_points),
        max(point[0] for point in placed_fit_points), max(point[1] for point in placed_fit_points),
    )
    return {
        "width": canvas_width,
        "height": canvas_height,
        "camera_yaw": float(yaw),
        "camera_pitch": float(pitch),
        "scale_px_per_projected_mm": scale,
        "bounds_2d": (min_x, min_y, max_x, max_y),
        "geometry_bounds_px": geometry_bounds,
        "content_bounds_px": content_bounds,
        "plot_bounds_px": (plot_left, plot_top, plot_right, plot_bottom),
        "axis_vectors_projected": global_preview_axis_vectors(yaw, pitch),
        "beam_path_projected": tuple(
            _global_preview_raw_point(_global_assembly_render_point(point, scene), yaw, pitch)
            for point in display_path
        ),
        "beam_path_px": tuple(place(_global_preview_raw_point(_global_assembly_render_point(point, scene), yaw, pitch)) for point in display_path),
        "components": components,
        "project": place,
    }


ORTHOGRAPHIC_VIEW_DEFINITIONS = {
    "xz": {"title": "XZ • longitudinal / beam", "axes": (0, 2), "x_label": "X (mm)", "y_label": "Z / beam (mm)", "axis_labels": ("X (mm)", "Z / beam (mm)")},
    "xy": {"title": "XY • transverse", "axes": (0, 1), "x_label": "X (mm)", "y_label": "Y (mm)", "axis_labels": ("X (mm)", "Y (mm)")},
    "yz": {"title": "YZ • terminal / transverse", "axes": (1, 2), "x_label": "Y (mm)", "y_label": "Z / beam (mm)", "axis_labels": ("Y (mm)", "Z / beam (mm)")},
}


def global_assembly_orthographic_layout(scene: Mapping[str, object], width: float, height: float) -> Dict[str, object]:
    """Fit the complete transformed global scene into XZ, XY, and YZ panes."""

    canvas_width = max(float(width), 420.0)
    canvas_height = max(float(height), 260.0)
    # Three columns are comfortable at normal widths.  A 2x2 grid keeps the
    # third pane and its labels visible on narrow notebook layouts.
    columns = 3 if canvas_width >= 720.0 else 2
    rows = 1 if columns == 3 else 2
    top = 30.0
    bottom = 44.0
    outer = 7.0
    gap = 6.0
    pane_width = (canvas_width - 2.0 * outer - gap * (columns - 1)) / columns
    pane_height = (canvas_height - top - bottom - gap * (rows - 1)) / rows
    if pane_width <= 40.0 or pane_height <= 40.0:
        raise ValueError("orthographic preview has no usable viewport")

    components = tuple(scene.get("render_components", ()))
    path = tuple(scene.get("particle_path_mm", scene.get("physical_particle_path_mm", ())))
    markers = tuple(scene.get("transformed_impact_markers_mm", scene.get("impact_markers_mm", ())))
    source_points = [list(point) for component in components for face in component.get("faces_mm", ()) for point in face]
    source_points += [list(point) for component in components for key in ("interface_faces_mm", "readout_faces_mm") for face in component.get(key, ()) for point in face]
    source_points += [list(point) for point in path] + [list(point) for point in markers] + [[0.0, 0.0, 0.0]]
    if not source_points:
        raise ValueError("global assembly has no renderable geometry")

    panes = {}
    for index, (view_name, definition) in enumerate(ORTHOGRAPHIC_VIEW_DEFINITIONS.items()):
        row, column = divmod(index, columns)
        left = outer + column * (pane_width + gap)
        top_edge = top + row * (pane_height + gap)
        right, bottom_edge = left + pane_width, top_edge + pane_height
        axis_a, axis_b = definition["axes"]
        values_a = [float(point[axis_a]) for point in source_points]
        values_b = [float(point[axis_b]) for point in source_points]
        low_a, high_a = min(values_a), max(values_a)
        low_b, high_b = min(values_b), max(values_b)
        span_a, span_b = max(high_a - low_a, 1.0), max(high_b - low_b, 1.0)
        pad_a, pad_b = max(span_a * 0.12, 1.0), max(span_b * 0.12, 1.0)
        low_a, high_a, low_b, high_b = low_a - pad_a, high_a + pad_a, low_b - pad_b, high_b + pad_b
        plot_left, plot_right = left + 34.0, right - 9.0
        plot_top, plot_bottom = top_edge + 27.0, bottom_edge - 24.0
        scale = min((plot_right - plot_left) / max(high_a - low_a, 1e-12), (plot_bottom - plot_top) / max(high_b - low_b, 1e-12))
        used_width = (high_a - low_a) * scale
        used_height = (high_b - low_b) * scale
        origin_x = plot_left + ((plot_right - plot_left) - used_width) / 2.0
        origin_y = plot_top + ((plot_bottom - plot_top) - used_height) / 2.0

        def project(point, _axis_a=axis_a, _axis_b=axis_b, _low_a=low_a, _low_b=low_b, _high_b=high_b, _scale=scale, _ox=origin_x, _oy=origin_y):
            return (_ox + (float(point[_axis_a]) - _low_a) * _scale, _oy + (_high_b - float(point[_axis_b])) * _scale)

        panes[view_name] = {
            "name": view_name, "title": definition["title"], "axes": definition["axes"],
            "axis_labels": definition["axis_labels"], "rect": (left, top_edge, right, bottom_edge),
            "plot_rect": (plot_left, plot_top, plot_right, plot_bottom),
            "bounds_mm": (low_a, low_b, high_a, high_b), "project": project,
            "scale_px_per_mm": scale, "origin_px": project((0.0, 0.0, 0.0)),
        }
    return {
        "width": canvas_width, "height": canvas_height, "columns": columns, "rows": rows,
        "views": panes, "view_order": tuple(ORTHOGRAPHIC_VIEW_DEFINITIONS),
        "legend_rect": (outer, canvas_height - bottom + 4.0, canvas_width - outer, canvas_height - 6.0),
        "path_mm": path, "markers_mm": markers, "origin_mm": (0.0, 0.0, 0.0),
        "responsive": "3 columns" if columns == 3 else "2 x 2 fallback",
    }


def _terminal_face_mm(solid: simulator.SolidGeometry, high: bool) -> list:
    """Return an explicit X-terminal polygon for interface/readout overlays."""

    low, upper = solid.bounds_min_m * 1e3, solid.bounds_max_m * 1e3
    x = float(upper[0] if high else low[0])
    if solid.shape == "tapered_rectangular_prism":
        width = solid.width_y_end_m if high else solid.width_y_start_m
        thickness = solid.thickness_z_end_m if high else solid.thickness_z_start_m
        y0, z0 = solid.center_y_m * 1e3, solid.center_z_m * 1e3
        hy, hz = width * 500.0, thickness * 500.0
        return [(x, y0 - hy, z0 - hz), (x, y0 + hy, z0 - hz), (x, y0 + hy, z0 + hz), (x, y0 - hy, z0 + hz)]
    if solid.shape == "box":
        y0, z0 = solid.center_y_m * 1e3, solid.center_z_m * 1e3
        hy, hz = solid.width_y_m * 500.0, solid.thickness_z_m * 500.0
        return [(x, y0 - hy, z0 - hz), (x, y0 + hy, z0 - hz), (x, y0 + hy, z0 + hz), (x, y0 - hy, z0 + hz)]
    rings = _surface_rings_mm(solid)
    return rings[1 if high else 0]


def _detector_render_components(detector: simulator.DetectorAssembly, detector_number: int) -> tuple:
    colors = {
        1: {"scintillator": "#32a89d", "guide_left": "#f0a35b", "guide_right": "#e98282"},
        2: {"scintillator": "#218c82", "guide_left": "#d17d32", "guide_right": "#c85b5b"},
    }[1 if detector_number == 1 else 2]
    entries = []
    for kind, solid in (("guide_left", detector.guide_left), ("scintillator", detector.scintillator), ("guide_right", detector.guide_right)):
        if solid is None:
            continue
        faces = solid_surface_faces_mm(solid)
        if detector.active_axis == "Z":
            if kind == "scintillator":
                interface_faces = []
                readout_faces = []
                excluded = set()
            elif kind == "guide_left":
                interface_faces = [_terminal_face_mm_axis(solid, 1, True, faces)]
                readout_faces = [_terminal_face_mm_axis(solid, 1, False, faces)]
                excluded = _axis_plane_face_indices(faces, 1, solid.bounds_min_m[1] * 1e3) | _axis_plane_face_indices(faces, 1, solid.bounds_max_m[1] * 1e3)
            else:
                interface_faces = [_terminal_face_mm_axis(solid, 1, False, faces)]
                readout_faces = [_terminal_face_mm_axis(solid, 1, True, faces)]
                excluded = _axis_plane_face_indices(faces, 1, solid.bounds_min_m[1] * 1e3) | _axis_plane_face_indices(faces, 1, solid.bounds_max_m[1] * 1e3)
            mylar_indices = sorted(set(range(len(faces))) - excluded)
        elif kind == "scintillator":
            interface_faces = [_terminal_face_mm(solid, False), _terminal_face_mm(solid, True)]
            readout_faces = []
            mylar_indices = sorted(set(range(len(faces))) - _x_plane_face_indices(faces, solid.bounds_min_m[0] * 1e3) - _x_plane_face_indices(faces, solid.bounds_max_m[0] * 1e3))
        elif kind == "guide_left":
            interface_faces = [_terminal_face_mm(solid, True)]
            readout_faces = [_terminal_face_mm(solid, False)]
            excluded = _x_plane_face_indices(faces, solid.bounds_min_m[0] * 1e3) | _x_plane_face_indices(faces, solid.bounds_max_m[0] * 1e3)
            mylar_indices = sorted(set(range(len(faces))) - excluded)
        else:
            interface_faces = [_terminal_face_mm(solid, False)]
            readout_faces = [_terminal_face_mm(solid, True)]
            excluded = _x_plane_face_indices(faces, solid.bounds_min_m[0] * 1e3) | _x_plane_face_indices(faces, solid.bounds_max_m[0] * 1e3)
            mylar_indices = sorted(set(range(len(faces))) - excluded)
        entries.append({
            "name": f"detector_{detector_number}_{kind}",
            "label": f"D{detector_number} {('SCINTILLATOR' if kind == 'scintillator' else ('LIGHT GUIDE A' if kind == 'guide_left' else 'LIGHT GUIDE B'))}",
            "faces_mm": faces,
            "color": colors[kind],
            "kind": kind,
            "detector": detector_number,
            "mylar_indices": mylar_indices,
            "interface_faces_mm": interface_faces,
            "readout_faces_mm": readout_faces,
        })
    return tuple(entries)


def _sensor_render_components(values: Mapping[str, object], assemblies: Optional[Mapping[int, object]] = None,
                              canonical_scenes: Optional[Mapping[str, Mapping[str, object]]] = None) -> tuple:
    """Create sensor solids from canonical local scenes for global rendering."""
    if canonical_scenes:
        entries = []
        for index, assembly_name in ((1, "cherenkov"), (2, "scintillator")):
            scene = canonical_scenes.get(assembly_name, {})
            sensors = [item for item in scene.get("components", ()) if item.get("kind") == "sensor"]
            for sensor in sensors:
                face = str(sensor.get("face", "terminal"))
                entries.append({
                    "name": f"sensor_{index}_{face}",
                    "label": sensor.get("label", f"SENSOR {index} {face}"),
                    "faces_mm": sensor["solid_faces_mm"], "color": sensor.get("color", "#9b59b6"),
                    "kind": "sensor", "detector": index, "mylar_indices": [],
                    "interface_faces_mm": sensor.get("coupling_face_mm", ()),
                    "readout_faces_mm": sensor.get("window_faces_mm", ()),
                    "solid_faces_mm": sensor["solid_faces_mm"], "window_faces_mm": sensor.get("window_faces_mm", ()),
                    "shape": sensor.get("shape"), "dimensions_mm": sensor.get("dimensions_mm"),
                    "position_mm": list(sensor["position_mm"]),
                    "coupling_position_mm": list(sensor["coupling_position_mm"]),
                    "coupled_terminal": sensor.get("coupled_terminal"),
                    "coupled_target": sensor.get("coupled_target"),
                    "transform": dict(sensor.get("transform", {})),
                    "readout_valid": bool(sensor.get("readout_valid", True)),
                    "readout_status": sensor.get("readout_status", "ok"),
                })
        return tuple(entries)
    sensor_profile_from_values(values)
    entries = []
    for index in (1, 2):
        component = _value_text(values, f"sensor_{index}_component", f"scintillator_{index}")
        coupling = _value_text(values, f"sensor_{index}_coupling", "direct")
        detector = assemblies.get(1 if index == 1 else 2) if assemblies else None
        detector_has_guide = bool(detector is not None and any(solid is not None for solid in (detector.guide_left, detector.guide_right)))
        if str(component).startswith("guide_") and not detector_has_guide:
            component = "quartz" if index == 1 else f"scintillator_{index}"
            coupling = "direct"
        side = _value_text(values, f"sensor_{index}_side", "Y+").upper()
        normal_axis = side[:1] if side[:1] in {"X", "Y", "Z"} else "Y"
        normal_sign = 1.0 if side.endswith("+") else -1.0
        try:
            x = float(_value_text(values, f"sensor_{index}_position_x_mm", "0"))
            y = float(_value_text(values, f"sensor_{index}_position_y_mm", "0"))
            z = float(_value_text(values, f"sensor_{index}_position_z_mm", "20" if index == 1 else "55"))
            width = float(_value_text(values, f"sensor_{index}_width_mm", "20"))
            height = float(_value_text(values, f"sensor_{index}_height_mm", "10"))
        except (TypeError, ValueError) as exc:
            raise ValueError(f"sensor {index} geometry must be numeric") from exc
        if min(width, height) <= 0.0:
            raise ValueError(f"sensor {index} dimensions must be positive")
        shape = _value_text(values, f"sensor_{index}_shape", "rectangular")
        depth = max(1.5, min(4.0, min(width, height) * 0.20))
        automatic_position = _bool_text(values, f"sensor_{index}_automatic_position")
        coupling_position = [x, y, z]
        if automatic_position and assemblies:
            target = assemblies.get(0) if component == "quartz" else (
                detector.guide_right if str(component).startswith("guide_") and detector.guide_right is not None
                else detector.guide_right if index == 2 and _scintillator_guide_configured(values) and detector.guide_right is not None
                else detector.scintillator
            )
            # BarGeometry is the legacy quartz object and exposes half-widths
            # rather than SolidGeometry bounds properties. Keep automatic
            # sensor placement compatible with both scene model families.
            bounds_low_m, bounds_high_m = simulator._geometry_bounds(target)
            bounds_min = bounds_low_m * 1e3
            bounds_max = bounds_high_m * 1e3
            axis_index = {"X": 0, "Y": 1, "Z": 2}[normal_axis]
            coupling_position = [
                float((bounds_min[0] + bounds_max[0]) * 0.5),
                float((bounds_min[1] + bounds_max[1]) * 0.5),
                float((bounds_min[2] + bounds_max[2]) * 0.5),
            ]
            coupling_position[axis_index] = float(bounds_max[axis_index] if normal_sign > 0 else bounds_min[axis_index])
        body_center = list(coupling_position)
        axis_index = {"X": 0, "Y": 1, "Z": 2}[normal_axis]
        if automatic_position and assemblies:
            body_center[axis_index] += normal_sign * depth / 2.0
        else:
            body_center = [x, y, z]
            coupling_position = list(body_center)
            coupling_position[axis_index] -= normal_sign * depth / 2.0
        sensor_solid = sensor_solid_geometry_mm(
            center_mm=body_center, width_mm=width, height_mm=height, depth_mm=depth,
            shape=shape, normal_axis=normal_axis, normal_sign=normal_sign,
        )
        faces = sensor_solid["faces_mm"]
        entries.append({
            "name": f"sensor_{index}", "label": f"SENSOR {index} • {component} • {coupling}",
            "faces_mm": faces, "color": "#9b59b6", "kind": "sensor", "detector": index,
            "mylar_indices": [], "interface_faces_mm": sensor_solid["coupling_face_mm"],
            "readout_faces_mm": sensor_solid["coupling_face_mm"],
            "solid_faces_mm": faces, "window_faces_mm": sensor_solid["coupling_face_mm"],
            "shape": sensor_solid["shape"], "dimensions_mm": sensor_solid["dimensions_mm"],
            "position_mm": list(body_center), "coupling_position_mm": list(coupling_position),
            "normal_axis": normal_axis, "normal_sign": normal_sign,
            "automatic_position": bool(automatic_position), "display_offset_mm": [0.0, 0.0, 0.0],
            "transform": {"position_mm": list(body_center), "coupling_position_mm": list(coupling_position),
                          "normal_axis": normal_axis, "normal_sign": normal_sign},
        })
    return tuple(entries)


def rigid_transform_point_mm(point: Sequence[float], transform: RigidAssemblyTransform, center_mm: Sequence[float] = (0.0, 0.0, 0.0)) -> tuple[float, float, float]:
    """Compatibility wrapper for the canonical model-layer transform."""

    if tuple(float(value) for value in center_mm) == transform.rotation_origin_mm:
        return transform.apply_point_mm(point)
    # Older callers can still provide a rotation center; keep that API while
    # delegating the matrix/convention to the canonical transform object.
    local = RigidAssemblyTransform(
        transform.translation_x_mm, transform.translation_y_mm,
        transform.translation_z_mm, transform.rotation_axis,
        transform.rotation_angle_deg, *tuple(float(value) for value in center_mm),
    )
    return local.apply_point_mm(point)


def _canonical_global_models(values: Mapping[str, object], quartz_geometry=None,
                             assembly_models: Optional[Mapping[str, AssemblyGeometryModel]] = None) -> tuple[dict[str, AssemblyGeometryModel], tuple[str, ...]]:
    """Select only the two live canonical editor models for the global view.

    The old ``detector_1``/``detector_2`` and ``scintillator_1``/``scintillator_2``
    stores are intentionally not accepted as geometry sources.  A headless
    caller without an editor instance receives small canonical fallback models
    so the public helper remains useful, but it still follows this same path.
    """

    selected: dict[str, AssemblyGeometryModel] = {}
    stale: list[str] = []
    if assembly_models is not None:
        for key, value in assembly_models.items():
            normalized = str(key).strip().lower()
            if normalized in {"cherenkov", "scintillator"} and isinstance(value, AssemblyGeometryModel):
                selected[normalized] = value
            else:
                stale.append(str(key))
    if "cherenkov" not in selected:
        quartz = quartz_geometry or simulator.BarGeometry()
        faces = geometry_surface_faces_mm(quartz)
        points = [point for face in faces for point in face]
        low = [min(point[index] for point in points) for index in range(3)]
        high = [max(point[index] for point in points) for index in range(3)]
        selected["cherenkov"] = AssemblyGeometryModel(
            "Cherenkov", [AssemblySection(
                name="Quartz radiator", material="fused silica", axis="X",
                length_x_mm=high[0] - low[0], width_y_mm=high[1] - low[1],
                thickness_z_mm=high[2] - low[2], center_x_mm=(low[0] + high[0]) / 2.0,
                center_y_mm=(low[1] + high[1]) / 2.0, center_z_mm=(low[2] + high[2]) / 2.0,
            )], terminal_plus_surface="Detector/Photosensor coupling",
            impact_x_mm=float(_value_text(values, "entry_x_mm", "0")),
            impact_y_mm=float(_value_text(values, "entry_y_mm", "0")))
        stale.append("missing:cherenkov")
    if "scintillator" not in selected:
        # Build the fallback from the same local-value adapter as transport;
        # the old unit placeholder was the source of the viewer/transport
        # separation and also dropped an enabled light guide.
        local = scintillator_local_geometry_from_values(values)
        def section_from_solid(solid, name):
            return AssemblySection(
                name=name, shape=solid.shape, material=solid.material,
                length_x_mm=solid.length_x_m * 1e3, width_y_mm=solid.width_y_m * 1e3,
                thickness_z_mm=solid.thickness_z_m * 1e3,
                width_y_start_mm=solid.width_y_start_m * 1e3, width_y_end_mm=solid.width_y_end_m * 1e3,
                thickness_z_start_mm=solid.thickness_z_start_m * 1e3, thickness_z_end_mm=solid.thickness_z_end_m * 1e3,
                axis=solid.axis, radius_mm=solid.radius_m * 1e3, length_mm=solid.length_m * 1e3,
                semi_axis_a_mm=solid.semi_axis_a_m * 1e3, semi_axis_b_mm=solid.semi_axis_b_m * 1e3,
                radius_start_mm=solid.radius_start_m * 1e3, radius_end_mm=solid.radius_end_m * 1e3,
                center_x_mm=solid.center_x_m * 1e3, center_y_mm=solid.center_y_m * 1e3,
                center_z_mm=solid.center_z_m * 1e3,
            )
        selected["scintillator"] = AssemblyGeometryModel(
            "Scintillator", [section_from_solid(local.scintillator, "Scintillator active")],
            [section_from_solid(local.light_guide, "Light guide")] if local.light_guide is not None else [],
            bool(local.light_guide is not None),
            _value_text(values, "scintillator_terminal_plus_surface", "Detector/Photosensor coupling"),
            _value_text(values, "scintillator_terminal_minus_surface", "Open"),
            _bool_text(values, "mylar_enabled"),
            float(_value_text(values, "mylar_reflectivity_percent", "95")) / 100.0,
            float(_value_text(values, "scintillator_impact_x_mm", "0")),
            float(_value_text(values, "scintillator_impact_y_mm", "0")),
        )
        stale.append("missing:scintillator")
    return {key: selected[key] for key in ("cherenkov", "scintillator")}, tuple(stale)


def _canonical_transform_from_values(values: Mapping[str, object], assembly_type: str) -> RigidAssemblyTransform:
    """Read one rigid transform per complete assembly, never per component."""

    if assembly_type == "cherenkov":
        prefix = "global"
        rotation_prefix = "assembly_1"
    else:
        prefix = "detector2"
        rotation_prefix = "assembly_2"
    try:
        transform = RigidAssemblyTransform(
            translation_x_mm=float(_value_text(values, f"{prefix}_offset_x_mm", "0")),
            translation_y_mm=float(_value_text(values, f"{prefix}_offset_y_mm", "0")),
            translation_z_mm=float(_value_text(values, f"{prefix}_offset_z_mm", "0")),
            rotation_axis=_value_text(values, f"{rotation_prefix}_rotation_axis", "None"),
            rotation_angle_deg=float(_value_text(values, f"{rotation_prefix}_rotation_angle_deg", "0")),
        )
        transform.validate()
        return transform
    except (TypeError, ValueError) as exc:
        raise ValueError(f"global rigid transform for {assembly_type} is invalid: {exc}") from exc


def _canonical_marker_face(component: Mapping[str, object]) -> list[list[float]]:
    """Turn a local YZ interface marker into a renderer-neutral 3-D face."""

    polygon = component.get("polygon", ())
    x = float(component.get("x_mm", 0.0))
    return [[x, float(point[0]), float(point[1])] for point in polygon if len(point) >= 2]


def _canonical_global_assembly(assembly_type: str, model: AssemblyGeometryModel,
                               values: Mapping[str, object], transform: RigidAssemblyTransform) -> dict[str, object]:
    """Build one complete global assembly from the exact local scene builder."""

    model.validate()
    model.to_simulation_geometry().validate()
    local = complete_local_assembly_scene(model, values, assembly=assembly_type)
    prefix = f"{assembly_type}_"
    section_components: list[dict[str, object]] = []
    sections = list(model.active_sections())
    guide_ids = {id(section) for section in model.valid_guide_sections}
    local_physical_by_name = {
        str(item.get("name")): item
        for item in local.get("physical_components", local.get("components", ()))
        if item.get("kind") in {"geometry", "guide"}
    }
    for section_index, section in enumerate(sections):
        canonical_component = local_physical_by_name.get(section.name, {})
        faces = deepcopy(list(canonical_component.get("faces_mm", canonical_section_faces_mm(section))))
        component_kind = "guide" if id(section) in guide_ids else "geometry"
        item = {
            "id": f"{prefix}section_{section_index}",
            "component_id": f"{prefix}section_{section_index}",
            "name": f"{prefix}{section.name}",
            "local_name": section.name,
            "assembly_type": assembly_type,
            "assembly_prefix": prefix,
            "section_index": section_index,
            "kind": component_kind,
            "material": section.material,
            "shape": section.shape,
            "dimensions_mm": section.to_dict(),
            "label": f"{assembly_type.title()} • {section.name}",
            "faces_mm": faces,
            "color": "#f0a35b" if component_kind == "guide" else ("#66c2a5" if assembly_type == "scintillator" else "#80b1d3"),
            "mylar_indices": list(canonical_component.get("mylar_face_indices", ())) if model.mylar_enabled else [],
            "mylar_overlay_faces_mm": deepcopy(list(canonical_component.get("mylar_overlay_faces_mm", ()))),
            "mylar_source_faces_mm": deepcopy(list(canonical_component.get("mylar_source_faces_mm", ()))),
            "interface_faces_mm": [], "readout_faces_mm": [],
        }
        section_components.append(item)

    # Terminal/mirror/readout markers and the visible sensor are taken from the
    # same local scene, so a terminal or sensor callback updates this assembly
    # without a second simplified geometry path.
    section_names = {section.name for section in sections}
    for local_component in local.get("physical_components", local.get("components", ())):
        kind = local_component.get("kind")
        if kind in {"geometry", "guide"} and local_component.get("name") in section_names:
            continue
        if kind == "sensor":
            faces = local_component.get("solid_faces_mm", local_component.get("faces_mm", ()))
            section_components.append({
                **{key: value for key, value in local_component.items() if key not in {"faces_mm", "solid_faces_mm", "window_faces_mm", "coupling_face_mm"}},
                "id": f"{prefix}sensor_{local_component.get('face', 'terminal')}",
                "component_id": f"{prefix}sensor_{local_component.get('face', 'terminal')}",
                "name": f"{prefix}sensor_{local_component.get('face', 'terminal')}",
                "local_name": local_component.get("name"), "assembly_type": assembly_type,
                "assembly_prefix": prefix, "faces_mm": faces,
                "solid_faces_mm": faces,
                "window_faces_mm": local_component.get("window_faces_mm", ()),
                "interface_faces_mm": local_component.get("coupling_face_mm", ()),
                "readout_faces_mm": local_component.get("window_faces_mm", ()),
                "mylar_indices": [],
            })
        elif kind in {"mirror", "readout_interface"}:
            face = _canonical_marker_face(local_component)
            if not face:
                continue
            section_components.append({
                "id": f"{prefix}{kind}_{local_component.get('face', 'terminal')}",
                "component_id": f"{prefix}{kind}_{local_component.get('face', 'terminal')}",
                "name": f"{prefix}{kind}_{local_component.get('face', 'terminal')}",
                "local_name": local_component.get("name"), "assembly_type": assembly_type,
                "assembly_prefix": prefix, "kind": kind, "faces_mm": [face],
                "color": "#d62728" if kind == "mirror" else "#6a1b9a",
                "mylar_indices": [], "interface_faces_mm": [face],
                "readout_faces_mm": [face] if kind == "readout_interface" else [],
                "label": local_component.get("label", kind.title()),
            })

    # Mylar is a display-only component family derived by the local canonical
    # scene.  It receives exactly the same one rigid transform as its source
    # section, including the explicit display epsilon, and never enters the
    # transport geometry.
    for local_overlay in local.get("overlay_components", ()):
        source_section = int(local_overlay.get("section_index", 0))
        source_face = int(local_overlay.get("source_face_index", 0))
        section_components.append({
            **dict(local_overlay),
            "id": f"{prefix}mylar_{source_section}_{source_face}",
            "component_id": f"{prefix}mylar_{source_section}_{source_face}",
            "name": f"{prefix}{local_overlay.get('name', 'Mylar')}",
            "local_name": local_overlay.get("name", "Mylar"),
            "assembly_type": assembly_type, "assembly_prefix": prefix,
            "kind": "mylar_overlay", "mylar_indices": [],
            "interface_faces_mm": [], "readout_faces_mm": [],
        })

    # This is deliberately the only component transform in the global scene.
    # It covers active/guide faces, sensor body and window, mirror/readout
    # plates, Mylar face metadata, and all component positions together.
    global_components = transform_assembly_components(section_components, transform)

    local_marker = [float(model.impact_x_mm), float(model.impact_y_mm),
                    float(np.mean([section.center_z_mm for section in model.active_sections()]))]
    terminal_states = dict(local.get("terminal_states", {}))
    global_terminals = transform_assembly_terminal_interfaces(local.get("terminal_interfaces", ()), transform)
    return {
        "assembly_type": assembly_type, "name": model.name, "prefix": prefix,
        "model": model, "transform": asdict(transform), "rigid_transform": transform,
        "local_scene": local, "scene": local, "components": global_components,
        "sections": [section.to_dict() for section in model.sections],
        "guide_sections": [section.to_dict() for section in model.valid_guide_sections],
        "guide_enabled": bool(model.guide_enabled), "guide_present": bool(model.has_light_guide),
        "terminal_states": terminal_states,
        "terminal_interfaces": global_terminals,
        "sensor_routing": dict(local.get("sensor_routing", {})),
        "optical_path": list(local.get("optical_path", ())),
        "mylar": dict(local.get("mylar", {})),
        "impact_marker_local_mm": local_marker,
        "impact_marker_global_mm": list(rigid_transform_point_mm(local_marker, transform)),
        "rigid_transform_status": {
            "convention": "global = rotation(local - origin) + origin + translation",
            "origin_mm": list(transform.rotation_origin_mm),
            "component_count": len(global_components),
            "sensor_count": sum(1 for item in global_components if item.get("kind") == "sensor"),
            "mylar_component_count": sum(1 for item in global_components if item.get("kind") == "mylar_overlay"),
            "single_transform_applied": True,
        },
        "component_ids": [str(item["component_id"]) for item in global_components],
    }


def global_assembly_scene(values: Mapping[str, object], quartz_geometry=None,
                          assembly_models: Optional[Mapping[str, AssemblyGeometryModel]] = None) -> Dict[str, object]:
    """Build exactly one current Cherenkov and one current Scintillator assembly.

    This is the sole scene source for the global preview and interactive
    viewer.  It deliberately does not call the legacy two-detector geometry
    builder; that builder remains available to the run/transport compatibility
    path only.
    """

    quartz_geometry = quartz_geometry or simulator.BarGeometry()
    models, stale_keys = _canonical_global_models(values, quartz_geometry, assembly_models)
    canonical = {
        assembly_type: _canonical_global_assembly(
            assembly_type, models[assembly_type], values,
            _canonical_transform_from_values(values, assembly_type),
        )
        for assembly_type in ("cherenkov", "scintillator")
    }
    # The simulator's scintillator transport is intentionally assembly-local.
    # Its local photon coordinates are transformed into the global scene with
    # the assembly transform below, exactly once.  The global assembly is the
    # placement authority: no transport-derived or gap-derived second
    # placement is allowed to move the event-view solid.
    path_transforms = {
        assembly_type: dict(canonical[assembly_type]["transform"])
        for assembly_type in ("cherenkov", "scintillator")
    }
    event_geometry_transforms = {
        assembly_type: dict(canonical[assembly_type]["transform"])
        for assembly_type in ("cherenkov", "scintillator")
    }
    scint_model = models["scintillator"]
    # Use the same complete canonical model that the worker adapts to
    # transport.  The value-field geometry is only a fallback editor input;
    # using it here would reintroduce a viewer/physics section mismatch.
    model_sections = tuple(scint_model.active_sections())
    model_active_solids = tuple(section.to_simulation_geometry() for section in model_sections)
    active_low = np.min(np.asarray([solid.bounds_min_m for solid in model_active_solids]), axis=0)
    active_high = np.max(np.asarray([solid.bounds_max_m for solid in model_active_solids]), axis=0)
    local_center_mm = ((active_low + active_high) / 2.0) * 1e3
    path_transform = canonical["scintillator"]["rigid_transform"]

    # The worker transports Scintillator photons in the local assembly frame.
    # Build the shared incident track in that same frame and transform it with
    # the exact path transform used for photon records.  v44 instead used the
    # quartz-origin global track here, so a translated Scintillator could show
    # a perfectly plausible photon count while the beam visibly missed it.
    impact_local_mm = np.asarray([
        float(scint_model.impact_x_mm),
        float(scint_model.impact_y_mm),
        float((active_low[2] + active_high[2]) * 500.0),
    ], dtype=float)
    track_margin = max(0.0, float(_value_text(values, "track_margin_mm", "20")))
    local_scintillator_track = [
        [float(active_low[0] * 1e3), float(active_low[1] * 1e3), float(active_low[2] * 1e3 - track_margin)],
        [float(active_low[0] * 1e3), float(active_low[1] * 1e3), float(active_high[2] * 1e3 + track_margin)],
    ]
    # Keep the X/Y impact coordinates authoritative; the active volume's
    # local bounds are used only for the two visible endpoints of the +Z beam.
    local_scintillator_track[0][0:2] = [float(impact_local_mm[0]), float(impact_local_mm[1])]
    local_scintillator_track[1][0:2] = [float(impact_local_mm[0]), float(impact_local_mm[1])]
    scintillator_particle_path = [
        list(path_transform.apply_point_mm(point)) for point in local_scintillator_track
    ]
    scintillator_impact_marker = list(path_transform.apply_point_mm(impact_local_mm))
    assemblies = [canonical["cherenkov"], canonical["scintillator"]]
    render_components = [component for assembly in assemblies for component in assembly["components"]]
    ids = [str(component["component_id"]) for component in render_components]
    if len(ids) != len(set(ids)):
        raise ValueError("canonical global component IDs are not unique")
    all_faces = [face for component in render_components for face in component.get("faces_mm", ())]
    if not all_faces:
        raise ValueError("canonical global scene has no physical component faces")
    all_points = [point for face in all_faces for point in face]
    lows = tuple(min(float(point[index]) for point in all_points) for index in range(3))
    highs = tuple(max(float(point[index]) for point in all_points) for index in range(3))
    beam_axis = _value_text(values, "global_beam_axis", "Z").upper()
    if beam_axis not in {"X", "Y", "Z"}:
        raise ValueError("global_beam_axis must be X, Y, or Z")
    axis_index = {"X": 0, "Y": 1, "Z": 2}[beam_axis]
    margin = float(_value_text(values, "track_margin_mm", "20"))
    if not math.isfinite(margin) or margin < 0.0:
        raise ValueError("track_margin_mm must be finite and non-negative")
    markers = [list(assembly["impact_marker_global_mm"]) for assembly in assemblies]
    # The Run inclination is a global +Z polar angle.  The GUI deliberately
    # keeps azimuth at zero: positive theta tilts +Z toward +X in the global
    # X-Z plane.  The selected single-detector mode anchors at that detector;
    # All/Cherenkov anchor at Cherenkov and use the same physical line.
    inclination_deg = float(_value_text(values, "inclination_deg", "0"))
    if not math.isfinite(inclination_deg) or not -89.0 <= inclination_deg <= 89.0:
        raise ValueError("inclination_deg must be in [-89, 89] degrees")
    theta = math.radians(inclination_deg)
    if beam_axis == "Z":
        selected_mode = detector_run_mode_from_values(values)
        anchor_index = 1 if selected_mode == "scintillator" else 0
        anchor = np.asarray(markers[anchor_index], dtype=float)
        global_track_direction = np.asarray([math.sin(theta), 0.0, math.cos(theta)], dtype=float)
        projections = [float(np.dot(np.asarray(point, dtype=float) - anchor, global_track_direction)) for point in all_points]
        entry_t = min(projections) - margin
        exit_t = max(projections) + margin
        entry = list(anchor + entry_t * global_track_direction)
        exit_point = list(anchor + exit_t * global_track_direction)
    else:
        # Retain the older engineering-axis preview for imported placement
        # profiles, but the Run particle convention is always +Z when the
        # normal Z beam is selected.
        cherenkov_marker = np.asarray(markers[0], dtype=float)
        entry = [float(cherenkov_marker[0]), float(cherenkov_marker[1]), lows[2] - margin]
        exit_point = [float(cherenkov_marker[0]), float(cherenkov_marker[1]), highs[2] + margin]
        if beam_axis == "X":
            entry = [lows[0] - margin, float(cherenkov_marker[1]), float(cherenkov_marker[2])]
            exit_point = [highs[0] + margin, float(cherenkov_marker[1]), entry[2]]
        elif beam_axis == "Y":
            entry = [float(cherenkov_marker[0]), lows[1] - margin, float(cherenkov_marker[2])]
            exit_point = [float(cherenkov_marker[0]), highs[1] + margin, entry[2]]
    particle_path = [entry, exit_point]
    path_start = np.asarray(entry, dtype=float)
    path_vector = np.asarray(exit_point, dtype=float) - path_start
    path_length_sq = float(np.dot(path_vector, path_vector))
    alignment_tolerance_mm = 1e-7
    beam_alignment = {}
    if path_length_sq <= 0.0:
        raise ValueError("canonical shared particle path has zero length")
    global_direction = path_vector / math.sqrt(path_length_sq)
    canonical_axis_vector = np.zeros(3, dtype=float)
    canonical_axis_vector[axis_index] = 1.0

    def inverse_point_mm(point, transform):
        origin = np.asarray(transform.rotation_origin_mm, dtype=float)
        translation = np.asarray((transform.translation_x_mm, transform.translation_y_mm, transform.translation_z_mm), dtype=float)
        return transform.rotation_matrix().T.dot(np.asarray(point, dtype=float) - origin - translation) + origin

    def inverse_vector(vector, transform):
        return transform.rotation_matrix().T.dot(np.asarray(vector, dtype=float))

    for assembly_type, marker in zip(("cherenkov", "scintillator"), markers):
        point = np.asarray(marker, dtype=float)
        fraction = float(np.dot(point - path_start, path_vector) / path_length_sq)
        projected = path_start + fraction * path_vector
        distance = float(np.linalg.norm(point - projected))
        transform = canonical[assembly_type]["rigid_transform"]
        local_impact_mm = inverse_point_mm(point, transform)
        active_solids = tuple(section.to_simulation_geometry() for section in models[assembly_type].active_sections())
        local_impact_m = local_impact_mm * 1e-3
        active_volume_contains = any(solid.contains(local_impact_m, tolerance=1e-10) for solid in active_solids)
        local_start_m = inverse_point_mm(path_start, transform) * 1e-3
        local_direction = inverse_vector(global_direction, transform)
        direction_norm = float(np.linalg.norm(local_direction))
        if direction_norm > 0.0:
            local_direction = local_direction / direction_norm
        segment_intersections = []
        segment_length_m = math.sqrt(path_length_sq) * 1e-3
        for section, solid in zip(models[assembly_type].active_sections(), active_solids):
            try:
                entry_t, exit_t = solid.ray_interval(local_start_m, local_direction)
            except (RuntimeError, ValueError):
                continue
            if exit_t >= -1e-10 and entry_t <= segment_length_m + 1e-10:
                segment_intersections.append({"section": section.name, "entry_m": float(entry_t), "exit_m": float(exit_t)})
        line_intersects_active = bool(segment_intersections)
        beam_alignment[assembly_type] = {
            "impact_point_mm": list(point),
            "track_interpolation_point_mm": list(projected),
            "track_fraction": fraction,
            "distance_mm": distance,
            "inside_track": bool(-alignment_tolerance_mm <= fraction <= 1.0 + alignment_tolerance_mm),
            "configured_impact_on_track": bool(distance <= alignment_tolerance_mm),
            "active_volume_contains": bool(active_volume_contains),
            "line_intersects_active_volume": bool(line_intersects_active),
            "local_impact_point_mm": [float(value) for value in local_impact_mm],
            "active_sections": [str(section.name) for section in models[assembly_type].active_sections()],
            "segment_intersections": segment_intersections,
            # A tilted track is allowed to cross the configured active solid
            # at a different local x/y than its saved nominal marker.  The
            # marker residual remains visible, while transport follows the
            # actual common line rather than silently moving the assembly.
            "aligned": bool(
                -alignment_tolerance_mm <= fraction <= 1.0 + alignment_tolerance_mm
                and line_intersects_active
                and (distance <= alignment_tolerance_mm or abs(inclination_deg) > alignment_tolerance_mm)
            ),
            "tolerance_mm": alignment_tolerance_mm,
        }
    alignment_mismatches = [
        f"{assembly_type.title()} impact is not on the shared {beam_axis} track / active-volume check failed: residual {item['distance_mm']:.9g} mm; "
        f"inside_track={item['inside_track']}; segment_intersects={item['line_intersects_active_volume']}"
        for assembly_type, item in beam_alignment.items() if not item["aligned"]
    ]
    cherenkov_faces = [face for component in assemblies[0]["components"] if component["kind"] == "geometry" for face in component["faces_mm"]]
    scintillator_faces = [face for component in assemblies[1]["components"] if component["kind"] == "geometry" for face in component["faces_mm"]]
    guide_faces = [face for component in render_components if component["kind"] == "guide" for face in component["faces_mm"]]
    ch_points = [point for face in cherenkov_faces for point in face]
    sc_points = [point for face in scintillator_faces for point in face]
    ch_low, ch_high = min(point[0] for point in ch_points), max(point[0] for point in ch_points)
    sc_low, sc_high = min(point[0] for point in sc_points), max(point[0] for point in sc_points)
    def component_face_x(assembly, *, kind=None, face=None, fallback=0.0):
        candidates = []
        for component in assembly["components"]:
            if kind is not None and component.get("kind") != kind:
                continue
            if face is not None and str(component.get("face", "")) != face:
                continue
            candidates.extend(point[0] for polygon in component.get("faces_mm", ()) for point in polygon)
        if candidates:
            return float(max(candidates) if face == "plus_x" else min(candidates))
        return float(fallback)

    # These are three independent engineering quantities.  The first is the
    # open/interface gap.  The second uses corresponding +X faces, so it is
    # intentionally not another spelling of the first.  The third is a
    # translation vector and is never folded into either distance.
    cherenkov_plus_terminal_x = float(ch_high)
    scintillator_minus_reference_x = float(sc_low)
    scintillator_plus_readout_x = component_face_x(
        canonical["scintillator"], kind="readout_interface", face="plus_x", fallback=sc_high,
    )
    cher_d1_gap = scintillator_minus_reference_x - cherenkov_plus_terminal_x
    d1_d2_separation = scintillator_plus_readout_x - cherenkov_plus_terminal_x
    gap = cher_d1_gap  # compatibility alias; the UI uses the named definition.
    span = max(ch_high - ch_low, sc_high - sc_low, 1.0)
    display_gap = max(gap, 0.0)
    ch_transform = asdict(canonical["cherenkov"]["rigid_transform"])
    sc_transform = asdict(canonical["scintillator"]["rigid_transform"])
    relative_offset = [
        float(sc_transform[key] - ch_transform[key])
        for key in ("translation_x_mm", "translation_y_mm", "translation_z_mm")
    ]
    relative = {
        "cherenkov_signed_offset_mm": [asdict(canonical["cherenkov"]["rigid_transform"])[key] for key in ("translation_x_mm", "translation_y_mm", "translation_z_mm")],
        "scintillator_signed_offset_mm": [asdict(canonical["scintillator"]["rigid_transform"])[key] for key in ("translation_x_mm", "translation_y_mm", "translation_z_mm")],
        "detector_1_signed_offset_mm": [asdict(canonical["cherenkov"]["rigid_transform"])[key] for key in ("translation_x_mm", "translation_y_mm", "translation_z_mm")],
        "detector_2_signed_offset_mm": [asdict(canonical["scintillator"]["rigid_transform"])[key] for key in ("translation_x_mm", "translation_y_mm", "translation_z_mm")],
        "relative_assembly_offset_mm": relative_offset,
    }
    measurement_definitions = {
        "reference_axis": "+X global assembly-longitudinal axis; positive is toward +X",
        "sign_convention": "signed_mm = coordinate(to_face) - coordinate(from_face) along +X",
        "cher_d1_gap": {
            "label": "Cherenkov +X terminal face -> Scintillator -X entrance/reference face",
            "from_face": "Cherenkov assembly / +X terminal face",
            "to_face": "Scintillator assembly / -X active entrance/reference face",
            "signed_mm": float(cher_d1_gap),
            "formula": "x(Scintillator -X active entrance) - x(Cherenkov +X terminal)",
        },
        "d1_d2_separation": {
            "label": "D1/D2 corresponding +X faces: Cherenkov terminal -> Scintillator readout terminal",
            "from_face": "D1 = Cherenkov assembly / +X terminal face",
            "to_face": "D2 = Scintillator assembly / +X terminal/readout face",
            "signed_mm": float(d1_d2_separation),
            "formula": "x(Scintillator +X terminal/readout) - x(Cherenkov +X terminal)",
        },
        "relative_global_assembly_offset": {
            "label": "Scintillator assembly - Cherenkov assembly translation",
            "vector_mm": relative_offset,
            "formula": "[dx, dy, dz] = T_scintillator - T_cherenkov",
        },
    }
    validation_lines = (
        "Cherenkov assembly: 1",
        "Scintillator assembly: 1",
        "Total detector assemblies: 2",
        f"Cherenkov sections: {len(assemblies[0]['sections'])}; guides: {len(assemblies[0]['guide_sections'])}; sensor(s): {assemblies[0]['sensor_routing'].get('sensor_count', 0)}; terminals: {assemblies[0]['terminal_states']}",
        f"Scintillator sections: {len(assemblies[1]['sections'])}; guides: {len(assemblies[1]['guide_sections'])}; sensor(s): {assemblies[1]['sensor_routing'].get('sensor_count', 0)}; terminals: {assemblies[1]['terminal_states']}",
        "Shared global particle path intersects both participating active volumes. Configured markers remain fixed."
        if not alignment_mismatches else
        "Shared global particle path mismatch: " + "; ".join(alignment_mismatches),
    )
    return {
        "model": "CanonicalTwoAssemblyGlobalScene", "two_detector": True,
        "assemblies": assemblies, "assembly_list": assemblies,
        "assembly_by_type": canonical, "canonical_assemblies": canonical,
        "assembly_count": len(assemblies), "assembly_types": ["cherenkov", "scintillator"],
        "render_components": tuple(render_components), "component_ids": ids,
        "mylar_overlay_components": tuple(item for item in render_components if item.get("kind") == "mylar_overlay"),
        "mylar_enabled": all(bool(value["mylar"].get("enabled")) for value in canonical.values()),
        "mylar_reflectivity": {key: float(value["mylar"].get("reflectivity", 0.0)) for key, value in canonical.items()},
        "canonical_local_scenes": {key: value["local_scene"] for key, value in canonical.items()},
        # This is a global-scene field.  v29 accidentally sourced it from the
        # local scene, which allowed a consumer to draw a stationary sensor
        # while render_components had already moved the assembly.
        "canonical_sensor_transforms": {
            key: [dict(item.get("transform", {})) for item in value["components"] if item.get("kind") == "sensor"]
            for key, value in canonical.items()
        },
        "quartz_faces_mm": cherenkov_faces, "scintillator_faces_mm": scintillator_faces, "guide_faces_mm": guide_faces,
        "terminal_interfaces": [
            {**dict(item), "assembly_type": assembly_type}
            for assembly_type in ("cherenkov", "scintillator")
            for item in canonical[assembly_type]["terminal_interfaces"]
        ],
        "terminal_states": {assembly_type: dict(canonical[assembly_type]["terminal_states"]) for assembly_type in ("cherenkov", "scintillator")},
        "bounds_mm": [*lows, *highs], "physical_bounds_mm": [*lows, *highs],
        "beam_axis": beam_axis, "beam_direction": f"+{beam_axis}",
        "global_beam_axis_vector": [float(value) for value in canonical_axis_vector],
        "particle_direction_global": [float(value) for value in global_direction],
        "global_axes": ("X", "Y", "Z"),
        "orthographic_views": {
            name: {"title": definition["title"], "axes": definition["axes"], "axis_labels": definition["axis_labels"], "projection": "global transformed vertices"}
            for name, definition in ORTHOGRAPHIC_VIEW_DEFINITIONS.items()
        },
        "particle_path_mm": particle_path, "physical_particle_path_mm": particle_path,
        "display_particle_path_mm": particle_path, "particle_path": {"points_mm": particle_path, "beam_axis": beam_axis},
        "transformed_impact_markers_mm": markers, "impact_markers_mm": markers,
        "particle_path_by_assembly": {
            "cherenkov": particle_path,
            "scintillator": scintillator_particle_path,
        },
        "impact_markers_by_assembly": {
            "cherenkov": list(markers[0]),
            "scintillator": scintillator_impact_marker,
        },
        "scintillator_particle_path_mm": scintillator_particle_path,
        "scintillator_impact_marker_mm": scintillator_impact_marker,
        "gap_mm": gap, "quartz_high_x_mm": ch_high, "assembly_low_x_mm": sc_low,
        "configured_scintillator_gap_mm": float(_value_text(values, "scintillator_gap_mm", "0")),
        "assembly_high_x_mm": sc_high, "quartz_span_x_mm": ch_high - ch_low,
        "assembly_span_x_mm": sc_high - ch_low, "assembly_separation_mm": d1_d2_separation,
        "detector_separation_mm": d1_d2_separation, "detector1_low_x_mm": ch_low, "detector1_high_x_mm": ch_high,
        "detector2_low_x_mm": sc_low, "detector2_high_x_mm": sc_high,
        "display_gap_mm": display_gap, "display_detector_separation_mm": abs(float(d1_d2_separation)),
        "display_mode": "true-scale", "guide_present": any(value["guide_present"] for value in canonical.values()),
        "readout_face": "canonical assembly terminal readouts", "quartz_label": "Cherenkov assembly",
        "readout_points_mm": [list(point) for component in render_components for face in component.get("readout_faces_mm", ()) for point in face],
        "readout_point_mm": tuple(np.mean(np.asarray([point for component in render_components for face in component.get("readout_faces_mm", ()) for point in face]), axis=0).tolist()) if any(component.get("readout_faces_mm") for component in render_components) else (sc_high, 0.0, 0.0),
        "gap_marker_y_mm": 0.0, "gap_marker_z_mm": 0.0, "separation_marker_y_mm": 0.0, "separation_marker_z_mm": 0.0,
        "cher_d1_gap_mm": float(cher_d1_gap), "d1_d2_separation_mm": float(d1_d2_separation),
        "cher_d1_from_x_mm": cherenkov_plus_terminal_x, "cher_d1_to_x_mm": scintillator_minus_reference_x,
        "d1_d2_from_x_mm": cherenkov_plus_terminal_x, "d1_d2_to_x_mm": scintillator_plus_readout_x,
        "measurement_definitions": measurement_definitions,
        "relative_global_assembly_offset_mm": relative_offset,
        "relative_placement": relative, "dimensions_local_unchanged": True,
        "rigid_transforms": {key: dict(value["transform"]) for key, value in canonical.items()},
        "rigid_transform_status": {
            key: dict(value["rigid_transform_status"]) for key, value in canonical.items()
        },
        "validation_summary": list(validation_lines), "status_summary": "\n".join(validation_lines),
        "stale_legacy_keys_ignored": list(stale_keys),
        "global_geometry_source": "canonical AssemblyGeometryModel -> build_complete_assembly_scene -> transport adapter",
        "path_transforms": path_transforms,
        "event_geometry_transforms": event_geometry_transforms,
        "beam_alignment": beam_alignment,
        "beam_alignment_tolerance_mm": alignment_tolerance_mm,
        "beam_alignment_mismatches": alignment_mismatches,
        "inclination_deg": float(inclination_deg),
        "azimuth_deg": 0.0,
        "beam_anchor_global_mm": list(anchor) if beam_axis == "Z" else list(markers[0]),
        "particle_axis_diagnostic": (
            f"particle direction: [{global_direction[0]:.9g}, {global_direction[1]:.9g}, {global_direction[2]:.9g}] global; "
            f"inclination {inclination_deg:g}° from +Z, positive toward +X"
        ),
        "assembly_frame_contract": "global assembly transform is authoritative; local event paths are transformed exactly once",
    }


def _event_view_json_safe(value):
    """Convert a canonical scene to a JSON-safe typed payload."""

    if is_dataclass(value):
        return _event_view_json_safe(asdict(value))
    if isinstance(value, Mapping):
        return {str(key): _event_view_json_safe(item) for key, item in value.items()}
    if isinstance(value, (tuple, list)):
        return [_event_view_json_safe(item) for item in value]
    if isinstance(value, np.ndarray):
        return [_event_view_json_safe(item) for item in value.tolist()]
    if isinstance(value, np.generic):
        return value.item()
    if isinstance(value, Path):
        return str(value)
    return value


def canonical_event_view_payload(scene: Mapping[str, object], detector_mode: str = "combined") -> Dict[str, object]:
    """Return the detector-separated canonical scene stored with a run.

    The payload deliberately keeps the complete two-assembly scene and adds
    the participation mode separately.  A single-detector viewer can then
    hide the non-participating assembly without falling back to legacy
    detector geometry or losing the current global transforms.
    """

    mode = detector_run_mode(detector_mode)
    if mode == "combined":
        mismatches = list(scene.get("beam_alignment_mismatches", ()))
        if mismatches:
            raise ValueError(
                "combined shared-beam contract is not representable by one global trajectory: "
                + "; ".join(str(item) for item in mismatches)
            )
    payload = _event_view_json_safe(scene)
    payload["event_detector_mode"] = mode
    payload["event_view_schema"] = "canonical-two-assembly-event-view-v50"
    payload["coordinate_system"] = "global_mm"
    payload["visible_assembly_types"] = (
        ["cherenkov"] if mode == "cherenkov" else
        ["scintillator"] if mode == "scintillator" else
        ["cherenkov", "scintillator"]
    )
    # Persist only participating detector assemblies in the event payload.
    # The global editor can still show both assemblies, but a detector-only
    # run must not leave a non-participating body/path for viewers to draw.
    visible = set(payload["visible_assembly_types"])
    for key in ("assemblies", "assembly_list"):
        if key in payload:
            payload[key] = [item for item in payload[key] if item.get("assembly_type") in visible]
    for key in ("assembly_by_type", "canonical_assemblies"):
        if isinstance(payload.get(key), Mapping):
            payload[key] = {name: item for name, item in payload[key].items() if name in visible}
    for key in ("render_components", "mylar_overlay_components"):
        if isinstance(payload.get(key), list):
            payload[key] = [item for item in payload[key] if item.get("assembly_type") in visible]
    if isinstance(payload.get("component_ids"), list):
        payload["component_ids"] = [
            str(component.get("component_id"))
            for assembly in payload.get("assemblies", ())
            for component in assembly.get("components", ())
            if component.get("component_id") is not None
        ]
    if isinstance(payload.get("impact_markers_mm"), list):
        indices = ["cherenkov", "scintillator"].index(next(iter(visible))) if len(visible) == 1 else None
        if indices is not None:
            payload["impact_markers_mm"] = [payload["impact_markers_mm"][indices]]
            payload["transformed_impact_markers_mm"] = list(payload["impact_markers_mm"])
    # The event viewer's shared beam must describe the participating detector
    # frame.  In Scintillator-only, the physical transport beam is the local
    # -Z path through the active scintillator, transformed once with the same
    # transform used for photon records and event geometry.
    paths_by_assembly = payload.get("particle_path_by_assembly", {})
    markers_by_assembly = payload.get("impact_markers_by_assembly", {})
    if mode == "scintillator" and isinstance(paths_by_assembly, Mapping):
        scint_path = paths_by_assembly.get("scintillator")
        if scint_path:
            payload["particle_path_mm"] = scint_path
            payload["physical_particle_path_mm"] = scint_path
            payload["display_particle_path_mm"] = scint_path
            payload["particle_path"] = {"points_mm": scint_path, "beam_axis": payload.get("beam_axis", "Z")}
    if mode == "scintillator" and isinstance(markers_by_assembly, Mapping):
        marker = markers_by_assembly.get("scintillator")
        if marker:
            payload["impact_markers_mm"] = [marker]
            payload["transformed_impact_markers_mm"] = [marker]
    payload["assembly_count"] = len(payload.get("assemblies", ()))
    payload["assembly_types"] = list(payload.get("visible_assembly_types", ()))
    return payload


def compact_event_view_payload(scene: Mapping[str, object]) -> Dict[str, object]:
    """Return the disk form of the event scene without model/raw-track copies.

    The full canonical scene remains attached to the in-memory worker payload
    for the interactive viewer.  The saved summary only needs transformed
    component faces and the global path/markers, so it deliberately drops
    ``model``, ``local_scene``, duplicate assembly aliases, and any other
    editor bookkeeping that used to make final JSON serialization needlessly
    expensive.
    """

    compact_assemblies = []
    for assembly in scene.get("assemblies", ()):
        if not isinstance(assembly, Mapping):
            continue
        components = []
        for component in assembly.get("components", ()):
            if not isinstance(component, Mapping):
                continue
            components.append({
                key: _json_safe(component.get(key))
                for key in (
                    "component_id", "assembly_type", "kind", "name", "label",
                    "color", "faces_mm", "mylar_indices", "interface_faces_mm",
                    "readout_faces_mm",
                )
                if key in component
            })
        compact_assemblies.append({
            "assembly_type": assembly.get("assembly_type"),
            "name": assembly.get("name"),
            "prefix": assembly.get("prefix"),
            "components": components,
            "transform": _json_safe(assembly.get("transform", {})),
            "terminal_states": _json_safe(assembly.get("terminal_states", {})),
            "sensor_routing": _json_safe(assembly.get("sensor_routing", {})),
            "mylar": _json_safe(assembly.get("mylar", {})),
            "rigid_transform_status": _json_safe(assembly.get("rigid_transform_status", {})),
        })
    keep = (
        "two_detector", "assembly_count", "assembly_types", "bounds_mm",
        "physical_bounds_mm", "beam_axis", "beam_direction", "global_axes",
        "particle_path_mm", "physical_particle_path_mm", "display_particle_path_mm",
        "transformed_impact_markers_mm", "impact_markers_mm", "readout_points_mm",
        "particle_path_by_assembly", "impact_markers_by_assembly",
        "scintillator_particle_path_mm", "scintillator_impact_marker_mm",
        "readout_point_mm", "gap_mm", "display_gap_mm", "display_mode",
        "cher_d1_gap_mm", "d1_d2_separation_mm", "measurement_definitions",
        "relative_global_assembly_offset_mm", "detector_separation_mm",
        "guide_present", "readout_face", "relative_placement", "dimensions_local_unchanged",
        "rigid_transforms", "rigid_transform_status", "validation_summary",
        "status_summary", "stale_legacy_keys_ignored", "global_geometry_source",
        "event_detector_mode", "event_view_schema", "coordinate_system",
        "visible_assembly_types",
    )
    compact = {key: _json_safe(scene[key]) for key in keep if key in scene}
    compact["assemblies"] = compact_assemblies
    compact["assembly_list"] = compact_assemblies
    compact["component_ids"] = [
        str(component.get("component_id"))
        for assembly in compact_assemblies
        for component in assembly.get("components", ())
    ]
    return compact


def _event_record_groups(result, mode: str):
    groups = []
    if mode in {"cherenkov", "combined"}:
        groups.append(("cherenkov", list(getattr(result, "photon_records", ()) or ())))
    scintillator = getattr(result, "scintillator", None)
    if mode in {"scintillator", "combined"} and scintillator is not None:
        groups.append(("scintillator", list(getattr(scintillator, "photon_records", ()) or ())))
    return groups


def _event_path_points_global_mm(record, detector_type: str, canonical_scene: Mapping[str, object]) -> list[list[float]]:
    """Convert a typed simulator path to the scene's global millimetres."""

    points = [list(point) for point in getattr(record, "path_points_m", ()) if len(point) >= 3]
    frame = str(getattr(record, "path_coordinate_frame", "geometry_m"))
    # Cherenkov transport uses the historical ``geometry_m`` frame, but the
    # current run geometry is assembly-local.  Its canonical assembly
    # transform is therefore required just like the explicit scintillator
    # local frame.  Applying it here is the one and only path transform; the
    # scene geometry and readout components are transformed by the same
    # canonical source.
    needs_canonical_transform = frame in {"cherenkov_local_m", "scintillator_local_m"}
    needs_canonical_transform = needs_canonical_transform or (
        detector_type == "cherenkov" and frame == "geometry_m"
    )
    if not needs_canonical_transform:
        return [[float(value) * 1e3 for value in point] for point in points]
    transforms = canonical_scene.get("path_transforms", {})
    transform_payload = transforms.get(detector_type) if isinstance(transforms, Mapping) else None
    if not isinstance(transform_payload, Mapping):
        raise ValueError(f"event scene is missing the {detector_type} local-to-global path transform")
    transform = RigidAssemblyTransform(**{key: transform_payload[key] for key in (
        "translation_x_mm", "translation_y_mm", "translation_z_mm", "rotation_axis",
        "rotation_angle_deg", "rotation_origin_x_mm", "rotation_origin_y_mm", "rotation_origin_z_mm",
    ) if key in transform_payload})
    return [list(transform.apply_point_mm(np.asarray(point, dtype=float) * 1e3)) for point in points]


def _rebase_scene_point(point: Sequence[float], source_payload: Mapping[str, object], target_payload: Mapping[str, object]) -> list[float]:
    """Move one already-global scene point between two rigid placements."""

    def transform(payload):
        return RigidAssemblyTransform(**{key: payload[key] for key in (
            "translation_x_mm", "translation_y_mm", "translation_z_mm", "rotation_axis",
            "rotation_angle_deg", "rotation_origin_x_mm", "rotation_origin_y_mm", "rotation_origin_z_mm",
        ) if key in payload})
    source = transform(source_payload)
    target = transform(target_payload)
    source_origin = np.asarray(source.rotation_origin_mm, dtype=float)
    source_translation = np.asarray([source.translation_x_mm, source.translation_y_mm, source.translation_z_mm], dtype=float)
    global_point = np.asarray(point, dtype=float)
    local_point = source.rotation_matrix().T @ (global_point - source_origin - source_translation) + source_origin
    return list(target.apply_point_mm(local_point))


def _rebase_scene_component(component: Mapping[str, object], canonical_scene: Mapping[str, object]) -> dict:
    """Use the transport-aligned event placement without mutating global editor data."""

    result = dict(component)
    assembly_type = str(component.get("assembly_type", ""))
    targets = canonical_scene.get("event_geometry_transforms", {})
    assemblies = {str(item.get("assembly_type")): item for item in canonical_scene.get("assemblies", ()) if isinstance(item, Mapping)}
    source = assemblies.get(assembly_type, {}).get("transform")
    target = targets.get(assembly_type) if isinstance(targets, Mapping) else None
    if not isinstance(source, Mapping) or not isinstance(target, Mapping):
        return result
    for key in ("faces_mm", "interface_faces_mm", "readout_faces_mm", "solid_faces_mm", "window_faces_mm", "mylar_overlay_faces_mm", "mylar_source_faces_mm"):
        if key in result and result[key]:
            result[key] = [[_rebase_scene_point(point, source, target) for point in face] for face in result[key]]
    return result


def build_interactive_event_scene(result, canonical_scene: Mapping[str, object], *,
                                  detector_mode: Optional[str] = None,
                                  selected_categories=None,
                                  display_photons: int = 40,
                                  track_margin_mm: float = 20.0) -> Dict[str, object]:
    """Build the complete global event-view scene for headless and Tk use."""

    inferred = detector_mode or getattr(result, "detector_mode", "auto")
    if inferred == "auto":
        inferred = "combined" if getattr(result, "scintillator", None) is not None else "cherenkov"
    mode = detector_run_mode(inferred)
    visible = ["cherenkov"] if mode == "cherenkov" else ["scintillator"] if mode == "scintillator" else ["cherenkov", "scintillator"]
    assemblies = {
        str(item.get("assembly_type")): item
        for item in canonical_scene.get("assemblies", canonical_scene.get("assembly_list", ()))
        if isinstance(item, Mapping)
    }
    # The interactive scene is a renderer-neutral view of the same frozen
    # canonical scene used by Global Assembly 3D.  In particular, do not
    # rebase components to a transport/gap placement here: local photon paths
    # are transformed once by ``path_transforms`` and assembly vertices remain
    # byte-for-byte sourced from the global frame (within float tolerance).
    components = [
        deepcopy(dict(component))
        for assembly_type in visible
        for component in assemblies.get(assembly_type, {}).get("components", ())
        if component.get("faces_mm")
    ]
    geometry_lines = []
    geometry_faces = []
    geometry_components = []
    all_points = []
    for component in components:
        faces = [list(face) for face in component.get("faces_mm", ()) if len(face) >= 3]
        geometry_components.append({
            "assembly_type": component.get("assembly_type"),
            "kind": component.get("kind", "geometry"),
            "label": component.get("label", component.get("name", "component")),
            "component_id": component.get("component_id", component.get("name", component.get("label", "component"))),
            "detail": component.get("label", component.get("name", "component")),
            "color": component.get("color", "#80b1d3"),
            "faces_mm": faces,
        })
        for face in faces:
            geometry_faces.append((face, component.get("color", "#80b1d3"), component.get("assembly_type"), component.get("kind", "geometry")))
            geometry_lines.append((face + [face[0]], component.get("assembly_type"), component.get("kind", "geometry")))
            all_points.extend(face)

    track = [list(point) for point in canonical_scene.get("particle_path_mm", canonical_scene.get("physical_particle_path_mm", ())) if len(point) >= 3]
    if len(track) < 2:
        track = []
    if track:
        all_points.extend(track)
    impact_markers = []
    markers_by_assembly = canonical_scene.get("impact_markers_by_assembly", {})
    canonical_markers = canonical_scene.get("impact_markers_mm", canonical_scene.get("transformed_impact_markers_mm", ()))
    for index, assembly_type in enumerate(("cherenkov", "scintillator")):
        if assembly_type not in visible:
            continue
        marker = None
        if isinstance(markers_by_assembly, Mapping) and markers_by_assembly.get(assembly_type) is not None:
            marker = list(markers_by_assembly[assembly_type])
        elif index < len(canonical_markers):
            marker = list(canonical_markers[index])
        elif len(visible) == 1 and canonical_markers:
            marker = list(canonical_markers[0])
        if marker is None:
            continue
        impact_markers.append({"assembly_type": assembly_type, "point_mm": marker})
        all_points.append(marker)

    if selected_categories is None:
        selected_categories = {"photon_paths"}
    photon_lines = []
    per_detector = {}
    colors = {
        "cherenkov": {"detected_plus_x": "#d95f02", "detected_minus_x": "#e7298a", "absorbed": "#b2182b", "max_bounces": "#7570b3", "escaped": "#1f78b4", "losses": "#ff7f00", "unknown": "#6b7280"},
        "scintillator": {"detected_plus_x": "#1b9e77", "detected_minus_x": "#66a61e", "absorbed": "#762a83", "max_bounces": "#386cb0", "escaped": "#4daf4a", "losses": "#a6761d", "unknown": "#008c95"},
    }
    for detector_type, records in _event_record_groups(result, mode):
        detector_lines = []
        filtered_records = filter_photon_records(records, selected_categories)
        pathless_count = 0
        for record in filtered_records[:max(0, int(display_photons))]:
            points = _event_path_points_global_mm(record, detector_type, canonical_scene)
            if len(points) < 2:
                pathless_count += 1
                continue
            category = classify_photon_record(record)
            line = {
                "detector_type": detector_type,
                "assembly_type": detector_type,
                "label": "Cherenkov photon" if detector_type == "cherenkov" else "Scintillator photon",
                "category": category,
                "outcome": str(getattr(record, "outcome", "unknown")),
                "outcome_reason": str(getattr(record, "outcome_reason", "")),
                "detected_marker": bool(getattr(record, "detected", False)),
                "photon_id": int(getattr(record, "photon_id", -1)),
                "endpoint_mm": list(points[-1]),
                "readout_endpoint_distance_mm": (
                    float(getattr(record, "readout_endpoint_distance_m", float("nan"))) * 1e3
                    if math.isfinite(float(getattr(record, "readout_endpoint_distance_m", float("nan")))) else None
                ),
                "points_mm": points,
                "color": colors[detector_type].get(category, colors[detector_type]["unknown"]),
            }
            photon_lines.append(line)
            detector_lines.append(line)
            all_points.extend(points)
        per_detector[detector_type] = {
            "record_count": len(records),
            "selected_record_count": len(filtered_records),
            "displayed_count": len(detector_lines),
            "pathless_count": pathless_count,
            "display_limit": max(0, int(display_photons)),
            "photon_lines": detector_lines,
            "readout_count": sum(1 for component in components if component.get("assembly_type") == detector_type and component.get("kind") in {"sensor", "readout_interface"}),
        }

    readout_points = []
    for component in components:
        if component.get("kind") not in {"sensor", "readout_interface"}:
            continue
        for face in component.get("faces_mm", ()):
            readout_points.extend(face)
    all_points.extend(readout_points)
    if not all_points:
        all_points = [[0.0, 0.0, 0.0]]
    lows = [min(float(point[index]) for point in all_points) for index in range(3)]
    highs = [max(float(point[index]) for point in all_points) for index in range(3)]
    margin = max(0.0, float(track_margin_mm))
    # The canonical path already has its configured margin. Expand by a small
    # display margin so line caps and event endpoints cannot touch the frame.
    lows = [value - margin * 0.10 for value in lows]
    highs = [value + margin * 0.10 for value in highs]
    center = [(lows[index] + highs[index]) / 2.0 for index in range(3)]
    extent = max(max(highs[index] - lows[index] for index in range(3)), 1.0)
    legend = list(compact_scene_legend_lines(geometry_components, detector_mode=mode))
    return {
        "schema": "canonical-two-assembly-event-view-v50",
        "assembly_frame_contract": "global assembly vertices/transforms are reused before camera projection",
        "canonical_assembly_transforms": {
            assembly_type: dict(assemblies[assembly_type].get("transform", {}))
            for assembly_type in visible if assembly_type in assemblies
        },
        "mode": mode,
        "visible_assembly_types": visible,
        "geometry_components": geometry_components,
        "geometry_faces": geometry_faces,
        "geometry_lines": geometry_lines,
        "track_mm": track,
        "impact_markers": impact_markers,
        "readout_points_mm": readout_points,
        "photon_lines": photon_lines,
        "per_detector": per_detector,
        # This is deliberately copied from the lines that survived category
        # filtering, path validation, and the independent per-detector budget.
        # Static PNG verification can therefore distinguish available records
        # from paths actually rendered without re-reading the SimulationResult.
        "static_rendered_count": {
            detector_type: int(details["displayed_count"])
            for detector_type, details in per_detector.items()
        },
        "bounds_mm": [*lows, *highs],
        "center_mm": center,
        "extent_mm": extent,
        "legend": legend,
    }


def complete_local_assembly_scene(model: AssemblyGeometryModel, values: Optional[Mapping[str, object]] = None,
                                  assembly: str = "cherenkov") -> Dict[str, object]:
    """Public headless helper used by both local editors and verification."""

    values = values or {}
    index = "1" if assembly == "cherenkov" else "2"
    dimensions = sensor_readout_dimensions_from_values(values, index)
    return build_complete_assembly_scene(
        model,
        sensor_label="Photosensor" if assembly == "cherenkov" else "Scintillator readout",
        # The old sensor coupling field is accepted only for profile
        # migration.  Scene routing is derived from Geometry Editor terminal
        # state and actual valid guide placement.
        sensor_coupling="derived",
        sensor_automatic_position=_bool_text(values, f"sensor_{index}_automatic_position"),
        sensor_shape=_value_text(values, f"sensor_{index}_shape", "rectangular"),
        sensor_width_mm=float(dimensions["width_mm"]),
        sensor_height_mm=float(dimensions["height_mm"]),
        sensor_radius_mm=float(dimensions["radius_mm"]),
        sensor_semi_axis_a_mm=float(dimensions["semi_axis_a_mm"]),
        sensor_semi_axis_b_mm=float(dimensions["semi_axis_b_mm"]),
    )


def global_assembly_viewer_binding_sequences() -> Tuple[str, ...]:
    """Headless-testable contract for the global interactive viewer."""

    return ("<Configure>", "<ButtonPress-1>", "<B1-Motion>", "<ButtonRelease-1>", "<MouseWheel>", "<Button-4>", "<Button-5>")


def global_assembly_view_preset_names() -> Tuple[str, ...]:
    """Headless-testable names exposed by the global viewer controls."""

    return (
        "Reset view", "Isometric", "View along X", "View along Y", "View along Z",
        "View XY", "View XZ", "View YZ", "Particle path through detectors", "Fit assembly",
    )


def local_geometry_view_binding_sequences() -> Tuple[str, ...]:
    """Headless-testable contract for the local interactive viewer."""

    return global_assembly_viewer_binding_sequences()


def default_values(output_directory: Optional[Path] = None) -> Dict[str, str]:
    """Return GUI defaults as strings, suitable for tkinter StringVars."""

    output = output_directory or (Path(__file__).resolve().parent / "quartz_gui_runs")
    return {
        "output_directory": str(output),
        "events": "10",
        "photons_per_event": "100",
        "sampling_mode": SAMPLING_MODE_LABELS["fixed"],
        "kinetic_mev": "2.0",
        "particle_type": "electron",
        "mass_mev_c2": str(PARTICLE_MASSES_MEV_C2["electron"]),
        "entry_x_mm": "0.0",
        "entry_y_mm": "0.0",
        "inclination_deg": "0.0",
        "azimuth_deg": "0.0",
        "seed": "12345",
        "detector_index": "1.0",
        "detector_acceptance_deg": "30.0",
        "detector_qe": "0.25",
        "collection_efficiency": "0.80",
        "opposite_end_mode": "open",
        "mirror_reflectivity": "1.0",
        "opposite_detector_index": "",
        "opposite_detector_acceptance_deg": "",
        "opposite_detector_qe": "",
        "opposite_collection_efficiency": "",
        "detector_shape": "rectangular",
        "detector_active_width_mm": "20.0",
        "detector_active_height_mm": "10.0",
        "detector_active_radius_mm": "5.0",
        "detector_active_semi_axis_a_mm": "10.0",
        "detector_active_semi_axis_b_mm": "5.0",
        "detector_active_offset_y_mm": "0.0",
        "detector_active_offset_z_mm": "0.0",
        "detector_active_rotation_deg": "0.0",
        "opposite_detector_shape": "rectangular",
        "opposite_detector_active_width_mm": "20.0",
        "opposite_detector_active_height_mm": "10.0",
        "opposite_detector_active_radius_mm": "5.0",
        "opposite_detector_active_semi_axis_a_mm": "10.0",
        "opposite_detector_active_semi_axis_b_mm": "5.0",
        "opposite_detector_active_offset_y_mm": "0.0",
        "opposite_detector_active_offset_z_mm": "0.0",
        "opposite_detector_active_rotation_deg": "0.0",
        "track_margin_mm": "20.0",
        "display_photons": "40",
        # One Run-tab authority for both participating transport models.  The
        # legacy scintillator/Cherenkov profile fields remain persisted for
        # compatibility but never override this value during a run.
        "maximum_photon_bounces": "1000",
        # One authoritative Run-tab participation choice.  The default keeps
        # the typical historical Cherenkov-only run.
        "detector_run_mode": DETECTOR_RUN_MODE_LABELS["cherenkov"],
        # Gamma transport is opt-in.  The approved profile is shipped beside
        # the GUI, while every discrete line remains explicitly unset.
        "gamma_source_profile_path": str(Path(__file__).resolve().parent / "gamma_source_profile_collimated_v1.json"),
        "gamma_line_weight_511_keV": "0.0", "gamma_line_weight_835_keV": "0.0",
        "gamma_line_weight_1173_keV": "0.0", "gamma_line_weight_1275_keV": "0.0",
        "gamma_line_weight_1332_keV": "0.0", "gamma_line_weight_1461_keV": "0.0",
        "gamma_max_interactions": "1000", "gamma_max_light_samples": "10000",
        # Local Cherenkov assembly controls.  The optional guide is appended
        # to the active ordered sections; it never replaces the radiator.
        "cherenkov_guide_enabled": "0",
        "cherenkov_mylar_enabled": "1", "cherenkov_mylar_reflectivity_percent": "95.0",
        "terminal_plus_surface": "Detector/Photosensor coupling",
        "terminal_minus_surface": "Open",
        "scintillator_terminal_plus_surface": "Detector/Photosensor coupling",
        "scintillator_terminal_minus_surface": "Open",
        "scintillator_guide_enabled": "0",
        "scintillator_impact_x_mm": "0.0",
        "scintillator_impact_y_mm": "0.0",
        # Beam-monitor scintillator/light-guide assembly.
        "scintillator_enabled": "0",
        "scintillator_shape": "box", "scintillator_length_x_mm": "20.0", "scintillator_width_y_mm": "20.0", "scintillator_thickness_z_mm": "10.0",
        "scintillator_width_y_start_mm": "20.0", "scintillator_width_y_end_mm": "20.0", "scintillator_thickness_z_start_mm": "10.0", "scintillator_thickness_z_end_mm": "10.0",
        "scintillator_axis": "X", "scintillator_radius_mm": "10.0", "scintillator_length_mm": "20.0", "scintillator_semi_axis_a_mm": "10.0", "scintillator_semi_axis_b_mm": "5.0", "scintillator_radius_start_mm": "10.0", "scintillator_radius_end_mm": "8.0",
        "scintillator_center_x_mm": "60.0", "scintillator_center_y_mm": "0.0", "scintillator_center_z_mm": "0.0", "scintillator_material": "scintillator", "scintillator_refractive_index": "1.58",
        "guide_shape": "box", "guide_length_x_mm": "40.0", "guide_width_y_mm": "20.0", "guide_thickness_z_mm": "10.0",
        "guide_width_y_start_mm": "20.0", "guide_width_y_end_mm": "20.0", "guide_thickness_z_start_mm": "10.0", "guide_thickness_z_end_mm": "10.0",
        "guide_axis": "X", "guide_radius_mm": "10.0", "guide_length_mm": "40.0", "guide_semi_axis_a_mm": "10.0", "guide_semi_axis_b_mm": "5.0", "guide_radius_start_mm": "10.0", "guide_radius_end_mm": "8.0",
        "guide_center_x_mm": "90.0", "guide_center_y_mm": "0.0", "guide_center_z_mm": "0.0", "guide_material": "light_guide", "guide_refractive_index": "1.49",
        "scintillator_gap_mm": "10.0", "scintillator_sampling_mode": SCINTILLATION_SAMPLING_MODE_LABELS["fixed"], "scintillator_photon_count": "100",
        "scintillator_yield_per_MeV": "10000.0", "scintillator_deposited_energy_MeV": "", "scintillator_dEdx_MeV_per_m": "200.0", "scintillator_decay_time_ns": str(SCINTILLATOR_DECAY_TIME_DEFAULT_NS),
        "scintillator_wavelength_min_nm": "350.0", "scintillator_wavelength_max_nm": "550.0", "scintillator_readout_index": "1.0", "scintillator_readout_acceptance_deg": "90.0", "scintillator_readout_qe": "1.0", "scintillator_readout_collection_efficiency": "1.0", "scintillator_absorption_length_mm": "", "scintillator_max_bounces": "100",
        "mylar_enabled": "1", "mylar_reflectivity_percent": "95.0",
        # Global placement is deliberately separate from local solid fields.
        "global_offset_x_mm": "0.0", "global_offset_y_mm": "0.0", "global_offset_z_mm": "0.0", "global_beam_axis": "Z",
        "assembly_1_rotation_axis": "None", "assembly_1_rotation_angle_deg": "0.0",
        "assembly_2_rotation_axis": "None", "assembly_2_rotation_angle_deg": "0.0",
        # Two-detector global beam-monitor placement.  The legacy global_* fields
        # are Detector 1's signed placement for profile compatibility.
        "detector_separation_mm": "25.0",
        "detector2_offset_x_mm": "0.0", "detector2_offset_y_mm": "0.0", "detector2_offset_z_mm": "0.0",
        "detector_1_z_mm": "20.0", "detector_2_z_mm": "55.0",
        "detector1_guide_left_enabled": "0", "detector1_guide_right_enabled": "0",
        "detector2_guide_left_enabled": "0", "detector2_guide_right_enabled": "0",
        "detector1_guide_left_shape": "box", "detector1_guide_right_shape": "box",
        "detector2_scintillator_shape": "box", "detector2_scintillator_material": "scintillator", "detector2_scintillator_refractive_index": "1.58",
        "detector2_guide_left_shape": "box", "detector2_guide_right_shape": "box",
        # Sensors are independent from active solids and may couple directly
        # or through an explicitly present light guide.
        "sensor_1_component": "quartz", "sensor_1_coupling": "direct", "sensor_1_side": "Y+",
        "sensor_1_position_x_mm": "0.0", "sensor_1_position_y_mm": "0.0", "sensor_1_position_z_mm": "20.0", "sensor_1_automatic_position": "1",
        "sensor_1_shape": "rectangular", "sensor_1_width_mm": "20.0", "sensor_1_height_mm": "10.0", "sensor_1_radius_mm": "5.0", "sensor_1_semi_axis_a_mm": "10.0", "sensor_1_semi_axis_b_mm": "5.0", "sensor_1_depth_mm": "2.0", "sensor_1_index": "1.0", "sensor_1_acceptance_deg": "30.0", "sensor_1_qe": "0.25", "sensor_1_collection": "0.80", "sensor_1_timing_sigma_ns": "0.0", "sensor_1_ideal_active": "0",
        "sensor_2_component": "scintillator_2", "sensor_2_coupling": "direct", "sensor_2_side": "Y+",
        "sensor_2_position_x_mm": "0.0", "sensor_2_position_y_mm": "0.0", "sensor_2_position_z_mm": "55.0", "sensor_2_automatic_position": "1",
        "sensor_2_shape": "rectangular", "sensor_2_width_mm": "20.0", "sensor_2_height_mm": "10.0", "sensor_2_radius_mm": "5.0", "sensor_2_semi_axis_a_mm": "10.0", "sensor_2_semi_axis_b_mm": "5.0", "sensor_2_depth_mm": "2.0", "sensor_2_index": "1.0", "sensor_2_acceptance_deg": "90.0", "sensor_2_qe": "1.0", "sensor_2_collection": "1.0", "sensor_2_timing_sigma_ns": "0.0", "sensor_2_ideal_active": "0",
    }


SCINTILLATOR_PROFILE_KEYS = tuple(
    key for key in default_values().keys()
    if key.startswith(("scintillator_", "guide_", "mylar_")) and key != "scintillator_enabled"
)

# This is intentionally a different envelope from the existing transform-only
# ``quartz-cherenkov-global-geometry`` profile and the per-detector v32 files.
# It is the one-file persistence boundary for the two complete local models,
# readouts, and beam-relative placement.
COMPLETE_GLOBAL_ASSEMBLY_PROFILE_FORMAT = "beam-monitor-global-assembly-profile-v1"
COMPLETE_GLOBAL_ASSEMBLY_PROFILE_VERSION = 1
COMPLETE_GLOBAL_ASSEMBLY_PROFILE_UNITS = {
    "length": "mm", "angle": "deg", "time": "ns", "efficiency": "fraction",
}


def _number(raw: Mapping[str, str], key: str, label: str, integer: bool = False, optional: bool = False) -> Optional[float]:
    text = str(raw.get(key, "")).strip()
    if optional and not text:
        return None
    try:
        value = int(text) if integer else float(text)
    except ValueError as exc:
        raise ValueError(f"{label} must be a number") from exc
    if not math.isfinite(value):
        raise ValueError(f"{label} must be finite")
    return value


def _value_text(values: Mapping[str, object], key: str, default: str = "") -> str:
    if key not in values:
        return default
    value = values[key]
    value = value.get() if hasattr(value, "get") and callable(value.get) else value
    return "" if value is None else str(value).strip()


def sensor_readout_dimensions_from_values(values: Mapping[str, object], index: str) -> dict[str, float | str]:
    """Parse the shape-specific readout controls for either detector tab."""

    shape = _value_text(values, f"sensor_{index}_shape", "rectangular")
    def number(suffix: str, default: str) -> float:
        text = _value_text(values, f"sensor_{index}_{suffix}", default)
        try:
            value = float(text)
        except (TypeError, ValueError) as exc:
            raise ValueError(f"sensor_{index}_{suffix} must be numeric") from exc
        if not math.isfinite(value) or value <= 0.0:
            raise ValueError(f"sensor_{index}_{suffix} must be finite and positive")
        return value

    radius_text = _value_text(values, f"sensor_{index}_radius_mm", "")
    axis_a_text = _value_text(values, f"sensor_{index}_semi_axis_a_mm", "")
    axis_b_text = _value_text(values, f"sensor_{index}_semi_axis_b_mm", "")
    shape_key = str(shape).lower()
    # Once explicit shape controls are present, hidden legacy width/height
    # values are derived compatibility data and cannot invalidate the active
    # shape or overwrite it.
    width_value = 20.0 if shape_key in {"circular", "circle"} and radius_text else number("width_mm", "20")
    height_value = 10.0 if shape_key in {"circular", "circle"} and radius_text else number("height_mm", "10")
    if shape_key in {"elliptical", "ellipse"} and axis_a_text and axis_b_text:
        width_value, height_value = 20.0, 10.0
    return resolve_sensor_readout_dimensions_mm(
        shape,
        width_mm=width_value, height_mm=height_value,
        radius_mm=number("radius_mm", "5") if radius_text and str(shape).lower() in {"circular", "circle"} else None,
        semi_axis_a_mm=number("semi_axis_a_mm", "10") if axis_a_text and str(shape).lower() in {"elliptical", "ellipse"} else None,
        semi_axis_b_mm=number("semi_axis_b_mm", "5") if axis_b_text and str(shape).lower() in {"elliptical", "ellipse"} else None,
    )


def detector_active_area_from_values(values: Mapping[str, object], prefix: str = "detector", label: str = "Detector") -> Dict[str, object]:
    """Parse and validate one serializable detector active-area specification."""

    if prefix not in {"detector", "opposite_detector"}:
        raise ValueError("detector area prefix must be detector or opposite_detector")
    shape = _value_text(values, f"{prefix}_shape", "rectangular").lower()
    aliases = {"rectangle": "rectangular", "box": "rectangular", "circle": "circular", "ellipse": "elliptical"}
    shape = aliases.get(shape, shape)
    if shape not in DETECTOR_ACTIVE_AREA_SHAPES:
        raise ValueError(f"{label} active-area shape must be rectangular, circular, or elliptical")

    def number(suffix: str, default: str = "0.0") -> float:
        key = f"{prefix}_{suffix}"
        text = _value_text(values, key, default)
        try:
            result = float(text)
        except (TypeError, ValueError) as exc:
            raise ValueError(f"{label} active-area {suffix} must be numeric") from exc
        if not math.isfinite(result):
            raise ValueError(f"{label} active-area {suffix} must be finite")
        return result

    spec = {
        "shape": shape,
        "active_width_mm": number("active_width_mm", "20.0"),
        "active_height_mm": number("active_height_mm", "10.0"),
        "active_radius_mm": number("active_radius_mm", "5.0"),
        "active_semi_axis_a_mm": number("active_semi_axis_a_mm", "10.0"),
        "active_semi_axis_b_mm": number("active_semi_axis_b_mm", "5.0"),
        "active_offset_y_mm": number("active_offset_y_mm"),
        "active_offset_z_mm": number("active_offset_z_mm"),
        "active_rotation_deg": number("active_rotation_deg"),
    }
    dimensions = (("active_width_mm", "width"), ("active_height_mm", "height"), ("active_radius_mm", "radius"), ("active_semi_axis_a_mm", "semi-axis a"), ("active_semi_axis_b_mm", "semi-axis b"))
    for key, name in dimensions:
        if spec[key] <= 0.0:
            raise ValueError(f"{label} active-area {name} must be positive")
    return spec


def detector_area_preview_geometry(values: Mapping[str, object], prefix: str = "detector", samples: int = 64):
    """Return a deterministic Y/Z-mm boundary for headless and GUI previews."""

    spec = detector_active_area_from_values(values, prefix)
    if spec["shape"] == "rectangular":
        points = [(-spec["active_width_mm"] / 2, -spec["active_height_mm"] / 2), (spec["active_width_mm"] / 2, -spec["active_height_mm"] / 2), (spec["active_width_mm"] / 2, spec["active_height_mm"] / 2), (-spec["active_width_mm"] / 2, spec["active_height_mm"] / 2)]
    else:
        count = max(16, int(samples))
        points = []
        for index in range(count):
            angle = 2.0 * math.pi * index / count
            if spec["shape"] == "circular":
                points.append((spec["active_radius_mm"] * math.cos(angle), spec["active_radius_mm"] * math.sin(angle)))
            else:
                points.append((spec["active_semi_axis_a_mm"] * math.cos(angle), spec["active_semi_axis_b_mm"] * math.sin(angle)))
    angle = math.radians(spec["active_rotation_deg"])
    c, s = math.cos(angle), math.sin(angle)
    return [(spec["active_offset_y_mm"] + c * y - s * z, spec["active_offset_z_mm"] + s * y + c * z) for y, z in points]


def optical_coupling_diagnostic(values: Mapping[str, object], assembly: str = "cherenkov") -> Dict[str, object]:
    """Compare current coupling with the established ideal reference."""

    if assembly == "cherenkov":
        keys = ("sensor_1_index", "sensor_1_acceptance_deg", "sensor_1_qe", "sensor_1_collection")
        ideal = ideal_sensor_reference(assembly)
    elif assembly == "scintillator":
        keys = ("sensor_2_index", "sensor_2_acceptance_deg", "sensor_2_qe", "sensor_2_collection")
        ideal = ideal_sensor_reference(assembly)
    else:
        raise ValueError("assembly must be cherenkov or scintillator")
    real = {
        "coupling_index": float(_value_text(values, keys[0], "1.0")),
        "acceptance_deg": float(_value_text(values, keys[1], "90.0")),
        "quantum_efficiency": float(_value_text(values, keys[2], "1.0")),
        "collection_efficiency": float(_value_text(values, keys[3], "1.0")),
        "timing_sigma_ns": float(_value_text(values, f"sensor_{'1' if assembly == 'cherenkov' else '2'}_timing_sigma_ns", "0.0")),
    }
    ideal_active = _bool_text(values, f"sensor_{'1' if assembly == 'cherenkov' else '2'}_ideal_active")
    return {
        "assembly": assembly, "real": real, "ideal": dict(ideal),
        "non_mutating": True,
        "ideal_active": ideal_active,
        "diagnosis": "ideal sensor settings active" if ideal_active else "ideal reference only; real settings are unchanged",
        "real_is_ideal": all(real.get(key) == ideal.get(key) for key in ("coupling_index", "acceptance_deg", "quantum_efficiency", "collection_efficiency")),
    }


def complete_sensor_readout_from_values(
    values: Mapping[str, object], assembly: str, model: Optional[AssemblyGeometryModel] = None,
) -> Dict[str, object]:
    """Build the sensor/readout block used by the complete v32 JSON envelope."""

    kind = str(assembly).strip().lower()
    if kind not in {"cherenkov", "scintillator"}:
        raise ValueError("assembly must be cherenkov or scintillator")
    index = "1" if kind == "cherenkov" else "2"
    def number(suffix: str, default: str) -> float:
        try:
            value = float(_value_text(values, f"sensor_{index}_{suffix}", default))
        except (TypeError, ValueError) as exc:
            raise ValueError(f"sensor_{index}_{suffix} must be numeric") from exc
        if not math.isfinite(value):
            raise ValueError(f"sensor_{index}_{suffix} must be finite")
        return value

    shape = _value_text(values, f"sensor_{index}_shape", "rectangular")
    dimensions = sensor_readout_dimensions_from_values(values, index)
    width, height, depth = float(dimensions["width_mm"]), float(dimensions["height_mm"]), number("depth_mm", "2")
    present = False
    canonical = None
    if model is not None:
        scene = complete_local_assembly_scene(model, values, assembly=kind)
        canonical_items = [item for item in scene.get("components", ()) if item.get("kind") == "sensor"]
        present = bool(canonical_items)
        canonical = canonical_items[0] if canonical_items else None
    return {
        "sensor_id": f"sensor_{index}",
        "present": present,
        "shape": shape,
        "dimensions_mm": {"width": width, "height": height, "depth": depth},
        "readout_radius_mm": float(dimensions["radius_mm"]),
        "readout_semiaxes_mm": {"a": float(dimensions["semi_axis_a_mm"]), "b": float(dimensions["semi_axis_b_mm"])},
        "quantum_efficiency": number("qe", "1"),
        "acceptance_deg": number("acceptance_deg", "90"),
        "optical_window_coupling_index": number("index", "1"),
        "collection_efficiency": number("collection", "1"),
        "timing_readout": {
            "sigma_ns": number("timing_sigma_ns", "0"),
            "mode": "Gaussian timing/readout sigma",
        },
        "automatic_position": _bool_text(values, f"sensor_{index}_automatic_position"),
        "coupling_mode": str((canonical or {}).get("sensor_coupling", "derived")),
        "coupled_terminal": (canonical or {}).get("coupled_terminal"),
        "coupled_target": (canonical or {}).get("coupled_target"),
        "position_mm": list((canonical or {}).get("position_mm", (
            number("position_x_mm", "0"), number("position_y_mm", "0"),
            number("position_z_mm", "20" if kind == "cherenkov" else "55"),
        ))),
        "detector_specific_response": (
            {
                "generation_model": "Cherenkov generation/transport response",
                "ideal_reference": ideal_sensor_reference(kind),
            } if kind == "cherenkov" else {
                "generation_model": "deposited-energy/scintillation generation/transport response",
                "sampling_mode": _value_text(values, "scintillator_sampling_mode", "fixed"),
                "yield_per_deposited_MeV": number_from_mapping(values, "scintillator_yield_per_MeV", "10000"),
                "decay_time_ns": number_from_mapping(values, "scintillator_decay_time_ns", "2"),
                "absorption_length_mm": _optional_float(values, "scintillator_absorption_length_mm"),
                "max_bounces": int(number_from_mapping(values, "scintillator_max_bounces", "100")),
            }
        ),
    }


def _complete_global_profile_value_keys(assembly: str) -> tuple[str, ...]:
    """Return the explicit value-store fields persisted for one assembly.

    The canonical section lists live in ``AssemblyGeometryModel``.  These
    fields are the response/readout and editor compatibility values that do
    not belong inside that model.  The Run-only bounce override is
    intentionally absent: it is persisted in neither geometry nor assembly
    profile state.
    """

    common = (
        "global_offset_x_mm", "global_offset_y_mm", "global_offset_z_mm",
        "global_beam_axis", "detector_separation_mm", "detector2_offset_x_mm",
        "detector2_offset_y_mm", "detector2_offset_z_mm", "detector_1_z_mm",
        "detector_2_z_mm", "scintillator_gap_mm",
        "assembly_1_rotation_axis", "assembly_1_rotation_angle_deg",
        "assembly_2_rotation_axis", "assembly_2_rotation_angle_deg",
        "detector1_guide_left_enabled", "detector1_guide_right_enabled",
        "detector2_guide_left_enabled", "detector2_guide_right_enabled",
        "detector1_guide_left_shape", "detector1_guide_right_shape",
        "detector2_scintillator_shape", "detector2_scintillator_material",
        "detector2_scintillator_refractive_index", "detector2_guide_left_shape",
        "detector2_guide_right_shape",
    )
    sensor_index = "1" if assembly == "cherenkov" else "2"
    sensor_keys = tuple(key for key in SENSOR_PROFILE_KEYS if key.startswith(f"sensor_{sensor_index}_"))
    if assembly == "cherenkov":
        own = DETECTOR_PROFILE_KEYS + sensor_keys + (
            "cherenkov_guide_enabled", "cherenkov_mylar_enabled",
            "cherenkov_mylar_reflectivity_percent", "terminal_plus_surface",
            "terminal_minus_surface", "entry_x_mm", "entry_y_mm",
        )
    elif assembly == "scintillator":
        own = SCINTILLATOR_LOCAL_PROFILE_KEYS + SCINTILLATOR_DETECTOR_PROFILE_KEYS + sensor_keys + (
            "scintillator_guide_enabled", "scintillator_terminal_plus_surface",
            "scintillator_terminal_minus_surface", "scintillator_impact_x_mm",
            "scintillator_impact_y_mm",
        )
    else:
        raise ValueError("assembly must be cherenkov or scintillator")
    return tuple(dict.fromkeys(common + own))


def _complete_global_profile_transform(values: Mapping[str, object], number: int) -> dict[str, object]:
    prefix = "assembly_1" if number == 1 else "assembly_2"
    offset_prefix = "global" if number == 1 else "detector2"
    transform = {
        "translation_mm": [
            float(_value_text(values, f"{offset_prefix}_offset_{axis}_mm", "0"))
            for axis in ("x", "y", "z")
        ],
        "rotation_axis": _value_text(values, f"{prefix}_rotation_axis", "None"),
        "rotation_angle_deg": float(_value_text(values, f"{prefix}_rotation_angle_deg", "0")),
    }
    if not all(math.isfinite(float(item)) for item in transform["translation_mm"] + [transform["rotation_angle_deg"]]):
        raise ValueError(f"assembly {number} transform contains a non-finite number")
    if transform["rotation_axis"] not in {"None", "X", "Y", "Z"}:
        raise ValueError(f"assembly {number} rotation axis is unsupported")
    RigidAssemblyTransform(
        translation_x_mm=transform["translation_mm"][0],
        translation_y_mm=transform["translation_mm"][1],
        translation_z_mm=transform["translation_mm"][2],
        rotation_axis=transform["rotation_axis"],
        rotation_angle_deg=transform["rotation_angle_deg"],
    ).validate()
    return transform


def complete_global_assembly_profile_from_state(
    values: Mapping[str, object], assembly_models: Mapping[str, AssemblyGeometryModel],
) -> Dict[str, object]:
    """Serialize both complete local models plus readouts and placement.

    This function validates before writing and never mutates either input.
    Lengths are canonical millimetres; the profile also records units so a
    future reader cannot silently interpret an inch/metre file as mm.
    """

    if set(assembly_models) != {"cherenkov", "scintillator"}:
        raise ValueError("complete profile requires exactly Cherenkov and Scintillator assemblies")
    assemblies: Dict[str, object] = {}
    for number, assembly in enumerate(("cherenkov", "scintillator"), 1):
        model = assembly_models[assembly]
        if not isinstance(model, AssemblyGeometryModel):
            raise ValueError(f"{assembly} assembly model is invalid")
        model_copy = deepcopy(model)
        model_copy.validate()
        sensor = complete_sensor_readout_from_values(values, assembly, model_copy)
        model_copy.validate_complete_assembly(
            sensor_shape=str(sensor["shape"]),
            sensor_width_mm=float(sensor["dimensions_mm"]["width"]),
            sensor_height_mm=float(sensor["dimensions_mm"]["height"]),
            sensor_depth_mm=float(sensor["dimensions_mm"]["depth"]),
            sensor_radius_mm=float(sensor["readout_radius_mm"]),
            sensor_semi_axis_a_mm=float(sensor["readout_semiaxes_mm"]["a"]),
            sensor_semi_axis_b_mm=float(sensor["readout_semiaxes_mm"]["b"]),
        )
        keys = _complete_global_profile_value_keys(assembly)
        value_snapshot = {key: _value_text(values, key, "") for key in keys if key in values}
        if assembly == "cherenkov":
            value_snapshot.update({
                "cherenkov_guide_enabled": "1" if model_copy.guide_enabled else "0",
                "cherenkov_mylar_enabled": "1" if model_copy.mylar_enabled else "0",
                "cherenkov_mylar_reflectivity_percent": str(float(model_copy.mylar_reflectivity) * 100.0),
                "terminal_plus_surface": terminal_surface_label(model_copy.terminal_plus_surface),
                "terminal_minus_surface": terminal_surface_label(model_copy.terminal_minus_surface),
                "entry_x_mm": str(float(model_copy.impact_x_mm)),
                "entry_y_mm": str(float(model_copy.impact_y_mm)),
            })
        else:
            value_snapshot.update({
                "scintillator_guide_enabled": "1" if model_copy.guide_enabled else "0",
                "mylar_enabled": "1" if model_copy.mylar_enabled else "0",
                "mylar_reflectivity_percent": str(float(model_copy.mylar_reflectivity) * 100.0),
                "scintillator_terminal_plus_surface": terminal_surface_label(model_copy.terminal_plus_surface),
                "scintillator_terminal_minus_surface": terminal_surface_label(model_copy.terminal_minus_surface),
                "scintillator_impact_x_mm": str(float(model_copy.impact_x_mm)),
                "scintillator_impact_y_mm": str(float(model_copy.impact_y_mm)),
            })
            for prefix, section in (("scintillator", model_copy.sections[0]), ("guide", model_copy.guide_sections[0] if model_copy.guide_sections else None)):
                if section is None:
                    continue
                for field_name in (
                    "shape", "length_x_mm", "width_y_mm", "thickness_z_mm", "width_y_start_mm",
                    "width_y_end_mm", "thickness_z_start_mm", "thickness_z_end_mm", "axis",
                    "radius_mm", "length_mm", "semi_axis_a_mm", "semi_axis_b_mm", "radius_start_mm",
                    "radius_end_mm", "center_x_mm", "center_y_mm", "center_z_mm", "material",
                ):
                    value_snapshot[f"{prefix}_{field_name}"] = str(getattr(section, field_name))
        assemblies[assembly] = {
            "identity": {"id": assembly, "name": model_copy.name},
            "geometry": complete_detector_assembly_geometry_to_dict(
                model_copy, assembly, sensor,
                local_settings={
                    "coordinate_system": "assembly-local millimetres",
                    "global_placement_offsets_excluded": True,
                    "run_particle_configuration_excluded": True,
                },
            ),
            "photosensor_readout": sensor,
            "values": value_snapshot,
            "transform": _complete_global_profile_transform(values, number),
        }
    placement_keys = tuple(dict.fromkeys(
        SCINTILLATOR_GLOBAL_PROFILE_KEYS + TWO_DETECTOR_GLOBAL_PROFILE_KEYS +
        ("assembly_1_rotation_axis", "assembly_1_rotation_angle_deg", "assembly_2_rotation_axis", "assembly_2_rotation_angle_deg")
    ))
    placement = {key: _value_text(values, key, "") for key in placement_keys if key in values}
    # Keep the two duplicate transform blocks and the shared placement block
    # explicitly consistent at write time.
    return {
        "format": COMPLETE_GLOBAL_ASSEMBLY_PROFILE_FORMAT,
        "schema_version": COMPLETE_GLOBAL_ASSEMBLY_PROFILE_VERSION,
        "units": dict(COMPLETE_GLOBAL_ASSEMBLY_PROFILE_UNITS),
        "assembly_order": ["cherenkov", "scintillator"],
        "assemblies": assemblies,
        "beam_relative_placement": {
            "coordinate_system": "global millimetres; rotations in degrees",
            "values": placement,
            "assembly_transforms": {
                "cherenkov": assemblies["cherenkov"]["transform"],
                "scintillator": assemblies["scintillator"]["transform"],
            },
        },
        "run_configuration": {
            "particle_configuration": "not included; remains in Run / Particle track",
            "maximum_photon_bounces": "not included; Run controls is the run-time authority",
        },
    }


def _complete_profile_candidate_values(profile: Mapping[str, object]) -> tuple[Dict[str, object], Dict[str, AssemblyGeometryModel]]:
    """Parse and validate every profile component without side effects."""

    if not isinstance(profile, Mapping) or profile.get("format") != COMPLETE_GLOBAL_ASSEMBLY_PROFILE_FORMAT:
        raise ValueError("Invalid complete Global Assembly profile format")
    if profile.get("schema_version") != COMPLETE_GLOBAL_ASSEMBLY_PROFILE_VERSION:
        raise ValueError("Unsupported complete Global Assembly profile schema")
    if profile.get("units") != COMPLETE_GLOBAL_ASSEMBLY_PROFILE_UNITS:
        raise ValueError("Complete Global Assembly profile units must be length=mm, angle=deg, time=ns, efficiency=fraction")
    if profile.get("assembly_order") != ["cherenkov", "scintillator"]:
        raise ValueError("Complete Global Assembly profile assembly order is invalid")
    blocks = profile.get("assemblies")
    if not isinstance(blocks, Mapping) or set(blocks) != {"cherenkov", "scintillator"}:
        raise ValueError("Complete Global Assembly profile must contain both assemblies")
    placement_block = profile.get("beam_relative_placement")
    if not isinstance(placement_block, Mapping) or placement_block.get("coordinate_system") != "global millimetres; rotations in degrees":
        raise ValueError("Complete Global Assembly profile beam placement block is invalid")
    placement = placement_block.get("values")
    required_placement = tuple(dict.fromkeys(
        SCINTILLATOR_GLOBAL_PROFILE_KEYS + TWO_DETECTOR_GLOBAL_PROFILE_KEYS +
        ("assembly_1_rotation_axis", "assembly_1_rotation_angle_deg", "assembly_2_rotation_axis", "assembly_2_rotation_angle_deg")
    ))
    if not isinstance(placement, Mapping) or any(key not in placement for key in required_placement):
        raise ValueError("Complete Global Assembly profile is missing beam placement values")
    candidate = default_values()
    candidate.update({key: placement[key] for key in required_placement})
    models: Dict[str, AssemblyGeometryModel] = {}
    for number, assembly in enumerate(("cherenkov", "scintillator"), 1):
        block = blocks[assembly]
        if not isinstance(block, Mapping):
            raise ValueError(f"{assembly} profile block must be an object")
        identity = block.get("identity")
        if not isinstance(identity, Mapping) or identity.get("id") != assembly or not str(identity.get("name", "")).strip():
            raise ValueError(f"{assembly} assembly identity is invalid")
        geometry_payload = block.get("geometry")
        readout = block.get("photosensor_readout")
        if not isinstance(geometry_payload, Mapping) or not isinstance(readout, Mapping):
            raise ValueError(f"{assembly} complete geometry/readout is missing")
        if geometry_payload.get("name") != identity.get("name"):
            raise ValueError(f"{assembly} geometry/name identity mismatch")
        model, geometry_sensor = complete_detector_assembly_geometry_from_dict(
            geometry_payload, assembly, name=str(identity["name"]),
        )
        legacy_readout = dict(readout)
        # v61 complete profiles lack the additive shape-specific keys.  The
        # geometry loader has already performed the diameter->radius and
        # extent->semiaxis migration, so compare against that normalized form.
        legacy_readout.setdefault("readout_radius_mm", geometry_sensor.get("readout_radius_mm"))
        legacy_readout.setdefault("readout_semiaxes_mm", geometry_sensor.get("readout_semiaxes_mm"))
        legacy_readout["dimensions_mm"] = geometry_sensor["dimensions_mm"]
        if legacy_readout != dict(geometry_sensor):
            raise ValueError(f"{assembly} photosensor readout disagrees with complete geometry")
        block_values = block.get("values")
        keys = _complete_global_profile_value_keys(assembly)
        if not isinstance(block_values, Mapping):
            raise ValueError(f"{assembly} profile values are missing")
        required_legacy = tuple(key for key in keys if not (key.endswith("_radius_mm") or key.endswith("_semi_axis_a_mm") or key.endswith("_semi_axis_b_mm")))
        if any(key not in block_values for key in required_legacy):
            missing = [key for key in required_legacy if key not in block_values]
            raise ValueError(f"{assembly} profile is missing values: {', '.join(missing)}")
        candidate.update({key: block_values[key] for key in required_legacy})
        for sensor_index in ("1", "2"):
            shape = _value_text(candidate, f"sensor_{sensor_index}_shape", "rectangular").lower()
            width = float(_value_text(candidate, f"sensor_{sensor_index}_width_mm", "20"))
            height = float(_value_text(candidate, f"sensor_{sensor_index}_height_mm", "10"))
            if f"sensor_{sensor_index}_radius_mm" not in block_values:
                candidate[f"sensor_{sensor_index}_radius_mm"] = min(width, height) / 2.0
            if f"sensor_{sensor_index}_semi_axis_a_mm" not in block_values:
                candidate[f"sensor_{sensor_index}_semi_axis_a_mm"] = width / 2.0
            if f"sensor_{sensor_index}_semi_axis_b_mm" not in block_values:
                candidate[f"sensor_{sensor_index}_semi_axis_b_mm"] = height / 2.0
        transform = block.get("transform")
        if not isinstance(transform, Mapping) or not isinstance(transform.get("translation_mm"), (list, tuple)) or len(transform["translation_mm"]) != 3:
            raise ValueError(f"{assembly} transform is invalid")
        transform_prefix = "global" if number == 1 else "detector2"
        rotation_prefix = "assembly_1" if number == 1 else "assembly_2"
        for axis, value in zip(("x", "y", "z"), transform["translation_mm"]):
            key = f"{transform_prefix}_offset_{axis}_mm"
            if key not in block_values:
                candidate[key] = value
            elif not math.isclose(float(block_values[key]), float(value), rel_tol=0.0, abs_tol=1e-12):
                raise ValueError(f"{assembly} transform disagrees with {key}")
        axis_key = f"{rotation_prefix}_rotation_axis"
        angle_key = f"{rotation_prefix}_rotation_angle_deg"
        if axis_key in block_values and str(block_values[axis_key]) != str(transform.get("rotation_axis")):
            raise ValueError(f"{assembly} transform disagrees with {axis_key}")
        if angle_key in block_values and not math.isclose(float(block_values[angle_key]), float(transform.get("rotation_angle_deg")), rel_tol=0.0, abs_tol=1e-12):
            raise ValueError(f"{assembly} transform disagrees with {angle_key}")
        if axis_key not in block_values:
            candidate[axis_key] = transform.get("rotation_axis")
        if angle_key not in block_values:
            candidate[angle_key] = transform.get("rotation_angle_deg")
        # Validate the typed transform and the readout against the detached model.
        _complete_global_profile_transform(candidate, number)
        model.validate_complete_assembly(
            sensor_shape=str(geometry_sensor["shape"]),
            sensor_width_mm=float(geometry_sensor["dimensions_mm"]["width"]),
            sensor_height_mm=float(geometry_sensor["dimensions_mm"]["height"]),
            sensor_depth_mm=float(geometry_sensor["dimensions_mm"]["depth"]),
            sensor_radius_mm=float(geometry_sensor["readout_radius_mm"]),
            sensor_semi_axis_a_mm=float(geometry_sensor["readout_semiaxes_mm"]["a"]),
            sensor_semi_axis_b_mm=float(geometry_sensor["readout_semiaxes_mm"]["b"]),
        )
        models[assembly] = model
        QuartzCherenkovGUI._apply_complete_sensor_candidate(candidate, assembly, geometry_sensor)
    # The common validators are run only on the detached candidate.
    validate_detector_profile({"format": DETECTOR_PROFILE_FORMAT, "version": DETECTOR_PROFILE_VERSION, **{key: candidate[key] for key in DETECTOR_PROFILE_KEYS}})
    validate_detector_characteristics_values(candidate)
    scintillator_global_geometry_from_values(candidate)
    two_detector_beam_monitor_from_values(candidate, simulator.BarGeometry())
    for assembly, model in models.items():
        guide_key = "cherenkov_guide_enabled" if assembly == "cherenkov" else "scintillator_guide_enabled"
        if _bool_text(candidate, guide_key) != bool(model.guide_enabled):
            raise ValueError(f"{assembly} guide_enabled disagrees with its canonical geometry")
    return candidate, models


def validate_complete_global_assembly_profile(profile: Mapping[str, object]) -> Dict[str, object]:
    """Validate a complete profile without mutating live models or values."""

    _complete_profile_candidate_values(profile)
    return dict(profile)


def apply_complete_global_assembly_profile(
    values: Mapping[str, object], assembly_models: Mapping[str, AssemblyGeometryModel], profile: Mapping[str, object],
) -> Dict[str, object]:
    """Atomically commit a fully detached two-assembly profile.

    All parsing, schema, enum, numeric, geometry, sensor, readout, and beam
    placement validation occurs before this function touches either mapping.
    """

    candidate, models = _complete_profile_candidate_values(profile)
    placement_keys = tuple(dict.fromkeys(
        SCINTILLATOR_GLOBAL_PROFILE_KEYS + TWO_DETECTOR_GLOBAL_PROFILE_KEYS +
        ("assembly_1_rotation_axis", "assembly_1_rotation_angle_deg", "assembly_2_rotation_axis", "assembly_2_rotation_angle_deg")
    ))
    controlled_keys = set(placement_keys)
    blocks = profile["assemblies"]
    for assembly in ("cherenkov", "scintillator"):
        controlled_keys.update(blocks[assembly]["values"].keys())
    commit_values = {key: candidate[key] for key in controlled_keys if key in candidate}
    old_values = {key: _value_text(values, key, "") for key in commit_values if key in values}
    old_models = {key: deepcopy(assembly_models[key]) for key in ("cherenkov", "scintillator") if key in assembly_models}
    try:
        for key, value in commit_values.items():
            if key in values:
                _set_value(values, key, value)
        assembly_models["cherenkov"] = models["cherenkov"]
        assembly_models["scintillator"] = models["scintillator"]
    except Exception:
        for key, value in old_values.items():
            _set_value(values, key, value)
        for key, model in old_models.items():
            assembly_models[key] = model
        raise
    return complete_global_assembly_profile_from_state(values, assembly_models)


def save_complete_global_assembly_profile(path: Path | str, values: Mapping[str, object], assembly_models: Mapping[str, AssemblyGeometryModel]) -> Dict[str, object]:
    profile = complete_global_assembly_profile_from_state(values, assembly_models)
    Path(path).expanduser().resolve().write_text(json.dumps(profile, indent=2, sort_keys=True) + "\n", encoding="utf-8")
    return profile


def load_complete_global_assembly_profile(path: Path | str, values: Mapping[str, object], assembly_models: Mapping[str, AssemblyGeometryModel]) -> Dict[str, object]:
    profile = json.loads(Path(path).expanduser().resolve().read_text(encoding="utf-8"))
    return apply_complete_global_assembly_profile(values, assembly_models, profile)


# Short public names make the Global Assembly feature discoverable to
# headless clients while retaining the explicit global-qualified names above.
save_complete_assembly_profile = save_complete_global_assembly_profile
load_complete_assembly_profile = load_complete_global_assembly_profile
apply_complete_assembly_profile = apply_complete_global_assembly_profile
validate_complete_assembly_profile = validate_complete_global_assembly_profile


def number_from_mapping(values: Mapping[str, object], key: str, default: str) -> float:
    """Parse one finite number from a StringVar/mapping store."""

    try:
        result = float(_value_text(values, key, default))
    except (TypeError, ValueError) as exc:
        raise ValueError(f"{key} must be numeric") from exc
    if not math.isfinite(result):
        raise ValueError(f"{key} must be finite")
    return result


def _diagnostic_section_summary(section: AssemblySection, role: str, index: int) -> Dict[str, object]:
    if section.shape == "box":
        dimensions = f"L {section.length_x_mm:g} × W {section.width_y_mm:g} × T {section.thickness_z_mm:g} mm"
    elif section.shape == "tapered_rectangular_prism":
        dimensions = f"L {section.length_x_mm:g}; W {section.width_y_start_mm:g}→{section.width_y_end_mm:g}; T {section.thickness_z_start_mm:g}→{section.thickness_z_end_mm:g} mm"
    elif section.shape == "cylinder":
        dimensions = f"R {section.radius_mm:g}; L {section.length_mm:g}; axis {section.axis}"
    elif section.shape == "elliptical_cylinder":
        dimensions = f"A {section.semi_axis_a_mm:g}; B {section.semi_axis_b_mm:g}; L {section.length_mm:g}; axis {section.axis}"
    else:
        dimensions = f"R {section.radius_start_mm:g}→{section.radius_end_mm:g}; L {section.length_mm:g}; axis {section.axis}"
    return {
        "role": role, "order": index, "name": section.name, "material": section.material,
        "shape": section.shape, "dimensions_mm": dimensions,
        "center_mm": [section.center_x_mm, section.center_y_mm, section.center_z_mm],
    }


def _diagnostic_bounds(section: AssemblySection) -> Tuple[float, float]:
    length = section.length_x_mm if section.shape in {"box", "tapered_rectangular_prism"} else section.length_mm
    return section.center_x_mm - length / 2.0, section.center_x_mm + length / 2.0


def detector_transport_diagnostic(
    values: Mapping[str, object], model: AssemblyGeometryModel, assembly: str,
) -> Dict[str, object]:
    """Return a detector-specific, read-only transport diagnostic report.

    The report deliberately separates local geometry, sensor/readout state,
    optical transport path, warnings/errors, and the physical generation
    model.  It is also the structured source for the GUI's copy/export text.
    """

    kind = str(assembly).strip().lower()
    if kind not in {"cherenkov", "scintillator"}:
        raise ValueError("assembly must be cherenkov or scintillator")
    warnings = []
    errors = []
    try:
        model.validate()
        assembly_valid = True
    except (TypeError, ValueError) as exc:
        assembly_valid = False
        errors.append(f"stale/invalid configuration: {exc}")
    active = [_diagnostic_section_summary(section, "active", index) for index, section in enumerate(model.sections)]
    guides = [_diagnostic_section_summary(section, "guide", index) for index, section in enumerate(model.guide_sections)]
    for label, items in (("active", model.sections), ("guide", model.guide_sections)):
        previous = None
        for section in items:
            try:
                low, high = _diagnostic_bounds(section)
            except (TypeError, ValueError):
                continue
            if previous is not None:
                delta = low - previous
                if delta > 1e-9:
                    warnings.append(f"{label} sections have a {delta:g} mm gap")
                elif delta < -1e-9:
                    warnings.append(f"{label} sections overlap by {-delta:g} mm")
            previous = high
    guide_present = bool(model.has_light_guide)
    guide_valid = bool(not model.guide_enabled or model.valid_guide_sections)
    if model.guide_enabled and not model.guide_sections:
        warnings.append("guide enabled but empty")
    if model.guide_enabled and model.guide_sections and not model.valid_guide_sections:
        warnings.append("guide enabled but has no valid sections")
    if guide_present and not model.guide_attachment_faces():
        warnings.append("light guide has no active-terminal attachment")
    readout = model.readout_status()
    if guide_present and not readout["valid"]:
        warnings.append("guide without a valid readout terminal")
    if any(model.terminal_state(face) == "detector" for face in ("plus_x", "minus_x")) and not readout["valid"]:
        warnings.append("detector terminal with no valid sensor coupling")
    try:
        scene = complete_local_assembly_scene(model, values, assembly=kind)
    except (TypeError, ValueError) as exc:
        scene = {"components": (), "invalid_couplings": (str(exc),), "mylar": {"overlay_surface_count": 0}}
        errors.append(f"invalid optical interface: {exc}")
    sensor_items = [item for item in scene.get("components", ()) if item.get("kind") == "sensor"]
    sensor = complete_sensor_readout_from_values(values, kind, model)
    if readout["valid"] and not sensor_items:
        warnings.append("photosensor readout is selected but no sensor placement was derived")
    if sensor_items:
        actual = sensor_items[0]
        sensor["present"] = True
        sensor["coupled_terminal"] = actual.get("coupled_terminal")
        sensor["coupled_target"] = actual.get("coupled_target")
        sensor["placement_mm"] = actual.get("position_mm")
        if not sensor.get("automatic_position"):
            configured = sensor.get("position_mm")
            if configured and any(abs(float(a) - float(b)) > 1e-6 for a, b in zip(configured, actual.get("position_mm", configured))):
                warnings.append("sensor placement mismatch: canonical terminal placement is active")
    mylar = dict(scene.get("mylar", {}))
    mylar.update({"enabled": bool(model.mylar_enabled), "reflectivity": float(model.mylar_reflectivity), "overlay_surface_count": len(scene.get("overlay_components", ())), "interface_conflicts": []})
    terminal_surfaces = {"+X": terminal_surface_label(model.terminal_plus_surface), "-X": terminal_surface_label(model.terminal_minus_surface)}
    path = " → ".join(readout.get("optical_path", ())) or "No valid photosensor transport path"
    generation = (
        {"model": "Cherenkov generation/transport response", "source": "charged-particle Cherenkov generation in active quartz radiator", "transport": "Cherenkov photon transport and terminal response"}
        if kind == "cherenkov" else
        {"model": "deposited-energy/scintillation generation/transport response", "source": "deposited energy converted by scintillation yield and decay", "transport": "scintillation photon transport and readout response"}
    )
    return {
        "detector_type": kind,
        "assembly_valid": assembly_valid and not errors,
        "local_geometry_state": {"active_section_count": len(active), "active_sections": active, "impact_xy_mm": {"x": model.impact_x_mm, "y": model.impact_y_mm}, "coordinate_system": "assembly-local millimetres"},
        "light_guide": {"enabled": bool(model.guide_enabled), "present": guide_present, "valid": guide_valid, "section_count": len(guides), "sections": guides, "summary": "not present" if not guide_present else f"{len(guides)} ordered section(s)"},
        "optical_transport_path": path,
        "terminal_surfaces": terminal_surfaces,
        "mylar": mylar,
        "sensor_readout_state": sensor,
        "physical_generation_model": generation,
        "warnings_errors": {"warnings": warnings, "errors": errors, "anomalies": list(scene.get("invalid_couplings", ()))},
        "diagnostic_actions": ["Refresh diagnostic", "Compare real vs ideal coupling", "Copy diagnostic text", "Export diagnostic text"],
        "read_only_note": "Status is read-only; refresh/compare/export actions do not edit geometry or sensor settings.",
    }


def detector_transport_diagnostic_text(report: Mapping[str, object]) -> str:
    """Render the structured detector report with explicit state headings."""

    geometry_state = report.get("local_geometry_state", {})
    guide = report.get("light_guide", {})
    sensor = report.get("sensor_readout_state", {})
    terminals = report.get("terminal_surfaces", {})
    mylar = report.get("mylar", {})
    warnings = report.get("warnings_errors", {})
    sections = geometry_state.get("active_sections", ()) if isinstance(geometry_state, Mapping) else ()
    section_text = "; ".join(f"{item.get('name')} [{item.get('material')}; {item.get('shape')}; {item.get('dimensions_mm')}]" for item in sections)
    guide_text = "not present" if not guide.get("present") else "; ".join(f"{item.get('name')} [{item.get('material')}; {item.get('shape')}; {item.get('dimensions_mm')}]" for item in guide.get("sections", ()))
    return "\n".join((
        f"{str(report.get('detector_type', '')).title()} Transport Diagnostic (read-only)",
        f"Assembly validity: {'VALID' if report.get('assembly_valid') else 'INVALID'}",
        f"Local geometry state: {geometry_state.get('active_section_count', 0)} active section(s): {section_text or 'none'}; impact XY {geometry_state.get('impact_xy_mm')}",
        f"Light guide: {'enabled' if guide.get('enabled') else 'disabled'}; {'valid/present' if guide.get('present') and guide.get('valid') else 'not present' if not guide.get('present') else 'invalid'}; {guide_text}",
        f"Optical transport path: {report.get('optical_transport_path')}",
        f"Terminal surfaces: +X {terminals.get('+X')}; -X {terminals.get('-X')}",
        f"Mylar: {'enabled' if mylar.get('enabled') else 'disabled'}; reflectivity {float(mylar.get('reflectivity', 0.0)) * 100:g}%; overlay surfaces {mylar.get('overlay_surface_count', 0)}",
        f"Mylar interface conflicts: {'; '.join(mylar.get('interface_conflicts', ())) or 'none; terminal optical/mirror interfaces are excluded from overlays'}",
        f"Sensor/readout state: {'present' if sensor.get('present') else 'absent'}; shape {sensor.get('shape')}; dimensions {sensor.get('dimensions_mm')}; QE {sensor.get('quantum_efficiency')}; acceptance {sensor.get('acceptance_deg')} deg; coupling index {sensor.get('optical_window_coupling_index')}; collection {sensor.get('collection_efficiency')}; timing {sensor.get('timing_readout', {}).get('sigma_ns')} ns; coupled terminal {sensor.get('coupled_terminal') or 'none'}; target {sensor.get('coupled_target') or 'none'}",
        f"Physical generation model: {report.get('physical_generation_model', {}).get('model')}; {report.get('physical_generation_model', {}).get('source')}",
        f"Warnings/errors: warnings={'; '.join(warnings.get('warnings', ())) or 'none'}; errors={'; '.join(warnings.get('errors', ())) or 'none'}; anomalies={'; '.join(warnings.get('anomalies', ())) or 'none'}",
        "Actions are temporary diagnostics only; geometry and sensor/readout fields are edited in their detector-specific tabs.",
    ))


def validate_settings(raw: Mapping[str, str]) -> RunSettings:
    """Validate GUI field values and return strongly typed simulation settings."""

    output_text = str(raw.get("output_directory", "")).strip()
    if not output_text:
        raise ValueError("Output directory is required")
    output_directory = Path(output_text).expanduser().resolve()
    events = int(_number(raw, "events", "Events", integer=True))
    photons = int(_number(raw, "photons_per_event", "Photons per event", integer=True))
    sampling_value = str(raw.get("sampling_mode", SAMPLING_MODE_LABELS["fixed"])).strip()
    sampling_mode = SAMPLING_MODE_BY_LABEL.get(sampling_value, sampling_value.lower())
    kinetic = float(_number(raw, "kinetic_mev", "Kinetic energy"))
    mass = float(_number(raw, "mass_mev_c2", "Particle mass"))
    entry_x = float(_number(raw, "entry_x_mm", "Entry x"))
    entry_y = float(_number(raw, "entry_y_mm", "Entry y"))
    inclination = float(_number(raw, "inclination_deg", "Inclination"))
    azimuth = float(_number(raw, "azimuth_deg", "Azimuth"))
    seed = int(_number(raw, "seed", "Seed", integer=True))
    detector_index = float(_number(raw, "sensor_1_index", "Cherenkov sensor coupling index"))
    acceptance = float(_number(raw, "sensor_1_acceptance_deg", "Cherenkov sensor acceptance"))
    detector_qe = float(_number(raw, "sensor_1_qe", "Cherenkov sensor QE"))
    collection = float(_number(raw, "sensor_1_collection", "Cherenkov sensor collection"))
    mirror_reflectivity = float(_number(raw, "mirror_reflectivity", "Mirror reflectivity"))
    track_margin = float(_number(raw, "track_margin_mm", "Track margin"))
    display_photons = int(_number(raw, "display_photons", "Display photons", integer=True))
    maximum_photon_bounces = int(_number(raw, "maximum_photon_bounces", "Maximum photon bounces", integer=True))
    opposite_index = _number(raw, "opposite_detector_index", "Opposite detector index", optional=True)
    opposite_acceptance = _number(raw, "opposite_detector_acceptance_deg", "Opposite detector acceptance", optional=True)
    opposite_qe = _number(raw, "opposite_detector_qe", "Opposite detector QE", optional=True)
    opposite_collection = _number(raw, "opposite_collection_efficiency", "Opposite collection efficiency", optional=True)
    plus_area = detector_active_area_from_values(raw, "detector", "First detector")
    minus_area = detector_active_area_from_values(raw, "opposite_detector", "Second detector")

    if events < 1 or photons < 0 or display_photons < 0:
        raise ValueError("Events must be at least 1; photon counts cannot be negative")
    if maximum_photon_bounces <= 0:
        raise ValueError("Maximum photon bounces must be a positive integer")
    if sampling_mode not in {"fixed", "poisson"}:
        raise ValueError("Photon sampling mode must be fixed or poisson")
    if kinetic < 0.0 or mass <= 0.0:
        raise ValueError("Kinetic energy must be non-negative and mass must be positive")
    if not -50.0 <= entry_x <= 50.0 or not -10.0 <= entry_y <= 10.0:
        raise ValueError("Entry coordinates must be within x +/-50 mm and y +/-10 mm")
    if not -89.0 <= inclination <= 89.0:
        raise ValueError("Inclination from +Z must be in [-89, 89] degrees")
    if detector_index <= 0.0 or not 0.0 <= acceptance <= 180.0:
        raise ValueError("Detector index must be positive and acceptance must be in [0, 180] degrees")
    for label, value in (("Detector QE", detector_qe), ("Collection efficiency", collection), ("Mirror reflectivity", mirror_reflectivity)):
        if not 0.0 <= value <= 1.0:
            raise ValueError(f"{label} must be between 0 and 1")
    if track_margin < 0.0:
        raise ValueError("Track margin must be non-negative")
    mode = str(raw.get("opposite_end_mode", "open")).strip().lower()
    if mode not in {"open", "mirror", "detector"}:
        raise ValueError("Opposite end mode must be open, mirror, or detector")
    if opposite_index is not None and opposite_index <= 0.0:
        raise ValueError("Opposite detector index must be positive")
    if opposite_acceptance is not None and not 0.0 <= opposite_acceptance <= 180.0:
        raise ValueError("Opposite detector acceptance must be in [0, 180] degrees")
    for label, value in (("Opposite detector QE", opposite_qe), ("Opposite collection efficiency", opposite_collection)):
        if value is not None and not 0.0 <= value <= 1.0:
            raise ValueError(f"{label} must be between 0 and 1")

    # Validate only the physical configuration for the selected detector.
    # The old checkbox may still be present in imported dictionaries, but it
    # is deliberately not consulted here.
    detector_run_mode = detector_run_mode_from_values(raw)
    assembly_values = dict(raw) if detector_run_mode in {"scintillator", "combined", GAMMA_RUN_MODE} else default_values()
    scintillator_local = scintillator_local_geometry_from_values(assembly_values)
    scintillator_global = scintillator_global_geometry_from_values(assembly_values)
    scintillator_config = scintillator_config_from_values(assembly_values)
    two_detector = two_detector_beam_monitor_from_values(assembly_values)
    guide = scintillator_local.light_guide
    scint_mode = scintillator_config.sampling_mode
    particle_type = str(raw.get("particle_type", "electron")).strip().lower()
    if particle_type not in set(PARTICLE_MASSES_MEV_C2) | {"custom"}:
        raise ValueError("Particle type is not supported")
    gamma_profile = str(raw.get("gamma_source_profile_path", "")).strip()
    gamma_line_keys = (
        "gamma_line_weight_511_keV", "gamma_line_weight_835_keV", "gamma_line_weight_1173_keV",
        "gamma_line_weight_1275_keV", "gamma_line_weight_1332_keV", "gamma_line_weight_1461_keV",
    )
    gamma_weights = {}
    for key in gamma_line_keys:
        try:
            value = float(raw.get(key, "0") or 0.0)
        except (TypeError, ValueError) as exc:
            raise ValueError(f"{key} must be numeric") from exc
        if not math.isfinite(value) or value < 0.0:
            raise ValueError(f"{key} must be finite and non-negative")
        gamma_weights[key] = value
    try:
        gamma_max_interactions = int(_number(raw, "gamma_max_interactions", "Gamma maximum interactions", integer=True))
        gamma_max_light_samples = int(_number(raw, "gamma_max_light_samples", "Gamma maximum light samples", integer=True))
    except (TypeError, ValueError) as exc:
        raise ValueError(f"Gamma limits are invalid: {exc}") from exc
    if not gamma_profile:
        gamma_profile = str(Path(__file__).resolve().parent / "gamma_source_profile_collimated_v1.json")
    if gamma_max_interactions < 1 or gamma_max_light_samples < 0:
        raise ValueError("Gamma limits must be positive/non-negative")
    if detector_run_mode == GAMMA_RUN_MODE:
        profile_path = Path(gamma_profile).expanduser().resolve()
        try:
            gamma_transport.GammaSource.from_profile(profile_path, {
                511.0: gamma_weights["gamma_line_weight_511_keV"], 835.0: gamma_weights["gamma_line_weight_835_keV"],
                1173.0: gamma_weights["gamma_line_weight_1173_keV"], 1275.0: gamma_weights["gamma_line_weight_1275_keV"],
                1332.0: gamma_weights["gamma_line_weight_1332_keV"], 1461.0: gamma_weights["gamma_line_weight_1461_keV"],
            })
        except (OSError, TypeError, ValueError, json.JSONDecodeError) as exc:
            raise ValueError(f"Gamma source profile is invalid: {exc}") from exc
    return RunSettings(
        output_directory=output_directory,
        events=events,
        photons_per_event=photons,
        sampling_mode=sampling_mode,
        kinetic_mev=kinetic,
        particle_type=particle_type,
        mass_mev_c2=mass,
        entry_x_mm=entry_x,
        entry_y_mm=entry_y,
        inclination_deg=inclination,
        azimuth_deg=azimuth,
        seed=seed,
        detector_index=detector_index,
        detector_acceptance_deg=acceptance,
        detector_qe=detector_qe,
        collection_efficiency=collection,
        opposite_end_mode=mode,
        mirror_reflectivity=mirror_reflectivity,
        opposite_detector_index=None if opposite_index is None else float(opposite_index),
        opposite_detector_acceptance_deg=None if opposite_acceptance is None else float(opposite_acceptance),
        opposite_detector_qe=None if opposite_qe is None else float(opposite_qe),
        opposite_collection_efficiency=None if opposite_collection is None else float(opposite_collection),
        detector_shape=str(plus_area["shape"]),
        detector_active_width_mm=float(plus_area["active_width_mm"]),
        detector_active_height_mm=float(plus_area["active_height_mm"]),
        detector_active_radius_mm=float(plus_area["active_radius_mm"]),
        detector_active_semi_axis_a_mm=float(plus_area["active_semi_axis_a_mm"]),
        detector_active_semi_axis_b_mm=float(plus_area["active_semi_axis_b_mm"]),
        detector_active_offset_y_mm=float(plus_area["active_offset_y_mm"]),
        detector_active_offset_z_mm=float(plus_area["active_offset_z_mm"]),
        detector_active_rotation_deg=float(plus_area["active_rotation_deg"]),
        opposite_detector_shape=str(minus_area["shape"]),
        opposite_detector_active_width_mm=float(minus_area["active_width_mm"]),
        opposite_detector_active_height_mm=float(minus_area["active_height_mm"]),
        opposite_detector_active_radius_mm=float(minus_area["active_radius_mm"]),
        opposite_detector_active_semi_axis_a_mm=float(minus_area["active_semi_axis_a_mm"]),
        opposite_detector_active_semi_axis_b_mm=float(minus_area["active_semi_axis_b_mm"]),
        opposite_detector_active_offset_y_mm=float(minus_area["active_offset_y_mm"]),
        opposite_detector_active_offset_z_mm=float(minus_area["active_offset_z_mm"]),
        opposite_detector_active_rotation_deg=float(minus_area["active_rotation_deg"]),
        track_margin_mm=track_margin,
        display_photons=display_photons,
        scintillator_shape=scintillator_local.scintillator.shape,
        scintillator_length_x_mm=scintillator_local.scintillator.length_x_m * 1e3,
        scintillator_width_y_mm=scintillator_local.scintillator.width_y_m * 1e3,
        scintillator_thickness_z_mm=scintillator_local.scintillator.thickness_z_m * 1e3,
        scintillator_width_y_start_mm=scintillator_local.scintillator.width_y_start_m * 1e3,
        scintillator_width_y_end_mm=scintillator_local.scintillator.width_y_end_m * 1e3,
        scintillator_thickness_z_start_mm=scintillator_local.scintillator.thickness_z_start_m * 1e3,
        scintillator_thickness_z_end_mm=scintillator_local.scintillator.thickness_z_end_m * 1e3,
        scintillator_axis=scintillator_local.scintillator.axis,
        scintillator_radius_mm=scintillator_local.scintillator.radius_m * 1e3,
        scintillator_length_mm=scintillator_local.scintillator.length_m * 1e3,
        scintillator_semi_axis_a_mm=scintillator_local.scintillator.semi_axis_a_m * 1e3,
        scintillator_semi_axis_b_mm=scintillator_local.scintillator.semi_axis_b_m * 1e3,
        scintillator_radius_start_mm=scintillator_local.scintillator.radius_start_m * 1e3,
        scintillator_radius_end_mm=scintillator_local.scintillator.radius_end_m * 1e3,
        scintillator_center_x_mm=scintillator_local.scintillator.center_x_m * 1e3,
        scintillator_center_y_mm=scintillator_local.scintillator.center_y_m * 1e3,
        scintillator_center_z_mm=scintillator_local.scintillator.center_z_m * 1e3,
        scintillator_material=scintillator_local.scintillator.material,
        scintillator_refractive_index=float(scintillator_local.scintillator.refractive_index_value),
        guide_shape="none" if guide is None else guide.shape,
        guide_length_x_mm=0.0 if guide is None else guide.length_x_m * 1e3,
        guide_width_y_mm=0.0 if guide is None else guide.width_y_m * 1e3,
        guide_thickness_z_mm=0.0 if guide is None else guide.thickness_z_m * 1e3,
        guide_width_y_start_mm=0.0 if guide is None else guide.width_y_start_m * 1e3,
        guide_width_y_end_mm=0.0 if guide is None else guide.width_y_end_m * 1e3,
        guide_thickness_z_start_mm=0.0 if guide is None else guide.thickness_z_start_m * 1e3,
        guide_thickness_z_end_mm=0.0 if guide is None else guide.thickness_z_end_m * 1e3,
        guide_axis="" if guide is None else guide.axis,
        guide_radius_mm=0.0 if guide is None else guide.radius_m * 1e3,
        guide_length_mm=0.0 if guide is None else guide.length_m * 1e3,
        guide_semi_axis_a_mm=0.0 if guide is None else guide.semi_axis_a_m * 1e3,
        guide_semi_axis_b_mm=0.0 if guide is None else guide.semi_axis_b_m * 1e3,
        guide_radius_start_mm=0.0 if guide is None else guide.radius_start_m * 1e3,
        guide_radius_end_mm=0.0 if guide is None else guide.radius_end_m * 1e3,
        guide_center_x_mm=0.0 if guide is None else guide.center_x_m * 1e3,
        guide_center_y_mm=0.0 if guide is None else guide.center_y_m * 1e3,
        guide_center_z_mm=0.0 if guide is None else guide.center_z_m * 1e3,
        guide_material="" if guide is None else guide.material,
        guide_refractive_index=0.0 if guide is None else float(guide.refractive_index_value),
        scintillator_gap_mm=scintillator_global.gap_from_quartz_m * 1e3,
        scintillator_sampling_mode=scint_mode,
        scintillator_photon_count=scintillator_config.photon_count,
        scintillator_yield_per_MeV=scintillator_config.photon_yield_per_MeV,
        scintillator_deposited_energy_MeV=scintillator_config.deposited_energy_MeV,
        scintillator_dEdx_MeV_per_m=scintillator_config.dEdx_MeV_per_m,
        scintillator_decay_time_ns=scintillator_config.decay_time_s * 1e9,
        scintillator_wavelength_min_nm=scintillator_config.emission_wavelength_min_m * 1e9,
        scintillator_wavelength_max_nm=scintillator_config.emission_wavelength_max_m * 1e9,
        scintillator_readout_index=scintillator_config.readout_index,
        scintillator_readout_acceptance_deg=math.degrees(scintillator_config.readout_acceptance_rad),
        scintillator_readout_qe=scintillator_config.readout_quantum_efficiency,
        scintillator_readout_collection_efficiency=scintillator_config.readout_collection_efficiency,
        scintillator_absorption_length_mm=None if scintillator_config.absorption_length_m is None else scintillator_config.absorption_length_m * 1e3,
        # Typed run settings expose the effective value, while the legacy
        # StringVar/profile field remains untouched for compatibility.
        scintillator_max_bounces=maximum_photon_bounces,
        maximum_photon_bounces=maximum_photon_bounces,
        mylar_enabled=_bool_text(assembly_values, "mylar_enabled"),
        mylar_reflectivity_percent=float(_value_text(assembly_values, "mylar_reflectivity_percent", "95")),
        global_offset_x_mm=scintillator_global.offset_x_m * 1e3,
        global_offset_y_mm=scintillator_global.offset_y_m * 1e3,
        global_offset_z_mm=scintillator_global.offset_z_m * 1e3,
        global_beam_axis=scintillator_global.beam_axis,
        detector_separation_mm=two_detector.detector_separation_m * 1e3,
        detector2_offset_x_mm=two_detector.placement_2.offset_x_m * 1e3,
        detector2_offset_y_mm=two_detector.placement_2.offset_y_m * 1e3,
        detector2_offset_z_mm=two_detector.placement_2.offset_z_m * 1e3,
        detector1_guide_left_shape="none" if two_detector.detector_1.guide_left is None else two_detector.detector_1.guide_left.shape,
        detector1_guide_right_shape="none" if two_detector.detector_1.guide_right is None else two_detector.detector_1.guide_right.shape,
        detector2_scintillator_shape=two_detector.detector_2.scintillator.shape,
        detector2_scintillator_material=two_detector.detector_2.scintillator.material,
        detector2_scintillator_refractive_index=float(two_detector.detector_2.scintillator.refractive_index_value),
        detector2_guide_left_shape="none" if two_detector.detector_2.guide_left is None else two_detector.detector_2.guide_left.shape,
        detector2_guide_right_shape="none" if two_detector.detector_2.guide_right is None else two_detector.detector_2.guide_right.shape,
        detector_1_z_mm=float(two_detector.detector_1_z_m or 0.0) * 1e3,
        detector_2_z_mm=float(two_detector.detector_2_z_m or 0.0) * 1e3,
        detector1_guide_left_enabled=two_detector.detector_1.guide_left is not None,
        detector1_guide_right_enabled=two_detector.detector_1.guide_right is not None,
        detector2_guide_left_enabled=two_detector.detector_2.guide_left is not None,
        detector2_guide_right_enabled=two_detector.detector_2.guide_right is not None,
        sensor_profile_values={key: _value_text(raw, key, "") for key in SENSOR_PROFILE_KEYS},
        detector_run_mode=detector_run_mode,
        assembly_1_rotation_axis=_value_text(raw, "assembly_1_rotation_axis", "None"),
        assembly_1_rotation_angle_deg=float(_value_text(raw, "assembly_1_rotation_angle_deg", "0")),
        assembly_2_rotation_axis=_value_text(raw, "assembly_2_rotation_axis", "None"),
        assembly_2_rotation_angle_deg=float(_value_text(raw, "assembly_2_rotation_angle_deg", "0")),
        gamma_source_profile_path=gamma_profile,
        gamma_line_weight_511_keV=gamma_weights["gamma_line_weight_511_keV"],
        gamma_line_weight_835_keV=gamma_weights["gamma_line_weight_835_keV"],
        gamma_line_weight_1173_keV=gamma_weights["gamma_line_weight_1173_keV"],
        gamma_line_weight_1275_keV=gamma_weights["gamma_line_weight_1275_keV"],
        gamma_line_weight_1332_keV=gamma_weights["gamma_line_weight_1332_keV"],
        gamma_line_weight_1461_keV=gamma_weights["gamma_line_weight_1461_keV"],
        gamma_max_interactions=gamma_max_interactions,
        gamma_max_light_samples=gamma_max_light_samples,
    )


DETECTOR_PROFILE_FORMAT = "quartz-cherenkov-detector-profile"
DETECTOR_PROFILE_VERSION = 2
DETECTOR_PROFILE_KEYS = (
    "detector_index",
    "detector_acceptance_deg",
    "detector_qe",
    "collection_efficiency",
    "opposite_end_mode",
    "mirror_reflectivity",
    "opposite_detector_index",
    "opposite_detector_acceptance_deg",
    "opposite_detector_qe",
    "opposite_collection_efficiency",
    "detector_shape",
    "detector_active_width_mm",
    "detector_active_height_mm",
    "detector_active_radius_mm",
    "detector_active_semi_axis_a_mm",
    "detector_active_semi_axis_b_mm",
    "detector_active_offset_y_mm",
    "detector_active_offset_z_mm",
    "detector_active_rotation_deg",
    "opposite_detector_shape",
    "opposite_detector_active_width_mm",
    "opposite_detector_active_height_mm",
    "opposite_detector_active_radius_mm",
    "opposite_detector_active_semi_axis_a_mm",
    "opposite_detector_active_semi_axis_b_mm",
    "opposite_detector_active_offset_y_mm",
    "opposite_detector_active_offset_z_mm",
    "opposite_detector_active_rotation_deg",
)


def _mapping_value(values: Mapping[str, object], key: str):
    if key not in values:
        raise ValueError(f"Detector profile is missing key: {key}")
    value = values[key]
    return value.get() if hasattr(value, "get") and callable(value.get) else value


def _profile_number(values: Mapping[str, object], key: str, label: str, optional: bool = False):
    value = _mapping_value(values, key)
    if value is None or (isinstance(value, str) and not value.strip()):
        if optional:
            return ""
        raise ValueError(f"{label} is required")
    if isinstance(value, bool):
        raise ValueError(f"{label} must be numeric")
    try:
        number = float(value)
    except (TypeError, ValueError) as exc:
        raise ValueError(f"{label} must be numeric") from exc
    if not math.isfinite(number):
        raise ValueError(f"{label} must be finite")
    return number


def validate_detector_profile(profile: Mapping[str, object]) -> Dict[str, object]:
    """Validate and normalize a standalone detector JSON profile."""

    if not isinstance(profile, Mapping):
        raise ValueError("Detector profile must be a JSON object")
    if not isinstance(profile.get("format"), str) or profile.get("format") != DETECTOR_PROFILE_FORMAT:
        raise ValueError("Invalid detector profile format")
    if isinstance(profile.get("version"), bool) or not isinstance(profile.get("version"), int) or profile.get("version") not in {1, DETECTOR_PROFILE_VERSION}:
        raise ValueError("Unsupported detector profile version")
    plus_area = detector_active_area_from_values(profile, "detector", "First detector")
    minus_area = detector_active_area_from_values(profile, "opposite_detector", "Second detector")
    normalized = {
        "format": DETECTOR_PROFILE_FORMAT,
        "version": DETECTOR_PROFILE_VERSION,
        "detector_index": _profile_number(profile, "detector_index", "Detector index"),
        "detector_acceptance_deg": _profile_number(profile, "detector_acceptance_deg", "Detector acceptance"),
        "detector_qe": _profile_number(profile, "detector_qe", "Detector QE"),
        "collection_efficiency": _profile_number(profile, "collection_efficiency", "Collection efficiency"),
        "opposite_end_mode": str(_mapping_value(profile, "opposite_end_mode")).strip().lower(),
        "mirror_reflectivity": _profile_number(profile, "mirror_reflectivity", "Mirror reflectivity"),
        "opposite_detector_index": _profile_number(profile, "opposite_detector_index", "Opposite detector index", optional=True),
        "opposite_detector_acceptance_deg": _profile_number(profile, "opposite_detector_acceptance_deg", "Opposite detector acceptance", optional=True),
        "opposite_detector_qe": _profile_number(profile, "opposite_detector_qe", "Opposite detector QE", optional=True),
        "opposite_collection_efficiency": _profile_number(profile, "opposite_collection_efficiency", "Opposite collection efficiency", optional=True),
        "detector_shape": plus_area["shape"],
        "detector_active_width_mm": plus_area["active_width_mm"],
        "detector_active_height_mm": plus_area["active_height_mm"],
        "detector_active_radius_mm": plus_area["active_radius_mm"],
        "detector_active_semi_axis_a_mm": plus_area["active_semi_axis_a_mm"],
        "detector_active_semi_axis_b_mm": plus_area["active_semi_axis_b_mm"],
        "detector_active_offset_y_mm": plus_area["active_offset_y_mm"],
        "detector_active_offset_z_mm": plus_area["active_offset_z_mm"],
        "detector_active_rotation_deg": plus_area["active_rotation_deg"],
        "opposite_detector_shape": minus_area["shape"],
        "opposite_detector_active_width_mm": minus_area["active_width_mm"],
        "opposite_detector_active_height_mm": minus_area["active_height_mm"],
        "opposite_detector_active_radius_mm": minus_area["active_radius_mm"],
        "opposite_detector_active_semi_axis_a_mm": minus_area["active_semi_axis_a_mm"],
        "opposite_detector_active_semi_axis_b_mm": minus_area["active_semi_axis_b_mm"],
        "opposite_detector_active_offset_y_mm": minus_area["active_offset_y_mm"],
        "opposite_detector_active_offset_z_mm": minus_area["active_offset_z_mm"],
        "opposite_detector_active_rotation_deg": minus_area["active_rotation_deg"],
    }
    normalized["version"] = DETECTOR_PROFILE_VERSION
    if normalized["opposite_end_mode"] not in {"open", "mirror", "detector"}:
        raise ValueError("Opposite end mode must be open, mirror, or detector")
    if normalized["detector_index"] <= 0.0:
        raise ValueError("Detector index must be positive")
    if not 0.0 <= normalized["detector_acceptance_deg"] <= 180.0:
        raise ValueError("Detector acceptance must be in [0, 180] degrees")
    for key, label in (
        ("detector_qe", "Detector QE"),
        ("collection_efficiency", "Collection efficiency"),
        ("mirror_reflectivity", "Mirror reflectivity"),
    ):
        if not 0.0 <= normalized[key] <= 1.0:
            raise ValueError(f"{label} must be between 0 and 1")
    optional_ranges = (
        ("opposite_detector_index", "Opposite detector index", lambda value: value > 0.0),
        ("opposite_detector_acceptance_deg", "Opposite detector acceptance", lambda value: 0.0 <= value <= 180.0),
        ("opposite_detector_qe", "Opposite detector QE", lambda value: 0.0 <= value <= 1.0),
        ("opposite_collection_efficiency", "Opposite collection efficiency", lambda value: 0.0 <= value <= 1.0),
    )
    for key, label, predicate in optional_ranges:
        value = normalized[key]
        if value != "" and not predicate(value):
            raise ValueError(f"{label} has an invalid range")
    return normalized


def detector_profile_from_values(values: Mapping[str, object]) -> Dict[str, object]:
    """Collect and validate detector fields from StringVars or plain values."""

    profile = {"format": DETECTOR_PROFILE_FORMAT, "version": DETECTOR_PROFILE_VERSION}
    plus_area = detector_active_area_from_values(values, "detector", "First detector")
    minus_area = detector_active_area_from_values(values, "opposite_detector", "Second detector")
    profile.update({
        "detector_index": _profile_number(values, "detector_index", "Detector index"),
        "detector_acceptance_deg": _profile_number(values, "detector_acceptance_deg", "Detector acceptance"),
        "detector_qe": _profile_number(values, "detector_qe", "Detector QE"),
        "collection_efficiency": _profile_number(values, "collection_efficiency", "Collection efficiency"),
        "opposite_end_mode": str(_mapping_value(values, "opposite_end_mode")).strip().lower(),
        "mirror_reflectivity": _profile_number(values, "mirror_reflectivity", "Mirror reflectivity"),
        "opposite_detector_index": _profile_number(values, "opposite_detector_index", "Opposite detector index", optional=True),
        "opposite_detector_acceptance_deg": _profile_number(values, "opposite_detector_acceptance_deg", "Opposite detector acceptance", optional=True),
        "opposite_detector_qe": _profile_number(values, "opposite_detector_qe", "Opposite detector QE", optional=True),
        "opposite_collection_efficiency": _profile_number(values, "opposite_collection_efficiency", "Opposite collection efficiency", optional=True),
        "detector_shape": plus_area["shape"],
        "detector_active_width_mm": plus_area["active_width_mm"],
        "detector_active_height_mm": plus_area["active_height_mm"],
        "detector_active_radius_mm": plus_area["active_radius_mm"],
        "detector_active_semi_axis_a_mm": plus_area["active_semi_axis_a_mm"],
        "detector_active_semi_axis_b_mm": plus_area["active_semi_axis_b_mm"],
        "detector_active_offset_y_mm": plus_area["active_offset_y_mm"],
        "detector_active_offset_z_mm": plus_area["active_offset_z_mm"],
        "detector_active_rotation_deg": plus_area["active_rotation_deg"],
        "opposite_detector_shape": minus_area["shape"],
        "opposite_detector_active_width_mm": minus_area["active_width_mm"],
        "opposite_detector_active_height_mm": minus_area["active_height_mm"],
        "opposite_detector_active_radius_mm": minus_area["active_radius_mm"],
        "opposite_detector_active_semi_axis_a_mm": minus_area["active_semi_axis_a_mm"],
        "opposite_detector_active_semi_axis_b_mm": minus_area["active_semi_axis_b_mm"],
        "opposite_detector_active_offset_y_mm": minus_area["active_offset_y_mm"],
        "opposite_detector_active_offset_z_mm": minus_area["active_offset_z_mm"],
        "opposite_detector_active_rotation_deg": minus_area["active_rotation_deg"],
    })
    return validate_detector_profile(profile)


def apply_detector_profile(values: Mapping[str, object], profile: Mapping[str, object]) -> Dict[str, object]:
    """Apply a fully validated profile atomically to StringVars or a dict."""

    normalized = validate_detector_profile(profile)
    old_values = {key: _mapping_value(values, key) for key in DETECTOR_PROFILE_KEYS if key in values}

    def set_value(key: str, value: object) -> None:
        target = values[key]
        display = "" if value == "" else str(value)
        if hasattr(target, "set") and callable(target.set):
            target.set(display)
        else:
            values[key] = display  # type: ignore[index]

    try:
        for key in DETECTOR_PROFILE_KEYS:
            if key not in values:
                continue
            if key == "opposite_end_mode":
                continue
            set_value(key, normalized[key])
        set_value("opposite_end_mode", normalized["opposite_end_mode"])
        # Keep the legacy run-settings aliases synchronized with the
        # detector-specific Cherenkov Sensor tab after profile import.
        for legacy, sensor_key in (("detector_index", "sensor_1_index"), ("detector_acceptance_deg", "sensor_1_acceptance_deg"), ("detector_qe", "sensor_1_qe"), ("collection_efficiency", "sensor_1_collection")):
            if sensor_key in values:
                set_value(sensor_key, normalized[legacy])
    except Exception:
        for key, value in old_values.items():
            set_value(key, value)
        raise
    return normalized


def save_detector_profile(path: Path | str, values: Mapping[str, object]) -> Dict[str, object]:
    """Validate and save detector fields; current values are never changed."""

    candidate = dict(values)
    for legacy, sensor_key in (("detector_index", "sensor_1_index"), ("detector_acceptance_deg", "sensor_1_acceptance_deg"), ("detector_qe", "sensor_1_qe"), ("collection_efficiency", "sensor_1_collection")):
        if sensor_key in values:
            candidate[legacy] = _value_text(values, sensor_key)
    profile = detector_profile_from_values(candidate)
    output_path = Path(path).expanduser().resolve()
    output_path.write_text(json.dumps(profile, indent=2, sort_keys=True) + "\n", encoding="utf-8")
    return profile


def load_detector_profile(path: Path | str, values: Mapping[str, object]) -> Dict[str, object]:
    """Load, validate, then apply a detector profile without partial mutation."""

    input_path = Path(path).expanduser().resolve()
    profile = json.loads(input_path.read_text(encoding="utf-8"))
    return apply_detector_profile(values, profile)


def create_run_directory(base_directory: Path) -> Path:
    """Create a unique timestamped output directory and return its path."""

    base = Path(base_directory).expanduser().resolve()
    base.mkdir(parents=True, exist_ok=True)
    stamp = datetime.now().strftime("%Y%m%d_%H%M%S")
    for suffix in range(1000):
        name = f"quartz_run_{stamp}" if suffix == 0 else f"quartz_run_{stamp}_{suffix:03d}"
        candidate = base / name
        try:
            candidate.mkdir()
            return candidate
        except FileExistsError:
            continue
    raise RuntimeError(f"Could not create a unique run directory under {base}")


def _particle_for_settings(settings: RunSettings) -> simulator.Particle:
    return simulator.Particle(
        mass_kg=settings.mass_mev_c2 * simulator.MEV / simulator.C_LIGHT**2,
        kinetic_energy_J=settings.kinetic_mev * simulator.MEV,
        entry_x_m=settings.entry_x_mm * 1e-3,
        entry_y_m=settings.entry_y_mm * 1e-3,
        # Run owns the incidence control.  The supported GUI convention is
        # the X-Z plane, so azimuth is fixed at zero even though the typed
        # field remains in RunSettings for old profiles/CLI compatibility.
        inclination_rad=math.radians(settings.inclination_deg),
        azimuth_rad=0.0,
        charge_e=1.0,
    )


def _align_scintillator_transport_to_global_track(
    assembly: simulator.ScintillatorAssembly,
    scene: Mapping[str, object],
) -> simulator.ScintillatorAssembly:
    """Use the common global line's local intersection as S transport input.

    The scintillator optical model is transported in its assembly-local frame.
    For an All run, convert the frozen global anchor/direction into that frame
    and use the active volume's +Z reference plane to derive the local x/y
    impact.  If the line misses the active solid, leave the configured impact
    untouched; the caller reports the preflight mismatch instead of moving a
    detector to make it fit.
    """

    if not scene.get("particle_direction_global") or not scene.get("beam_anchor_global_mm"):
        return assembly
    payload = scene.get("rigid_transforms", {}).get("scintillator")
    if not isinstance(payload, Mapping):
        return assembly
    transform = RigidAssemblyTransform(**{
        key: payload[key] for key in (
            "translation_x_mm", "translation_y_mm", "translation_z_mm",
            "rotation_axis", "rotation_angle_deg", "rotation_origin_x_mm",
            "rotation_origin_y_mm", "rotation_origin_z_mm",
        ) if key in payload
    })
    origin = np.asarray(transform.rotation_origin_mm, dtype=float)
    anchor_global = np.asarray(scene["beam_anchor_global_mm"], dtype=float)
    direction_global = np.asarray(scene["particle_direction_global"], dtype=float)
    local_anchor = transform.rotation_matrix().T.dot(
        anchor_global - origin - np.asarray((transform.translation_x_mm, transform.translation_y_mm, transform.translation_z_mm), dtype=float)
    ) + origin
    local_direction = transform.rotation_matrix().T.dot(direction_global)
    norm = float(np.linalg.norm(local_direction))
    if norm <= 0.0 or abs(float(local_direction[2])) < 1e-12:
        return assembly
    local_direction /= norm
    active_low = min(float(section.bounds_min_m[2]) for section in assembly.active_sections)
    active_high = max(float(section.bounds_max_m[2]) for section in assembly.active_sections)
    reference_z = active_high if local_direction[2] > 0.0 else active_low
    t = (reference_z - float(local_anchor[2])) / float(local_direction[2])
    reference = local_anchor + t * local_direction
    try:
        start = reference - local_direction * max(4.0 * assembly.transport_geometry.characteristic_length_m, 1e-6)
        intervals = assembly.transport_geometry._particle_intervals(start, local_direction)
    except (AttributeError, RuntimeError, ValueError):
        return assembly
    if not intervals:
        return assembly
    return replace(
        assembly,
        particle_impact_x_m=float(reference[0]),
        particle_impact_y_m=float(reference[1]),
    )


def _optical_for_settings(settings: RunSettings) -> simulator.OpticalConfig:
    return simulator.OpticalConfig(
        photon_count=settings.photons_per_event,
        max_bounces=settings.maximum_photon_bounces,
        sampling_mode=settings.sampling_mode,
        detector_index=settings.detector_index,
        detector_acceptance_rad=math.radians(settings.detector_acceptance_deg),
        detector_quantum_efficiency=settings.detector_qe,
        detector_collection_efficiency=settings.collection_efficiency,
        opposite_end_mode=settings.opposite_end_mode,
        opposite_mirror_reflectivity=settings.mirror_reflectivity,
        opposite_detector_index=settings.opposite_detector_index,
        opposite_detector_acceptance_rad=None if settings.opposite_detector_acceptance_deg is None else math.radians(settings.opposite_detector_acceptance_deg),
        opposite_detector_quantum_efficiency=settings.opposite_detector_qe,
        opposite_detector_collection_efficiency=settings.opposite_collection_efficiency,
        detector_shape=settings.detector_shape,
        detector_active_width_m=settings.detector_active_width_mm * 1e-3,
        detector_active_height_m=settings.detector_active_height_mm * 1e-3,
        detector_active_radius_m=settings.detector_active_radius_mm * 1e-3,
        detector_active_semi_axis_a_m=settings.detector_active_semi_axis_a_mm * 1e-3,
        detector_active_semi_axis_b_m=settings.detector_active_semi_axis_b_mm * 1e-3,
        detector_active_offset_y_m=settings.detector_active_offset_y_mm * 1e-3,
        detector_active_offset_z_m=settings.detector_active_offset_z_mm * 1e-3,
        detector_active_rotation_rad=math.radians(settings.detector_active_rotation_deg),
        opposite_detector_shape=settings.opposite_detector_shape,
        opposite_detector_active_width_m=settings.opposite_detector_active_width_mm * 1e-3,
        opposite_detector_active_height_m=settings.opposite_detector_active_height_mm * 1e-3,
        opposite_detector_active_radius_m=settings.opposite_detector_active_radius_mm * 1e-3,
        opposite_detector_active_semi_axis_a_m=settings.opposite_detector_active_semi_axis_a_mm * 1e-3,
        opposite_detector_active_semi_axis_b_m=settings.opposite_detector_active_semi_axis_b_mm * 1e-3,
        opposite_detector_active_offset_y_m=settings.opposite_detector_active_offset_y_mm * 1e-3,
        opposite_detector_active_offset_z_m=settings.opposite_detector_active_offset_z_mm * 1e-3,
        opposite_detector_active_rotation_rad=math.radians(settings.opposite_detector_active_rotation_deg),
    )


def validate_canonical_run_assemblies(
    values: Mapping[str, object], assembly_models: Mapping[str, AssemblyGeometryModel],
    detector_run_mode: str,
) -> Dict[str, object]:
    """Validate the participating complete assemblies for a current run.

    The editor models are the authority here.  In particular, no generic
    legacy detector-area fields are consulted; the photosensor window is
    checked against the actual active-terminal or guide-output dimensions.
    """

    mode = detector_run_mode_from_values({"detector_run_mode": detector_run_mode})
    selected = ("cherenkov",) if mode == "cherenkov" else ("scintillator",) if mode == "scintillator" else ("cherenkov", "scintillator")
    results: Dict[str, object] = {}
    for assembly in ("cherenkov", "scintillator"):
        model = assembly_models.get(assembly)
        if assembly not in selected:
            results[assembly] = {"assembly": assembly.title(), "participates": False, "validated": False}
            continue
        if not isinstance(model, AssemblyGeometryModel):
            raise ValueError(f"{assembly.title()} canonical assembly is missing")
        index = "1" if assembly == "cherenkov" else "2"
        sensor = complete_sensor_readout_from_values(values, assembly, model)
        try:
            results[assembly] = model.validate_complete_assembly(
                sensor_shape=str(sensor["shape"]),
                sensor_width_mm=float(sensor["dimensions_mm"]["width"]),
                sensor_height_mm=float(sensor["dimensions_mm"]["height"]),
                sensor_depth_mm=float(sensor["dimensions_mm"]["depth"]),
                sensor_radius_mm=float(sensor["readout_radius_mm"]),
                sensor_semi_axis_a_mm=float(sensor["readout_semiaxes_mm"]["a"]),
                sensor_semi_axis_b_mm=float(sensor["readout_semiaxes_mm"]["b"]),
            )
        except (TypeError, ValueError) as exc:
            raise ValueError(f"{model.name} run validation failed: {exc}") from exc
        results[assembly]["participates"] = True
    return results


def _safe_plot(log: Callable[[str], None], filename: Path, function: Callable[..., None], *args) -> Optional[str]:
    try:
        function(args[0], str(filename), *args[1:])
        log(f"Saved {filename.name}")
        return None
    except Exception as exc:
        message = f"Could not save {filename.name}: {exc}"
        log(message)
        return message


def _gamma_line_weights_from_settings(settings: RunSettings) -> Dict[float, float]:
    """Return only explicitly positive user line weights; zero stays disabled."""

    return {
        511.0: float(settings.gamma_line_weight_511_keV),
        835.0: float(settings.gamma_line_weight_835_keV),
        1173.0: float(settings.gamma_line_weight_1173_keV),
        1275.0: float(settings.gamma_line_weight_1275_keV),
        1332.0: float(settings.gamma_line_weight_1332_keV),
        1461.0: float(settings.gamma_line_weight_1461_keV),
    }


def _gamma_rigid_transform_from_values(values: Mapping[str, object], assembly: str) -> RigidAssemblyTransform:
    """Return the same canonical assembly transform used by the global view."""

    return _canonical_transform_from_values(values, assembly)


def _gamma_snapshot_volumes(
    assembly_models: Mapping[str, AssemblyGeometryModel],
    values: Mapping[str, object],
) -> list[gamma_transport.TransportVolume]:
    """Build Gamma transport volumes from the two complete run assemblies.

    These are the same detached Cherenkov and Scintillator models and rigid
    transforms consumed by ``global_assembly_scene``.  Missing run snapshots
    are errors: reconstructing detector solids from legacy dimensions would
    make the event viewer disagree with the geometry that was transported.
    """

    required = ("cherenkov", "scintillator")
    if not isinstance(assembly_models, Mapping) or not all(
        isinstance(assembly_models.get(key), AssemblyGeometryModel) for key in required
    ):
        raise ValueError("Gamma transport requires both complete assembly models from the run snapshot")

    volumes: list[gamma_transport.TransportVolume] = []
    for assembly_type in required:
        model = assembly_models[assembly_type]
        model.validate()
        guide_ids = {id(section) for section in model.valid_guide_sections}
        entries = [
            (section, "guide" if id(section) in guide_ids else "active")
            for section in model.active_sections()
        ]
        if not any(role == "active" for _section, role in entries):
            raise ValueError(f"{assembly_type.title()} run snapshot has no active Gamma transport section")
        transform = _gamma_rigid_transform_from_values(values, assembly_type)
        rotation = transform.rotation_matrix()
        translation = np.asarray(
            (transform.translation_x_mm, transform.translation_y_mm, transform.translation_z_mm), dtype=float,
        ) * 1e-3
        origin = np.asarray(transform.rotation_origin_mm, dtype=float) * 1e-3
        for index, (section, role) in enumerate(entries, 1):
            solid = section.to_simulation_geometry()
            solid.validate()
            placed = gamma_transport.rigid_transform_solid(solid, rotation, translation, origin)
            component = (
                "quartz" if assembly_type == "cherenkov" and role == "active"
                else "scintillator" if assembly_type == "scintillator" and role == "active"
                else "guide"
            )
            volumes.append(gamma_transport.TransportVolume(
                f"{assembly_type.title()} {component} section {index}",
                assembly_type, component, placed, solid.material,
                10000.0 if component == "scintillator" else 0.0,
            ))
    return volumes


def gamma_run_geometry_snapshot(
    settings: RunSettings,
    geometry,
    assembly_models: Mapping[str, AssemblyGeometryModel],
) -> Dict[str, object]:
    """Freeze the exact complete Global Assembly scene used by a Gamma run."""

    if not isinstance(assembly_models, Mapping) or any(
        not isinstance(assembly_models.get(key), AssemblyGeometryModel)
        for key in ("cherenkov", "scintillator")
    ):
        raise ValueError("Gamma geometry snapshot requires both complete run-snapshot assembly models")
    models = {key: assembly_models[key] for key in ("cherenkov", "scintillator")}
    for model in models.values():
        model.validate()
    values = _canonical_scene_values_from_settings(settings)
    scene = global_assembly_scene(values, geometry, models)
    compact_scene = compact_event_view_payload(scene)
    if compact_scene.get("assembly_types") != ["cherenkov", "scintillator"]:
        raise ValueError("Gamma Global Assembly snapshot must contain Cherenkov and Scintillator")
    if any(not item.get("components") for item in compact_scene.get("assemblies", ())):
        raise ValueError("Gamma Global Assembly snapshot contains an assembly with no renderable components")
    profiles = {
        key: _json_safe(assembly_geometry_to_dict(deepcopy(models[key]), key))
        for key in ("cherenkov", "scintillator")
    }
    return {
        "schema": "gamma-global-assembly-geometry-v1",
        "coordinate_system": "global_mm",
        "source": "detached RunSnapshot -> global_assembly_scene",
        "assembly_count": 2,
        "assembly_types": ["cherenkov", "scintillator"],
        "profiles": profiles,
        "sensor_profile": _json_safe(sensor_profile_from_values(values, models)),
        "assembly_transforms": {
            key: _json_safe(scene["assembly_by_type"][key]["transform"])
            for key in ("cherenkov", "scintillator")
        },
        "scene": _json_safe(compact_scene),
    }


def run_gamma_transport_to_directory(
    settings: RunSettings,
    run_directory: Path,
    geometry,
    *,
    log: Callable[[str], None],
    progress_callback=None,
    status_callback=None,
    stop_event=None,
    assembly_models: Optional[Mapping[str, AssemblyGeometryModel]] = None,
) -> Dict[str, object]:
    """Run gamma transport from the detached GUI snapshot and persist its JSON.

    Gamma solids are reconstructed only from the detached complete assembly
    models. Their rigid transforms and the saved renderer scene come from the
    same canonical Global Assembly builder; no fallback/default geometry is
    substituted when the run snapshot is missing.
    """

    values = _canonical_scene_values_from_settings(settings)
    geometry_snapshot = gamma_run_geometry_snapshot(settings, geometry, assembly_models or {})
    profile_path = Path(settings.gamma_source_profile_path).expanduser().resolve()
    source = gamma_transport.GammaSource.from_profile(profile_path, _gamma_line_weights_from_settings(settings))
    volumes = _gamma_snapshot_volumes(assembly_models or {}, values)
    log("Gamma source: thin centered collimated +Z; continuum uses approved relative profile")
    active_lines = [(energy, weight) for energy, weight in _gamma_line_weights_from_settings(settings).items() if weight > 0.0]
    log(f"Gamma active user line weights: {active_lines or 'none (continuum-only default)'}")
    log("Gamma semantic boundary: local expected light only; no optical transport or sensor detection")
    result = gamma_transport.simulate(
        source, volumes, events=settings.events, seed=settings.seed,
        max_interactions=settings.gamma_max_interactions,
        max_light_samples=settings.gamma_max_light_samples,
        progress_callback=progress_callback,
        status_callback=status_callback,
        stop_event=stop_event,
    )
    gamma_payload = result.to_dict()
    gamma_summary = dict(result.summary())
    gamma_summary.update({
        "detector_mode": GAMMA_RUN_MODE,
        "detector_mode_label": GAMMA_RUN_MODE_LABEL,
        "geometry_source": geometry_snapshot["source"],
        "geometry_snapshot_schema": geometry_snapshot["schema"],
        "geometry_assemblies": list(geometry_snapshot["assembly_types"]),
    })
    gamma_payload["summary"] = gamma_summary
    gamma_payload["geometry_snapshot"] = geometry_snapshot
    gamma_path = run_directory / "gamma_results.json"
    gamma_path.write_text(json.dumps(_json_safe(gamma_payload), indent=2, sort_keys=True) + "\n", encoding="utf-8")
    log(f"Saved gamma_results.json ({gamma_path.stat().st_size} bytes)")
    try:
        gamma_transport._write_plot(result, run_directory / "gamma_transport.png")
        log("Saved gamma_transport.png")
    except Exception as exc:
        log(f"Could not save gamma_transport.png: {type(exc).__name__}: {exc}")
    return {
        "run_directory": str(run_directory),
        "settings": asdict(settings),
        "geometry": asdict(geometry) if hasattr(geometry, "__dataclass_fields__") else {"shape": "unknown"},
        "summary": gamma_summary,
        "gamma_results_path": str(gamma_path),
        "gamma_model_status": gamma_summary.get("model_status"),
        "stopped": bool(result.stopped),
        "root_path": None,
        "root_error": None,
        "readout_arrivals_path": None,
        "readout_arrivals_error": None,
        "plot_errors": [],
        "_gamma_result": gamma_payload,
        "_settings_object": settings,
    }


def event_table_rows(results) -> list[Dict[str, object]]:
    """Return concise per-event rows for the GUI table or headless tests."""

    rows = []
    for result in results:
        mode = detector_run_mode(getattr(result, "detector_mode", "cherenkov"))
        records = result.photon_records
        crossings = [record for record in records if record.detector_transmitted]
        detected = [record for record in records if record.detected]
        scint = result.scintillator.summary() if result.scintillator is not None else {}
        arrival_times = [record.arrival_time_s * 1e9 for record in detected if not math.isnan(record.arrival_time_s)]
        rows.append({
            "event_id": result.event_id,
            "sampling_mode": result.optical.sampling_mode,
            "detector_mode": mode,
            "requested_photons": result.optical.photon_count if mode in {"cherenkov", "combined"} else 0,
            "generated_photons": result.generated_photons,
            "transported_photons": len(records),
            "detector_face_crossings": len(crossings),
            "valid_sensor_hits": sum(record.detector_valid_sensor_hit for record in records),
            "detected_photons": len(detected),
            "detected_plus_x": sum(record.detector_face == "plus_x" for record in detected),
            "detected_minus_x": sum(record.detector_face == "minus_x" for record in detected),
            "mean_detected_arrival_ns": (sum(arrival_times) / len(arrival_times) if arrival_times else None),
            "physical_expected_photons": result.physical_expected_photons if mode in {"cherenkov", "combined"} else 0.0,
            "poisson_mean": result.physical_expected_photons if mode in {"cherenkov", "combined"} and result.optical.sampling_mode == "poisson" else None,
            "generation_fraction": result.summary()["generation_fraction"],
            "detector_attempts": sum(record.detector_attempts for record in records),
            "tir_plus_x": sum(record.tir_plus_x for record in records),
            "tir_minus_x": sum(record.tir_minus_x for record in records),
            "detector_interface_reflections": sum(record.detector_interface_reflections for record in records),
            "detector_angle_rejected": sum(record.detector_angle_rejected for record in records),
            "detector_active_area_rejected": sum(record.detector_active_area_rejected for record in records),
            "detector_collection_loss": sum(record.detector_collection_loss for record in records),
            "detector_qe_loss": sum(record.detector_qe_loss for record in records),
            "detector_losses": sum(record.detector_angle_rejected + record.detector_active_area_rejected + record.detector_collection_loss + record.detector_qe_loss for record in records),
            "outcome_reasons": dict(sorted(simulator.detector_outcome_reason_counts(records).items())),
            **scint,
        })
    return rows


def event_table_rows_from_summary(summary: Mapping[str, object]) -> list[Dict[str, object]]:
    """Convert a disk summary into table rows when no in-memory run exists."""

    rows = []
    for item in summary.get("event_summaries", []) if isinstance(summary, Mapping) else []:
        if not isinstance(item, Mapping):
            continue
        rows.append({
            "event_id": item.get("event_id", len(rows)),
            "detector_mode": item.get("detector_mode", summary.get("detector_mode", "cherenkov")),
            "sampling_mode": item.get("sampling_mode", "fixed"),
            "requested_photons": item.get("requested_photons", 0),
            "generated_photons": item.get("generated_photons", 0),
            "transported_photons": item.get("transported_photons", 0),
            "detector_face_crossings": item.get("detector_face_crossings", 0),
            "valid_sensor_hits": item.get("valid_sensor_hits", 0),
            "detected_photons": item.get("detected_photons", 0),
            "detected_plus_x": item.get("detected_plus_x", 0),
            "detected_minus_x": item.get("detected_minus_x", 0),
            "mean_detected_arrival_ns": item.get("mean_detected_time_ns"),
            "physical_expected_photons": item.get("physical_expected_photons", 0.0),
            "poisson_mean": item.get("poisson_mean"),
            "generation_fraction": item.get("generation_fraction", 0.0),
            "detector_attempts": item.get("detector_attempts", 0),
            "tir_plus_x": item.get("tir_plus_x", 0),
            "tir_minus_x": item.get("tir_minus_x", 0),
            "detector_interface_reflections": item.get("detector_interface_reflections", 0),
            "detector_angle_rejected": item.get("detector_angle_rejected", 0),
            "detector_active_area_rejected": item.get("detector_active_area_rejected", 0),
            "detector_collection_loss": item.get("detector_collection_loss", 0),
            "detector_qe_loss": item.get("detector_qe_loss", 0),
            "detector_losses": item.get("detector_losses", 0),
            "outcome_reasons": item.get("outcome_reasons", {}),
            **{key: item.get(key, 0) for key in ("scintillator_requested_photons", "scintillator_generated_photons", "scintillator_transported_photons", "scintillator_source_rows", "scintillator_accepted_records", "scintillator_excluded_records", "scintillator_photons_reaching_guide", "scintillator_photons_reaching_readout", "scintillator_detected_photons", "scintillator_mylar_reflections", "scintillator_mylar_losses", "scintillator_readout_alignment_mismatches", "scintillator_active_to_guide_crossings", "scintillator_guide_to_readout_crossings", "scintillator_interface_transmissions", "scintillator_interface_reflections", "scintillator_tir_reflections", "scintillator_readout_endpoint_tolerance_um", "scintillator_mean_arrival_ns", "scintillator_earliest_arrival_ns")},
        })
    return rows


def format_event_details(result, display_photons: Optional[int] = None) -> str:
    """Format a selected event's geometry, counts, and timing details."""

    mode = detector_run_mode(getattr(result, "detector_mode", "cherenkov"))
    geometry = result.geometry
    entry = [result.particle.entry_x_m, result.particle.entry_y_m, getattr(geometry, "half_z", 0.0)]
    exit_point = [entry[index] + result.track_length_m * result.particle.direction[index] for index in range(3)]
    records = result.photon_records
    crossings = [record for record in records if record.detector_transmitted]
    detected = [record for record in records if record.detected]
    arrivals = [record.arrival_time_s * 1e9 for record in detected if not math.isnan(record.arrival_time_s)]
    tir_plus_x = sum(record.tir_plus_x for record in records)
    tir_minus_x = sum(record.tir_minus_x for record in records)
    interface_reflections = sum(record.detector_interface_reflections for record in records)
    angle_losses = sum(record.detector_angle_rejected for record in records)
    collection_losses = sum(record.detector_collection_loss for record in records)
    qe_losses = sum(record.detector_qe_loss for record in records)
    particle_time_ns = result.track_length_m / max(result.beta_start * simulator.C_LIGHT, 1e-30) * 1e9
    lines = [
        f"Event {result.event_id} — detector mode: {detector_run_mode_label(mode)} ({mode})",
        f"Entry point: ({entry[0] * 1e3:.3f}, {entry[1] * 1e3:.3f}, {entry[2] * 1e3:.3f}) mm",
        f"Exit point: ({exit_point[0] * 1e3:.3f}, {exit_point[1] * 1e3:.3f}, {exit_point[2] * 1e3:.3f}) mm",
        f"Particle track length: {result.track_length_m * 1e3:.3f} mm",
        f"Particle path time: {particle_time_ns:.4f} ns",
        f"Beta start/end: {result.beta_start:.6f} / {result.beta_end:.6f}",
    ]
    if mode in {"cherenkov", "combined"}:
        lines.extend([
            f"+x coupling index / acceptance: {result.optical.detector_index:.3f} / {math.degrees(result.optical.detector_acceptance_rad):.1f} deg",
            f"Cherenkov sampling mode: {SAMPLING_MODE_LABELS.get(result.optical.sampling_mode, result.optical.sampling_mode)}",
            f"Cherenkov requested/generated/transported: {result.optical.photon_count} / {result.generated_photons} / {len(records)}",
            f"Cherenkov physical expected yield: {result.physical_expected_photons:.3g} photons/event",
            f"Detector attempts: {sum(record.detector_attempts for record in records)}; face crossings: {len(crossings)}; valid sensor hits: {sum(record.detector_valid_sensor_hit for record in records)}; detected: {len(detected)}",
            f"TIR +x / -x: {tir_plus_x} / {tir_minus_x}; interface reflections: {interface_reflections}",
            f"Detector losses angle / collection / QE: {angle_losses} / {collection_losses} / {qe_losses}",
            f"Outcome reasons: {dict(sorted(simulator.detector_outcome_reason_counts(records).items()))}",
            f"Effective sensor response: {result.summary().get('effective_detector_response', {})}",
            f"Detected +x / -x: {sum(record.detector_face == 'plus_x' for record in detected)} / {sum(record.detector_face == 'minus_x' for record in detected)}",
            f"Mean detected arrival: {(sum(arrivals) / len(arrivals)) if arrivals else 'n/a'} ns",
        ])
    if result.scintillator is not None:
        scint = result.scintillator.summary()
        lines.extend([
            "",
            "Scintillator / light guide",
            f"Requested/generated/transported: {scint['scintillator_requested_photons']} / {scint['scintillator_generated_photons']} / {scint['scintillator_transported_photons']}",
            f"Source rows/accepted/excluded: {scint.get('scintillator_source_rows', 0)} / {scint.get('scintillator_accepted_records', 0)} / {scint.get('scintillator_excluded_records', 0)}",
            f"Guide/readout/detected: {scint['scintillator_photons_reaching_guide']} / {scint['scintillator_photons_reaching_readout']} / {scint['scintillator_detected_photons']}",
            f"Active→guide / guide→readout crossings: {scint.get('scintillator_active_to_guide_crossings', 0)} / {scint.get('scintillator_guide_to_readout_crossings', 0)}",
            f"Interface transmit/reflection/TIR: {scint.get('scintillator_interface_transmissions', 0)} / {scint.get('scintillator_interface_reflections', 0)} / {scint.get('scintillator_tir_reflections', 0)}",
            f"Readout endpoint mismatches: {scint.get('scintillator_readout_alignment_mismatches', 0)} (tolerance {scint.get('scintillator_readout_endpoint_tolerance_um', 0.0):g} µm)",
            f"Mylar reflections/losses: {scint['scintillator_mylar_reflections']} / {scint['scintillator_mylar_losses']}",
            f"Absorbed/escaped/max-bounces: {scint['scintillator_absorbed']} / {scint['scintillator_escaped']} / {scint['scintillator_max_bounces']}",
            f"Mean/earliest arrival: {scint['scintillator_mean_arrival_ns'] or 'n/a'} / {scint['scintillator_earliest_arrival_ns'] or 'n/a'} ns",
        ])
    if tir_plus_x > 0 and not crossings and result.optical.detector_index <= 1.0:
        lines.append("Diagnostic: +x photons are undergoing quartz-air total internal reflection; QE is applied only after transmission.")
    if display_photons is not None:
        lines.append(f"Photon paths displayed: {min(max(display_photons, 0), len(records))}")
    for face in ("plus_x", "minus_x"):
        hits = [record for record in detected if record.detector_face == face]
        if hits:
            mean_y = sum(record.impact_position_m[1] for record in hits) / len(hits) * 1e3
            mean_z = sum(record.impact_position_m[2] for record in hits) / len(hits) * 1e3
            lines.append(f"{face} mean hit y/z: {mean_y:.3f} / {mean_z:.3f} mm ({len(hits)} hits)")
    return "\n".join(lines)


def latest_run_directory(base_directory: Path) -> Optional[Path]:
    """Find the newest GUI run directory, if one exists."""

    base = Path(base_directory).expanduser().resolve()
    candidates = [path for path in base.glob("quartz_run_*") if path.is_dir()]
    return max(candidates, key=lambda path: path.stat().st_mtime) if candidates else None


LATEST_GAMMA_RESULT_INDEX = "latest_gamma_result.json"
LATEST_GAMMA_RESULT_SCHEMA = "quartz-latest-gamma-result-v1"


def retain_latest_completed_gamma_payload(existing: Optional[Mapping[str, object]], run_payload: Mapping[str, object]) -> Optional[Mapping[str, object]]:
    """Replace cached Gamma data only after a complete Gamma transport run."""

    candidate = run_payload.get("_gamma_result")
    if isinstance(candidate, Mapping) and not bool(run_payload.get("stopped")):
        return candidate
    return existing


def persist_latest_gamma_result_index(output_directory: Path, run_directory: Path, gamma_results_path: Path) -> Path:
    """Atomically point to the last completed Gamma result, independent of optical runs."""

    root = Path(output_directory).expanduser().resolve()
    run = Path(run_directory).expanduser().resolve()
    artifact = Path(gamma_results_path).expanduser().resolve()
    if not artifact.is_file():
        raise FileNotFoundError(f"Gamma result artifact does not exist: {artifact}")
    root.mkdir(parents=True, exist_ok=True)
    try:
        stored_artifact = str(artifact.relative_to(root))
    except ValueError:
        stored_artifact = str(artifact)
    try:
        stored_run = str(run.relative_to(root))
    except ValueError:
        stored_run = str(run)
    index_path = root / LATEST_GAMMA_RESULT_INDEX
    temporary = index_path.with_name(index_path.name + ".tmp")
    temporary.write_text(json.dumps({
        "schema": LATEST_GAMMA_RESULT_SCHEMA,
        "gamma_results_path": stored_artifact,
        "run_directory": stored_run,
        "updated_utc": datetime.now().astimezone().isoformat(),
    }, indent=2, sort_keys=True) + "\n", encoding="utf-8")
    os.replace(temporary, index_path)
    return index_path


def _read_gamma_result_artifact(path: Path) -> Optional[Mapping[str, object]]:
    try:
        payload = json.loads(Path(path).read_text(encoding="utf-8"))
    except (OSError, json.JSONDecodeError, TypeError, ValueError):
        return None
    if isinstance(payload, Mapping) and payload.get("model_status") and isinstance(payload.get("events", ()), list):
        return payload
    return None


def load_latest_gamma_result(output_directory: Path) -> tuple[Optional[Mapping[str, object]], Optional[Path], Optional[Path]]:
    """Restore the latest Gamma payload without coupling it to the newest run folder.

    The index is authoritative for v71 and later.  A reverse-chronological scan
    migrates complete v70 Gamma run folders that predate the independent index.
    """

    root = Path(output_directory).expanduser().resolve()
    index_path = root / LATEST_GAMMA_RESULT_INDEX
    try:
        index = json.loads(index_path.read_text(encoding="utf-8"))
    except (OSError, json.JSONDecodeError, TypeError, ValueError):
        index = {}
    if isinstance(index, Mapping) and index.get("schema") == LATEST_GAMMA_RESULT_SCHEMA:
        artifact = Path(str(index.get("gamma_results_path", ""))).expanduser()
        if not artifact.is_absolute():
            artifact = root / artifact
        payload = _read_gamma_result_artifact(artifact)
        if payload is not None:
            run = Path(str(index.get("run_directory", artifact.parent))).expanduser()
            if not run.is_absolute():
                run = root / run
            return payload, run.resolve(), artifact.resolve()

    candidates = [path for path in root.glob("quartz_run_*") if path.is_dir()]
    for run in sorted(candidates, key=lambda path: path.stat().st_mtime, reverse=True):
        summary_path = run / "summary.json"
        try:
            summary = json.loads(summary_path.read_text(encoding="utf-8"))
        except (OSError, json.JSONDecodeError, TypeError, ValueError):
            continue
        if not isinstance(summary, Mapping):
            continue
        raw_path = summary.get("gamma_results_path")
        if not raw_path and isinstance(summary.get("summary"), Mapping):
            if summary["summary"].get("detector_mode") == GAMMA_RUN_MODE:
                raw_path = "gamma_results.json"
        if not raw_path:
            continue
        artifact = Path(str(raw_path)).expanduser()
        if not artifact.is_absolute():
            artifact = run / artifact
        payload = _read_gamma_result_artifact(artifact)
        if payload is not None:
            return payload, run.resolve(), artifact.resolve()
    return None, None, None


def _gamma_number(value, default: float = 0.0) -> float:
    try:
        result = float(value)
    except (TypeError, ValueError):
        return float(default)
    return result if math.isfinite(result) else float(default)


def gamma_event_table_rows(events: Sequence[object]) -> list[Dict[str, object]]:
    """Summarize persisted Gamma events for the selectable per-event table."""

    rows = []
    for index, event in enumerate(events or ()):
        if not isinstance(event, Mapping):
            continue
        deposits = [item for item in (event.get("deposits", ()) or ()) if isinstance(item, Mapping)]
        tracks = [item for item in (event.get("charged_tracks", ()) or ()) if isinstance(item, Mapping)]
        raw_counts = event.get("interaction_counts", {})
        counts: Dict[str, int] = {}
        if isinstance(raw_counts, Mapping):
            for name, count in raw_counts.items():
                parsed = max(0, int(_gamma_number(count)))
                if parsed:
                    counts[str(name)] = parsed
        if not counts and deposits:
            for item in deposits:
                name = str(item.get("interaction", "interaction"))
                counts[name] = counts.get(name, 0) + 1
        interaction_count = sum(counts.values())
        processes = ", ".join(
            f"{name} ×{count}" for name, count in sorted(counts.items())
        ) or "none"
        deposit_total = sum(_gamma_number(item.get("energy_keV", 0.0)) for item in deposits)
        has_light = isinstance(event.get("light"), (list, tuple))
        light_rows = [item for item in (event.get("light", ()) or ()) if isinstance(item, Mapping)]
        light_by_mechanism: Dict[str, float] = {}
        for item in light_rows:
            mechanism = str(item.get("mechanism", "local"))
            light_by_mechanism[mechanism] = light_by_mechanism.get(mechanism, 0.0) + _gamma_number(item.get("expected_photons", 0.0))
        light_total = sum(light_by_mechanism.values()) if has_light else None
        light_summary = (
            "; ".join(f"{name} {value:.4g}" for name, value in sorted(light_by_mechanism.items()))
            if has_light and light_by_mechanism else "0 expected" if has_light else "n/a"
        )
        rows.append({
            "event_id": int(_gamma_number(event.get("event_id", index), index)),
            "source_energy_keV": _gamma_number(event.get("source_energy_keV", 0.0)),
            "source_component": str(event.get("source_component", "unknown")),
            "interaction_count": interaction_count,
            "interaction_process_summary": processes,
            "deposit_total_keV": deposit_total,
            "deposit_record_count": len(deposits),
            "charged_track_count": len(tracks),
            "local_light_expected_photons": light_total,
            "local_light_summary": light_summary,
            "status": "INTERACTION" if interaction_count or deposits or tracks else "EMPTY",
        })
    return rows


def gamma_geometry_view_data(geometry_snapshot: Optional[Mapping[str, object]]) -> Dict[str, object]:
    """Read renderable solids only from a saved v72 run geometry snapshot."""

    empty = {"assembly_types": [], "components": []}
    if not isinstance(geometry_snapshot, Mapping):
        return empty
    if geometry_snapshot.get("schema") != "gamma-global-assembly-geometry-v1":
        return empty
    scene = geometry_snapshot.get("scene")
    if not isinstance(scene, Mapping):
        return empty
    assemblies = scene.get("assemblies", ())
    components = []
    assembly_types = []
    for assembly in assemblies if isinstance(assemblies, (list, tuple)) else ():
        if not isinstance(assembly, Mapping):
            continue
        assembly_type = str(assembly.get("assembly_type", ""))
        if assembly_type not in {"cherenkov", "scintillator"}:
            continue
        assembly_types.append(assembly_type)
        for item in assembly.get("components", ()) or ():
            if not isinstance(item, Mapping):
                continue
            faces = []
            for face in item.get("faces_mm", ()) or ():
                try:
                    points = [tuple(float(coordinate) for coordinate in point[:3]) for point in face]
                except (TypeError, ValueError):
                    continue
                if len(points) >= 3 and all(len(point) == 3 for point in points):
                    faces.append(points)
            if faces:
                components.append({
                    "component_id": str(item.get("component_id", item.get("name", "component"))),
                    "assembly_type": assembly_type,
                    "kind": str(item.get("kind", "geometry")),
                    "name": str(item.get("name", "")),
                    "color": str(item.get("color", "#8db6ce")),
                    "faces_mm": faces,
                })
    if sorted(set(assembly_types)) != ["cherenkov", "scintillator"]:
        return empty
    if {str(item["assembly_type"]) for item in components} != {"cherenkov", "scintillator"}:
        return empty
    return {"assembly_types": ["cherenkov", "scintillator"], "components": components}


def gamma_event_view_data(event: Mapping[str, object]) -> Dict[str, list]:
    """Convert one persisted Gamma event into global-mm drawable segments/points."""

    def vector(values, default=(0.0, 0.0, 0.0)):
        try:
            result = tuple(float(value) for value in values)
        except (TypeError, ValueError):
            return tuple(default)
        return result if len(result) == 3 else tuple(default)

    gamma_segments = []
    for item in event.get("gamma_tracks", ()) or ():
        if not isinstance(item, Mapping):
            continue
        start = tuple(value * 1000.0 for value in vector(item.get("start_position_m")))
        end = tuple(value * 1000.0 for value in vector(item.get("end_position_m")))
        gamma_segments.append({"start_mm": start, "end_mm": end, "volume": str(item.get("volume", "")), "interaction": str(item.get("interaction", "transmitted")), "energy_keV": float(item.get("energy_keV", 0.0))})
    charged_segments = []
    for item in event.get("charged_tracks", ()) or ():
        if not isinstance(item, Mapping):
            continue
        detailed_segments = item.get("segments", ())
        if isinstance(detailed_segments, (list, tuple)) and detailed_segments:
            for segment in detailed_segments:
                if not isinstance(segment, Mapping):
                    continue
                start = tuple(value * 1000.0 for value in vector(segment.get("start_position_m")))
                end = tuple(value * 1000.0 for value in vector(segment.get("end_position_m")))
                charged_segments.append({
                    "start_mm": start,
                    "end_mm": end,
                    "kind": str(item.get("kind", "charged")),
                    "volume": str(segment.get("volume", item.get("volume", ""))),
                    "path_mm": max(0.0, _gamma_number(segment.get("path_mm", 0.0))),
                    "deposited_energy_keV": max(0.0, _gamma_number(segment.get("deposited_energy_keV", 0.0))),
                })
            continue
        start_m = vector(item.get("position_m"))
        direction = vector(item.get("direction"))
        length_m = max(0.0, float(item.get("path_mm", 0.0))) / 1000.0
        end_m = tuple(start_m[index] + direction[index] * length_m for index in range(3))
        charged_segments.append({"start_mm": tuple(value * 1000.0 for value in start_m), "end_mm": tuple(value * 1000.0 for value in end_m), "kind": str(item.get("kind", "charged")), "volume": str(item.get("volume", "")), "path_mm": length_m * 1000.0})
    deposits = []
    for item in event.get("deposits", ()) or ():
        if not isinstance(item, Mapping):
            continue
        deposits.append({"position_mm": tuple(value * 1000.0 for value in vector(item.get("position_m"))), "interaction": str(item.get("interaction", "interaction")), "energy_keV": float(item.get("energy_keV", 0.0)), "volume": str(item.get("volume", ""))})
    return {"gamma_segments": gamma_segments, "charged_segments": charged_segments, "deposits": deposits}


def gamma_global_view_binding_sequences() -> Tuple[str, ...]:
    """Headless contract for the detached Gamma viewer's camera controls."""

    return ("<Configure>", "<ButtonPress-1>", "<B1-Motion>", "<MouseWheel>", "<Button-4>", "<Button-5>")


def gamma_global_view_control_names() -> Tuple[str, ...]:
    """Visible controls exposed by the independent global Gamma window."""

    return ("Event", "Reset view", "Show geometry", "Show source", "Show tracks / deposits")


def build_gamma_global_view_scene(
    event: Optional[Mapping[str, object]],
    geometry_snapshot: Optional[Mapping[str, object]],
    source_direction_global=(0.0, 0.0, 1.0),
) -> Dict[str, object]:
    """Build renderer-neutral global-mm data for the detached Gamma viewer."""

    geometry = gamma_geometry_view_data(geometry_snapshot)
    event = event if isinstance(event, Mapping) else {}
    tracks = gamma_event_view_data(event)
    components = geometry["components"]
    bounds_points = [
        point
        for component in components
        for face in component["faces_mm"]
        for point in face
    ]
    for key in ("gamma_segments", "charged_segments"):
        for segment in tracks[key]:
            bounds_points.extend((segment["start_mm"], segment["end_mm"]))
    bounds_points.extend(item["position_mm"] for item in tracks["deposits"])

    try:
        direction = [float(value) for value in source_direction_global[:3]]
    except (TypeError, ValueError, IndexError):
        direction = [0.0, 0.0, 1.0]
    if len(direction) != 3 or not all(math.isfinite(value) for value in direction):
        direction = [0.0, 0.0, 1.0]
    magnitude = math.sqrt(sum(value * value for value in direction))
    direction = [value / magnitude for value in direction] if magnitude > 1e-12 else [0.0, 0.0, 1.0]

    source_point = None
    for key, scale in (("source_position_global_mm", 1.0), ("source_position_mm", 1.0),
                       ("source_position_global_m", 1000.0), ("source_position_m", 1000.0)):
        value = event.get(key)
        if isinstance(value, (list, tuple)) and len(value) >= 3:
            try:
                candidate = [float(item) * scale for item in value[:3]]
            except (TypeError, ValueError):
                continue
            if all(math.isfinite(item) for item in candidate):
                source_point = candidate
                break
    source_inferred = False
    if source_point is None and tracks["gamma_segments"]:
        source_point = list(tracks["gamma_segments"][0]["start_mm"])
    if source_point is None and bounds_points and components:
        # The v76 Gamma transport starts its +Z beam 1 mm before the lowest
        # active/guide volume. Those component kinds match the saved transport
        # geometry while excluding sensors and other non-transport solids.
        transport_points = [
            point for component in components if component["kind"] in {"geometry", "guide"}
            for face in component["faces_mm"] for point in face
        ] or bounds_points
        low_z = min(float(point[2]) for point in transport_points)
        source_point = [0.0, 0.0, low_z - 1.0]
        source_inferred = True
    if source_point is not None:
        bounds_points.append(source_point)

    if bounds_points:
        low = [min(float(point[axis]) for point in bounds_points) for axis in range(3)]
        high = [max(float(point[axis]) for point in bounds_points) for axis in range(3)]
        spans = [max(high[axis] - low[axis], 1.0) for axis in range(3)]
        center = [(low[axis] + high[axis]) / 2.0 for axis in range(3)]
        extent = max(spans)
    else:
        low, high, center, extent = [-1.0] * 3, [1.0] * 3, [0.0] * 3, 2.0
    ray_length = max(extent * 0.12, 1.0)
    source_ray_end = (
        [source_point[axis] + direction[axis] * ray_length for axis in range(3)]
        if source_point is not None else None
    )
    if source_ray_end is not None:
        bounds_points.append(source_ray_end)
        low = [min(low[axis], source_ray_end[axis]) for axis in range(3)]
        high = [max(high[axis], source_ray_end[axis]) for axis in range(3)]
        spans = [max(high[axis] - low[axis], 1.0) for axis in range(3)]
        center = [(low[axis] + high[axis]) / 2.0 for axis in range(3)]
        extent = max(spans)
    rows = gamma_event_table_rows((event,)) if event else []
    return {
        "event_id": int(_gamma_number(event.get("event_id", 0), 0)),
        "status": rows[0]["status"] if rows else "NO EVENT",
        "source_energy_keV": _gamma_number(event.get("source_energy_keV", 0.0)),
        "source_component": str(event.get("source_component", "unknown")),
        "source_position_mm": source_point,
        "source_position_inferred": source_inferred,
        "source_direction_global": direction,
        "source_ray_end_mm": source_ray_end,
        "geometry_components": components,
        "geometry_available": bool(components),
        **tracks,
        "bounds_mm": [*low, *high],
        "center_mm": center,
        "extent_mm": max(extent, 1.0),
    }


def draw_gamma_event_view(
    figure, event: Optional[Mapping[str, object]],
    geometry_snapshot: Optional[Mapping[str, object]] = None,
):
    """Draw the saved complete detector scene, then event data in global mm."""

    figure.clear()
    axis = figure.add_subplot(111, projection="3d")
    axis.set_xlabel("Global X (mm)")
    axis.set_ylabel("Global Y (mm)")
    axis.set_zlabel("Global Z (mm)")
    geometry_scene = gamma_geometry_view_data(geometry_snapshot)
    bounds = []
    if geometry_scene["components"]:
        from mpl_toolkits.mplot3d.art3d import Poly3DCollection
        for component in geometry_scene["components"]:
            faces = component["faces_mm"]
            bounds.extend(point for face in faces for point in face)
            collection = Poly3DCollection(
                faces, facecolors=component["color"], edgecolors="#45515a",
                linewidths=0.45, alpha=0.24,
            )
            collection.set_label("_nolegend_")
            axis.add_collection3d(collection)
    else:
        axis.text2D(
            0.02, 0.02,
            "Detector geometry snapshot unavailable for this saved Gamma result",
            transform=axis.transAxes, ha="left", va="bottom", fontsize=8,
            color="#8a3a20",
        )
    if not isinstance(event, Mapping):
        axis.text2D(0.5, 0.5, "No completed Gamma event is selected", transform=axis.transAxes, ha="center", va="center")
        axis.set_title("Gamma event viewer")
        if bounds:
            points = np.asarray(bounds, dtype=float)
            low, high = points.min(axis=0), points.max(axis=0)
            spans = np.maximum(high - low, 1.0)
            center = (high + low) / 2.0
            radius = max(float(spans.max()) * 0.55, 1.0)
            axis.set_xlim(center[0] - radius, center[0] + radius)
            axis.set_ylim(center[1] - radius, center[1] + radius)
            axis.set_zlim(center[2] - radius, center[2] + radius)
        return axis
    scene = gamma_event_view_data(event)
    for segment in scene["gamma_segments"]:
        start, end = segment["start_mm"], segment["end_mm"]
        bounds.extend((start, end))
        color = "#414141" if segment["interaction"] == "transmitted" else "#111111"
        axis.plot((start[0], end[0]), (start[1], end[1]), (start[2], end[2]), color=color, linewidth=1.5, alpha=0.85, linestyle="--" if segment["interaction"] == "transmitted" else "-")
    kind_colors = {"electron": "#147d64", "positron": "#8246af"}
    for segment in scene["charged_segments"]:
        start, end = segment["start_mm"], segment["end_mm"]
        bounds.extend((start, end))
        axis.plot((start[0], end[0]), (start[1], end[1]), (start[2], end[2]), color=kind_colors.get(segment["kind"], "#147d64"), linewidth=1.7, alpha=0.9)
    interaction_colors = {"compton": "#e69f00", "photoelectric": "#d62728", "pair": "#8246af"}
    for interaction in sorted({item["interaction"] for item in scene["deposits"]}):
        points = [item["position_mm"] for item in scene["deposits"] if item["interaction"] == interaction]
        bounds.extend(points)
        axis.scatter([point[0] for point in points], [point[1] for point in points], [point[2] for point in points], color=interaction_colors.get(interaction, "#e69f00"), s=34, label=f"{interaction} deposit")
    event_id = event.get("event_id", "?")
    event_rows = gamma_event_table_rows((event,))
    event_status = event_rows[0]["status"] if event_rows else "UNKNOWN"
    axis.set_title(f"Gamma event {event_id}: {event_status} — photon paths, charged secondaries, deposits")
    if bounds:
        points = np.asarray(bounds, dtype=float)
        low, high = points.min(axis=0), points.max(axis=0)
        spans = np.maximum(high - low, 1.0)
        center = (high + low) / 2.0
        radius = max(float(spans.max()) * 0.55, 1.0)
        axis.set_xlim(center[0] - radius, center[0] + radius)
        axis.set_ylim(center[1] - radius, center[1] + radius)
        axis.set_zlim(center[2] - radius, center[2] + radius)
        axis.set_box_aspect(spans)
    elif scene["gamma_segments"] or scene["charged_segments"] or scene["deposits"]:
        axis.text2D(0.5, 0.5, "No interactions or visible tracks in this event", transform=axis.transAxes, ha="center", va="center")
    axis.view_init(elev=22, azim=-55)
    axis.legend(loc="upper left", fontsize=8) if scene["deposits"] else None
    return axis


def summary_with_diagnostic_defaults(summary: Mapping[str, object]) -> Dict[str, object]:
    """Keep older summary.json files readable after diagnostic fields were added."""

    result = dict(summary)
    result.setdefault("sampling_mode", "fixed")
    result.setdefault("physical_expected_photons", 0.0)
    result.setdefault("poisson_mean", None)
    for key in (
        "detector_attempts", "tir_plus_x", "tir_minus_x", "detector_interface_reflections",
        "detector_angle_rejected", "detector_collection_loss", "detector_qe_loss", "detector_losses",
        "detector_active_area_rejected", "valid_sensor_hits",
    ):
        result.setdefault(key, 0)
    return result


def create_selected_event_plot(
    results,
    event_index: int,
    output_directory: Path,
    display_photons: int = 40,
    track_margin_mm: float = 20.0,
    selected_categories=None,
    canonical_scene: Optional[Mapping[str, object]] = None,
) -> Path:
    """Headless-friendly helper used by the GUI and automated checks."""

    if event_index < 0 or event_index >= len(results):
        raise IndexError(f"event index must be in [0, {len(results) - 1}]")
    output_directory = Path(output_directory).expanduser().resolve()
    output_directory.mkdir(parents=True, exist_ok=True)
    output_path = output_directory / f"event_display_event_{event_index}_selected.png"
    result = results[event_index]
    if selected_categories is None:
        selected_categories = {"photon_paths"}
    filtered = filter_photon_records(result.photon_records, selected_categories)
    filtered_scintillator = None
    if getattr(result, "scintillator", None) is not None:
        filtered_scintillator = replace(
            result.scintillator,
            photon_records=filter_photon_records(result.scintillator.photon_records, selected_categories),
        )
    # SimulationResult is a dataclass.  Replacing both detector-owned record
    # collections preserves mode metadata while preventing a static plot from
    # drawing a disabled category.
    display_result = replace(result, photon_records=filtered, scintillator=filtered_scintillator) if is_dataclass(result) else copy(result)
    if not is_dataclass(result):
        display_result.photon_records = filtered
        if getattr(display_result, "scintillator", None) is not None:
            display_result.scintillator.photon_records = filtered_scintillator.photon_records
    if not filtered and not (filtered_scintillator and filtered_scintillator.photon_records):
        if selected_categories:
            message = "No photon paths match the selected viewer categories."
        else:
            message = "No photon categories selected; solid and particle track remain visible."
    else:
        message = None
    if canonical_scene is not None:
        # The static image consumes the exact renderer-neutral scene used by
        # the interactive event viewer.  In particular this retains the
        # complete guide section faces, sensor/readout components, global
        # placement, detector-owned paths, and mode filtering.
        event_scene = build_interactive_event_scene(
            display_result,
            canonical_scene,
            detector_mode=getattr(result, "detector_mode", "auto"),
            selected_categories=selected_categories,
            display_photons=max(0, int(display_photons)),
            track_margin_mm=max(0.0, float(track_margin_mm)),
        )
        simulator.plot_event_display(
            display_result, str(output_path), max(0, int(display_photons)), False,
            max(0.0, float(track_margin_mm)) * 1e-3, message,
            event_scene=event_scene,
        )
    else:
        # Compatibility for callers that only have a bare SimulationResult.
        # GUI runs always provide canonical_scene, so this fallback is not
        # used for saved run products or the selected-event button.
        simulator.plot_event_display(
            display_result, str(output_path), max(0, int(display_photons)), False,
            max(0.0, float(track_margin_mm)) * 1e-3, message,
        )
    return output_path


def _plot_canonical_event_display(
    result, output_path: str, canonical_scene: Mapping[str, object],
    display_photons: int, track_margin_mm: float,
) -> None:
    """Save the run's event PNG from the same scene as its interactive view."""

    event_scene = build_interactive_event_scene(
        result, canonical_scene,
        detector_mode=getattr(result, "detector_mode", "auto"),
        selected_categories={"photon_paths"},
        display_photons=max(0, int(display_photons)),
        track_margin_mm=max(0.0, float(track_margin_mm)),
    )
    simulator.plot_event_display(
        result, output_path, max(0, int(display_photons)), False,
        max(0.0, float(track_margin_mm)) * 1e-3,
        event_scene=event_scene,
    )


def _canonical_scene_values_from_settings(settings: RunSettings) -> Dict[str, str]:
    """Reconstruct the scene-value map needed by the canonical editor model."""

    values = default_values()
    typed = asdict(settings)
    for key in values:
        if key in typed:
            value = typed[key]
            values[key] = "1" if value is True else "0" if value is False else str(value)
    values.update({str(key): str(value) for key, value in (settings.sensor_profile_values or {}).items()})
    return values


def run_simulation_to_directory(
    settings: RunSettings,
    run_directory: Path,
    log: Optional[Callable[[str], None]] = None,
    geometry=None,
    progress_callback: Optional[Callable[[int, int, float], None]] = None,
    status_callback: Optional[Callable[[str], None]] = None,
    stop_event=None,
    assembly_models: Optional[Mapping[str, AssemblyGeometryModel]] = None,
) -> Dict[str, object]:
    """Run simulation and save ROOT, JSON, and all available PNG products."""

    external_log = log or (lambda message: None)
    status_callback = status_callback or (lambda message: None)
    run_directory = Path(run_directory).expanduser().resolve()
    run_directory.mkdir(parents=True, exist_ok=True)
    run_log_path = run_directory / "run.log"

    def log(message: object) -> None:
        """Fan one message to both the flushed disk log and the GUI queue."""

        text = str(message)
        _append_persistent_run_log(run_log_path, text)
        try:
            external_log(text)
        except Exception:
            # A closed Tk queue/widget must not prevent saved output or the
            # worker's terminal error route.
            pass

    def phase(name: str) -> None:
        """Publish and persist a major phase before doing synchronous work."""

        text = str(name)
        log(f"PHASE: {text}")
        try:
            status_callback(text)
        except Exception as exc:
            # A GUI/status producer is advisory.  Keep the simulation alive,
            # but leave a diagnostic in the persistent run log so callback
            # failures cannot masquerade as a frozen worker.
            log(f"STATUS CALLBACK ERROR: {type(exc).__name__}: {exc}")

    mode = detector_run_mode(settings.detector_run_mode)
    phase(f"Building {detector_run_mode_label(mode)} assemblies")
    geometry = geometry or simulator.BarGeometry()
    if mode == GAMMA_RUN_MODE:
        phase("Starting Gamma transport")
        gamma_payload = run_gamma_transport_to_directory(
            settings, run_directory, geometry, log=log,
            progress_callback=progress_callback, status_callback=status_callback,
            stop_event=stop_event, assembly_models=assembly_models,
        )
        summary_payload = {key: value for key, value in gamma_payload.items() if not key.startswith("_")}
        summary_text = json.dumps(_json_safe(summary_payload), allow_nan=False, indent=2, sort_keys=True) + "\n"
        (run_directory / "summary.json").write_text(summary_text, encoding="utf-8")
        log(f"Saved summary.json ({len(summary_text)} bytes)")
        if not bool(gamma_payload.get("stopped")):
            index_path = persist_latest_gamma_result_index(
                Path(settings.output_directory), run_directory, Path(str(gamma_payload["gamma_results_path"])),
            )
            log(f"Updated independent latest Gamma result index: {index_path}")
        phase("Updating GUI")
        return gamma_payload
    particle = _particle_for_settings(settings)
    optical = _optical_for_settings(settings)
    terminal_metadata = {item.get("face"): item.get("state") for item in getattr(geometry, "terminal_interfaces", ())}
    if terminal_metadata:
        optical = replace(
            optical,
            positive_end_mode=terminal_metadata.get("plus_x", "detector"),
            opposite_end_mode=terminal_metadata.get("minus_x", optical.opposite_end_mode),
        )
    effective_bounces = {
        "cherenkov": settings.maximum_photon_bounces if mode in {"cherenkov", "combined"} else None,
        "scintillator": settings.maximum_photon_bounces if mode in {"scintillator", "combined"} else None,
    }
    log(f"Effective maximum photon bounces: Cherenkov={effective_bounces['cherenkov']}, Scintillator={effective_bounces['scintillator']} (Run controls override)")
    settings_values = asdict(settings)
    scintillator_assembly = None
    scintillator_output_assembly = None
    scintillator_config = None
    two_detector_monitor = None
    if mode in {"scintillator", "combined"}:
        # The complete editor model is the optical authority.  The prior
        # value-field adapter rebuilt only one optional guide and therefore
        # could render more sections than the simulator transported.
        canonical_models, _ = _canonical_global_models(settings_values, geometry, assembly_models)
        scintillator_assembly = scintillator_assembly_from_model(
            canonical_models["scintillator"], values=settings_values, particle=particle,
        )
        scintillator_output_assembly = scintillator_assembly
        # The Run-tab control is the sole run-time bounce authority.  Legacy
        # scintillator_max_bounces remains in profiles for compatibility, but
        # cannot silently create a second transport limit.
        scintillator_config = scintillator_config_from_values(
            settings_values, max_bounces_override=settings.maximum_photon_bounces,
        )
        if mode == "combined":
            two_detector_monitor = two_detector_beam_monitor_from_values(settings_values, geometry)
    phase("Building global geometry and event view")
    scene_values = _canonical_scene_values_from_settings(settings)
    canonical_scene = global_assembly_scene(scene_values, geometry, assembly_models)
    if mode in {"scintillator", "combined"}:
        # The scintillator transport stays local, but its charged-particle
        # source is the intersection of the same frozen global line used by
        # the event viewer.  No assembly offset or rotation is mutated.
        scintillator_assembly = _align_scintillator_transport_to_global_track(
            scintillator_assembly, canonical_scene,
        )
        if canonical_scene.get("beam_alignment_mismatches"):
            message = "Run preflight: " + "; ".join(canonical_scene["beam_alignment_mismatches"])
            if mode == "combined":
                raise ValueError(message)
            log(message)
    event_view_scene = canonical_event_view_payload(canonical_scene, settings.detector_run_mode)
    if isinstance(geometry, simulator.CompoundGeometry):
        geometry_description = "CompoundGeometry(" + ", ".join(
            f"{index}:{section.shape}/{section.normalized_material()}"
            for index, section in enumerate(geometry.sections)
        ) + ")"
        log(f"Running {geometry_description} with internal Fresnel/Snell interfaces")
    else:
        log(f"Running {getattr(geometry, 'shape', 'box')} geometry with material {getattr(geometry, 'material', 'fused_silica')}")
    phase(f"Starting {detector_run_mode_label(mode)} event 1/{settings.events}")
    results = simulator.simulate_events(
        particle,
        settings.events,
        geometry=geometry,
        optical=optical,
        seed=settings.seed,
        progress_callback=progress_callback,
        status_callback=status_callback,
        stop_event=stop_event,
        scintillator_assembly=scintillator_assembly,
        scintillator_config=scintillator_config,
        detector_mode=mode,
    )
    stopped = bool(stop_event is not None and stop_event.is_set()) or any(getattr(result, "stopped", False) for result in results)
    if stopped:
        log(f"Stop requested; retained {len(results)} completed/partial event(s)")
    # ROOT is deliberately completed before the remaining finalization work.
    # This makes the line printed by write_root a hard boundary in the run
    # log and exposes every potentially slow post-ROOT operation to the GUI.
    root_path: Optional[str] = None
    root_error: Optional[str] = None
    if results:
        phase("Writing ROOT output")
        try:
            root_path = simulator.write_root(results, str(run_directory / "quartz_cherenkov.root"))
            log(f"ROOT output written to {root_path}")
            # write_root returns only after uproot's recreate/open contexts
            # have closed and the atomic replace has succeeded.
            phase("Closing ROOT output")
        except RuntimeError as exc:
            root_error = str(exc)
            log(f"ROOT output unavailable: {root_error}")
            phase("Closing ROOT output (ROOT unavailable)")

    phase("Aggregating results")
    summary = simulator.aggregate_summary(results, detector_mode=mode)
    phase("Writing readout arrival dataset")
    readout_path: Optional[str] = None
    readout_error: Optional[str] = None
    try:
        readout_path = simulator.write_readout_arrivals(run_directory / "readout_arrivals.json", results, detector_mode=mode)
        log(f"Readout arrival dataset written to {readout_path}")
    except Exception as exc:
        # The completed in-memory run and its summary remain valid if this
        # optional compact persistence artifact cannot be replaced.
        readout_error = f"{type(exc).__name__}: {exc}"
        log(f"Readout arrival dataset unavailable: {readout_error}")
    plot_errors = []
    plot_specs = []
    if results:
        plot_specs = [
            ("event_display_event_0.png", _plot_canonical_event_display, (results[0], event_view_scene, settings.display_photons, settings.track_margin_mm)),
            ("detector_impacts_cumulative.png", simulator.plot_cumulative, (results, False)),
            ("photons_per_event.png", simulator.plot_photons_per_event, (results, False)),
            ("arrival_times.png", simulator.plot_arrival_times, (results, False)),
            ("wavelengths.png", simulator.plot_wavelengths, (results, False)),
            ("outcome_face_counts.png", simulator.plot_outcome_face_counts, (results, False)),
            ("transport_diagnostic.png", simulator.plot_result, (results[0],)),
        ]
    phase("Writing plots")
    for filename, function, args in plot_specs:
        if stop_event is not None and stop_event.is_set() and filename != "event_display_event_0.png":
            break
        error = _safe_plot(log, run_directory / filename, function, *args)
        if error:
            plot_errors.append(error)
    stopped = stopped or bool(stop_event is not None and stop_event.is_set())

    phase("Saving summary")
    settings_payload = asdict(settings)
    settings_payload["output_directory"] = str(settings.output_directory)
    sensor_values = dict(default_values())
    sensor_values.update(getattr(settings, "sensor_profile_values", {}) or {})
    geometry_payload = asdict(geometry) if hasattr(geometry, "__dataclass_fields__") else {"shape": "box"}
    payload = {
        "run_directory": str(run_directory),
        "settings": settings_payload,
        "geometry": geometry_payload,
        "scintillator_assembly": None if scintillator_output_assembly is None else scintillator_output_assembly.preview(),
        "two_detector_beam_monitor": None if two_detector_monitor is None else two_detector_monitor.preview(geometry),
        "sensors": sensor_profile_from_values(sensor_values, assembly_models),
        "two_detector_note": "The existing scintillation transport result remains the single-assembly compatibility path; the two-detector object supplies explicit global placement, readout guides, and the shared engineering particle path.",
        # The public field is the compact persisted scene.  The private field
        # is intentionally excluded from summary.json and retains the full
        # transformed scene for the in-process event viewer.
        "event_view_scene": compact_event_view_payload(event_view_scene),
        "event_view_payload": {
            "schema": "canonical-two-assembly-event-view-v50",
            "coordinate_system": "global_mm",
            "detector_mode": mode,
            "detector_mode_label": detector_run_mode_label(mode),
            "visible_assembly_types": list(event_view_scene.get("visible_assembly_types", ())),
            "detectors": {
                "cherenkov": {"participates": mode in {"cherenkov", "combined"}, "photon_kind": "cherenkov", "readout_owner": "cherenkov"},
                "scintillator": {"participates": mode in {"scintillator", "combined"}, "photon_kind": "scintillator", "readout_owner": "scintillator"},
            },
            "scene_source": event_view_scene.get("global_geometry_source"),
        },
        "summary": summary,
        "run_controls": {
            "maximum_photon_bounces": settings.maximum_photon_bounces,
            "effective_max_bounces_by_detector": effective_bounces,
            "precedence": "Run controls > Maximum photon bounces overrides legacy sensor-tab max-bounce fields at run start",
        },
        "root_path": root_path,
        "root_error": root_error,
        "readout_arrivals_path": readout_path,
        "readout_arrivals_error": readout_error,
        "plot_errors": plot_errors,
        "stopped": stopped,
        "storage_note": "Photon and event scalars are stored in ROOT TTrees when available. The compact readout_arrivals.json stores only detector-owned in-window arrival times and response outcomes; PNG plots are separate files in this directory.",
        "_results": results,
        "_settings_object": settings,
        "_event_view_scene": event_view_scene,
    }
    summary_path = run_directory / "summary.json"
    disk_payload = {key: value for key, value in payload.items() if not key.startswith("_")}
    safe_disk_payload = _json_safe(disk_payload)
    summary_text = json.dumps(safe_disk_payload, allow_nan=False, indent=2, sort_keys=True) + "\n"
    summary_path.write_text(summary_text, encoding="utf-8")
    log(f"Saved {summary_path.name} ({len(summary_text)} bytes)")
    phase("Updating GUI")
    return payload


class QuartzCherenkovGUI(tk.Tk):
    """Tkinter front end with separate run, input, detector, and geometry workspaces."""

    GEOMETRY_MATERIAL_COLORS = {
        "fused silica": "#9fd3e6",
        "quartz": "#9fd3e6",
        "air": "#d9d9d9",
        "detector glass": "#f4b183",
        "custom": "#c5b0d5",
    }

    def __init__(self) -> None:
        super().__init__()
        self.title(f"Quartz Cherenkov Simulator | GUI {APP_VERSION} {APP_BUILD_ID}")
        self.geometry("1100x760")
        self.minsize(900, 620)
        self.values: Dict[str, tk.StringVar] = {
            key: tk.StringVar(value=value) for key, value in default_values().items()
        }
        self.selected_event = tk.StringVar(value="0")
        self.gamma_selected_event = tk.StringVar(value="0")
        self.latest_run_text = tk.StringVar(value="Latest run: none")
        self.total_summary_text = tk.StringVar(value="No simulation results loaded")
        self.status = tk.StringVar(value="Ready")
        self.queue: queue.Queue = queue.Queue()
        self.worker: Optional[threading.Thread] = None
        self._worker_terminal_seen = True
        self.stop_requested = threading.Event()
        self.run_state = "idle"
        self.last_run_directory: Optional[Path] = None
        self.latest_results = []
        self.latest_settings: Optional[RunSettings] = None
        self.latest_readout_payload: Optional[Mapping[str, object]] = None
        self.latest_summary: Optional[Mapping[str, object]] = None
        self.latest_gamma_payload: Optional[Mapping[str, object]] = None
        self.results_event_selector = None
        self.results_selected_figure = None
        self.results_aggregate_figure = None
        self.results_selected_canvas = None
        self.results_aggregate_canvas = None
        self.last_geometry_io: Dict[str, str] = {"operation": "none", "assembly": "", "path": "", "format": ""}
        self._visible_event_rows = []
        self._field_widgets: Dict[str, tk.Widget] = {}
        self._field_labels: Dict[str, ttk.Label] = {}
        self.assembly_models = {
            "cherenkov": AssemblyGeometryModel(
                "Cherenkov", sections=[AssemblySection(name="Quartz radiator", material="fused silica")],
            ),
            "scintillator": AssemblyGeometryModel(
                "Scintillator", sections=[AssemblySection(name="Scintillator active", material="scintillator", length_x_mm=20.0, width_y_mm=20.0, thickness_z_mm=10.0, axis="X")],
            ),
        }
        self.geometry_sections = self.assembly_models["cherenkov"].sections
        self.scintillator_geometry_sections = self.assembly_models["scintillator"].sections
        self.selected_geometry_index = 0
        self.geometry_property_vars: Dict[str, tk.StringVar] = {}
        self.geometry_property_widgets: Dict[str, tk.Widget] = {}
        self._loading_geometry_properties = False
        self._suppress_geometry_selection_event = False
        self.geometry_canvases: Dict[str, tk.Canvas] = {}
        self.geometry_viewer: Optional[tk.Toplevel] = None
        self.geometry_viewer_canvas: Optional[tk.Canvas] = None
        self.geometry_viewer_bindings = []
        self.geometry_viewer_yaw = 0.75
        self.geometry_viewer_pitch = 0.45
        self.geometry_viewer_zoom = 1.0
        self.geometry_viewer_drag = None
        self.progress_started_at = 0.0
        self.progress_total_events = 0
        self.progress_event_index = 0
        self.progress_event_fraction = 0.0
        self.interactive_event_viewer: Optional[tk.Toplevel] = None
        self.interactive_event_canvas: Optional[tk.Canvas] = None
        self.interactive_event_result = None
        self.interactive_event_yaw = 0.75
        self.interactive_event_pitch = 0.45
        self.interactive_event_zoom = 1.0
        self.interactive_event_drag = None
        # Gamma's detached global viewer owns a separate window, camera, event
        # selector, and visibility controls. It never borrows the optical
        # photon viewer's state.
        self.gamma_global_3d_viewer: Optional[tk.Toplevel] = None
        self.gamma_global_3d_canvas: Optional[tk.Canvas] = None
        self.gamma_global_3d_event_choice = tk.StringVar(value="0")
        self.gamma_global_3d_yaw = -0.96
        self.gamma_global_3d_pitch = 0.38
        self.gamma_global_3d_zoom = 1.0
        self.gamma_global_3d_drag = None
        self.gamma_global_3d_scene = None
        self.gamma_global_3d_show_geometry = tk.BooleanVar(value=True)
        self.gamma_global_3d_show_source = tk.BooleanVar(value=True)
        self.gamma_global_3d_show_tracks = tk.BooleanVar(value=True)
        self.show_3d_labels = tk.BooleanVar(value=True)
        self.global_viewer_label_layout = None
        self.global_viewer = None
        self.global_viewer_canvas = None
        self.global_viewer_info = None
        self.global_viewer_info_text = None
        self.global_viewer_buttons = {}
        self._global_viewer_draw_after = None
        self.global_viewer_render_error = ""
        self.interactive_event_label_layout = None
        # These BooleanVars live for the whole GUI lifetime.  They are not
        # recreated when the selected event changes, so filter choices remain
        # persistent across event selections.
        self.event_view_category_vars = {
            key: tk.BooleanVar(value=(key == "photon_paths"))
            for key in EVENT_VIEW_CATEGORY_LABELS
        }
        self.interactive_event_filter_message = ""
        self.diagnostic_text = tk.StringVar(value=source_diagnostic_text())
        self._scrollable_regions = []
        self._scrollable_mousewheel_passthrough = []
        self._scroll_bindings_installed = False
        self.local_preview_view = "isometric"
        self.local_viewer: Optional[tk.Toplevel] = None
        self.local_viewer_canvas: Optional[tk.Canvas] = None
        self.local_viewer_bindings = []
        # The local and Cherenkov viewers keep both their window/canvas handle
        # and the last canonical scene snapshot.  The snapshot is diagnostic
        # and headless-testable; every refresh rebuilds it from the live model.
        self.local_viewer_scene = None
        self.geometry_viewer_scene = None
        # Atomic complete-assembly loads update all Tk variables in one
        # transaction.  Traces are held until the canonical model and sensor
        # snapshot are both committed, then one consolidated redraw is queued.
        self._suppress_assembly_refresh = False
        self._pending_assembly_refresh = {}
        self.local_viewer_yaw = 0.75
        self.local_viewer_pitch = 0.45
        self.local_viewer_zoom = 1.0
        self.local_viewer_drag = None
        self._build_widgets()
        # Queue polling belongs to the Tk thread and must be scheduled once
        # when the window is ready.  v41 accidentally left this call after a
        # return in _three_d_labels_enabled, so worker messages had no
        # consumer and the visible status could remain stale forever.
        self._queue_poll_scheduled = False
        self._load_latest_summary()
        self._schedule_queue_poll()

    def _three_d_labels_enabled(self) -> bool:
        """Read the compact-label toggle safely for Tk and headless probes."""

        variable = _get_instance_attr(self, "show_3d_labels", None)
        if variable is None:
            return True
        try:
            return bool(variable.get())
        except (AttributeError, tk.TclError):
            return True

    def _schedule_queue_poll(self) -> None:
        """Schedule exactly one Tk-side queue poll, including after errors."""

        if _get_instance_attr(self, "_queue_poll_scheduled", False):
            return
        try:
            self._queue_poll_scheduled = True
            self.after(100, self._poll_queue)
        except (RuntimeError, tk.TclError):
            self._queue_poll_scheduled = False

    def _build_widgets(self) -> None:
        self.columnconfigure(0, weight=1)
        self.rowconfigure(1, weight=1)
        header = ttk.Frame(self, padding=(6, 5, 6, 3))
        header.grid(row=0, column=0, sticky="ew")
        header.columnconfigure(0, weight=1)
        ttk.Label(header, text="Quartz Cherenkov Simulator", font=("TkDefaultFont", 12, "bold")).grid(row=0, column=0, sticky="w")
        ttk.Label(header, text="Run results are isolated from the geometry editor; optical transport uses the active supported solid geometry.").grid(row=1, column=0, sticky="w", pady=(1, 0))
        ttk.Label(header, textvariable=self.latest_run_text, foreground="#245a35").grid(row=2, column=0, sticky="w", pady=(2, 0))
        diagnostic_frame = ttk.LabelFrame(header, text="Loaded source diagnostic", padding=(6, 3))
        diagnostic_frame.grid(row=3, column=0, sticky="ew", pady=(4, 0))
        diagnostic_frame.columnconfigure(0, weight=1)
        ttk.Label(
            diagnostic_frame,
            textvariable=self.diagnostic_text,
            justify="left",
            anchor="w",
            wraplength=1040,
        ).grid(row=0, column=0, sticky="ew")
        diagnostic_buttons = ttk.Frame(diagnostic_frame)
        diagnostic_buttons.grid(row=0, column=1, sticky="ne", padx=(10, 0))
        ttk.Button(diagnostic_buttons, text="Copy diagnostic", command=self._copy_diagnostic).pack(side="left")
        ttk.Button(diagnostic_buttons, text="Refresh diagnostic", command=self._refresh_diagnostic).pack(side="left", padx=(5, 0))

        self.notebook = ttk.Notebook(self)
        self.notebook.grid(row=1, column=0, sticky="nsew", padx=5, pady=(0, 5))
        # v13 deliberately has four root tabs.  Run retains the existing
        # workflow in a small nested notebook; assembly-owned controls live
        # only below their corresponding assembly root tab.
        run_tab = ttk.Frame(self.notebook)
        results_tab = ttk.Frame(self.notebook)
        gamma_results_tab = ttk.Frame(self.notebook)
        cherenkov_tab = ttk.Frame(self.notebook)
        scintillator_tab = ttk.Frame(self.notebook)
        global_geometry_tab = ttk.Frame(self.notebook)
        for tab, label in zip((run_tab, results_tab, gamma_results_tab, cherenkov_tab, scintillator_tab, global_geometry_tab), ROOT_NOTEBOOK_TAB_LABELS):
            self.notebook.add(tab, text=label)

        run_notebook = ttk.Notebook(run_tab)
        run_notebook.pack(fill="both", expand=True)
        run_controls_tab = ttk.Frame(run_notebook)
        particle_tab = ttk.Frame(run_notebook)
        detector_tab = ttk.Frame(run_notebook)
        run_notebook.add(run_controls_tab, text="Run controls")
        run_notebook.add(particle_tab, text="Particle track")
        run_notebook.add(detector_tab, text="Transport diagnostics")
        self._build_run_tab(self._make_scrollable_tab(
            run_controls_tab, "RUN", "Run controls, progress, event results, and event viewers."
        ))
        self._build_particle_tab(self._make_scrollable_tab(
            particle_tab, "PARTICLE TRACK", "Shared particle and beam configuration used by all assembly modes."
        ))
        self._build_detector_tab(self._make_scrollable_tab(
            detector_tab, "TRANSPORT DIAGNOSTICS", "Legacy transport validation remains available from Run. Assembly sensors are configured in their own tabs."
        ))

        self._build_results_tab(results_tab)
        self._build_gamma_results_tab(self._make_scrollable_tab(
            gamma_results_tab, "GAMMA RESULTS", "Scroll through the Gamma run summary, event list, and selected event viewer and details."
        ))

        cherenkov_notebook = ttk.Notebook(cherenkov_tab)
        cherenkov_notebook.pack(fill="both", expand=True)
        self.cherenkov_notebook = cherenkov_notebook
        ch_geo_tab = ttk.Frame(cherenkov_notebook)
        ch_sensor_tab = ttk.Frame(cherenkov_notebook)
        cherenkov_notebook.add(ch_geo_tab, text=ASSEMBLY_SUBTAB_LABELS["Cherenkov"][0])
        cherenkov_notebook.add(ch_sensor_tab, text=ASSEMBLY_SUBTAB_LABELS["Cherenkov"][1])
        self._build_cherenkov_geometry_tab(self._make_scrollable_tab(
            ch_geo_tab, "CHERENKOV GEOMETRY EDITOR", "Edit the ordered local Cherenkov assembly, optional appended light guide, terminal interfaces, Mylar, and local XY particle impact."
        ))
        self._build_cherenkov_sensor_tab(self._make_scrollable_tab(
            ch_sensor_tab, "CHERENKOV SENSOR", "Configure the Cherenkov photosensor/readout and inspect real versus ideal optical coupling."
        ))

        scintillator_notebook = ttk.Notebook(scintillator_tab)
        scintillator_notebook.pack(fill="both", expand=True)
        self.scintillator_notebook = scintillator_notebook
        sc_geo_tab = ttk.Frame(scintillator_notebook)
        sc_sensor_tab = ttk.Frame(scintillator_notebook)
        scintillator_notebook.add(sc_geo_tab, text=ASSEMBLY_SUBTAB_LABELS["Scintillator"][0])
        scintillator_notebook.add(sc_sensor_tab, text=ASSEMBLY_SUBTAB_LABELS["Scintillator"][1])
        self._build_scintillator_geometry_tab_v13(self._make_scrollable_tab(
            sc_geo_tab, "SCINTILLATOR GEOMETRY EDITOR", "Use the same ordered-section model as Cherenkov. Quartz geometry is intentionally not shown here."
        ))
        self._build_scintillator_sensor_tab_v13(self._make_scrollable_tab(
            sc_sensor_tab, "SCINTILLATOR SENSOR", "Configure scintillator generation, readout physics, timing, Mylar, and coupling diagnostics."
        ))
        self._build_global_assembly_tab_v13(self._make_scrollable_tab(
            global_geometry_tab, "GLOBAL ASSEMBLY 3D", "Display and position complete assemblies only. Internal geometry, guide, sensor, terminal, Mylar, and local impact edits stay in assembly tabs."
        ))

        for notebook in (self.notebook, self.cherenkov_notebook, self.scintillator_notebook):
            notebook.bind("<<NotebookTabChanged>>", self._refresh_visible_assembly_previews, add="+")

        self.values["particle_type"].trace_add("write", self._particle_type_changed)
        self.values["opposite_end_mode"].trace_add("write", self._mode_changed)
        self.values["sampling_mode"].trace_add("write", self._sampling_mode_changed)
        self._sampling_mode_changed()
        self._mode_changed()
        self._refresh_geometry_views()
        self._assembly_draw_preview("cherenkov")
        self._assembly_draw_preview("scintillator")

    def _refresh_diagnostic(self) -> None:
        """Refresh the visible provenance information for the loaded files."""

        self.diagnostic_text.set(source_diagnostic_text())
        self.status.set("Source diagnostic refreshed")

    def _copy_diagnostic(self) -> None:
        """Copy the full provenance diagnostic when a Tk clipboard exists."""

        try:
            self.clipboard_clear()
            self.clipboard_append(self.diagnostic_text.get())
            self.update()
            self.status.set("Source diagnostic copied to clipboard")
        except tk.TclError:
            # Clipboard ownership can be unavailable on minimal/headless X
            # sessions; the diagnostic remains fully visible in the header.
            self.status.set("Clipboard unavailable")

    def _make_scrollable_tab(self, parent: ttk.Frame, title: str, instruction: str) -> ttk.Frame:
        """Return a vertically scrollable content frame for a root tab."""

        parent.columnconfigure(0, weight=1)
        parent.rowconfigure(1, weight=1)
        banner = ttk.LabelFrame(parent, text=title, padding=(8, 5))
        banner.grid(row=0, column=0, sticky="ew", padx=5, pady=(5, 3))
        ttk.Label(banner, text=instruction, foreground="#245a35", wraplength=1050,
                  justify="left").grid(row=0, column=0, sticky="w")
        viewport = ttk.Frame(parent)
        viewport.grid(row=1, column=0, sticky="nsew", padx=5, pady=(0, 5))
        viewport.columnconfigure(0, weight=1)
        viewport.rowconfigure(0, weight=1)
        canvas = tk.Canvas(viewport, borderwidth=0, highlightthickness=0)
        scrollbar = ttk.Scrollbar(viewport, orient="vertical", command=canvas.yview)
        canvas.grid(row=0, column=0, sticky="nsew")
        scrollbar.grid(row=0, column=1, sticky="ns")
        canvas.configure(yscrollcommand=scrollbar.set)
        inner = ttk.Frame(canvas, padding=5)
        window = canvas.create_window((0, 0), window=inner, anchor="nw")
        inner.bind("<Configure>", lambda _event: canvas.configure(scrollregion=canvas.bbox("all")))
        canvas.bind("<Configure>", lambda event: canvas.itemconfigure(window, width=event.width))
        self._scrollable_regions.append((canvas, inner))
        if not self._scroll_bindings_installed:
            self.bind_all("<MouseWheel>", self._scrollable_mousewheel, add="+")
            self.bind_all("<Button-4>", self._scrollable_mousewheel, add="+")
            self.bind_all("<Button-5>", self._scrollable_mousewheel, add="+")
            self._scroll_bindings_installed = True
        return inner

    def _scrollable_mousewheel(self, event):
        """Scroll the root tab whose content contains the pointer."""

        try:
            widget = self.winfo_containing(event.x_root, event.y_root)
        except tk.TclError:
            widget = None
        target = None
        if widget is not None:
            widget_path = str(widget)
            for nested_widget in _get_instance_attr(self, "_scrollable_mousewheel_passthrough", ()):
                nested_path = str(nested_widget)
                if widget_path == nested_path or widget_path.startswith(nested_path + "."):
                    # Keep wheel input on Gamma's independently scrollable
                    # event list and text panes; all other points in the tab
                    # scroll the complete results page.
                    return
            for canvas, inner in self._scrollable_regions:
                if widget_path == str(inner) or widget_path.startswith(str(inner) + "."):
                    target = canvas
                    break
        if target is None:
            return
        delta = getattr(event, "delta", 0)
        if delta:
            units = -int(delta / 120) or (-1 if delta > 0 else 1)
        else:
            units = -1 if getattr(event, "num", 4) == 4 else 1
        target.yview_scroll(units, "units")
        return "break"

    # v13 assembly editors -------------------------------------------------

    def _assembly_model(self, assembly: str) -> AssemblyGeometryModel:
        """Return the one canonical local model used by editor and preview."""

        return self.assembly_models[assembly]

    @staticmethod
    def _assembly_section_dimension_summary(section: AssemblySection) -> str:
        if section.shape == "box":
            return f"{section.length_x_mm:g} × {section.width_y_mm:g} × {section.thickness_z_mm:g}"
        if section.shape == "tapered_rectangular_prism":
            return f"X {section.length_x_mm:g}; Y {section.width_y_start_mm:g}→{section.width_y_end_mm:g}; Z {section.thickness_z_start_mm:g}→{section.thickness_z_end_mm:g}"
        if section.shape == "cylinder":
            return f"R {section.radius_mm:g}; L {section.length_mm:g}; {section.axis}"
        if section.shape == "elliptical_cylinder":
            return f"A {section.semi_axis_a_mm:g}; B {section.semi_axis_b_mm:g}; L {section.length_mm:g}; {section.axis}"
        return f"R {section.radius_start_mm:g}→{section.radius_end_mm:g}; L {section.length_mm:g}; {section.axis}"

    def _build_shared_assembly_geometry_editor(self, parent: ttk.Frame, assembly: str) -> None:
        """Build the same ordered-section editor for Cherenkov and scintillator."""

        model = self._assembly_model(assembly)
        state = {"selected": ("main", 0), "tree": None, "property_vars": {}, "property_widgets": {}, "preview": None}
        if not _has_instance_attr(self, "_assembly_editor_state"):
            self._assembly_editor_state = {}
        self._assembly_editor_state[assembly] = state
        guide_key = "cherenkov_guide_enabled" if assembly == "cherenkov" else "scintillator_guide_enabled"

        parent.columnconfigure(0, weight=1)
        parent.rowconfigure(1, weight=1)
        editor = ttk.LabelFrame(parent, text=f"{assembly.title()} local assembly — canonical ordered sections", padding=6)
        editor.grid(row=0, column=0, sticky="ew", pady=(0, 6))
        editor.columnconfigure(0, weight=1)
        editor.columnconfigure(1, weight=2)
        list_frame = ttk.Frame(editor)
        list_frame.grid(row=0, column=0, rowspan=2, sticky="nsew", padx=(0, 8))
        tree = ttk.Treeview(list_frame, columns=("role", "name", "shape", "dimensions"), show="headings", height=8, selectmode="browse")
        for column, heading, width in (("role", "Role", 70), ("name", "Name", 150), ("shape", "Shape", 130), ("dimensions", "Dimensions mm", 240)):
            tree.heading(column, text=heading)
            tree.column(column, width=width, anchor="w")
        tree.pack(side="left", fill="both", expand=True)
        ttk.Scrollbar(list_frame, orient="vertical", command=tree.yview).pack(side="right", fill="y")
        tree.configure(yscrollcommand=lambda first, last: None)
        state["tree"] = tree

        button_row = ttk.Frame(editor)
        button_row.grid(row=2, column=0, sticky="w", pady=(5, 0))
        ttk.Button(button_row, text="Add section", command=lambda: self._assembly_add_section(assembly, "main")).pack(side="left")
        ttk.Button(button_row, text="Remove", command=lambda: self._assembly_remove_section(assembly)).pack(side="left", padx=3)
        ttk.Button(button_row, text="Move up", command=lambda: self._assembly_move_section(assembly, -1)).pack(side="left", padx=3)
        ttk.Button(button_row, text="Move down", command=lambda: self._assembly_move_section(assembly, 1)).pack(side="left", padx=3)
        ttk.Button(button_row, text="Add light-guide section", command=lambda: self._assembly_add_section(assembly, "guide")).pack(side="left", padx=3)
        ttk.Checkbutton(button_row, text="Enable physically separate light guide (append to active geometry)", variable=self.values[guide_key], onvalue="1", offvalue="0", command=lambda: self._assembly_model_changed(assembly)).pack(side="left", padx=(8, 0))

        props = ttk.LabelFrame(editor, text="Selected section properties", padding=5)
        props.grid(row=0, column=1, sticky="nsew")
        props.columnconfigure(1, weight=1)
        state["properties"] = props
        tree.bind("<<TreeviewSelect>>", lambda _event: self._assembly_selection_changed(assembly))
        self._assembly_refresh_tree(assembly)
        self._assembly_build_properties(assembly)

        preview_frame = ttk.LabelFrame(parent, text=f"Complete local {assembly} preview — active sections, guide, interfaces, sensor/readout, Mylar, axes, beam +Z", padding=5)
        preview_frame.grid(row=1, column=0, sticky="nsew")
        preview_frame.columnconfigure(0, weight=1)
        preview_frame.rowconfigure(1, weight=1)
        state["preview_text"] = tk.StringVar(value="Validating local assembly…")
        ttk.Label(preview_frame, textvariable=state["preview_text"], justify="left", wraplength=1000).grid(row=0, column=0, sticky="w")
        canvas = tk.Canvas(preview_frame, background="white", height=420, highlightthickness=1, highlightbackground="#b0b0b0")
        canvas.grid(row=1, column=0, sticky="nsew", pady=(5, 0))
        state["preview"] = canvas
        state["preview_redraw_after"] = None
        canvas.bind("<Configure>", lambda _event, a=assembly: self._schedule_assembly_preview_redraw(a), add="+")
        canvas.bind("<Map>", lambda _event, a=assembly: self._schedule_assembly_preview_redraw(a), add="+")
        canvas.bind("<Visibility>", lambda _event, a=assembly: self._schedule_assembly_preview_redraw(a), add="+")
        preview_frame.bind("<Configure>", lambda _event, a=assembly: self._schedule_assembly_preview_redraw(a), add="+")
        if assembly == "cherenkov":
            self.geometry_canvases = {"iso": canvas}
            ttk.Button(preview_frame, text="Open interactive complete local 3D", command=self._open_geometry_3d_viewer).grid(row=2, column=0, sticky="e", pady=(5, 0))
        else:
            self.scintillator_preview_canvas = canvas
            ttk.Button(preview_frame, text="Open interactive complete local 3D", command=self._open_local_3d_viewer).grid(row=2, column=0, sticky="e", pady=(5, 0))
        model.guide_enabled = _bool_text(self.values, guide_key)
        self._assembly_draw_preview(assembly)

    def _schedule_assembly_preview_redraw(self, assembly: str) -> None:
        """Coalesce resize/tab redraws while preserving one canvas per preview."""

        state = _get_instance_attr(self, "_assembly_editor_state", {}).get(assembly)
        if not state or state.get("preview") is None:
            return
        canvas = state["preview"]
        pending = state.get("preview_redraw_after")
        if pending is not None:
            try:
                canvas.after_cancel(pending)
            except tk.TclError:
                pass

        def redraw() -> None:
            state["preview_redraw_after"] = None
            self._assembly_draw_preview(assembly)

        try:
            state["preview_redraw_after"] = canvas.after_idle(redraw)
        except tk.TclError:
            state["preview_redraw_after"] = None

    def _refresh_loaded_assembly_views(self, assembly: str) -> None:
        """Redraw every local representation from the just-committed snapshot.

        The immediate call updates the renderer-neutral scene payload and any
        mapped canvas.  The idle call runs after Tk has completed notebook and
        scroll-frame layout, so an initially unmapped 1×1 canvas cannot remain
        blank until the user changes tabs.
        """

        if _get_instance_attr(self, "_suppress_assembly_refresh", False):
            self._pending_assembly_refresh[assembly] = True
            return
        self._assembly_draw_preview(assembly)
        self._refresh_live_assembly_viewer(assembly)
        self._refresh_global_assembly_preview()
        self._schedule_assembly_preview_redraw(assembly)

    def _refresh_visible_assembly_previews(self, *_args) -> None:
        """Refresh both complete local previews after notebook/tab changes."""

        for assembly in ("cherenkov", "scintillator"):
            self._schedule_assembly_preview_redraw(assembly)

    def _assembly_refresh_tree(self, assembly: str) -> None:
        state = self._assembly_editor_state[assembly]
        tree = state["tree"]
        tree.delete(*tree.get_children())
        model = self._assembly_model(assembly)
        rows = [("main", i, section) for i, section in enumerate(model.sections)] + [("guide", i, section) for i, section in enumerate(model.guide_sections)]
        for role, index, section in rows:
            tree.insert("", "end", iid=f"{role}:{index}", values=("active" if role == "main" else "guide", section.name, GEOMETRY_SHAPE_LABELS.get(section.shape, section.shape), self._assembly_section_dimension_summary(section)))
        role, index = state["selected"]
        iid = f"{role}:{index}"
        if tree.exists(iid):
            tree.selection_set(iid)
            tree.focus(iid)

    def _assembly_selection_changed(self, assembly: str) -> None:
        selection = self._assembly_editor_state[assembly]["tree"].selection()
        if not selection:
            return
        role, raw_index = selection[0].split(":", 1)
        self._assembly_editor_state[assembly]["selected"] = (role, int(raw_index))
        self._assembly_build_properties(assembly)

    def _assembly_selected_section(self, assembly: str):
        role, index = self._assembly_editor_state[assembly]["selected"]
        items = self._assembly_model(assembly).sections if role == "main" else self._assembly_model(assembly).guide_sections
        if not items:
            return role, index, None
        index = max(0, min(index, len(items) - 1))
        self._assembly_editor_state[assembly]["selected"] = (role, index)
        return role, index, items[index]

    def _assembly_build_properties(self, assembly: str) -> None:
        state = self._assembly_editor_state[assembly]
        props = state["properties"]
        for child in props.winfo_children():
            child.destroy()
        role, index, section = self._assembly_selected_section(assembly)
        if section is None:
            return
        state["property_vars"] = {}
        base_fields = [("Name", "name", "entry"), ("Shape", "shape", "shape"), ("Material", "material", "material")]
        if section.shape == "box":
            fields = base_fields + [("Length X (mm)", "length_x_mm", "entry"), ("Width Y (mm)", "width_y_mm", "entry"), ("Thickness Z (mm)", "thickness_z_mm", "entry")]
        elif section.shape == "tapered_rectangular_prism":
            fields = base_fields + [("Length X (mm)", "length_x_mm", "entry"), ("Width Y start (mm)", "width_y_start_mm", "entry"), ("Width Y end (mm)", "width_y_end_mm", "entry"), ("Thickness Z start (mm)", "thickness_z_start_mm", "entry"), ("Thickness Z end (mm)", "thickness_z_end_mm", "entry")]
        elif section.shape == "cylinder":
            fields = base_fields + [("Axis", "axis", "axis"), ("Radius (mm)", "radius_mm", "entry"), ("Length (mm)", "length_mm", "entry")]
        elif section.shape == "elliptical_cylinder":
            fields = base_fields + [("Axis", "axis", "axis"), ("Semi-axis A (mm)", "semi_axis_a_mm", "entry"), ("Semi-axis B (mm)", "semi_axis_b_mm", "entry"), ("Length (mm)", "length_mm", "entry")]
        else:
            fields = base_fields + [("Axis", "axis", "axis"), ("Start radius (mm)", "radius_start_mm", "entry"), ("End radius (mm)", "radius_end_mm", "entry"), ("Length (mm)", "length_mm", "entry")]
        fields += [("Center X (mm)", "center_x_mm", "entry"), ("Center Y (mm)", "center_y_mm", "entry"), ("Center Z (mm)", "center_z_mm", "entry")]
        for row, (label, key, kind) in enumerate(fields):
            ttk.Label(props, text=label).grid(row=row, column=0, sticky="w", padx=(0, 8), pady=2)
            value = GEOMETRY_SHAPE_LABELS.get(section.shape, section.shape) if key == "shape" else getattr(section, key)
            variable = tk.StringVar(value=str(value))
            state["property_vars"][key] = variable
            if kind == "shape":
                widget = ttk.Combobox(props, textvariable=variable, values=tuple(GEOMETRY_SHAPE_LABELS.values()), state="readonly")
            elif kind == "material":
                widget = ttk.Combobox(props, textvariable=variable, values=("fused silica", "quartz", "air", "scintillator", "light guide", "custom"), state="readonly")
            elif kind == "axis":
                widget = ttk.Combobox(props, textvariable=variable, values=GEOMETRY_AXES, state="readonly")
            else:
                widget = ttk.Entry(props, textvariable=variable)
            widget.grid(row=row, column=1, sticky="ew", pady=2)
            if kind == "entry":
                widget.bind("<Return>", lambda _event, a=assembly: self._assembly_apply_properties(a))
                widget.bind("<FocusOut>", lambda _event, a=assembly: self._assembly_apply_properties(a))
            else:
                widget.bind("<<ComboboxSelected>>", lambda _event, a=assembly: self._assembly_apply_properties(a))
        ttk.Button(props, text="Apply section", command=lambda: self._assembly_apply_properties(assembly)).grid(row=len(fields), column=1, sticky="e", pady=(5, 0))

    def _assembly_apply_properties(self, assembly: str) -> None:
        state = self._assembly_editor_state[assembly]
        role, index, section = self._assembly_selected_section(assembly)
        if section is None:
            return
        old = section.to_dict()
        numeric = {"length_x_mm", "width_y_mm", "thickness_z_mm", "width_y_start_mm", "width_y_end_mm", "thickness_z_start_mm", "thickness_z_end_mm", "radius_mm", "length_mm", "semi_axis_a_mm", "semi_axis_b_mm", "radius_start_mm", "radius_end_mm", "center_x_mm", "center_y_mm", "center_z_mm"}
        try:
            for key, variable in state["property_vars"].items():
                value = variable.get()
                if key == "shape":
                    value = GEOMETRY_SHAPE_BY_LABEL.get(value, value)
                elif key in numeric:
                    value = float(value)
                elif key == "axis":
                    value = value.upper()
                setattr(section, key, value)
            section.validate()
        except (TypeError, ValueError) as exc:
            for key, value in old.items():
                setattr(section, key, value)
            self._assembly_editor_state[assembly]["preview_text"].set(f"Invalid section; previous value restored: {exc}")
        self._assembly_refresh_tree(assembly)
        self._assembly_build_properties(assembly)
        self._assembly_model_changed(assembly)

    def _assembly_add_section(self, assembly: str, role: str) -> None:
        model = self._assembly_model(assembly)
        is_guide = role == "guide"
        index = len(model.guide_sections) + 1 if is_guide else len(model.sections) + 1
        material = "light guide" if is_guide else ("scintillator" if assembly == "scintillator" else "fused silica")
        existing = list(model.sections) + list(model.guide_sections if is_guide else ())
        current_high = max((section.center_x_mm + (section.length_mm if section.shape in {"cylinder", "elliptical_cylinder", "truncated_cone"} and section.axis == "X" else section.length_x_mm) / 2.0 for section in existing), default=0.0)
        center = current_high + 10.0
        section = AssemblySection(name=f"{'Guide' if is_guide else 'Section'} {index}", material=material, length_x_mm=20.0, width_y_mm=20.0, thickness_z_mm=10.0, axis="X", center_x_mm=center)
        (model.guide_sections if is_guide else model.sections).append(section)
        self._assembly_editor_state[assembly]["selected"] = (role, index - 1)
        self._assembly_refresh_tree(assembly)
        self._assembly_build_properties(assembly)
        self._assembly_model_changed(assembly)

    def _assembly_remove_section(self, assembly: str) -> None:
        model = self._assembly_model(assembly)
        role, index, _section = self._assembly_selected_section(assembly)
        items = model.sections if role == "main" else model.guide_sections
        if role == "main" and len(items) <= 1:
            self._assembly_editor_state[assembly]["preview_text"].set("At least one active section is required.")
            return
        if items:
            del items[index]
        self._assembly_editor_state[assembly]["selected"] = (role, max(0, index - 1))
        self._assembly_refresh_tree(assembly)
        self._assembly_build_properties(assembly)
        self._assembly_model_changed(assembly)

    def _assembly_move_section(self, assembly: str, delta: int) -> None:
        model = self._assembly_model(assembly)
        role, index, _section = self._assembly_selected_section(assembly)
        items = model.sections if role == "main" else model.guide_sections
        target = index + delta
        if not (0 <= target < len(items)):
            return
        items[index], items[target] = items[target], items[index]
        self._assembly_editor_state[assembly]["selected"] = (role, target)
        self._assembly_refresh_tree(assembly)
        self._assembly_model_changed(assembly)

    def _assembly_model_changed(self, assembly: str) -> None:
        if _get_instance_attr(self, "_suppress_assembly_refresh", False):
            self._pending_assembly_refresh[assembly] = True
            return
        model = self._assembly_model(assembly)
        guide_key = "cherenkov_guide_enabled" if assembly == "cherenkov" else "scintillator_guide_enabled"
        model.guide_enabled = _bool_text(self.values, guide_key)
        if assembly == "cherenkov":
            try:
                model.impact_x_mm = float(_value_text(self.values, "entry_x_mm", "0"))
                model.impact_y_mm = float(_value_text(self.values, "entry_y_mm", "0"))
            except (TypeError, ValueError):
                self._assembly_editor_state[assembly]["preview_text"].set("Impact X/Y is being edited; enter finite numbers to validate the assembly.")
                return
            self._refresh_geometry_views()
        else:
            try:
                model.impact_x_mm = float(_value_text(self.values, "scintillator_impact_x_mm", "0"))
                model.impact_y_mm = float(_value_text(self.values, "scintillator_impact_y_mm", "0"))
                model.mylar_enabled = _bool_text(self.values, "mylar_enabled")
                model.mylar_reflectivity = float(_value_text(self.values, "mylar_reflectivity_percent", "95")) / 100.0
                if not all(math.isfinite(value) for value in (model.impact_x_mm, model.impact_y_mm, model.mylar_reflectivity)):
                    raise ValueError("impact and Mylar values must be finite")
            except (TypeError, ValueError) as exc:
                self._assembly_editor_state[assembly]["preview_text"].set(f"Scintillator assembly values are being edited: {exc}")
                return
            self._sync_scintillator_model_to_values()
            self._refresh_scintillator_preview()
        self._assembly_draw_preview(assembly)
        self._refresh_global_assembly_preview()

    def _refresh_live_assembly_viewer(self, assembly: str) -> None:
        """Refresh an already-open local viewer from the live canonical model.

        The three orthographic preview and the interactive window deliberately
        share ``complete_local_assembly_scene``.  A destroyed Toplevel can
        still be referenced by Tk/Python briefly, so existence is checked
        defensively and stale handles are cleared here rather than being used
        on the next model callback.
        """

        if assembly == "scintillator":
            viewer_attr, canvas_attr, draw = "local_viewer", "local_viewer_canvas", self._draw_local_3d_view
            scene_attr = "local_viewer_scene"
        elif assembly == "cherenkov":
            viewer_attr, canvas_attr, draw = "geometry_viewer", "geometry_viewer_canvas", self._draw_interactive_geometry_view
            scene_attr = "geometry_viewer_scene"
        else:
            raise ValueError("assembly must be cherenkov or scintillator")
        viewer = _get_instance_attr(self, viewer_attr, None)
        canvas = _get_instance_attr(self, canvas_attr, None)
        if viewer is None or canvas is None:
            return
        try:
            alive = bool(viewer.winfo_exists()) and bool(canvas.winfo_exists())
        except (AttributeError, tk.TclError):
            alive = False
        if not alive:
            setattr(self, viewer_attr, None)
            setattr(self, canvas_attr, None)
            setattr(self, scene_attr, None)
            if assembly == "cherenkov":
                self.geometry_viewer_bindings = []
            return
        draw()

    def _refresh_scintillator_viewer(self) -> None:
        """Live-refresh hook used by all Scintillator model/GUI callbacks."""

        self._refresh_live_assembly_viewer("scintillator")

    def _refresh_cherenkov_viewer(self) -> None:
        """Live-refresh hook used by all Cherenkov model/GUI callbacks."""

        self._refresh_live_assembly_viewer("cherenkov")

    def _assembly_draw_preview(self, assembly: str) -> None:
        state = _get_instance_attr(self, "_assembly_editor_state", {}).get(assembly)
        if not state or state.get("preview") is None:
            return
        canvas = state["preview"]
        if not canvas.winfo_exists():
            return
        canvas.delete("all")
        model = self._assembly_model(assembly)
        try:
            model.validate()
            sensor_index = "1" if assembly == "cherenkov" else "2"
            # Use the public canonical scene helper so the preview receives
            # the same configured shape/dimensions and terminal transform as
            # both interactive 3D viewers and transport metadata.
            scene = complete_local_assembly_scene(model, self.values, assembly=assembly)
            state["scene"] = scene
            self._draw_assembly_scene_canvas(canvas, scene, assembly)
            state["preview_text"].set(
                scene["summary"] +
                f" • {len(scene['physical_components'])} physical components • "
                f"Mylar overlay: {'enabled' if scene['mylar']['enabled'] else 'disabled'} "
                f"({scene['mylar']['overlay_surface_count']} eligible surfaces, "
                f"{scene['mylar']['reflectivity'] * 100:g}%) • three synchronized views"
            )
            # This is the single callback fan-out for section, guide, terminal,
            # Mylar, impact, and sensor changes made by the Geometry Editor.
            self._refresh_live_assembly_viewer(assembly)
        except (TypeError, ValueError) as exc:
            state["preview_text"].set(f"Invalid local assembly: {exc}")
            canvas.create_text(max(float(canvas.winfo_width()), 1.0) / 2,
                               max(float(canvas.winfo_height()), 1.0) / 2,
                               text="Invalid local assembly", fill="#b00020")

    @staticmethod
    def _draw_assembly_scene_canvas(canvas: tk.Canvas, scene: Mapping[str, object], assembly: str) -> None:
        """Draw XZ, XY, and YZ plus visible physical axes from one scene.

        Axis lines are intentionally drawn on the Tk canvas after component
        polygons. The scene payload supplies their physical-coordinate
        endpoints and camera bounds, so this path stays synchronized with
        geometry, guide state, terminal state, sensors, and Mylar.
        """

        width = max(int(canvas.winfo_width()), 1)
        height = max(int(canvas.winfo_height()), 1)
        panel_bounds = responsive_three_view_panel_bounds(width, height)
        if width < 20 or height < 20:
            return
        for panel_index, view_name in enumerate(("xz", "xy", "yz")):
            view = scene["views"][view_name]
            left, top, right, bottom = panel_bounds[panel_index]
            panel_width = right - left
            canvas.create_rectangle(left, top, right, bottom, outline="#b0bec5")
            canvas.create_text(left + 8, 8, text=view["title"], anchor="nw", fill="#17324d", font=("TkDefaultFont", 9, "bold"))
            polygons = [item.get("polygon", ()) for item in view["components"]]
            points = [point for polygon in polygons for point in polygon]
            if not points:
                continue
            axes = view.get("axes", {})
            axis_bounds = axes.get("bounds_mm", ())
            if len(axis_bounds) == 4:
                min_a = min(min(point[0] for point in points), float(axis_bounds[0]))
                max_a = max(max(point[0] for point in points), float(axis_bounds[2]))
                min_b = min(min(point[1] for point in points), float(axis_bounds[1]))
                max_b = max(max(point[1] for point in points), float(axis_bounds[3]))
            else:
                min_a, max_a = min(point[0] for point in points), max(point[0] for point in points)
                min_b, max_b = min(point[1] for point in points), max(point[1] for point in points)
            span_a, span_b = max(max_a - min_a, 1.0), max(max_b - min_b, 1.0)
            inner_left, inner_right = left + 16, right - 16
            inner_top, inner_bottom = top + 34, bottom - 42
            inner_width = max(inner_right - inner_left, 1.0)
            inner_height = max(inner_bottom - inner_top, 1.0)
            scale = min(inner_width / span_a, inner_height / span_b) * 0.92
            center_a, center_b = (min_a + max_a) / 2, (min_b + max_b) / 2
            def project(point):
                return (left + panel_width / 2 + (point[0] - center_a) * scale,
                        (inner_top + inner_bottom) / 2 - (point[1] - center_b) * scale)
            for item in view["components"]:
                polygon = [project(point) for point in item.get("polygon", ())]
                if len(polygon) >= 3:
                    overlay = item.get("kind") == "mylar_overlay"
                    canvas.create_polygon(
                        polygon, fill=item.get("color", "#c5b0d5"),
                        outline=item.get("edge_color", "#263238") if overlay else "#263238",
                        stipple="gray25" if overlay else "gray50",
                    )
                    if overlay:
                        canvas.create_line(
                            polygon + [polygon[0]], fill=item.get("edge_color", MYLAR_DISPLAY_EDGE),
                            dash=(5, 3), width=2,
                        )
                elif len(polygon) == 2:
                    line_options = {"fill": item.get("edge_color", item.get("color", "#263238")), "width": 4}
                    if item.get("kind") == "mylar_overlay":
                        line_options["dash"] = (5, 3)
                    canvas.create_line(polygon, **line_options)
                if polygon:
                    canvas.create_text(polygon[0][0] + 4, polygon[0][1] - 4, text=item.get("label", item.get("name", "")), anchor="sw", fill=item.get("color", "#17324d"), font=("TkDefaultFont", 7))

            # Metadata alone cannot draw on Tk. Render the canonical physical
            # line/arrow/label/origin primitives after the geometry instead.
            for primitive in axes.get("primitives", ()):
                kind = primitive.get("kind")
                if kind == "line":
                    start = project(primitive["start"])
                    end = project(primitive["end"])
                    canvas.create_line(start[0], start[1], end[0], end[1],
                                       fill=primitive.get("color", "#455a64"),
                                       width=2, arrow=TK_ARROW_LAST)
                    canvas.create_text(end[0] + 5, end[1] - 3,
                                       text=primitive.get("label", ""), anchor="sw",
                                       fill=primitive.get("color", "#263238"),
                                       font=("TkDefaultFont", 9, "bold"))
                elif kind == "origin":
                    point = project(primitive["point"])
                    radius = 3
                    canvas.create_oval(point[0] - radius, point[1] - radius,
                                       point[0] + radius, point[1] + radius,
                                       fill="white", outline=primitive.get("color", "#263238"), width=2)
                    canvas.create_text(point[0] + 5, point[1] + 5,
                                       text=primitive.get("label", "O"), anchor="nw",
                                       fill=primitive.get("color", "#263238"),
                                       font=("TkDefaultFont", 8, "bold"))
            canvas.create_text(right - 8, top + 20, text=axes.get("direction", "+X"), anchor="ne",
                               fill="#455a64", font=("TkDefaultFont", 7))
            canvas.create_text(left + 8, bottom - 26,
                               text=axes.get("legend", "O = origin • arrows = positive axes"),
                               anchor="sw", fill="#455a64", font=("TkDefaultFont", 7))
            if view_name != "yz":
                canvas.create_text(left + 10, bottom - 10, text=f"{view['x_label']} / {view['y_label']}", anchor="sw", fill="#455a64", font=("TkDefaultFont", 7))
            else:
                canvas.create_text(left + 10, bottom - 10, text=axes.get("terminal_direction", "+X / −X terminal direction"), anchor="sw", fill="#455a64", font=("TkDefaultFont", 7))
        legend = "  ".join(f"{item['label']}" for item in scene["legend"])
        canvas.create_text(8, height - 2, text=f"{assembly.upper()} • {legend} • Mylar external only", anchor="sw", fill="#17324d", font=("TkDefaultFont", 7))

    def _sync_scintillator_model_to_values(self) -> None:
        """Keep the legacy transport adapter fed from the canonical editor model."""

        model = self._assembly_model("scintillator")
        if model.sections:
            section = model.sections[0]
            prefix = "scintillator"
            for suffix in ("shape", "length_x_mm", "width_y_mm", "thickness_z_mm", "width_y_start_mm", "width_y_end_mm", "thickness_z_start_mm", "thickness_z_end_mm", "axis", "radius_mm", "length_mm", "semi_axis_a_mm", "semi_axis_b_mm", "radius_start_mm", "radius_end_mm", "center_x_mm", "center_y_mm", "center_z_mm", "material"):
                key = f"{prefix}_{suffix}"
                if key in self.values:
                    self.values[key].set(str(getattr(section, suffix)))
        guide = model.valid_guide_sections[0] if model.has_light_guide else None
        if guide is None:
            # Clear the legacy adapter's optional volume and all stale side
            # flags whenever the canonical editor has no effective guide.
            if "guide_shape" in self.values:
                self.values["guide_shape"].set("none")
            for key in ("detector1_guide_left_enabled", "detector1_guide_right_enabled", "detector2_guide_left_enabled", "detector2_guide_right_enabled"):
                if key in self.values:
                    self.values[key].set("0")
            return
        for suffix in ("shape", "length_x_mm", "width_y_mm", "thickness_z_mm", "width_y_start_mm", "width_y_end_mm", "thickness_z_start_mm", "thickness_z_end_mm", "axis", "radius_mm", "length_mm", "semi_axis_a_mm", "semi_axis_b_mm", "radius_start_mm", "radius_end_mm", "center_x_mm", "center_y_mm", "center_z_mm", "material"):
            key = f"guide_{suffix}"
            if key in self.values:
                self.values[key].set(str(getattr(guide, suffix)))
        for key in ("detector1_guide_left_enabled", "detector1_guide_right_enabled", "detector2_guide_left_enabled", "detector2_guide_right_enabled"):
            if key in self.values:
                self.values[key].set("1")

    def _build_cherenkov_geometry_tab(self, parent: ttk.Frame) -> None:
        self._build_shared_assembly_geometry_editor(parent, "cherenkov")
        terminal = ttk.LabelFrame(parent, text="Cherenkov terminal YZ surfaces and local particle impact", padding=6)
        terminal.grid(row=2, column=0, sticky="ew", pady=(6, 0))
        terminal.columnconfigure(1, weight=1)
        rows = (
            ("+X terminal surface (YZ)", "terminal_plus_surface"),
            ("-X terminal surface (YZ)", "terminal_minus_surface"),
        )
        for row, (label, key) in enumerate(rows):
            ttk.Label(terminal, text=label).grid(row=row, column=0, sticky="w", padx=(0, 8), pady=2)
            ttk.Combobox(terminal, textvariable=self.values[key], values=TERMINAL_SURFACE_CHOICES, state="readonly").grid(row=row, column=1, sticky="ew", pady=2)
            self.values[key].trace_add("write", lambda *_args: self._cherenkov_terminal_changed())
        ttk.Label(terminal, text="Open, Mirror, and Detector/Photosensor coupling are terminal surface choices. A valid configured guide is appended physically and the optical path is derived.").grid(row=2, column=0, columnspan=2, sticky="w", pady=(5, 2))
        ttk.Label(terminal, text="Mylar enabled").grid(row=3, column=0, sticky="w", padx=(0, 8), pady=2)
        ttk.Checkbutton(terminal, variable=self.values["cherenkov_mylar_enabled"], onvalue="1", offvalue="0", command=self._cherenkov_terminal_changed).grid(row=3, column=1, sticky="w")
        ttk.Label(terminal, text="Mylar reflectivity (%)").grid(row=4, column=0, sticky="w", padx=(0, 8), pady=2)
        ttk.Entry(terminal, textvariable=self.values["cherenkov_mylar_reflectivity_percent"]).grid(row=4, column=1, sticky="ew", pady=2)
        ttk.Label(terminal, text="Impact X (mm) — editable").grid(row=5, column=0, sticky="w", padx=(0, 8), pady=2)
        ttk.Entry(terminal, textvariable=self.values["entry_x_mm"]).grid(row=5, column=1, sticky="ew", pady=2)
        ttk.Label(terminal, text="Impact Y (mm) — editable").grid(row=6, column=0, sticky="w", padx=(0, 8), pady=2)
        ttk.Entry(terminal, textvariable=self.values["entry_y_mm"]).grid(row=6, column=1, sticky="ew", pady=2)
        ttk.Label(terminal, text="The beam crosses the local assembly along +Z. Z is derived from the entry surface; there is no independent local impact-Z input.", foreground="#245a35", wraplength=900, justify="left").grid(row=7, column=0, columnspan=2, sticky="w", pady=(5, 0))
        self.values["entry_x_mm"].trace_add("write", lambda *_args: self._assembly_model_changed("cherenkov"))
        self.values["entry_y_mm"].trace_add("write", lambda *_args: self._assembly_model_changed("cherenkov"))
        self.values["cherenkov_mylar_enabled"].trace_add("write", lambda *_args: self._cherenkov_terminal_changed())
        self.values["cherenkov_mylar_reflectivity_percent"].trace_add("write", lambda *_args: self._cherenkov_terminal_changed())
        profile = ttk.Frame(terminal)
        profile.grid(row=8, column=0, columnspan=2, sticky="w", pady=(6, 0))
        ttk.Button(profile, text="Save Cherenkov 3D Geometry JSON", command=lambda: self._save_assembly_model("cherenkov")).pack(side="left")
        ttk.Button(profile, text="Load Cherenkov 3D Geometry JSON", command=lambda: self._load_assembly_model("cherenkov")).pack(side="left", padx=5)
        self._cherenkov_terminal_changed()

    def _cherenkov_terminal_changed(self, *_args) -> None:
        self._terminal_surface_changed("cherenkov")

    def _terminal_surface_changed(self, assembly: str) -> None:
        if _get_instance_attr(self, "_suppress_assembly_refresh", False):
            self._pending_assembly_refresh[assembly] = True
            return
        model = self._assembly_model(assembly)
        prefix = "cherenkov" if assembly == "cherenkov" else "scintillator"
        plus_key = "terminal_plus_surface" if assembly == "cherenkov" else "scintillator_terminal_plus_surface"
        minus_key = "terminal_minus_surface" if assembly == "cherenkov" else "scintillator_terminal_minus_surface"
        sensor_index = "1" if assembly == "cherenkov" else "2"
        model.terminal_plus_surface = _value_text(self.values, plus_key, "Detector/Photosensor coupling")
        model.terminal_minus_surface = _value_text(self.values, minus_key, "Open")
        plus_mode = normalize_terminal_surface(model.terminal_plus_surface)
        minus_mode = normalize_terminal_surface(model.terminal_minus_surface)
        if assembly == "cherenkov" and _value_text(self.values, "opposite_end_mode", "open") != minus_mode:
            self.values["opposite_end_mode"].set(minus_mode)
        if assembly == "cherenkov":
            model.mylar_enabled = _bool_text(self.values, "cherenkov_mylar_enabled")
            try:
                model.mylar_reflectivity = float(_value_text(self.values, "cherenkov_mylar_reflectivity_percent", "95")) / 100.0
            except (TypeError, ValueError):
                model.mylar_reflectivity = 0.95
        else:
            model.mylar_enabled = _bool_text(self.values, "mylar_enabled")
            try:
                model.mylar_reflectivity = float(_value_text(self.values, "mylar_reflectivity_percent", "95")) / 100.0
            except (TypeError, ValueError):
                model.mylar_reflectivity = 0.95
            for key, attribute in (("scintillator_impact_x_mm", "impact_x_mm"), ("scintillator_impact_y_mm", "impact_y_mm")):
                try:
                    setattr(model, attribute, float(_value_text(self.values, key, "0")))
                except (TypeError, ValueError):
                    pass
        if _has_instance_attr(self, "geometry_support_text") and assembly == "cherenkov":
            self.geometry_support_text.set("Terminal state is propagated to the canonical scene/transport; external Mylar excludes optical/readout and mirror interfaces.")
        self._assembly_draw_preview(assembly)
        self._refresh_global_assembly_preview()

    def _build_sensor_common(self, parent: ttk.Frame, assembly: str, fields: Sequence[tuple[str, str]]) -> None:
        frame = ttk.LabelFrame(parent, text=f"{assembly.title()} sensor / readout", padding=7)
        frame.grid(row=0, column=0, sticky="new", padx=(0, 8))
        frame.columnconfigure(1, weight=1)
        dimension_widgets = {}
        index = "1" if assembly == "cherenkov" else "2"
        for row, (label, key) in enumerate(fields):
            label_widget = ttk.Label(frame, text=label)
            label_widget.grid(row=row, column=0, sticky="w", padx=(0, 8), pady=2)
            if key.endswith("_shape"):
                widget = ttk.Combobox(frame, textvariable=self.values[key], values=("rectangular", "circular", "elliptical"), state="readonly")
                widget.bind("<<ComboboxSelected>>", lambda _event, a=assembly: self._sensor_shape_changed(a))
            else:
                widget = ttk.Entry(frame, textvariable=self.values[key])
            widget.grid(row=row, column=1, sticky="ew", pady=2)
            if key in {f"sensor_{index}_width_mm", f"sensor_{index}_height_mm"}:
                dimension_widgets[key] = (label_widget, widget)
            self.values[key].trace_add("write", lambda *_args, a=assembly: self._sensor_status_changed(a))
        next_row = len(fields)
        for label, key in (("Readout radius (mm)", f"sensor_{index}_radius_mm"),
                           ("Semi-axis A (mm)", f"sensor_{index}_semi_axis_a_mm"),
                           ("Semi-axis B (mm)", f"sensor_{index}_semi_axis_b_mm")):
            label_widget = ttk.Label(frame, text=label)
            label_widget.grid(row=next_row, column=0, sticky="w", padx=(0, 8), pady=2)
            widget = ttk.Entry(frame, textvariable=self.values[key])
            widget.grid(row=next_row, column=1, sticky="ew", pady=2)
            dimension_widgets[key] = (label_widget, widget)
            self.values[key].trace_add("write", lambda *_args, a=assembly: self._sensor_status_changed(a))
            next_row += 1
        if not _has_instance_attr(self, "_sensor_dimension_widgets"):
            self._sensor_dimension_widgets = {}
        self._sensor_dimension_widgets[assembly] = dimension_widgets
        status = tk.StringVar(value="Not validated")
        setattr(self, f"{assembly}_sensor_status", status)
        ttk.Label(frame, textvariable=status, foreground="#245a35", wraplength=620, justify="left").grid(row=next_row, column=0, columnspan=2, sticky="w", pady=(6, 0))
        diagnostic = tk.StringVar(value="")
        setattr(self, f"{assembly}_coupling_diagnostic", diagnostic)
        ttk.Button(frame, text="Refresh real-versus-ideal coupling diagnostic", command=lambda a=assembly: self._sensor_status_changed(a)).grid(row=next_row + 1, column=0, columnspan=2, sticky="w", pady=(6, 0))
        ttk.Label(frame, textvariable=diagnostic, foreground="#555555", wraplength=620, justify="left").grid(row=next_row + 2, column=0, columnspan=2, sticky="w")
        ttk.Button(frame, text="Apply ideal sensor settings", command=lambda a=assembly: self._apply_ideal_sensor_settings(a)).grid(row=next_row + 3, column=0, columnspan=2, sticky="w", pady=(7, 0))
        profile_row = ttk.Frame(frame)
        profile_row.grid(row=next_row + 4, column=0, columnspan=2, sticky="w", pady=(6, 0))
        ttk.Button(profile_row, text="Load sensor/readout JSON", command=self._load_sensor_profile).pack(side="left")
        ttk.Button(profile_row, text="Save sensor/readout JSON", command=self._save_sensor_profile).pack(side="left", padx=5)
        self._sensor_shape_changed(assembly)

    def _sensor_shape_changed(self, assembly: str) -> None:
        """Show only the dimensional convention that matches the selected shape."""

        index = "1" if assembly == "cherenkov" else "2"
        shape = _value_text(self.values, f"sensor_{index}_shape", "rectangular").lower()
        if not _has_instance_attr(self, "_sensor_last_shape"):
            self._sensor_last_shape = {}
        previous_shape = self._sensor_last_shape.get(assembly)
        # A shape switch migrates the previous full extents exactly once. This
        # is important for a v61 screen/profile showing circular 2 x 2 mm:
        # its first explicit radius must be 1 mm, not the generic 5 mm default.
        if previous_shape != shape and previous_shape is not None:
            try:
                legacy_width = float(_value_text(self.values, f"sensor_{index}_width_mm", "20"))
                legacy_height = float(_value_text(self.values, f"sensor_{index}_height_mm", "10"))
                if shape in {"circular", "circle"}:
                    self.values[f"sensor_{index}_radius_mm"].set(f"{min(legacy_width, legacy_height) / 2.0:g}")
                elif shape in {"elliptical", "ellipse"}:
                    self.values[f"sensor_{index}_semi_axis_a_mm"].set(f"{legacy_width / 2.0:g}")
                    self.values[f"sensor_{index}_semi_axis_b_mm"].set(f"{legacy_height / 2.0:g}")
            except (KeyError, TypeError, ValueError):
                pass
        self._sensor_last_shape[assembly] = shape
        widgets = _get_instance_attr(self, "_sensor_dimension_widgets", {}).get(assembly, {})
        visible = {
            f"sensor_{index}_width_mm": shape in {"rectangular", "rectangle"},
            f"sensor_{index}_height_mm": shape in {"rectangular", "rectangle"},
            f"sensor_{index}_radius_mm": shape in {"circular", "circle"},
            f"sensor_{index}_semi_axis_a_mm": shape in {"elliptical", "ellipse"},
            f"sensor_{index}_semi_axis_b_mm": shape in {"elliptical", "ellipse"},
        }
        for key, (label_widget, entry) in widgets.items():
            if visible.get(key, False):
                label_widget.grid()
                entry.grid()
            else:
                label_widget.grid_remove()
                entry.grid_remove()
        # Keep legacy full extents synchronized as derived compatibility data;
        # they are not editable when a radius or semiaxes control is active.
        try:
            dimensions = sensor_readout_dimensions_from_values(self.values, index)
            for key, value in ((f"sensor_{index}_width_mm", dimensions["width_mm"]), (f"sensor_{index}_height_mm", dimensions["height_mm"])):
                if key in self.values:
                    self.values[key].set(f"{float(value):g}")
        except (TypeError, ValueError):
            pass
        self._sensor_status_changed(assembly)

    def _build_cherenkov_sensor_tab(self, parent: ttk.Frame) -> None:
        parent.columnconfigure(0, weight=1)
        self._build_sensor_common(parent, "cherenkov", (
            ("Readout shape", "sensor_1_shape"), ("Readout width (mm)", "sensor_1_width_mm"), ("Readout height (mm)", "sensor_1_height_mm"), ("Readout depth (mm)", "sensor_1_depth_mm"),
            ("Quantum efficiency", "sensor_1_qe"), ("Acceptance angle (deg)", "sensor_1_acceptance_deg"), ("Optical window / coupling index", "sensor_1_index"),
            ("Collection efficiency", "sensor_1_collection"), ("Timing/readout sigma (ns)", "sensor_1_timing_sigma_ns"),
        ))
        ttk.Label(parent, text="Terminal, mirror, and opposite-end settings are edited only in Cherenkov Geometry Editor; this Sensor tab owns response and routing.", foreground="#245a35", wraplength=900, justify="left").grid(row=1, column=0, sticky="w", pady=(7, 0))
        sensor = ttk.LabelFrame(parent, text="Cherenkov derived optical readout", padding=7)
        sensor.grid(row=2, column=0, sticky="ew", pady=(7, 0))
        sensor.columnconfigure(1, weight=1)
        ttk.Label(sensor, text="Coupled component (derived status)").grid(row=0, column=0, sticky="w", padx=(0, 8), pady=2)
        self.cherenkov_component_status = tk.StringVar()
        ttk.Label(sensor, textvariable=self.cherenkov_component_status, foreground="#555555").grid(row=0, column=1, sticky="w", pady=2)
        ttk.Label(sensor, text="Optical path (derived, read-only)").grid(row=1, column=0, sticky="w", padx=(0, 8), pady=2)
        self.cherenkov_optical_path = tk.StringVar()
        ttk.Label(sensor, textvariable=self.cherenkov_optical_path, foreground="#555555").grid(row=1, column=1, sticky="w", pady=2)
        ttk.Label(sensor, text="Readout terminal / guide status (read-only)").grid(row=2, column=0, sticky="w", padx=(0, 8), pady=2)
        self.cherenkov_sensor_position_text = tk.StringVar()
        ttk.Label(sensor, textvariable=self.cherenkov_sensor_position_text, foreground="#245a35", wraplength=900, justify="left").grid(row=2, column=1, sticky="w", pady=(5, 0))
        ttk.Label(sensor, text="Terminal surface choices and optional ordered guide sections are edited only in Geometry Editor. This tab derives the complete path and keeps it read-only.", foreground="#245a35", wraplength=900, justify="left").grid(row=3, column=0, columnspan=2, sticky="w", pady=(5, 0))
        self._build_transport_diagnostic_panel(parent, "cherenkov", 4)
        self._sensor_status_changed("cherenkov")

    def _build_scintillator_geometry_tab_v13(self, parent: ttk.Frame) -> None:
        self._build_shared_assembly_geometry_editor(parent, "scintillator")
        self._build_scintillator_terminal_controls(parent)
        ttk.Label(parent, text="Scintillator local assembly only: active ordered sections plus optional appended light-guide sections. This panel contains no other assembly geometry.", foreground="#245a35", wraplength=1000, justify="left").grid(row=3, column=0, sticky="w", pady=(5, 0))
        self._sync_scintillator_model_to_values()
        self._terminal_surface_changed("scintillator")

    def _build_scintillator_terminal_controls(self, parent: ttk.Frame) -> None:
        terminal = ttk.LabelFrame(parent, text="Scintillator terminal surfaces, local impact, and Mylar", padding=6)
        terminal.grid(row=2, column=0, sticky="ew", pady=(6, 0))
        terminal.columnconfigure(1, weight=1)
        for row, (label, key) in enumerate((("+X terminal surface (YZ)", "scintillator_terminal_plus_surface"), ("-X terminal surface (YZ)", "scintillator_terminal_minus_surface"))):
            ttk.Label(terminal, text=label).grid(row=row, column=0, sticky="w", padx=(0, 8), pady=2)
            ttk.Combobox(terminal, textvariable=self.values[key], values=TERMINAL_SURFACE_CHOICES, state="readonly").grid(row=row, column=1, sticky="ew", pady=2)
            self.values[key].trace_add("write", lambda *_args: self._terminal_surface_changed("scintillator"))
        ttk.Label(terminal, text="Open has no marker; Mirror adds a red YZ plate; Detector/Photosensor coupling adds the derived readout at the active terminal or complete guide output. Mylar excludes terminal interfaces.", foreground="#245a35", wraplength=900, justify="left").grid(row=2, column=0, columnspan=2, sticky="w", pady=(5, 0))
        for row, (label, key) in enumerate((("Impact X (mm) — editable", "scintillator_impact_x_mm"), ("Impact Y (mm) — editable", "scintillator_impact_y_mm")), start=3):
            ttk.Label(terminal, text=label).grid(row=row, column=0, sticky="w", padx=(0, 8), pady=2)
            ttk.Entry(terminal, textvariable=self.values[key]).grid(row=row, column=1, sticky="ew", pady=2)
            self.values[key].trace_add("write", lambda *_args: self._assembly_model_changed("scintillator"))
        ttk.Label(terminal, text="Impact Z is derived from the assembly/beam configuration; there is no independent local impact-Z input.", foreground="#245a35", wraplength=900, justify="left").grid(row=5, column=0, columnspan=2, sticky="w", pady=(5, 0))
        ttk.Label(terminal, text="Mylar enabled").grid(row=6, column=0, sticky="w", padx=(0, 8), pady=2)
        ttk.Checkbutton(terminal, variable=self.values["mylar_enabled"], onvalue="1", offvalue="0", command=lambda: self._terminal_surface_changed("scintillator")).grid(row=6, column=1, sticky="w")
        ttk.Label(terminal, text="Mylar reflectivity (%)").grid(row=7, column=0, sticky="w", padx=(0, 8), pady=2)
        ttk.Entry(terminal, textvariable=self.values["mylar_reflectivity_percent"]).grid(row=7, column=1, sticky="ew", pady=2)
        self.values["mylar_enabled"].trace_add("write", lambda *_args: self._terminal_surface_changed("scintillator"))
        self.values["mylar_reflectivity_percent"].trace_add("write", lambda *_args: self._terminal_surface_changed("scintillator"))
        profile = ttk.Frame(terminal)
        profile.grid(row=8, column=0, columnspan=2, sticky="w", pady=(6, 0))
        # Legacy labels "Save Scintillator Geometry JSON" / "Load Scintillator Geometry JSON" remain recognized by v31 files; the visible commands
        # below are the complete 3D assembly commands for v32.
        ttk.Button(profile, text="Save Scintillator 3D Geometry JSON", command=lambda: self._save_assembly_model("scintillator")).pack(side="left")
        ttk.Button(profile, text="Load Scintillator 3D Geometry JSON", command=lambda: self._load_assembly_model("scintillator")).pack(side="left", padx=5)

    def _build_scintillator_sensor_tab_v13(self, parent: ttk.Frame) -> None:
        parent.columnconfigure(0, weight=1)
        self._build_sensor_common(parent, "scintillator", (
            ("Readout shape", "sensor_2_shape"), ("Readout width (mm)", "sensor_2_width_mm"), ("Readout height (mm)", "sensor_2_height_mm"), ("Readout depth (mm)", "sensor_2_depth_mm"),
            ("Quantum efficiency", "scintillator_readout_qe"), ("Acceptance angle (deg)", "scintillator_readout_acceptance_deg"), ("Optical window / coupling index", "scintillator_readout_index"),
            ("Collection efficiency", "scintillator_readout_collection_efficiency"), ("Timing/readout sigma (ns)", "sensor_2_timing_sigma_ns"),
        ))
        physics = ttk.LabelFrame(parent, text="Scintillator-specific generation and optical/readout physics", padding=7)
        physics.grid(row=1, column=0, sticky="ew", pady=(7, 0))
        physics.columnconfigure(1, weight=1)
        for row, (label, key) in enumerate((("Sampling mode", "scintillator_sampling_mode"), ("Photon yield / deposited MeV", "scintillator_yield_per_MeV"), ("Scintillation decay time (ns)", "scintillator_decay_time_ns"), ("Absorption length (mm)", "scintillator_absorption_length_mm"), ("Legacy maximum bounces (Run control overrides)", "scintillator_max_bounces"), ("Mylar enabled", "mylar_enabled"), ("Mylar reflectivity (%)", "mylar_reflectivity_percent"))):
            ttk.Label(physics, text=label).grid(row=row, column=0, sticky="w", padx=(0, 8), pady=2)
            if key in {"mylar_enabled"}:
                ttk.Checkbutton(physics, variable=self.values[key], onvalue="1", offvalue="0", command=lambda a="scintillator": self._sensor_status_changed(a)).grid(row=row, column=1, sticky="w")
            elif key == "scintillator_sampling_mode":
                ttk.Combobox(physics, textvariable=self.values[key], values=tuple(SCINTILLATION_SAMPLING_MODE_LABELS.values()), state="readonly").grid(row=row, column=1, sticky="ew")
            else:
                widget = ttk.Entry(physics, textvariable=self.values[key], state="readonly" if key == "scintillator_max_bounces" else "normal")
                widget.grid(row=row, column=1, sticky="ew")
        ttk.Label(physics, text="The legacy bounce value is retained for old detector profiles only. Runs use Run → Maximum photon bounces for Cherenkov, Scintillator, or All.", foreground="#245a35", wraplength=900, justify="left").grid(row=7, column=0, columnspan=2, sticky="w", pady=(4, 0))
        ttk.Label(physics, text="This is the exponential emission-delay mean. Units are ns; 0 means instantaneous scintillation emission with no added delay.", foreground="#245a35", wraplength=900, justify="left").grid(row=8, column=0, columnspan=2, sticky="w", pady=(2, 0))
        coupling = ttk.LabelFrame(parent, text="Scintillator derived optical readout", padding=7)
        coupling.grid(row=2, column=0, sticky="ew", pady=(7, 0))
        coupling.columnconfigure(1, weight=1)
        ttk.Label(coupling, text="Coupled component (derived status)").grid(row=0, column=0, sticky="w", padx=(0, 8), pady=2)
        self.scintillator_component_status = tk.StringVar()
        ttk.Label(coupling, textvariable=self.scintillator_component_status, foreground="#555555").grid(row=0, column=1, sticky="w")
        ttk.Label(coupling, text="Optical path (derived, read-only)").grid(row=1, column=0, sticky="w", padx=(0, 8), pady=2)
        self.scintillator_optical_path = tk.StringVar()
        ttk.Label(coupling, textvariable=self.scintillator_optical_path, foreground="#555555").grid(row=1, column=1, sticky="w")
        ttk.Label(coupling, text="Readout terminal / guide status (read-only)").grid(row=2, column=0, sticky="w", padx=(0, 8), pady=2)
        self.scintillator_sensor_position_text = tk.StringVar()
        ttk.Label(coupling, textvariable=self.scintillator_sensor_position_text, foreground="#245a35", wraplength=900, justify="left").grid(row=2, column=1, sticky="w", pady=(5, 0))
        ttk.Label(coupling, text="Terminal surface choices and optional ordered guide sections are edited only in Geometry Editor. This tab derives the complete path and keeps it read-only.", foreground="#245a35", wraplength=900, justify="left").grid(row=3, column=0, columnspan=2, sticky="w", pady=(5, 0))
        self._build_transport_diagnostic_panel(parent, "scintillator", 4)
        self._sensor_status_changed("scintillator")

    def _apply_ideal_sensor_settings(self, assembly: str) -> None:
        apply_ideal_sensor_settings(self.values, assembly)
        self._sensor_status_changed(assembly)
        self.status.set(f"{assembly.title()}: ideal sensor settings active")

    def _sensor_scene_status(self, assembly: str) -> tuple[dict, str]:
        model = self._assembly_model(assembly)
        scene = complete_local_assembly_scene(model, self.values, assembly=assembly)
        if scene.get("invalid_couplings"):
            return scene, str(scene.get("sensor_routing", {}).get("status", "Invalid derived optical readout configuration."))
        return scene, ""

    def _sensor_status_changed(self, assembly: str) -> None:
        if _get_instance_attr(self, "_suppress_assembly_refresh", False):
            self._pending_assembly_refresh[assembly] = True
            return
        try:
            diagnostic = optical_coupling_diagnostic(self.values, assembly)
            index = "1" if assembly == "cherenkov" else "2"
            scene, routing_error = self._sensor_scene_status(assembly)
            text = _get_instance_attr(self, f"{assembly}_coupling_diagnostic", None)
            status = _get_instance_attr(self, f"{assembly}_sensor_status", None)
            component_status = _get_instance_attr(self, f"{assembly}_component_status", None)
            position_text = _get_instance_attr(self, f"{assembly}_sensor_position_text", None)
            if text is not None:
                real, ideal = diagnostic["real"], diagnostic["ideal"]
                active = "IDEAL SENSOR SETTINGS ACTIVE" if diagnostic["ideal_active"] else "ideal reference only; real settings are unchanged"
                text.set(f"REAL {real}  |  IDEAL REFERENCE {ideal}  |  {active}")
            if component_status is not None:
                component_status.set(scene.get("sensor_routing", {}).get("derived_component", "No photosensor readout"))
            optical_path = _get_instance_attr(self, f"{assembly}_optical_path", None)
            if optical_path is not None:
                optical_path.set(" → ".join(scene.get("optical_path", ())) or "No valid photosensor readout")
            if position_text is not None:
                readout = scene.get("sensor_routing", {})
                terminal = readout.get("readout_terminal") or "none"
                position_text.set(f"Selected readout terminal: {terminal}; guide configuration: {readout.get('guide_configuration', 'absent')}")
            if status is not None:
                if routing_error:
                    status.set(routing_error)
                else:
                    defaults_note = "; configured sensor defaults used" if scene.get("sensor_routing", {}).get("sensor_defaults_used") else ""
                    status.set(f"{assembly.title()} sensor/readout metadata valid; canonical placement refreshed" + defaults_note + ("; ideal sensor settings active" if diagnostic["ideal_active"] else ""))
            self._assembly_draw_preview(assembly)
            self._refresh_geometry_views()
            if assembly == "scintillator":
                self._refresh_scintillator_preview()
            self._refresh_global_assembly_preview()
            # Keep an immediate renderer/transport-facing snapshot in sync;
            # the same profile is copied into RunSettings on execution.
            self.sensor_transport_metadata = sensor_profile_from_values(self.values, self.assembly_models)
            self._refresh_transport_diagnostic(assembly)
        except (TypeError, ValueError) as exc:
            status = _get_instance_attr(self, f"{assembly}_sensor_status", None)
            if status is not None:
                status.set(f"Invalid {assembly} sensor settings: {exc}")

    def _refresh_transport_diagnostic(self, assembly: str) -> None:
        """Refresh only the selected detector's read-only transport report."""

        target = _get_instance_attr(self, f"{assembly}_transport_diagnostic_text", None)
        if target is None:
            return
        try:
            report = detector_transport_diagnostic(self.values, self._assembly_model(assembly), assembly)
            text = detector_transport_diagnostic_text(report)
            target.set(text)
            setattr(self, f"{assembly}_transport_diagnostic_report", report)
            run_target = _get_instance_attr(self, f"run_{assembly}_transport_diagnostic_text", None)
            if run_target is not None:
                run_target.set(text)
        except (TypeError, ValueError) as exc:
            target.set(f"{assembly.title()} Transport Diagnostic (read-only)\nWarnings/errors: invalid/stale configuration: {exc}")

    def _compare_transport_diagnostic(self, assembly: str) -> None:
        """Run the non-mutating real-versus-ideal comparison for one detector."""

        self._sensor_status_changed(assembly)
        self.status.set(f"{assembly.title()} real vs ideal coupling compared (temporary diagnostic; settings unchanged)")

    def _copy_transport_diagnostic(self, assembly: str) -> None:
        text_var = _get_instance_attr(self, f"{assembly}_transport_diagnostic_text", None)
        if text_var is None:
            return
        try:
            self.clipboard_clear()
            self.clipboard_append(text_var.get())
            self.update()
            self.status.set(f"{assembly.title()} transport diagnostic copied")
        except tk.TclError:
            self.status.set("Clipboard unavailable")

    def _export_transport_diagnostic(self, assembly: str) -> None:
        text_var = _get_instance_attr(self, f"{assembly}_transport_diagnostic_text", None)
        if text_var is None:
            return
        path = filedialog.asksaveasfilename(title=f"Export {assembly.title()} transport diagnostic", defaultextension=".txt", filetypes=(("Diagnostic text", "*.txt"), ("All files", "*.*")))
        if path:
            try:
                Path(path).write_text(text_var.get() + "\n", encoding="utf-8")
                self.status.set(f"{assembly.title()} transport diagnostic exported: {Path(path).name}")
            except OSError as exc:
                messagebox.showerror("Export diagnostic", str(exc), parent=self)

    def _build_transport_diagnostic_panel(self, parent: ttk.Frame, assembly: str, row: int) -> None:
        """Build a detector-specific read-only diagnostic section."""

        title = "Cherenkov Transport Diagnostic" if assembly == "cherenkov" else "Scintillator Transport Diagnostic"
        frame = ttk.LabelFrame(parent, text=title, padding=7)
        frame.grid(row=row, column=0, sticky="ew", pady=(8, 0))
        frame.columnconfigure(0, weight=1)
        text = tk.StringVar(value="Preparing detector-specific transport diagnostic…")
        setattr(self, f"{assembly}_transport_diagnostic_text", text)
        ttk.Label(frame, textvariable=text, justify="left", wraplength=980, foreground="#17324d").grid(row=0, column=0, sticky="ew")
        buttons = ttk.Frame(frame)
        buttons.grid(row=1, column=0, sticky="w", pady=(6, 0))
        ttk.Button(buttons, text=f"Refresh {assembly.title()} Diagnostic", command=lambda a=assembly: self._refresh_transport_diagnostic(a)).pack(side="left")
        ttk.Button(buttons, text="Compare real vs ideal coupling", command=lambda a=assembly: self._compare_transport_diagnostic(a)).pack(side="left", padx=5)
        ttk.Button(buttons, text="Copy diagnostic text", command=lambda a=assembly: self._copy_transport_diagnostic(a)).pack(side="left", padx=5)
        ttk.Button(buttons, text="Export diagnostic text", command=lambda a=assembly: self._export_transport_diagnostic(a)).pack(side="left")
        ttk.Label(frame, text="Read-only local geometry, sensor/readout, optical path, terminal/Mylar status, warnings/errors, and detector-specific physical generation model. Actions are temporary diagnostics.", foreground="#555555", wraplength=980, justify="left").grid(row=2, column=0, sticky="w", pady=(5, 0))
        self._refresh_transport_diagnostic(assembly)

    def _save_assembly_model(self, assembly: str) -> None:
        path = filedialog.asksaveasfilename(title=f"Save {assembly.title()} 3D geometry", defaultextension=".json", filetypes=(("Complete 3D Geometry JSON", "*.json"), ("All files", "*.*")))
        if path:
            try:
                model = self._assembly_model(assembly)
                sensor = complete_sensor_readout_from_values(self.values, assembly, model)
                report = detector_transport_diagnostic(self.values, model, assembly)
                payload = complete_detector_assembly_geometry_to_dict(
                    model, assembly, sensor,
                    local_settings={
                        "coordinate_system": "assembly-local millimetres",
                        "global_placement_offsets_excluded": True,
                        "terminal_derived_sensor_placement": True,
                        "sensor_transform_source": "canonical terminal/guide output; recompute on load",
                    },
                    derived_diagnostics={"optical_transport_path": report["optical_transport_path"], "assembly_valid": report["assembly_valid"]},
                )
                Path(path).write_text(json.dumps(payload, indent=2, sort_keys=True) + "\n", encoding="utf-8")
                self.last_geometry_io = {"operation": "save", "assembly": assembly, "path": str(Path(path).resolve()), "format": payload["format"]}
                self.status.set(f"{assembly.title()} 3D geometry saved: {Path(path).name}")
            except (OSError, TypeError, ValueError) as exc:
                self.status.set(f"Could not save {assembly.title()} geometry: {exc}")
                messagebox.showerror("Save assembly", str(exc), parent=self)

    def _load_assembly_model(self, assembly: str) -> None:
        path = filedialog.askopenfilename(title=f"Load {assembly.title()} 3D geometry", filetypes=(("Complete 3D Geometry JSON", "*.json"), ("Assembly JSON", "*.json"), ("All files", "*.*")))
        if path:
            try:
                if not _has_instance_attr(self, "_pending_assembly_refresh"):
                    self._pending_assembly_refresh = {}
                # Parse, validate, and derive the complete candidate before
                # touching the live model or any Tk variable.  This is the
                # transaction boundary for malformed/incompatible files.
                payload = json.loads(Path(path).read_text(encoding="utf-8"))
                complete_sensor = None
                if isinstance(payload, Mapping) and payload.get("format") == "complete-detector-assembly-geometry":
                    model, complete_sensor = complete_detector_assembly_geometry_from_dict(payload, assembly, name=assembly.title())
                else:
                    # Existing v27 typed files remain readable.  They do not
                    # carry sensor fields, so current detector-specific
                    # readout settings remain unchanged on this compatibility
                    # path.
                    model = assembly_geometry_from_dict(payload, assembly, name=assembly.title())
                candidate_values = {key: _value_text(self.values, key, "") for key in self.values}
                if complete_sensor is not None:
                    self._apply_complete_sensor_candidate(candidate_values, assembly, complete_sensor)
                model.complete_scene(
                    sensor_label="Photosensor" if assembly == "cherenkov" else "Scintillator readout",
                    sensor_automatic_position=_bool_text(candidate_values, f"sensor_{'1' if assembly == 'cherenkov' else '2'}_automatic_position"),
                    sensor_shape=_value_text(candidate_values, f"sensor_{'1' if assembly == 'cherenkov' else '2'}_shape", "rectangular"),
                    sensor_width_mm=float(_value_text(candidate_values, f"sensor_{'1' if assembly == 'cherenkov' else '2'}_width_mm", "20")),
                    sensor_height_mm=float(_value_text(candidate_values, f"sensor_{'1' if assembly == 'cherenkov' else '2'}_height_mm", "10")),
                )
                model.to_simulation_geometry().validate()
                if complete_sensor is not None:
                    # Validate the full candidate readout before the first
                    # live Tk variable is touched.  The returned placement is
                    # derived again below and is never trusted from JSON.
                    complete_sensor_readout_from_values(candidate_values, assembly, model)
                # Commit only after every validation above succeeds.  Hold
                # all StringVar traces while the model and complete sensor
                # snapshot cross the transaction boundary together.
                self._suppress_assembly_refresh = True
                try:
                    self.assembly_models[assembly] = model
                    if complete_sensor is not None:
                        self._commit_candidate_values(candidate_values)
                    guide_key = "cherenkov_guide_enabled" if assembly == "cherenkov" else "scintillator_guide_enabled"
                    if guide_key in self.values:
                        self.values[guide_key].set("1" if model.guide_enabled else "0")
                    if assembly == "cherenkov":
                        self.geometry_sections = model.sections
                        self.values["terminal_plus_surface"].set(terminal_surface_label(model.terminal_plus_surface))
                        self.values["terminal_minus_surface"].set(terminal_surface_label(model.terminal_minus_surface))
                        # The canonical editor owns both endpoints; the older
                        # Run control is a compatibility view of -X.  Traces
                        # are suppressed during this transaction, so refresh
                        # the alias explicitly from the validated candidate.
                        if "opposite_end_mode" in self.values:
                            self.values["opposite_end_mode"].set(
                                normalize_terminal_surface(model.terminal_minus_surface)
                            )
                        self.values["cherenkov_mylar_enabled"].set("1" if model.mylar_enabled else "0")
                        self.values["cherenkov_mylar_reflectivity_percent"].set(str(model.mylar_reflectivity * 100.0))
                        self.values["entry_x_mm"].set(str(model.impact_x_mm))
                        self.values["entry_y_mm"].set(str(model.impact_y_mm))
                    else:
                        self.scintillator_geometry_sections = model.sections
                        self.values["scintillator_terminal_plus_surface"].set(terminal_surface_label(model.terminal_plus_surface))
                        self.values["scintillator_terminal_minus_surface"].set(terminal_surface_label(model.terminal_minus_surface))
                        self.values["mylar_enabled"].set("1" if model.mylar_enabled else "0")
                        self.values["mylar_reflectivity_percent"].set(str(model.mylar_reflectivity * 100.0))
                        self.values["scintillator_impact_x_mm"].set(str(model.impact_x_mm))
                        self.values["scintillator_impact_y_mm"].set(str(model.impact_y_mm))
                        self._sync_scintillator_model_to_values()
                finally:
                    self._suppress_assembly_refresh = False
                self._pending_assembly_refresh[assembly] = False
                self._assembly_refresh_tree(assembly)
                self._assembly_build_properties(assembly)
                # Refresh the legacy -X Run-panel labels/enable states after
                # the atomic commit.  The canonical local/global viewers are
                # rebuilt below from the same model snapshot.
                if assembly == "cherenkov" and _has_instance_attr(self, "detector_mode_text"):
                    self._mode_changed()
                # One post-commit fan-out rebuilds the local three-view scene,
                # interactive/local 3D, and global preview from the same
                # canonical snapshot.  Its idle Configure/Map fallback runs
                # again after the notebook/scroll frame has a real size.
                self._sensor_status_changed(assembly)
                self._refresh_loaded_assembly_views(assembly)
                self.last_geometry_io = {"operation": "load", "assembly": assembly, "path": str(Path(path).resolve()), "format": payload.get("format", "quartz-assembly-geometry")}
                self.status.set(f"{assembly.title()} 3D geometry loaded: {Path(path).name} • canonical placement recomputed")
            except (OSError, json.JSONDecodeError, TypeError, ValueError) as exc:
                self.status.set(f"Could not load {assembly.title()} geometry; current values preserved: {exc}")
                messagebox.showerror("Load assembly", str(exc), parent=self)

    @staticmethod
    def _apply_complete_sensor_candidate(candidate: Dict[str, object], assembly: str, sensor: Mapping[str, object]) -> None:
        """Copy complete-file readout fields into an isolated candidate map."""

        index = "1" if assembly == "cherenkov" else "2"
        dimensions = sensor["dimensions_mm"]
        mapping = {
            f"sensor_{index}_shape": sensor["shape"],
            f"sensor_{index}_width_mm": dimensions["width"],
            f"sensor_{index}_height_mm": dimensions["height"],
            f"sensor_{index}_depth_mm": dimensions["depth"],
            f"sensor_{index}_radius_mm": sensor.get("readout_radius_mm", min(float(dimensions["width"]), float(dimensions["height"])) / 2.0),
            f"sensor_{index}_qe": sensor["quantum_efficiency"],
            f"sensor_{index}_acceptance_deg": sensor["acceptance_deg"],
            f"sensor_{index}_index": sensor["optical_window_coupling_index"],
            f"sensor_{index}_collection": sensor["collection_efficiency"],
            f"sensor_{index}_timing_sigma_ns": sensor["timing_readout"]["sigma_ns"],
            f"sensor_{index}_automatic_position": "1" if sensor["automatic_position"] else "0",
        }
        axes = sensor.get("readout_semiaxes_mm", {})
        if isinstance(axes, Mapping):
            mapping[f"sensor_{index}_semi_axis_a_mm"] = axes.get("a", float(dimensions["width"]) / 2.0)
            mapping[f"sensor_{index}_semi_axis_b_mm"] = axes.get("b", float(dimensions["height"]) / 2.0)
        if sensor.get("position_mm") is not None:
            for suffix, value in zip(("position_x_mm", "position_y_mm", "position_z_mm"), sensor["position_mm"]):
                mapping[f"sensor_{index}_{suffix}"] = value
        for key, value in mapping.items():
            candidate[key] = str(value)
        specific = sensor.get("detector_specific_response", {})
        if assembly == "scintillator" and isinstance(specific, Mapping):
            for target, source in (("scintillator_sampling_mode", "sampling_mode"), ("scintillator_yield_per_MeV", "yield_per_deposited_MeV"), ("scintillator_decay_time_ns", "decay_time_ns"), ("scintillator_max_bounces", "max_bounces")):
                if source in specific:
                    candidate[target] = str(specific[source])

    def _commit_candidate_values(self, candidate: Mapping[str, object]) -> None:
        """Commit an already validated candidate value map atomically enough for Tk."""

        for key, value in candidate.items():
            if key in self.values:
                self.values[key].set(str(value))

    def _build_global_assembly_tab_v13(self, parent: ttk.Frame) -> None:
        """Global scene: transforms and inherited summaries only."""

        parent.columnconfigure(0, weight=0)
        parent.columnconfigure(1, weight=1)
        parent.rowconfigure(0, weight=1)
        controls = ttk.LabelFrame(parent, text="Rigid transforms for complete assemblies", padding=7)
        controls.grid(row=0, column=0, sticky="new", padx=(0, 8))
        controls.columnconfigure(1, weight=1)
        ttk.Label(controls, text="Manual fields here edit only rigid transforms and beam-relative placement. The placement-only JSON below is retained for compatibility; the complete-profile buttons restore both local assemblies, guides, readouts, terminals, Mylar, impacts, and transforms.", foreground="#245a35", wraplength=390, justify="left").grid(row=0, column=0, columnspan=2, sticky="w", pady=(0, 7))
        fields = (
            ("Cherenkov assembly translation X (mm)", "global_offset_x_mm"), ("Cherenkov assembly translation Y (mm)", "global_offset_y_mm"), ("Cherenkov assembly translation Z (mm)", "global_offset_z_mm"),
            ("Cherenkov rotation axis", "assembly_1_rotation_axis"), ("Cherenkov rotation angle (deg)", "assembly_1_rotation_angle_deg"),
            ("Scintillator assembly translation X (mm)", "detector2_offset_x_mm"), ("Scintillator assembly translation Y (mm)", "detector2_offset_y_mm"), ("Scintillator assembly translation Z (mm)", "detector2_offset_z_mm"),
            ("Scintillator rotation axis", "assembly_2_rotation_axis"), ("Scintillator rotation angle (deg)", "assembly_2_rotation_angle_deg"),
        )
        for row, (label, key) in enumerate(fields, start=1):
            ttk.Label(controls, text=label).grid(row=row, column=0, sticky="w", padx=(0, 8), pady=2)
            if key.endswith("rotation_axis"):
                widget = ttk.Combobox(controls, textvariable=self.values[key], values=("None", "X", "Y", "Z"), state="readonly")
            else:
                widget = ttk.Entry(controls, textvariable=self.values[key])
            widget.grid(row=row, column=1, sticky="ew", pady=2)
            self.values[key].trace_add("write", lambda *_args: self._refresh_global_assembly_preview())
        summary = tk.StringVar(value="")
        self.global_assembly_inherited_summary = summary
        ttk.Label(controls, text="Inherited read-only assembly summary", foreground="#555555").grid(row=11, column=0, columnspan=2, sticky="w", pady=(8, 2))
        ttk.Label(controls, textvariable=summary, justify="left", wraplength=370).grid(row=12, column=0, columnspan=2, sticky="w")
        profile_row = ttk.LabelFrame(controls, text="Complete assembly profile files", padding=5)
        profile_row.grid(row=13, column=0, columnspan=2, sticky="ew", pady=(8, 0))
        profile_row.columnconfigure(0, weight=1)
        profile_row.columnconfigure(1, weight=1)
        ttk.Button(profile_row, text="Save complete assembly profile", command=self._save_complete_global_assembly_profile).grid(row=0, column=0, sticky="ew", padx=(0, 4))
        ttk.Button(profile_row, text="Load complete assembly profile", command=self._load_complete_global_assembly_profile).grid(row=0, column=1, sticky="ew", padx=(4, 0))
        ttk.Label(
            profile_row,
            text="Saves/restores both canonical assemblies, guides, terminal/Mylar faces, sensors, impacts, and transforms.",
            foreground="#245a35", wraplength=380, justify="left",
        ).grid(row=1, column=0, columnspan=2, sticky="w", pady=(4, 0))
        ttk.Button(controls, text="Open interactive global 3D", command=self._open_global_assembly_3d_viewer).grid(row=14, column=0, columnspan=2, sticky="w", pady=(8, 0))
        preview = ttk.LabelFrame(parent, text="Global interactive assembly scene", padding=5)
        preview.grid(row=0, column=1, sticky="nsew")
        preview.columnconfigure(0, weight=1)
        preview.columnconfigure(1, weight=0)
        preview.rowconfigure(1, weight=1)
        self.global_preview_text = tk.StringVar(value="Preparing global assembly…")
        self.global_geometry_status = tk.StringVar(value="Not validated")
        ttk.Label(preview, textvariable=self.global_preview_text, justify="left", wraplength=760).grid(row=0, column=0, columnspan=2, sticky="w")
        self.global_preview_canvas = tk.Canvas(preview, background="white", width=700, height=520, highlightthickness=1, highlightbackground="#b0b0b0")
        self.global_preview_canvas.grid(row=1, column=0, sticky="nsew", pady=(5, 0))
        # Re-fit the same camera on tab/window resize so the two assemblies,
        # sensors, and both beam endpoints remain in frame.
        self.global_preview_canvas.bind("<Configure>", self._refresh_global_assembly_preview)
        details = ttk.LabelFrame(preview, text="Scene details / legend", padding=6)
        details.grid(row=1, column=1, sticky="ns", padx=(8, 0), pady=(5, 0))
        details.columnconfigure(0, weight=1)
        details.rowconfigure(0, weight=1)
        try:
            details_text = tk.Text(details, width=36, height=28, wrap="word", relief="flat", borderwidth=0,
                                   background=_ttk_background(details), takefocus=True)
            details_text.grid(row=0, column=0, sticky="nsew")
            details_scroll = ttk.Scrollbar(details, orient="vertical", command=details_text.yview)
            details_scroll.grid(row=0, column=1, sticky="ns")
            details_text.configure(yscrollcommand=details_scroll.set, state="disabled")
            self.global_preview_info_text = details_text
        except Exception:
            # Keep the geometry canvas usable on unusual Tk theme/Text builds.
            self.global_preview_info_text = None
        self._refresh_global_assembly_preview()

    def _assembly_transform(self, number: int) -> RigidAssemblyTransform:
        prefix = "assembly_1" if number == 1 else "assembly_2"
        transform = RigidAssemblyTransform(
            translation_x_mm=float(_value_text(self.values, "global_offset_x_mm" if number == 1 else "detector2_offset_x_mm", "0")),
            translation_y_mm=float(_value_text(self.values, "global_offset_y_mm" if number == 1 else "detector2_offset_y_mm", "0")),
            translation_z_mm=float(_value_text(self.values, "global_offset_z_mm" if number == 1 else "detector2_offset_z_mm", "0")),
            rotation_axis=_value_text(self.values, f"{prefix}_rotation_axis", "None"),
            rotation_angle_deg=float(_value_text(self.values, f"{prefix}_rotation_angle_deg", "0")),
        )
        transform.validate()
        return transform
        if model.guide_sections:
            section = model.guide_sections[0]
            for suffix in ("shape", "length_x_mm", "width_y_mm", "thickness_z_mm", "width_y_start_mm", "width_y_end_mm", "thickness_z_start_mm", "thickness_z_end_mm", "axis", "radius_mm", "length_mm", "semi_axis_a_mm", "semi_axis_b_mm", "radius_start_mm", "radius_end_mm", "center_x_mm", "center_y_mm", "center_z_mm", "material"):
                key = f"guide_{suffix}"
                if key in self.values:
                    self.values[key].set(str(getattr(section, suffix)))

    def _build_scintillator_geometry_tab(self, parent: ttk.Frame) -> None:
        """Edit only local scintillator/guide solids and their local preview."""

        parent.columnconfigure(0, weight=3)
        parent.columnconfigure(1, weight=2)
        left = ttk.Frame(parent)
        left.grid(row=0, column=0, sticky="nsew", padx=(0, 10))
        left.columnconfigure(1, weight=1)
        controls = ttk.LabelFrame(left, text="Local scintillator and light-guide geometry", padding=6)
        controls.grid(row=0, column=0, columnspan=2, sticky="new")
        controls.columnconfigure(1, weight=1)
        rows = [
            ("Scintillator shape", "scintillator_shape", "shape"),
            ("Scintillator material", "scintillator_material", "material"),
            ("Scintillator refractive index", "scintillator_refractive_index", "entry"),
            ("Scintillator length X / width Y / thickness Z (mm)", "scintillator_length_x_mm", "entry"),
            ("Scintillator radius / length / semi-axes (mm)", "scintillator_radius_mm", "entry"),
            ("Scintillator cylinder length (mm)", "scintillator_length_mm", "entry"), ("Scintillator semi-axis A (mm)", "scintillator_semi_axis_a_mm", "entry"), ("Scintillator semi-axis B (mm)", "scintillator_semi_axis_b_mm", "entry"),
            ("Guide shape", "guide_shape", "shape"),
            ("Guide material", "guide_material", "material"),
            ("Guide refractive index", "guide_refractive_index", "entry"),
            ("Guide length X / width Y / thickness Z (mm)", "guide_length_x_mm", "entry"),
            ("Guide radius / length / semi-axes (mm)", "guide_radius_mm", "entry"),
            ("Guide cylinder length (mm)", "guide_length_mm", "entry"), ("Guide semi-axis A (mm)", "guide_semi_axis_a_mm", "entry"), ("Guide semi-axis B (mm)", "guide_semi_axis_b_mm", "entry"),
        ]
        rows += [
            ("Scintillator width Y (mm)", "scintillator_width_y_mm", "entry"), ("Scintillator thickness Z (mm)", "scintillator_thickness_z_mm", "entry"),
            ("Scintillator taper Y start (mm)", "scintillator_width_y_start_mm", "entry"), ("Scintillator taper Y end (mm)", "scintillator_width_y_end_mm", "entry"),
            ("Scintillator taper Z start (mm)", "scintillator_thickness_z_start_mm", "entry"), ("Scintillator taper Z end (mm)", "scintillator_thickness_z_end_mm", "entry"),
            ("Scintillator axis", "scintillator_axis", "axis"),
            ("Scintillator cone start radius (mm)", "scintillator_radius_start_mm", "entry"), ("Scintillator cone end radius (mm)", "scintillator_radius_end_mm", "entry"),
            ("Guide width Y (mm)", "guide_width_y_mm", "entry"), ("Guide thickness Z (mm)", "guide_thickness_z_mm", "entry"),
            ("Guide taper Y start (mm)", "guide_width_y_start_mm", "entry"), ("Guide taper Y end (mm)", "guide_width_y_end_mm", "entry"),
            ("Guide taper Z start (mm)", "guide_thickness_z_start_mm", "entry"), ("Guide taper Z end (mm)", "guide_thickness_z_end_mm", "entry"),
            ("Guide axis", "guide_axis", "axis"),
            ("Guide cone start radius (mm)", "guide_radius_start_mm", "entry"), ("Guide cone end radius (mm)", "guide_radius_end_mm", "entry"),
        ]
        for row, (label, key, kind) in enumerate(rows):
            ttk.Label(controls, text=label).grid(row=row, column=0, sticky="w", padx=(0, 8), pady=2)
            if kind == "bool":
                widget = ttk.Checkbutton(controls, variable=self.values[key], command=self._refresh_scintillator_preview)
            elif kind == "shape":
                widget = ttk.Combobox(controls, textvariable=self.values[key], values=tuple(GEOMETRY_SHAPE_LABELS.keys()), state="readonly")
            elif kind == "material":
                widget = ttk.Combobox(controls, textvariable=self.values[key], values=("scintillator", "light_guide", "fused silica", "air"), state="readonly")
            elif kind == "scintmode":
                widget = ttk.Combobox(controls, textvariable=self.values[key], values=tuple(SCINTILLATION_SAMPLING_MODE_LABELS.values()), state="readonly")
            elif kind == "axis":
                widget = ttk.Combobox(controls, textvariable=self.values[key], values=GEOMETRY_AXES, state="readonly")
            else:
                widget = ttk.Entry(controls, textvariable=self.values[key])
            widget.grid(row=row, column=1, sticky="ew", pady=2)
            if kind != "bool":
                self.values[key].trace_add("write", lambda *_args: self._scintillator_local_value_changed())
        ttk.Label(controls, text="LOCAL GEOMETRY only: boxes/tapers use X/Y/Z dimensions; cylinders, elliptical cylinders, and truncated cones use radius/axis fields. The two solids are coupled locally along +X. There is deliberately no quartz gap, global offset, or beam-monitor placement control on this page. Mylar response is configured in Scintillator Detector.", foreground="#245a35", wraplength=620, justify="left").grid(row=len(rows), column=0, columnspan=2, sticky="w", pady=(7, 2))
        profile_row = ttk.Frame(controls)
        profile_row.grid(row=len(rows) + 1, column=0, columnspan=2, sticky="w", pady=(5, 0))
        ttk.Button(profile_row, text="Load Local Geometry JSON", command=self._load_local_geometry_profile).pack(side="left")
        ttk.Button(profile_row, text="Save Local Geometry JSON", command=self._save_local_geometry_profile).pack(side="left", padx=(5, 0))
        preview = ttk.LabelFrame(parent, text="Local views: scintillator / guide / interface / Mylar", padding=6)
        preview.grid(row=0, column=1, sticky="nsew")
        preview.columnconfigure(0, weight=1)
        preview.rowconfigure(2, weight=1)
        self.scintillator_preview_text = tk.StringVar(value="Local geometry not validated yet")
        ttk.Label(preview, textvariable=self.scintillator_preview_text, wraplength=520, justify="left").grid(row=0, column=0, sticky="nw")
        view_buttons = ttk.Frame(preview)
        view_buttons.grid(row=1, column=0, sticky="w", pady=(7, 3))
        for key, label in (("xy", "View XY"), ("xz", "View XZ"), ("yz", "View YZ"), ("isometric", "View isometric")):
            ttk.Button(view_buttons, text=label, command=lambda value=key: self._set_local_preview_view(value)).pack(side="left", padx=(0, 4))
        ttk.Button(view_buttons, text="Open Local 3D Viewer", command=self._open_local_3d_viewer).pack(side="left", padx=(8, 0))
        self.scintillator_preview_canvas = tk.Canvas(preview, height=420, background="white", highlightthickness=1, highlightbackground="#b0b0b0")
        self.scintillator_preview_canvas.grid(row=2, column=0, sticky="nsew", pady=(5, 0))

    def _refresh_scintillator_preview(self, *_args) -> None:
        # The active v24 Scintillator Geometry tab uses the complete shared
        # three-view renderer. Keep legacy callers routed to that same canvas
        # so a refresh cannot paint a stale single-view layout over it.
        shared_state = _get_instance_attr(self, "_assembly_editor_state", {}).get("scintillator")
        if shared_state and shared_state.get("preview") is _get_instance_attr(self, "scintillator_preview_canvas", None):
            self._assembly_draw_preview("scintillator")
            return
        if not _has_instance_attr(self, "scintillator_preview_canvas"):
            return
        canvas = self.scintillator_preview_canvas
        canvas.delete("all")
        try:
            scene = local_geometry_scene(self.values)
            local = scintillator_local_geometry_from_values(self.values)
            faces_by_solid = (scene["scintillator_faces_mm"],) + ((scene["guide_faces_mm"],) if scene.get("guide_present") else ())
            points = [point for solid_faces in faces_by_solid for face in solid_faces for point in face]
            lows = [min(point[index] for point in points) for index in range(3)]
            highs = [max(point[index] for point in points) for index in range(3)]
            width, height = max(canvas.winfo_width(), 560), max(canvas.winfo_height(), 360)
            view = self.local_preview_view
            if view == "xy":
                axes, labels = (0, 1), ("X (mm)", "Y (mm)")
            elif view == "xz":
                axes, labels = (0, 2), ("X (mm)", "Z (mm)")
            elif view == "yz":
                axes, labels = (1, 2), ("Y (mm)", "Z (mm)")
            else:
                axes, labels = None, ("isometric X/Y/Z", "")
            if axes is not None:
                span_x = max(highs[axes[0]] - lows[axes[0]], 1.0)
                span_y = max(highs[axes[1]] - lows[axes[1]], 1.0)
                scale = min(width - 90, height - 90) / max(span_x, span_y)
                def project(point):
                    return (45 + (point[axes[0]] - lows[axes[0]]) * scale,
                            height - 45 - (point[axes[1]] - lows[axes[1]]) * scale)
            else:
                center = tuple((lows[index] + highs[index]) / 2.0 for index in range(3))
                extent = max(max(highs[index] - lows[index] for index in range(3)), 1.0)
                scale = min(width - 100, height - 100) * 0.72 / extent
                def project(point):
                    x, y, z = (point[index] - center[index] for index in range(3))
                    return (width / 2 + (x - 0.65 * y) * scale,
                            height / 2 + (0.55 * z - 0.35 * y) * scale)
            colors = ("#80cbc4", "#90caf9")
            names = ("SCINTILLATOR",) + (("LIGHT GUIDE",) if scene.get("guide_present") else ())
            mylar_indices = (set(scene["scintillator_mylar_face_indices"]),) + ((set(scene["guide_mylar_face_indices"]),) if scene.get("guide_present") else ())
            for index, solid_faces in enumerate(faces_by_solid):
                for face_index, face in enumerate(solid_faces):
                    projected = [project(point) for point in face]
                    canvas.create_polygon([coord for point in projected for coord in point],
                                          fill=colors[index], outline="#1f4e79", stipple="gray50")
                    canvas.create_line([coord for point in projected + [projected[0]] for coord in point],
                                       fill="#263238", width=2)
                    if _bool_text(self.values, "mylar_enabled") and face_index in mylar_indices[index]:
                        canvas.create_line([coord for point in projected + [projected[0]] for coord in point],
                                           fill="#b26a00", dash=(5, 3), width=2)
                if solid_faces:
                    label_point = project(solid_faces[0][0])
                    canvas.create_text(label_point[0], label_point[1], text=names[index],
                                       fill="#17324d", anchor="sw")
            interface_x = scene["coupling_interface_x_mm"]
            if axes is not None and axes[0] == 0:
                ix = 45 + (interface_x - lows[0]) * scale
                canvas.create_line(ix, 30, ix, height - 45, fill="#7b1fa2", dash=(3, 2), width=2)
                canvas.create_text(ix + 4, 32, text="coupling interface (unwrapped)", anchor="nw", fill="#7b1fa2")
            canvas.create_text(
                45, 18, anchor="w",
                text=f"LOCAL {view.upper()} view • gold striped Mylar ({scene['mylar_reflectivity'] * 100:g}%) • external faces only",
                fill="#245a35",
            )
            canvas.create_text(width - 12, height - 28, anchor="e", text=labels[0], fill="#333333")
            self.scintillator_preview_text.set(
                f"Valid local geometry ({view.upper()}). "
                + ("Scintillator and guide are separate volumes; coupling and readout/output faces are interfaces and remain unwrapped. " if scene.get("guide_present") else "Active scintillator only; the terminal surface is the readout/output face. ")
                + "Global gap/offset is intentionally excluded."
            )
        except (TypeError, ValueError) as exc:
            self.scintillator_preview_text.set(f"Invalid local geometry: {exc}")
            canvas.create_text(max(canvas.winfo_width(), 560) / 2, 180, text="Invalid local geometry", fill="#b00020")
        self._refresh_scintillator_viewer()

    def _set_local_preview_view(self, view: str) -> None:
        if view not in {"xy", "xz", "yz", "isometric"}:
            raise ValueError("unknown local preview view")
        self.local_preview_view = view
        self._refresh_scintillator_preview()

    def _local_3d_projection(self, point, center, scale, width, height):
        x, y, z = (float(point[index]) - center[index] for index in range(3))
        yaw, pitch = self.local_viewer_yaw, self.local_viewer_pitch
        cy, sy, cp, sp = math.cos(yaw), math.sin(yaw), math.cos(pitch), math.sin(pitch)
        rx, ry = cy * x - sy * y, sy * x + cy * y
        screen_y = cp * ry - sp * z
        depth = sp * ry + cp * z
        return width / 2 + rx * scale, height / 2 - screen_y * scale, depth

    def _draw_local_3d_view(self) -> None:
        canvas = self.local_viewer_canvas
        if canvas is None or not canvas.winfo_exists():
            return
        canvas.delete("all")
        try:
            data = self._canonical_local_view_data("scintillator")
            assembly_scene = data["scene"]
            faces = data["faces"]
            edges = data["edges"]
            centers = data["centers"]
            sensor_labels = data["sensor_labels"]
            sensor_windows = data["sensor_windows"]
            mirror_labels = data["mirror_labels"]
            center, extent = data["center"], data["extent"]
            self.local_viewer_scene = assembly_scene
        except (TypeError, ValueError) as exc:
            canvas.create_text(300, 200, text=f"Invalid local geometry: {exc}", fill="#b00020")
            return
        width, height = max(canvas.winfo_width(), 600), max(canvas.winfo_height(), 420)
        scale = min(width, height) * 0.40 / extent * self.local_viewer_zoom
        projected = []
        for face_index, (_average, face, component_index, color, shade) in enumerate(faces):
            points_2d = [self._local_3d_projection(point, center, scale, width, height) for point in face]
            projected.append((sum(item[2] for item in points_2d) / len(points_2d),
                              face_index, face, component_index, color, shade, points_2d))
        for _depth, face_index, face, _component_index, color, shade, polygon in sorted(projected, key=lambda item: item[0]):
            canvas.create_polygon([coord for point in polygon for coord in point[:2]],
                                  fill=self._geometry_viewer_color(color, shade), outline="#263238", stipple="gray50")
            canvas.create_line([coord for point in polygon + [polygon[0]] for coord in point[:2]],
                               fill="#263238", width=2)
            if _bool_text(self.values, "mylar_enabled") and id(face) in data["mylar_face_ids"]:
                canvas.create_line([coord for point in polygon + [polygon[0]] for coord in point[:2]],
                                   fill="#b26a00", dash=(5, 3), width=2)
        for mylar_face, mylar_label in data.get("mylar_faces", ()):
            polygon = [self._local_3d_projection(point, center, scale, width, height) for point in mylar_face]
            canvas.create_polygon([coord for point in polygon for coord in point[:2]],
                                  fill=MYLAR_DISPLAY_COLOR, outline=MYLAR_DISPLAY_EDGE,
                                  stipple="gray25")
            canvas.create_line([coord for point in polygon + [polygon[0]] for coord in point[:2]],
                               fill=MYLAR_DISPLAY_EDGE, dash=(5, 3), width=2)
            if polygon:
                canvas.create_text(polygon[0][0] + 3, polygon[0][1] - 3,
                                   text=mylar_label, fill=MYLAR_DISPLAY_EDGE,
                                   anchor="sw", font=("TkDefaultFont", 8, "bold"))
        for edge, _component_index in edges:
            polygon = [self._local_3d_projection(point, center, scale, width, height) for point in edge]
            canvas.create_line([coord for point in polygon for coord in point[:2]], fill="#263238", width=2)
        for label, sensor_center in sensor_labels:
            point = self._local_3d_projection(sensor_center, center, scale, width, height)
            canvas.create_text(point[0], point[1], text=label, fill="#4a235a", font=("TkDefaultFont", 9, "bold"))
        for window in sensor_windows:
            projected_window = [self._local_3d_projection(point, center, scale, width, height) for point in window]
            canvas.create_polygon([coord for point in projected_window for coord in point[:2]],
                                  fill="#e1bee7", outline="#4a235a", width=3)
        # Keep a small physical XYZ triad visible in the refreshed scene.
        origin = self._local_3d_projection(center, center, scale, width, height)
        for label, color, delta in (("X", "#c62828", (extent * .55, 0, 0)), ("Y", "#2e7d32", (0, extent * .55, 0)), ("Z", "#1565c0", (0, 0, extent * .55))):
            endpoint = self._local_3d_projection(tuple(center[i] + delta[i] for i in range(3)), center, scale, width, height)
            canvas.create_line(origin[0], origin[1], endpoint[0], endpoint[1], fill=color, width=2, arrow=TK_ARROW_LAST)
            canvas.create_text(endpoint[0] + 5, endpoint[1], text=label, fill=color, anchor="w")
        if projected:
            names = [str(item.get("name", "")) for item in assembly_scene.get("components", ()) if item.get("kind") in {"geometry", "guide", "mirror", "sensor"}]
            canvas.create_text(10, 12, text="LOCAL 3D • " + " / ".join(names) +
                               " • dashed gold = external Mylar", anchor="nw", fill="#245a35")
            if mirror_labels:
                legend_x, legend_y = 10, 30
                legend_width = 172
                legend_height = 22 + 18 * len(mirror_labels)
                canvas.create_rectangle(legend_x, legend_y, legend_x + legend_width, legend_y + legend_height,
                                        fill="#ffffff", outline="#b0b0b0")
                canvas.create_text(legend_x + 7, legend_y + 4, text="Terminal mirrors", anchor="nw",
                                   fill="#333333", font=("TkDefaultFont", 9, "bold"))
                for row, (label, _mirror_center) in enumerate(mirror_labels):
                    y = legend_y + 20 + row * 18
                    canvas.create_rectangle(legend_x + 7, y, legend_x + 18, y + 11,
                                            fill="#d62728", outline="#7f1717")
                    canvas.create_text(legend_x + 24, y + 5, text=label, anchor="w",
                                       fill="#7f1717", font=("TkDefaultFont", 9, "bold"))
            canvas.create_text(10, height - 10, text="Drag: rotate • Wheel: zoom • Reset view",
                               anchor="sw", fill="#555555")

    def _reset_local_view(self) -> None:
        self.local_viewer_yaw, self.local_viewer_pitch, self.local_viewer_zoom = .75, .45, 1.0
        self.local_viewer_drag = None
        self._draw_local_3d_view()

    def _local_viewer_press(self, event) -> None:
        self.local_viewer_drag = (event.x, event.y, self.local_viewer_yaw, self.local_viewer_pitch)
        event.widget.focus_set()

    def _local_viewer_drag(self, event) -> None:
        if self.local_viewer_drag is None:
            return
        x0, y0, yaw, pitch = self.local_viewer_drag
        self.local_viewer_yaw = yaw + (event.x - x0) * .012
        self.local_viewer_pitch = max(-1.45, min(1.45, pitch + (event.y - y0) * .012))
        self._draw_local_3d_view()

    def _local_viewer_zoom(self, event) -> None:
        delta = getattr(event, "delta", 0) or (120 if getattr(event, "num", 4) == 4 else -120)
        self.local_viewer_zoom = max(.2, min(8.0, self.local_viewer_zoom * (1.12 if delta > 0 else 1 / 1.12)))
        self._draw_local_3d_view()

    def _close_local_3d_viewer(self) -> None:
        """Close the Scintillator viewer and clear all live references."""

        viewer = _get_instance_attr(self, "local_viewer", None)
        self.local_viewer = None
        self.local_viewer_canvas = None
        self.local_viewer_scene = None
        self.local_viewer_drag = None
        try:
            if viewer is not None and viewer.winfo_exists():
                viewer.destroy()
        except (AttributeError, tk.TclError):
            pass

    def _open_local_3d_viewer(self) -> None:
        if self.local_viewer is not None and self.local_viewer.winfo_exists():
            self.local_viewer.deiconify()
            self.local_viewer.lift()
            self._draw_local_3d_view()
            return
        viewer = tk.Toplevel(self)
        viewer.title("Local Scintillator / Light-guide Viewer")
        viewer.geometry("900x650")
        viewer.rowconfigure(0, weight=1)
        viewer.columnconfigure(0, weight=1)
        canvas = tk.Canvas(viewer, background="white", highlightthickness=1,
                           highlightbackground="#909090", takefocus=True)
        canvas.grid(row=0, column=0, sticky="nsew", padx=6, pady=6)
        controls = ttk.Frame(viewer)
        controls.grid(row=1, column=0, sticky="ew", padx=6, pady=(0, 6))
        ttk.Button(controls, text="Reset view", command=self._reset_local_view).pack(side="left")
        ttk.Label(controls, text="Local geometry only; quartz/gap are not shown",
                  foreground="#555555").pack(side="left", padx=12)
        ttk.Button(controls, text="Close", command=self._close_local_3d_viewer).pack(side="right")
        self.local_viewer, self.local_viewer_canvas = viewer, canvas
        # Keep a live callback on the window object as well as on the GUI.
        # This makes the refresh contract explicit for integrations that hold
        # the viewer reference directly.
        viewer.refresh = self._refresh_scintillator_viewer
        self.local_viewer_yaw, self.local_viewer_pitch, self.local_viewer_zoom, self.local_viewer_drag = .75, .45, 1.0, None
        bindings = local_geometry_view_binding_sequences()
        self.local_viewer_bindings = list(bindings)
        canvas.bind("<Configure>", lambda _event: self._draw_local_3d_view())
        canvas.bind("<ButtonPress-1>", self._local_viewer_press)
        canvas.bind("<B1-Motion>", self._local_viewer_drag)
        canvas.bind("<ButtonRelease-1>", lambda _event: setattr(self, "local_viewer_drag", None))
        canvas.bind("<MouseWheel>", self._local_viewer_zoom)
        canvas.bind("<Button-4>", self._local_viewer_zoom)
        canvas.bind("<Button-5>", self._local_viewer_zoom)
        viewer.protocol("WM_DELETE_WINDOW", self._close_local_3d_viewer)
        self._draw_local_3d_view()

    def _scintillator_local_value_changed(self) -> None:
        self._refresh_scintillator_preview()
        self._refresh_scintillator_viewer()
        self._refresh_global_assembly_preview()

    def _save_scintillator_profile(self) -> None:
        path = filedialog.asksaveasfilename(title="Save scintillator assembly", defaultextension=".json", filetypes=(("Assembly JSON", "*.json"), ("All files", "*.*")))
        if path:
            try:
                save_scintillator_profile(path, self.values)
                self.status.set(f"Scintillator assembly saved: {Path(path).name}")
            except (OSError, TypeError, ValueError) as exc:
                messagebox.showerror("Save scintillator assembly", str(exc), parent=self)

    def _load_scintillator_profile(self) -> None:
        path = filedialog.askopenfilename(title="Load scintillator assembly", filetypes=(("Assembly JSON", "*.json"), ("All files", "*.*")))
        if path:
            try:
                load_scintillator_profile(path, self.values)
                self._refresh_scintillator_preview()
                self.status.set(f"Scintillator assembly loaded: {Path(path).name}")
            except (OSError, json.JSONDecodeError, TypeError, ValueError) as exc:
                messagebox.showerror("Load scintillator assembly", str(exc), parent=self)

    def _build_scintillator_detector_tab(self, parent: ttk.Frame) -> None:
        """Edit detector response, generation, timing, and Mylar once."""

        parent.columnconfigure(0, weight=1)
        parent.columnconfigure(1, weight=1)
        response = ttk.LabelFrame(parent, text="Scintillator detector characteristics", padding=8)
        response.grid(row=0, column=0, sticky="new", padx=(0, 8))
        response.columnconfigure(1, weight=1)
        fields = [
            ("Sampling mode", "scintillator_sampling_mode", "mode"),
            ("Fixed sample photons", "scintillator_photon_count", "entry"),
            ("Yield (photons / deposited MeV)", "scintillator_yield_per_MeV", "entry"),
            ("Deposited energy override (MeV)", "scintillator_deposited_energy_MeV", "entry"),
            ("Constant dE/dx (MeV/m)", "scintillator_dEdx_MeV_per_m", "entry"),
            ("Decay time (ns)", "scintillator_decay_time_ns", "entry"),
            ("Emission wavelength min (nm)", "scintillator_wavelength_min_nm", "entry"),
            ("Emission wavelength max (nm)", "scintillator_wavelength_max_nm", "entry"),
            ("Readout index", "scintillator_readout_index", "entry"),
            ("Readout acceptance (deg)", "scintillator_readout_acceptance_deg", "entry"),
            ("Readout QE", "scintillator_readout_qe", "entry"),
            ("Readout collection efficiency", "scintillator_readout_collection_efficiency", "entry"),
            ("Absorption length (mm, blank none)", "scintillator_absorption_length_mm", "entry"),
            ("Legacy maximum bounces (Run control overrides)", "scintillator_max_bounces", "legacy_entry"),
            ("Mylar enabled", "mylar_enabled", "bool"),
            ("Mylar reflectivity (%) [80–100]", "mylar_reflectivity_percent", "entry"),
        ]
        for row, (label, key, kind) in enumerate(fields):
            ttk.Label(response, text=label).grid(row=row, column=0, sticky="w", padx=(0, 8), pady=2)
            if kind == "bool":
                widget = ttk.Checkbutton(response, variable=self.values[key], command=self._refresh_detector_characteristics_status)
            elif kind == "mode":
                widget = ttk.Combobox(response, textvariable=self.values[key], values=tuple(SCINTILLATION_SAMPLING_MODE_LABELS.values()), state="readonly")
            else:
                widget = ttk.Entry(response, textvariable=self.values[key], state="readonly" if kind == "legacy_entry" else "normal")
            widget.grid(row=row, column=1, sticky="ew", pady=2)
            if kind != "bool":
                self.values[key].trace_add("write", lambda *_args: self._detector_characteristics_value_changed())
        buttons = ttk.Frame(response)
        buttons.grid(row=len(fields), column=0, columnspan=2, sticky="w", pady=(8, 0))
        ttk.Button(buttons, text="Load Detector Characteristics JSON", command=self._load_detector_characteristics_profile).pack(side="left")
        ttk.Button(buttons, text="Save Detector Characteristics JSON", command=self._save_detector_characteristics_profile).pack(side="left", padx=(5, 0))
        self.scintillator_detector_status = tk.StringVar(value="Not validated")
        ttk.Label(response, textvariable=self.scintillator_detector_status, foreground="#245a35", wraplength=600, justify="left").grid(row=len(fields) + 1, column=0, columnspan=2, sticky="w", pady=(6, 0))
        ttk.Label(parent, text="These values are the single detector-characteristics source used by the run and by the detector JSON profile. Mylar applies to external assembly surfaces only; the scintillator-guide interface and readout face remain unwrapped.", foreground="#245a35", wraplength=420, justify="left").grid(row=0, column=1, sticky="nw")
        self._refresh_detector_characteristics_status()

    def _detector_characteristics_value_changed(self) -> None:
        self._refresh_detector_characteristics_status()
        self._refresh_global_assembly_preview()

    def _refresh_detector_characteristics_status(self, *_args) -> None:
        if not _has_instance_attr(self, "scintillator_detector_status"):
            return
        try:
            config = scintillator_config_from_values(self.values)
            reflectivity = float(_value_text(self.values, "mylar_reflectivity_percent", "95"))
            if not 80.0 <= reflectivity <= 100.0:
                raise ValueError("Mylar reflectivity must be between 80 and 100 percent")
            self.scintillator_detector_status.set(f"Valid detector characteristics: {config.sampling_mode} mode, yield {config.photon_yield_per_MeV:g}/MeV, Mylar {reflectivity:g}%")
        except (TypeError, ValueError) as exc:
            self.scintillator_detector_status.set(f"Invalid detector characteristics: {exc}")

    def _save_local_geometry_profile(self) -> None:
        path = filedialog.asksaveasfilename(title="Save local scintillator geometry", defaultextension=".json", filetypes=(("Local geometry JSON", "*.json"), ("All files", "*.*")))
        if path:
            try:
                save_component_profile(path, local_geometry_profile_from_values(self.values))
                self.status.set(f"Local geometry saved: {Path(path).name}")
            except (OSError, TypeError, ValueError) as exc:
                messagebox.showerror("Save local geometry", str(exc), parent=self)

    def _load_local_geometry_profile(self) -> None:
        path = filedialog.askopenfilename(title="Load local scintillator geometry", filetypes=(("Local geometry JSON", "*.json"), ("All files", "*.*")))
        if path:
            try:
                load_component_profile(path, self.values, apply_local_geometry_profile)
                self._refresh_scintillator_preview()
                self._refresh_global_assembly_preview()
                self.status.set(f"Local geometry loaded: {Path(path).name}")
            except (OSError, json.JSONDecodeError, TypeError, ValueError) as exc:
                messagebox.showerror("Load local geometry", str(exc), parent=self)

    def _save_detector_characteristics_profile(self) -> None:
        path = filedialog.asksaveasfilename(title="Save scintillator detector characteristics", defaultextension=".json", filetypes=(("Detector JSON", "*.json"), ("All files", "*.*")))
        if path:
            try:
                save_component_profile(path, detector_characteristics_profile_from_values(self.values))
                self.status.set(f"Detector characteristics saved: {Path(path).name}")
            except (OSError, TypeError, ValueError) as exc:
                messagebox.showerror("Save detector characteristics", str(exc), parent=self)

    def _load_detector_characteristics_profile(self) -> None:
        path = filedialog.askopenfilename(title="Load scintillator detector characteristics", filetypes=(("Detector JSON", "*.json"), ("All files", "*.*")))
        if path:
            try:
                load_component_profile(path, self.values, apply_detector_characteristics_profile)
                self._refresh_detector_characteristics_status()
                self._refresh_scintillator_preview()
                self._refresh_global_assembly_preview()
                self.status.set(f"Detector characteristics loaded: {Path(path).name}")
            except (OSError, json.JSONDecodeError, TypeError, ValueError) as exc:
                messagebox.showerror("Load detector characteristics", str(exc), parent=self)

    def _build_global_geometry_tab(self, parent: ttk.Frame) -> None:
        """Edit only placement/gap and show the global assembly viewers."""

        parent.columnconfigure(0, weight=1)
        parent.columnconfigure(1, weight=2)
        controls = ttk.LabelFrame(parent, text="Canonical Cherenkov + Scintillator assemblies • GLOBAL PLACEMENT", padding=8)
        controls.grid(row=0, column=0, sticky="new", padx=(0, 8))
        controls.columnconfigure(1, weight=1)
        ttk.Label(
            controls,
            text=f"Loaded GUI: {APP_VERSION} / {APP_BUILD_ID} • Cherenkov assembly → Scintillator assembly",
            foreground="#b00020",
        ).grid(row=0, column=0, columnspan=2, sticky="w", pady=(0, 4))
        ttk.Label(
            controls,
            text="Signed offsets accept negative and positive values; only gap/separation are constrained non-negative.",
            foreground="#245a35", wraplength=360, justify="left",
        ).grid(row=1, column=0, columnspan=2, sticky="w", pady=(0, 6))
        fields = [
            ("Quartz → Detector 1 gap (mm)", "scintillator_gap_mm"),
            ("Detector 1 active-plane Z (mm, nominal)", "detector_1_z_mm"),
            ("Detector 1 signed offset X (mm)", "global_offset_x_mm"),
            ("Detector 1 signed offset Y (mm)", "global_offset_y_mm"),
            ("Detector 1 signed offset Z (mm)", "global_offset_z_mm"),
            ("Detector 1 → Detector 2 separation (mm)", "detector_separation_mm"),
            ("Detector 2 active-plane Z (mm, nominal)", "detector_2_z_mm"),
            ("Detector 2 signed offset X (mm)", "detector2_offset_x_mm"),
            ("Detector 2 signed offset Y (mm)", "detector2_offset_y_mm"),
            ("Detector 2 signed offset Z (mm)", "detector2_offset_z_mm"),
        ]
        for row, (label, key) in enumerate(fields, start=2):
            ttk.Label(controls, text=label).grid(row=row, column=0, sticky="w", padx=(0, 8), pady=3)
            ttk.Entry(controls, textvariable=self.values[key]).grid(row=row, column=1, sticky="ew", pady=3)
            self.values[key].trace_add("write", lambda *_args: self._refresh_global_assembly_preview())
        next_row = 2 + len(fields)
        ttk.Label(controls, text="Beam/assembly axis",).grid(row=next_row, column=0, sticky="w", padx=(0, 8), pady=3)
        ttk.Combobox(controls, textvariable=self.values["global_beam_axis"], values=GEOMETRY_AXES, state="readonly").grid(row=next_row, column=1, sticky="ew", pady=3)
        self.values["global_beam_axis"].trace_add("write", lambda *_args: self._refresh_global_assembly_preview())
        guide_row = next_row + 1
        ttk.Label(controls, text="Detector 1 Y-/Y+ guides (optional)",).grid(row=guide_row, column=0, sticky="w", padx=(0, 8), pady=3)
        guide_frame = ttk.Frame(controls)
        guide_frame.grid(row=guide_row, column=1, sticky="ew", pady=3)
        ttk.Combobox(guide_frame, textvariable=self.values["detector1_guide_left_shape"], values=tuple(GEOMETRY_SHAPE_LABELS.keys()), state="readonly", width=18).pack(side="left")
        ttk.Combobox(guide_frame, textvariable=self.values["detector1_guide_right_shape"], values=tuple(GEOMETRY_SHAPE_LABELS.keys()), state="readonly", width=18).pack(side="left", padx=(4, 0))
        ttk.Checkbutton(guide_frame, text="Y-", variable=self.values["detector1_guide_left_enabled"], onvalue="1", offvalue="0", command=self._refresh_global_assembly_preview).pack(side="left", padx=(6, 0))
        ttk.Checkbutton(guide_frame, text="Y+", variable=self.values["detector1_guide_right_enabled"], onvalue="1", offvalue="0", command=self._refresh_global_assembly_preview).pack(side="left", padx=(3, 0))
        detector2_row = guide_row + 1
        ttk.Label(controls, text="Detector 2 scintillator shape",).grid(row=detector2_row, column=0, sticky="w", padx=(0, 8), pady=3)
        ttk.Combobox(controls, textvariable=self.values["detector2_scintillator_shape"], values=tuple(GEOMETRY_SHAPE_LABELS.keys()), state="readonly").grid(row=detector2_row, column=1, sticky="ew", pady=3)
        detector2_guide_row = detector2_row + 1
        ttk.Label(controls, text="Detector 2 Y-/Y+ guides (optional)",).grid(row=detector2_guide_row, column=0, sticky="w", padx=(0, 8), pady=3)
        guide2_frame = ttk.Frame(controls)
        guide2_frame.grid(row=detector2_guide_row, column=1, sticky="ew", pady=3)
        ttk.Combobox(guide2_frame, textvariable=self.values["detector2_guide_left_shape"], values=tuple(GEOMETRY_SHAPE_LABELS.keys()), state="readonly", width=18).pack(side="left")
        ttk.Combobox(guide2_frame, textvariable=self.values["detector2_guide_right_shape"], values=tuple(GEOMETRY_SHAPE_LABELS.keys()), state="readonly", width=18).pack(side="left", padx=(4, 0))
        ttk.Checkbutton(guide2_frame, text="Y-", variable=self.values["detector2_guide_left_enabled"], onvalue="1", offvalue="0", command=self._refresh_global_assembly_preview).pack(side="left", padx=(6, 0))
        ttk.Checkbutton(guide2_frame, text="Y+", variable=self.values["detector2_guide_right_enabled"], onvalue="1", offvalue="0", command=self._refresh_global_assembly_preview).pack(side="left", padx=(3, 0))
        for key in ("detector1_guide_left_shape", "detector1_guide_right_shape", "detector2_scintillator_shape", "detector2_guide_left_shape", "detector2_guide_right_shape", "detector_1_z_mm", "detector_2_z_mm"):
            self.values[key].trace_add("write", lambda *_args: self._refresh_global_assembly_preview())
        profile_row = ttk.Frame(controls)
        profile_row.grid(row=detector2_guide_row + 1, column=0, columnspan=2, sticky="w", pady=(8, 0))
        ttk.Button(profile_row, text="Load placement-only Global Geometry JSON", command=self._load_global_geometry_profile).pack(side="left")
        ttk.Button(profile_row, text="Save placement-only Global Geometry JSON", command=self._save_global_geometry_profile).pack(side="left", padx=(5, 0))
        ttk.Button(profile_row, text="Open Assembly 3D Viewer", command=self._open_global_assembly_3d_viewer).pack(side="left", padx=(12, 0))
        complete_profile_row = ttk.Frame(controls)
        complete_profile_row.grid(row=detector2_guide_row + 2, column=0, columnspan=2, sticky="w", pady=(5, 0))
        ttk.Button(complete_profile_row, text="Save complete assembly profile", command=self._save_complete_global_assembly_profile).pack(side="left")
        ttk.Button(complete_profile_row, text="Load complete assembly profile", command=self._load_complete_global_assembly_profile).pack(side="left", padx=(5, 0))
        ttk.Label(
            controls,
            text="Complete profile = both local canonical assemblies + ordered guides, terminals/Mylar, sensors/readouts, impacts, transforms, and beam gaps. The existing Global Geometry JSON above restores placement values only; Run particle settings and Maximum photon bounces stay outside this profile.",
            foreground="#245a35", wraplength=430, justify="left",
        ).grid(row=detector2_guide_row + 3, column=0, columnspan=2, sticky="w", pady=(4, 0))
        self.global_geometry_status = tk.StringVar(value="Not validated")
        ttk.Label(controls, textvariable=self.global_geometry_status, foreground="#245a35", wraplength=360, justify="left").grid(row=detector2_guide_row + 4, column=0, columnspan=2, sticky="w", pady=(8, 0))
        preview = ttk.LabelFrame(parent, text="Canonical complete assemblies • orthographic XZ / XY / YZ", padding=6)
        preview.grid(row=0, column=1, sticky="nsew")
        preview.columnconfigure(0, weight=1)
        preview.rowconfigure(1, weight=1)
        self.global_preview_text = tk.StringVar(value="Not validated")
        ttk.Label(preview, textvariable=self.global_preview_text, wraplength=600, justify="left").grid(row=0, column=0, sticky="nw")
        self.global_preview_canvas = tk.Canvas(preview, width=700, height=430, background="white", highlightthickness=1, highlightbackground="#b0b0b0")
        self.global_preview_canvas.grid(row=1, column=0, sticky="nsew")
        self.global_preview_canvas.bind("<Configure>", self._refresh_global_assembly_preview)
        self._refresh_global_assembly_preview()

    def _draw_global_orthographic_preview(self, canvas, scene, width, height):
        """Draw the complete global scene in three physical-coordinate panes."""

        layout = global_assembly_orthographic_layout(scene, width, height)
        canvas.delete("all")
        canvas.create_rectangle(0, 0, layout["width"], layout["height"], fill="white", outline="")
        axis_colors = {"X": "#c62828", "Y": "#6a1b9a", "Z": "#1565c0"}
        component_list = tuple(scene.get("render_components", ()))
        for view_name in layout["view_order"]:
            pane = layout["views"][view_name]
            left, top_edge, right, bottom_edge = pane["rect"]
            canvas.create_rectangle(left, top_edge, right, bottom_edge, fill="#fbfdff", outline="#90a4ae", width=1)
            canvas.create_text(left + 7, top_edge + 6, text=pane["title"], anchor="nw", fill="#17324d", font=("TkDefaultFont", 9, "bold"))
            project = pane["project"]
            axis_a, axis_b = pane["axes"]
            for component in component_list:
                faces = list(component.get("faces_mm", ()))
                is_mylar = component.get("kind") == "mylar_overlay"
                for face in faces:
                    points = [project(point) for point in face]
                    if len(points) < 2:
                        continue
                    canvas.create_polygon([coord for point in points for coord in point],
                                          fill=MYLAR_DISPLAY_COLOR if is_mylar else component.get("color", "#80b1d3"),
                                          outline=MYLAR_DISPLAY_EDGE if is_mylar else "#263238",
                                          stipple="gray25" if is_mylar else "gray50")
                    line = {"fill": MYLAR_DISPLAY_EDGE if is_mylar else "#263238", "width": 1 if is_mylar else 2}
                    if is_mylar:
                        line["dash"] = (5, 3)
                    canvas.create_line([coord for point in points + [points[0]] for coord in point], **line)
                for key, color, dash, fill in (
                    ("interface_faces_mm", "#00695c", (2, 2), None),
                    ("readout_faces_mm", "#4a235a", None, "#9b59b6"),
                ):
                    for face in component.get(key, ()):
                        points = [project(point) for point in face]
                        if fill:
                            canvas.create_polygon([coord for point in points for coord in point], fill=fill, outline=color, stipple="gray50")
                        options = {"fill": color, "width": 2}
                        if dash:
                            options["dash"] = dash
                        canvas.create_line([coord for point in points + [points[0]] for coord in point], **options)
            path = [project(point) for point in layout["path_mm"]]
            for first, second in zip(path, path[1:]):
                canvas.create_line(first[0], first[1], second[0], second[1], fill=GLOBAL_PARTICLE_PATH_COLOR, dash=(7, 3), width=2)
            for point in path:
                canvas.create_oval(point[0] - 2.5, point[1] - 2.5, point[0] + 2.5, point[1] + 2.5, fill=GLOBAL_PARTICLE_PATH_COLOR, outline="white")
            for marker in layout["markers_mm"]:
                point = project(marker)
                canvas.create_oval(point[0] - 4, point[1] - 4, point[0] + 4, point[1] + 4, fill="#d32f2f", outline="white", width=1)
            origin = project(layout["origin_mm"])
            canvas.create_oval(origin[0] - 3, origin[1] - 3, origin[0] + 3, origin[1] + 3, fill="#263238", outline="white")
            canvas.create_text(origin[0] + 5, origin[1] - 5, text="O", anchor="sw", fill="#263238", font=("TkDefaultFont", 8, "bold"))
            # Physical reference arrows are drawn from projected global origin;
            # the two labels are outside the component silhouettes.
            for axis_index, axis in ((axis_a, pane["axis_labels"][0][0]), (axis_b, pane["axis_labels"][1][0])):
                vector = [0.0, 0.0, 0.0]
                vector[axis_index] = max(pane["bounds_mm"][2] - pane["bounds_mm"][0], pane["bounds_mm"][3] - pane["bounds_mm"][1]) * 0.16
                endpoint = project(vector)
                start = project((0.0, 0.0, 0.0))
                canvas.create_line(start[0], start[1], endpoint[0], endpoint[1], fill=axis_colors[axis], width=2, arrow=TK_ARROW_LAST)
                canvas.create_text(endpoint[0] + 4, endpoint[1], text=pane["axis_labels"][0] if axis_index == axis_a else pane["axis_labels"][1], anchor="w", fill=axis_colors[axis], font=("TkDefaultFont", 7, "bold"))
            canvas.create_text((left + right) / 2.0, bottom_edge - 7, text=pane["axis_labels"][0] + "  →   " + pane["axis_labels"][1] + " ↑", anchor="s", fill="#455a64", font=("TkDefaultFont", 7))
        legend_left, legend_top, legend_right, legend_bottom = layout["legend_rect"]
        canvas.create_text(legend_left, legend_top, anchor="nw", text="Legend: filled outlines = detector / guide surfaces   purple = readout face   teal dashed = interface   gold dashed = Mylar   green dashed = one shared particle path   red = configured impact", fill="#37474f", font=("TkDefaultFont", 7))
        canvas.create_text(legend_left, legend_top + 16, anchor="nw", text=f"Global +Z beam convention • inclination {scene.get('inclination_deg', 0.0):g}° from +Z, positive toward +X • origin O shown in every projection • {layout['responsive']}", fill="#0b6b43", font=("TkDefaultFont", 7, "bold"))
        return layout


    def _refresh_global_assembly_preview(self, *_args) -> None:
        canvas = _get_instance_attr(self, "global_preview_canvas", None)
        if canvas is None:
            # The interactive viewer is a live view of the same canonical
            # scene.  Refresh it in place even when the global tab has not
            # been visited or its preview canvas is not available.
            viewer = _get_instance_attr(self, "global_viewer", None)
            if viewer is not None and _get_instance_attr(self, "global_viewer_canvas", None) is not None:
                self._draw_global_assembly_3d_view()
            return
        canvas.delete("all")
        try:
            scene = global_assembly_scene(self.values, self._geometry_for_simulation(), self.assembly_models)
            width, height = max(int(canvas.winfo_width()), 420), max(int(canvas.winfo_height()), 260)
            layout = global_assembly_orthographic_layout(scene, width, height)
            width, height = layout["width"], layout["height"]
            self._draw_global_orthographic_preview(canvas, scene, width, height)
            alignments = scene["beam_alignment"]
            self.global_preview_text.set(
                f"Three orthographic views • shared global direction {scene['particle_direction_global']} • "
                f"configured-marker residuals C={alignments['cherenkov']['distance_mm']:.3g} mm, "
                f"S={alignments['scintillator']['distance_mm']:.3g} mm"
            )
            self.global_geometry_status.set(scene["status_summary"] + "\n" + scene["particle_axis_diagnostic"])
            details_lines = [
                "THREE ORTHOGRAPHIC VIEWS",
                "XZ = longitudinal / beam; XY = transverse; YZ = terminal/transverse.",
                "All panes use the same transformed global component vertices.",
                "Green dashed = one shared particle path; red = saved impact marker.",
                f"Direction = {scene['particle_direction_global']}; anchor = {scene.get('beam_anchor_global_mm')}",
                f"Cherenkov line/active check: {alignments['cherenkov']['line_intersects_active_volume']}",
                f"Scintillator line/active check: {alignments['scintillator']['line_intersects_active_volume']}",
                "",
                "DISTANCE SUMMARY",
                f"C +X terminal → S −X entrance/reference: {scene['cher_d1_gap_mm']:g} mm",
                f"D1(C) +X face → D2(S) +X face: {scene['d1_d2_separation_mm']:g} mm",
                f"S−C translation: {scene['relative_global_assembly_offset_mm']} mm",
                "Gap and assembly offset are independent.",
            ]
            self._set_viewer_details_text(_get_instance_attr(self, "global_preview_info_text", None), "\n".join(details_lines))
            # Keep the detached interactive viewer as an optional legacy tool;
            # its scene refresh still uses this same canonical source.
            self._draw_global_assembly_3d_view()
            return
            project = layout["project"]
            # The embedded canvas is geometry-only.  Details, measurements,
            # and the component legend live in the docked scrollable pane so
            # no annotation can hide a sensitive volume or readout.
            for component, placed in zip(scene["render_components"], layout["components"]):
                faces = placed["faces"]
                for face_index, (face, points) in enumerate(zip(faces, placed["faces"])):
                    is_mylar = component.get("kind") == "mylar_overlay"
                    canvas.create_polygon(
                        [coord for point in points for coord in point],
                        fill=MYLAR_DISPLAY_COLOR if is_mylar else component["color"],
                        outline=MYLAR_DISPLAY_EDGE if is_mylar else "#263238",
                        stipple="gray25",
                    )
                    line_options = {"fill": MYLAR_DISPLAY_EDGE if is_mylar else "#263238", "width": 2}
                    if is_mylar:
                        line_options["dash"] = (6, 3)
                    canvas.create_line(
                        [coord for point in points + [points[0]] for coord in point], **line_options,
                    )
                    if face_index in set(component.get("mylar_indices", ())):
                        canvas.create_line(
                            [coord for point in points + [points[0]] for coord in point],
                            fill="#b26a00", dash=(6, 3), width=2,
                        )
                for interface in component.get("interface_faces_mm", ()):
                    points = [project(_global_assembly_render_point(point, scene)) for point in interface]
                    canvas.create_line([coord for point in points + [points[0]] for coord in point],
                                       fill="#00695c", dash=(2, 2), width=3)
                for readout in component.get("readout_faces_mm", ()):
                    points = [project(_global_assembly_render_point(point, scene)) for point in readout]
                    canvas.create_polygon([coord for point in points for coord in point],
                                          fill="#9b59b6", outline="#4a235a", stipple="gray25")
                    canvas.create_line([coord for point in points + [points[0]] for coord in point],
                                       fill="#4a235a", width=3)

            # One physical shared path is drawn from the canonical scene.  At
            # the default no-rotation placement this is exactly global +Z;
            # offsets do not move either assembly to make the line fit.
            path_points = [project(_global_assembly_render_point(point, scene)) for point in scene.get("display_particle_path_mm", scene["particle_path_mm"])]
            for first, second in zip(path_points, path_points[1:]):
                canvas.create_line(first[0], first[1], second[0], second[1],
                                   fill=GLOBAL_PARTICLE_PATH_COLOR, dash=(7, 3), width=3)
            for point in path_points:
                canvas.create_oval(point[0] - 3, point[1] - 3, point[0] + 3, point[1] + 3,
                                   fill=GLOBAL_PARTICLE_PATH_COLOR, outline="#ffffff")
            for marker_index, marker in enumerate(scene.get("transformed_impact_markers_mm", ())):
                marker_point = project(_global_assembly_render_point(marker, scene))
                canvas.create_oval(marker_point[0] - 5, marker_point[1] - 5, marker_point[0] + 5, marker_point[1] + 5, fill="#d32f2f", outline="#ffffff", width=2)

            # Title and XYZ triad occupy reserved strips outside the physical
            # drawing bounds.  The triad directions are the same projected
            # vectors used by the geometry, so its Z arrow cannot contradict
            # the green path.
            canvas.create_rectangle(0, 0, width, 42, fill="#f5f8fb", outline="")
            canvas.create_text(10, 7, anchor="nw", text=f"GLOBAL ASSEMBLY • {scene['display_mode']}",
                               fill="#17324d", font=("TkDefaultFont", 11, "bold"))
            canvas.create_text(width - 10, 9, anchor="ne", text="GLOBAL +Z • charged-particle direction",
                               fill="#0b6b43", font=("TkDefaultFont", 9, "bold"))
            canvas.create_text(width - 10, 27, anchor="ne", text="C Cherenkov → S scintillator • true geometry",
                               fill="#8a5200", font=("TkDefaultFont", 8))
            triad_origin = (38.0, height - 17.0)
            axis_colors = {"X": "#c62828", "Y": "#6a1b9a", "Z": "#1565c0"}
            for axis in ("X", "Y", "Z"):
                vector = layout["axis_vectors_projected"][axis]
                vector_length = max(math.hypot(vector[0], vector[1]), 1e-12)
                endpoint = (
                    triad_origin[0] + vector[0] / vector_length * 25.0,
                    triad_origin[1] - vector[1] / vector_length * 25.0,
                )
                canvas.create_line(triad_origin[0], triad_origin[1], endpoint[0], endpoint[1],
                                   fill=axis_colors[axis], width=2, arrow=TK_ARROW_LAST)
                canvas.create_text(endpoint[0] + 4, endpoint[1], text=f"+{axis}",
                                   fill=axis_colors[axis], anchor="w", font=("TkDefaultFont", 8, "bold"))
            canvas.create_text(82, height - 14, anchor="w", text="global XYZ • camera-coupled",
                               fill="#455a64", font=("TkDefaultFont", 8))
            canvas.create_text(width - 10, height - 12, anchor="e", text="green dashed = shared charged-particle path",
                               fill=GLOBAL_PARTICLE_PATH_COLOR, font=("TkDefaultFont", 8, "bold"))
            mode = scene["display_mode"]
            definitions = scene["measurement_definitions"]
            alignments = scene["beam_alignment"]
            self.global_preview_text.set(
                f"{scene['status_summary'].splitlines()[0]} • particle axis: global +{scene['beam_axis']} • "
                f"shared path residuals C={alignments['cherenkov']['distance_mm']:.3g} mm, "
                f"S={alignments['scintillator']['distance_mm']:.3g} mm • "
                f"S−C offset {scene['relative_global_assembly_offset_mm']} mm"
            )
            details_lines = [
                "PARTICLE / BEAM",
                f"Particle axis: global +{scene['beam_axis']}",
                "Green dashed line: one shared charged-particle path",
                f"Track endpoints: {scene['particle_path_mm'][0]} → {scene['particle_path_mm'][-1]} mm",
                "",
                "IMPACTS (local → global)",
            ]
            for assembly_type in ("cherenkov", "scintillator"):
                alignment = alignments[assembly_type]
                details_lines.append(
                    f"{assembly_type.title()}: {alignment['local_impact_point_mm']} → "
                    f"{alignment['impact_point_mm']} mm; residual {alignment['distance_mm']:.3g} mm"
                )
            details_lines.extend([
                "",
                "MEASUREMENTS (+X assembly axis)",
                f"C +X terminal → S −X entrance/reference: {definitions['cher_d1_gap']['signed_mm']:g} mm",
                f"D1(C) +X face → D2(S) +X face: {definitions['d1_d2_separation']['signed_mm']:g} mm",
                f"Relative translation S−C: {definitions['relative_global_assembly_offset']['vector_mm']} mm",
                "Gap and Z offset are independent quantities.",
                "",
                "LEGEND",
                "C quartz / guide; S scintillator / guide",
                "Purple: C/S readout bodies; teal: interfaces",
                "Gold dashed: Mylar; red markers: configured impacts",
                "Geometry dimensions and transforms are unchanged by projection.",
            ])
            self._set_viewer_details_text(_get_instance_attr(self, "global_preview_info_text", None), "\n".join(details_lines))
            if _has_instance_attr(self, "global_assembly_inherited_summary"):
                self.global_assembly_inherited_summary.set(
                    f"Cherenkov: {len(self._assembly_model('cherenkov').sections)} active section(s); "
                    f"guide {'on' if self._assembly_model('cherenkov').guide_enabled else 'off'}; "
                    f"terminal {self._assembly_model('cherenkov').terminal_plus_surface}.\n"
                    f"Scintillator: {len(self._assembly_model('scintillator').sections)} active section(s); "
                    f"guide {'on' if self._assembly_model('scintillator').guide_enabled else 'off'}.\n"
                    "Summary is inherited read-only; this tab controls rigid transforms only."
                )
            self.global_geometry_status.set(scene["status_summary"] + "\nParticle axis: global +" + scene["beam_axis"] + "; signed offsets accepted; current shared path shown")
            # Do not open a second Toplevel.  If it already exists, redraw
            # that same object from the current canonical scene immediately.
            self._draw_global_assembly_3d_view()
        except (TypeError, ValueError) as exc:
            self.global_preview_text.set(f"Invalid global assembly: {exc}")
            self.global_geometry_status.set(f"Invalid: {exc}")
            canvas.create_text(max(canvas.winfo_width(), 420) / 2, 180, text="Invalid global geometry", fill="#b00020")

    def _save_global_geometry_profile(self) -> None:
        path = filedialog.asksaveasfilename(title="Save placement-only Global Geometry", defaultextension=".json", filetypes=(("Placement-only Global geometry JSON", "*.json"), ("All files", "*.*")))
        if path:
            try:
                save_component_profile(path, global_geometry_profile_from_values(self.values))
                self.status.set(f"Global geometry saved: {Path(path).name}")
            except (OSError, TypeError, ValueError) as exc:
                messagebox.showerror("Save global geometry", str(exc), parent=self)

    def _load_global_geometry_profile(self) -> None:
        path = filedialog.askopenfilename(title="Load placement-only Global Geometry", filetypes=(("Placement-only Global geometry JSON", "*.json"), ("All files", "*.*")))
        if path:
            try:
                load_component_profile(path, self.values, apply_global_geometry_profile)
                self._refresh_global_assembly_preview()
                self.status.set(f"Global geometry loaded: {Path(path).name}")
            except (OSError, json.JSONDecodeError, TypeError, ValueError) as exc:
                messagebox.showerror("Load global geometry", str(exc), parent=self)

    def _save_complete_global_assembly_profile(self) -> None:
        path = filedialog.asksaveasfilename(
            title="Save complete Global Assembly profile", defaultextension=".json",
            filetypes=(("Complete Global Assembly profile", "*.json"), ("All files", "*.*")),
        )
        if not path:
            return
        try:
            save_complete_global_assembly_profile(path, self.values, self.assembly_models)
            self.status.set(f"Complete Global Assembly profile saved: {Path(path).name}")
        except (OSError, TypeError, ValueError) as exc:
            messagebox.showerror("Save complete assembly profile", str(exc), parent=self)

    def _load_complete_global_assembly_profile(self) -> None:
        path = filedialog.askopenfilename(
            title="Load complete Global Assembly profile",
            filetypes=(("Complete Global Assembly profile", "*.json"), ("All files", "*.*")),
        )
        if not path:
            return
        try:
            profile = json.loads(Path(path).expanduser().resolve().read_text(encoding="utf-8"))
            # The detached parse/validation is deliberately outside the Tk
            # suppression window.  Invalid files therefore cannot fire a
            # trace, replace a model, or alter any visible StringVar.
            validate_complete_global_assembly_profile(profile)
            self._suppress_assembly_refresh = True
            try:
                apply_complete_global_assembly_profile(self.values, self.assembly_models, profile)
                self.geometry_sections = self.assembly_models["cherenkov"].sections
                self.scintillator_geometry_sections = self.assembly_models["scintillator"].sections
            finally:
                self._suppress_assembly_refresh = False
            for assembly in ("cherenkov", "scintillator"):
                self._assembly_refresh_tree(assembly)
                self._assembly_build_properties(assembly)
            self._sensor_status_changed("cherenkov")
            self._sensor_status_changed("scintillator")
            self._refresh_loaded_assembly_views("cherenkov")
            self._refresh_loaded_assembly_views("scintillator")
            self._refresh_global_assembly_preview()
            self.status.set(f"Complete Global Assembly profile loaded: {Path(path).name}")
        except (OSError, json.JSONDecodeError, TypeError, ValueError) as exc:
            # Keep this error path free of status/StringVar changes; the
            # profile API has already guaranteed no partial state mutation.
            messagebox.showerror("Load complete assembly profile", str(exc), parent=self)

    def _global_assembly_viewer_projection(self, point, center, scale, width, height):
        yaw = _get_instance_attr(self, "global_viewer_yaw", GLOBAL_PREVIEW_DEFAULT_YAW)
        pitch = _get_instance_attr(self, "global_viewer_pitch", GLOBAL_PREVIEW_DEFAULT_PITCH)
        return global_assembly_camera_project(point, center, scale, width, height, yaw, pitch)

    @staticmethod
    def _set_viewer_details_text(widget, text: str) -> None:
        """Refresh a read-only, scrollable viewer details pane."""

        if widget is None:
            return
        try:
            widget.configure(state="normal")
            widget.delete("1.0", "end")
            widget.insert("1.0", text)
            widget.configure(state="disabled")
        except tk.TclError:
            # A viewer may be closing while a canvas Configure event fires.
            pass

    def _global_viewer_render_error(self, canvas, error: Exception) -> None:
        """Leave a visible diagnostic when a Tk callback cannot render."""

        message = f"Global assembly render error: {type(error).__name__}: {error}"
        self.global_viewer_render_error = message
        try:
            self.status.set(message)
        except Exception:
            pass
        try:
            self._append_log("ERROR: " + message)
        except Exception:
            pass
        try:
            canvas.delete("all")
            width = max(int(canvas.winfo_width()), 500)
            height = max(int(canvas.winfo_height()), 350)
            canvas.create_rectangle(0, 0, width, height, fill="white", outline="")
            canvas.create_text(
                width / 2, height / 2,
                text=message + "\nSee the main-window status/log for details.",
                fill="#b00020", width=max(width - 48, 240), justify="center",
            )
        except Exception:
            # A window can be destroyed between the callback and this
            # diagnostic.  The status/log attempt above is still useful.
            pass

    def _draw_global_assembly_3d_view(self, width=None, height=None, *, viewer=None, canvas=None) -> None:
        """Render the current global scene, with no silent Tk callback failure."""

        canvas = canvas or _get_instance_attr(self, "global_viewer_canvas", None)
        current_viewer = _get_instance_attr(self, "global_viewer", None)
        if viewer is not None and (viewer is not current_viewer or canvas is not _get_instance_attr(self, "global_viewer_canvas", None)):
            return
        try:
            if canvas is None or not canvas.winfo_exists():
                return
            canvas.delete("all")
            scene = global_assembly_scene(self.values, self._geometry_for_simulation(), self.assembly_models)
            if width is None:
                width = canvas.winfo_width()
            if height is None:
                height = canvas.winfo_height()
            width, height = max(int(width), 500), max(int(height), 350)
            self._draw_global_assembly_3d_view_impl(canvas, scene, width, height)
            self.global_viewer_render_error = ""
        except Exception as exc:
            self._global_viewer_render_error(canvas, exc)

    def _draw_global_assembly_3d_view_impl(self, canvas, scene, width, height) -> None:
        """Canvas primitive implementation kept behind the safe callback wrapper."""

        fit = global_assembly_fit(scene, width, height)
        center = fit["center_mm"]
        scale = fit["scale_px_per_mm"] * _get_instance_attr(self, "global_viewer_zoom", 1.0)
        label_font = tkfont.Font(root=canvas, family="TkDefaultFont", size=9, weight="bold") if tkfont else ("TkDefaultFont", 9, "bold")
        label_measure = (lambda text: (label_font.measure(text) + 8, label_font.metrics("linespace") + 2)) if tkfont else None
        callout_items = []
        for component in scene["render_components"]:
            faces, color, label = component["faces_mm"], component["color"], component["label"]
            projected = []
            for face_index, face in enumerate(faces):
                display_face = [_global_assembly_render_point(point, scene) for point in face]
                points = [self._global_assembly_viewer_projection(point, center, scale, width, height) for point in display_face]
                projected.append((sum(point[2] for point in points) / len(points), face_index, points))
            for _depth, original_face_index, points in sorted(projected, key=lambda item: item[0]):
                is_mylar = component.get("kind") == "mylar_overlay"
                canvas.create_polygon([coord for point in points for coord in point[:2]],
                                      fill=MYLAR_DISPLAY_COLOR if is_mylar else color,
                                      outline=MYLAR_DISPLAY_EDGE if is_mylar else "#263238",
                                      stipple="gray25" if is_mylar else "gray50")
                line_options = {"fill": MYLAR_DISPLAY_EDGE if is_mylar else "#263238", "width": 2}
                if is_mylar:
                    line_options["dash"] = (5, 3)
                canvas.create_line([coord for point in points + [points[0]] for coord in point[:2]], **line_options)
                if original_face_index in set(component.get("mylar_indices", ())):
                    canvas.create_line([coord for point in points + [points[0]] for coord in point[:2]],
                                       fill="#b26a00", dash=(5, 3), width=2)
            if projected:
                visible_points = [point[:2] for _depth, _index, polygon in projected for point in polygon]
                target = (
                    sum(point[0] for point in visible_points) / len(visible_points),
                    sum(point[1] for point in visible_points) / len(visible_points),
                )
                callout_items.append({
                    "component_id": component.get("component_id", component.get("name", label)),
                    "assembly_type": component.get("assembly_type"),
                    "kind": component.get("kind", "geometry"),
                    "label": detector_callout_label(component),
                    "detail": label,
                    "target": target,
                })
            for interface in component.get("interface_faces_mm", ()):
                points = [self._global_assembly_viewer_projection(_global_assembly_render_point(point, scene), center, scale, width, height) for point in interface]
                canvas.create_line([coord for point in points + [points[0]] for coord in point[:2]], fill="#00695c", dash=(2, 2), width=3)
            for readout in component.get("readout_faces_mm", ()):
                points = [self._global_assembly_viewer_projection(_global_assembly_render_point(point, scene), center, scale, width, height) for point in readout]
                canvas.create_polygon([coord for point in points for coord in point[:2]], fill="#9b59b6", outline="#4a235a", stipple="gray50")
                canvas.create_line([coord for point in points + [points[0]] for coord in point[:2]], fill="#4a235a", width=3)
        label_layout = collision_aware_callout_layout(
            callout_items, width, height,
            occupied_rectangles=scene_reserved_rectangles(width, height),
            show_labels=self._three_d_labels_enabled(), text_measure=label_measure,
        )
        self.global_viewer_label_layout = label_layout
        for callout in label_layout["callouts"]:
            start_x, start_y = callout["leader_start"]
            target_x, target_y = callout["target"]
            canvas.create_line(start_x, start_y, target_x, target_y,
                               fill="#455a64", dash=(3, 2), width=1)
            canvas.create_text(callout["position"][0], callout["position"][1],
                               text=callout["label"], anchor="nw", fill="#17324d",
                               font=label_font)
        legend_lines = compact_scene_legend_lines(callout_items, label_layout["omitted"], scene.get("mode", "combined"))
        legend_x1, legend_y1, legend_x2, legend_y2 = scene_legend_rect(width, height)
        canvas.create_rectangle(legend_x1, legend_y1, legend_x2, legend_y2,
                                fill="#ffffff", outline="#b0bec5")
        canvas_legend_lines = (
            f"MODE: {scene.get('mode', 'combined').upper()}",
            "C quartz + guide | S active + guide",
            "C/S sensors + terminals",
            "green dashed = shared particle path",
            "teal = mirrors | gold = Mylar",
            "photon paths = detector outcomes",
        )
        for row, line in enumerate(canvas_legend_lines):
            canvas.create_text(legend_x1 + 7, legend_y1 + 7 + row * 14,
                               text=line, anchor="nw", fill="#37474f",
                               font=("TkDefaultFont", 8, "bold" if row == 0 else "normal"))
        measurement_x1, measurement_y1, measurement_x2, measurement_y2 = scene_measurement_rect(width, height)
        canvas.create_rectangle(measurement_x1, measurement_y1, measurement_x2, measurement_y2,
                                fill="#fffdf5", outline="#c79a37")
        measurement_lines = (
            "MEASUREMENTS • +X longitudinal",
            f"C +X → S −X gap: {scene['cher_d1_gap_mm']:g} mm",
            f"D1(C) → D2(S) separation: {scene['d1_d2_separation_mm']:g} mm",
            f"Δ assembly offset S−C: {scene['relative_global_assembly_offset_mm']} mm",
            "purple = C/S readout bodies",
        )
        for row, line in enumerate(measurement_lines):
            canvas.create_text(measurement_x1 + 7, measurement_y1 + 7 + row * 21,
                               text=line, anchor="nw", fill="#604a12",
                               font=("TkDefaultFont", 8, "bold" if row == 0 else "normal"))
        gap_start = _global_assembly_render_point(
            (scene["quartz_high_x_mm"], scene["gap_marker_y_mm"], scene["gap_marker_z_mm"]), scene
        )
        gap_end = _global_assembly_render_point(
            (scene["assembly_low_x_mm"], scene["gap_marker_y_mm"], scene["gap_marker_z_mm"]), scene
        )
        gap_start_screen = self._global_assembly_viewer_projection(gap_start, center, scale, width, height)
        gap_end_screen = self._global_assembly_viewer_projection(gap_end, center, scale, width, height)
        canvas.create_line(gap_start_screen[0], gap_start_screen[1], gap_end_screen[0], gap_end_screen[1],
                           fill="#b26a00", dash=(5, 3), width=2, arrow="both")
        # Numeric gap/separation details live in the fixed sidebar; the canvas
        # retains only the color-coded measurement strokes.
        sep_start = self._global_assembly_viewer_projection(_global_assembly_render_point((scene["d1_d2_from_x_mm"], scene["separation_marker_y_mm"], scene["separation_marker_z_mm"]), scene), center, scale, width, height)
        sep_end = self._global_assembly_viewer_projection(_global_assembly_render_point((scene["d1_d2_to_x_mm"], scene["separation_marker_y_mm"], scene["separation_marker_z_mm"]), scene), center, scale, width, height)
        canvas.create_line(sep_start[0], sep_start[1], sep_end[0], sep_end[1], fill="#d97706", dash=(4, 2), width=2, arrow="both")
        # Keep the canvas free of long geometry-adjacent diagnostics.
        display_path = scene.get("display_particle_path_mm", scene["particle_path_mm"])
        for point_index, point in enumerate(display_path):
            projected_point = self._global_assembly_viewer_projection(_global_assembly_render_point(point, scene), center, scale, width, height)
            if point_index:
                previous = self._global_assembly_viewer_projection(_global_assembly_render_point(display_path[point_index - 1], scene), center, scale, width, height)
                canvas.create_line(previous[0], previous[1], projected_point[0], projected_point[1], fill=GLOBAL_PARTICLE_PATH_COLOR, dash=(7, 3), width=2)
            canvas.create_oval(projected_point[0] - 3, projected_point[1] - 3, projected_point[0] + 3, projected_point[1] + 3, fill=GLOBAL_PARTICLE_PATH_COLOR, outline="#ffffff")
        # Orientation gizmo: the origin is screen-anchored for readability,
        # while every endpoint is the current camera projection of the
        # canonical global basis.  It therefore rotates with the orbiting
        # figure and never implies that the assemblies themselves were
        # physically rotated.
        triad_x, triad_y = 42, height - 58
        triad_vectors = global_assembly_orientation_triad(
            _get_instance_attr(self, "global_viewer_yaw", GLOBAL_PREVIEW_DEFAULT_YAW),
            _get_instance_attr(self, "global_viewer_pitch", GLOBAL_PREVIEW_DEFAULT_PITCH),
        )
        axis_colors = {"X": "#c62828", "Y": "#6a1b9a", "Z": "#1565c0"}
        for axis in ("X", "Y", "Z"):
            vector = triad_vectors[axis]
            vector_length = max(math.hypot(vector[0], vector[1]), 1e-12)
            endpoint = (
                triad_x + vector[0] / vector_length * 42.0,
                triad_y - vector[1] / vector_length * 42.0,
            )
            canvas.create_line(triad_x, triad_y, endpoint[0], endpoint[1],
                               fill=axis_colors[axis], width=2, arrow=TK_ARROW_LAST)
            canvas.create_text(endpoint[0] + 5, endpoint[1], text=f"+{axis}",
                               fill=axis_colors[axis], anchor="w", font=("TkDefaultFont", 9, "bold"))
        canvas.create_text(triad_x, triad_y + 13, text="global XYZ • camera-coupled",
                           fill="#455a64", anchor="n", font=("TkDefaultFont", 8))
        mode = scene["display_mode"]
        canvas.create_rectangle(0, 0, width, 52, fill="#f5f8fb", outline="")
        canvas.create_text(width / 2, 8, text=f"GLOBAL ASSEMBLY 3D • {mode}",
                           font=("TkDefaultFont", 12, "bold"), fill="#17324d", anchor="n")
        canvas.create_text(width / 2, 31, text="C Cherenkov quartz → S scintillator • true geometry • shared particle path",
                           anchor="n", fill="#8a5200", font=("TkDefaultFont", 9))
        canvas.create_text(8, height - 10, text="Drag: rotate • Wheel: zoom • preset buttons change camera • Reset restores view",
                           anchor="sw", fill="#555555")
        info = _get_instance_attr(self, "global_viewer_info", None)
        if info is not None:
            detail_lines = compact_scene_legend_lines(callout_items, label_layout["omitted"])
            self._set_viewer_details_text(
                _get_instance_attr(self, "global_viewer_info_text", None),
                f"{scene['status_summary']}\n\n"
                f"MODE\n{mode}\n\n"
                f"MEASUREMENTS (+X global assembly axis; signed to−from)\n"
                f"{scene['measurement_definitions']['cher_d1_gap']['label']}\n{scene['cher_d1_gap_mm']:g} mm\n"
                f"{scene['measurement_definitions']['d1_d2_separation']['label']}\n{scene['d1_d2_separation_mm']:g} mm\n"
                f"{scene['measurement_definitions']['relative_global_assembly_offset']['label']}\n{scene['relative_global_assembly_offset_mm']} mm\n\n"
                "LEGEND\n" + "\n".join(f"  {line}" for line in detail_lines) + "\n\n"
                f"Callouts shown: {len(label_layout['callouts'])}; in legend: {len(label_layout['omitted'])}"
            )

    def _reset_global_assembly_view(self) -> None:
        self.global_viewer_yaw, self.global_viewer_pitch, self.global_viewer_zoom = .65, .55, 1.0
        self.global_view_preset = "Isometric"
        self.global_viewer_drag = None
        self._draw_global_assembly_3d_view()

    def _set_global_assembly_view(self, preset: str) -> None:
        presets = {
            "View along X": (math.pi / 2.0, 1.50),
            "View along Y": (0.0, 1.50),
            "View along Z": (0.0, 0.0),
            "View XY": (0.0, 0.0),
            "View XZ": (0.0, 1.50),
            "View YZ": (math.pi / 2.0, 1.50),
            "Isometric": (.65, .55),
            "Particle path through detectors": (.65, .55),
        }
        if preset not in presets:
            raise ValueError(f"unknown global assembly view: {preset}")
        self.global_viewer_yaw, self.global_viewer_pitch = presets[preset]
        self.global_viewer_zoom = 1.0
        self.global_view_preset = preset
        self._draw_global_assembly_3d_view()

    def _fit_global_assembly(self) -> None:
        self.global_viewer_zoom = 1.0
        self._draw_global_assembly_3d_view()

    def _global_view_along_x(self) -> None:
        self._set_global_assembly_view("View along X")

    def _global_view_along_y(self) -> None:
        self._set_global_assembly_view("View along Y")

    def _global_view_along_z(self) -> None:
        self._set_global_assembly_view("View along Z")

    def _global_view_isometric(self) -> None:
        self._set_global_assembly_view("Isometric")

    def _global_view_xy(self) -> None:
        self._set_global_assembly_view("View XY")

    def _global_view_xz(self) -> None:
        self._set_global_assembly_view("View XZ")

    def _global_view_yz(self) -> None:
        self._set_global_assembly_view("View YZ")

    def _global_view_particle_path(self) -> None:
        self._set_global_assembly_view("Particle path through detectors")

    def _global_assembly_viewer_press(self, event) -> None:
        self.global_viewer_drag = (event.x, event.y, self.global_viewer_yaw, self.global_viewer_pitch)
        event.widget.focus_set()

    def _global_assembly_viewer_drag(self, event) -> None:
        if self.global_viewer_drag is None:
            return
        x0, y0, yaw, pitch = self.global_viewer_drag
        self.global_viewer_yaw = yaw + (event.x - x0) * .012
        self.global_viewer_pitch = max(-1.45, min(1.45, pitch + (event.y - y0) * .012))
        self._draw_global_assembly_3d_view()

    def _global_assembly_viewer_zoom(self, event) -> None:
        delta = getattr(event, "delta", 0) or (120 if getattr(event, "num", 4) == 4 else -120)
        self.global_viewer_zoom = max(.2, min(8.0, self.global_viewer_zoom * (1.12 if delta > 0 else 1 / 1.12)))
        self._draw_global_assembly_3d_view()

    def _schedule_global_assembly_3d_draw(self, viewer=None, canvas=None, _event=None) -> None:
        """Debounce mapped/resize redraws and reject callbacks from old popups."""

        viewer = viewer or _get_instance_attr(self, "global_viewer", None)
        canvas = canvas or _get_instance_attr(self, "global_viewer_canvas", None)
        if viewer is not _get_instance_attr(self, "global_viewer", None) or canvas is not _get_instance_attr(self, "global_viewer_canvas", None):
            return
        try:
            if not viewer.winfo_exists() or not canvas.winfo_exists():
                return
            width, height = int(canvas.winfo_width()), int(canvas.winfo_height())
        except (AttributeError, tk.TclError, TypeError, ValueError):
            return
        if _get_instance_attr(self, "_global_viewer_draw_after", None) is not None:
            return

        def redraw_when_ready():
            self._global_viewer_draw_after = None
            if viewer is not _get_instance_attr(self, "global_viewer", None) or canvas is not _get_instance_attr(self, "global_viewer_canvas", None):
                return
            try:
                current_width, current_height = int(canvas.winfo_width()), int(canvas.winfo_height())
            except (AttributeError, tk.TclError, TypeError, ValueError):
                return
            if current_width <= 1 or current_height <= 1:
                try:
                    self._global_viewer_draw_after = viewer.after(30, redraw_when_ready)
                except (AttributeError, tk.TclError, RuntimeError):
                    pass
                return
            self._draw_global_assembly_3d_view(
                current_width, current_height, viewer=viewer, canvas=canvas,
            )

        try:
            # after_idle lets grid/map settle; the callback then reads fresh
            # dimensions instead of trusting a Configure event's old size.
            self._global_viewer_draw_after = viewer.after_idle(redraw_when_ready)
        except (AttributeError, tk.TclError, RuntimeError):
            self._global_viewer_draw_after = None
            if width > 1 and height > 1:
                self._draw_global_assembly_3d_view(width, height, viewer=viewer, canvas=canvas)

    def _close_global_assembly_3d_viewer(self, viewer=None) -> None:
        """Destroy one popup without allowing its queued callbacks to draw in a new one."""

        viewer = viewer or _get_instance_attr(self, "global_viewer", None)
        if viewer is _get_instance_attr(self, "global_viewer", None):
            pending = _get_instance_attr(self, "_global_viewer_draw_after", None)
            if pending is not None:
                try:
                    viewer.after_cancel(pending)
                except (AttributeError, tk.TclError, RuntimeError):
                    pass
            self._global_viewer_draw_after = None
            self.global_viewer = None
            self.global_viewer_canvas = None
            self.global_viewer_info_text = None
            self.global_viewer_buttons = {}
        try:
            if viewer is not None and viewer.winfo_exists():
                viewer.destroy()
        except (AttributeError, tk.TclError, RuntimeError):
            pass

    def _open_global_assembly_3d_viewer(self) -> None:
        existing = _get_instance_attr(self, "global_viewer", None)
        try:
            existing_alive = existing is not None and bool(existing.winfo_exists())
        except (AttributeError, tk.TclError, RuntimeError):
            existing_alive = False
        if existing_alive:
            required_buttons = global_assembly_view_preset_names()
            current_buttons = _get_instance_attr(self, "global_viewer_buttons", {})
            current_banner = _get_instance_attr(self, "global_viewer_banner", None)
            try:
                complete = (
                    current_banner is not None and current_banner.winfo_exists()
                    and all(label in current_buttons and current_buttons[label].winfo_exists() for label in required_buttons)
                    and _get_instance_attr(self, "global_viewer_canvas", None) is not None
                )
            except (AttributeError, tk.TclError, RuntimeError):
                complete = False
            if complete:
                existing.title("Global Beam Monitor Assembly - 3D Viewer")
                current_banner.configure(text="NEW GLOBAL ASSEMBLY 3D VIEWER")
                for button in current_buttons.values():
                    button.configure(state="normal")
                existing.deiconify()
                existing.lift()
                existing.update_idletasks()
                self._schedule_global_assembly_3d_draw(existing, self.global_viewer_canvas)
                return
            self._close_global_assembly_3d_viewer(existing)
        viewer = tk.Toplevel(self)
        viewer.title("Global Beam Monitor Assembly - 3D Viewer")
        viewer.geometry("1120x760")
        viewer.minsize(900, 600)
        viewer.rowconfigure(0, weight=0)
        viewer.rowconfigure(1, weight=1)
        viewer.rowconfigure(2, weight=0)
        viewer.columnconfigure(0, weight=1)
        viewer.columnconfigure(1, weight=0)
        self.global_viewer, self.global_viewer_canvas = viewer, None
        self.global_viewer_yaw, self.global_viewer_pitch, self.global_viewer_zoom, self.global_viewer_drag = .65, .55, 1.0, None
        self.global_view_preset = "Isometric"
        title_frame = ttk.Frame(viewer, padding=(8, 6, 8, 2))
        title_frame.grid(row=0, column=0, columnspan=2, sticky="ew")
        self.global_viewer_banner = ttk.Label(
            title_frame, text="NEW GLOBAL ASSEMBLY 3D VIEWER",
            font=("TkDefaultFont", 12, "bold"), foreground="#b00020",
        )
        self.global_viewer_banner.pack(side="left")
        ttk.Label(title_frame, text="Cherenkov assembly → separation → Scintillator assembly → current readouts",
                  foreground="#555555").pack(side="left", padx=14)
        canvas = tk.Canvas(viewer, background="white", highlightthickness=1, highlightbackground="#909090", takefocus=True)
        # Store both references before grid/map can emit a Configure event.
        self.global_viewer_canvas = canvas
        canvas.grid(row=1, column=0, sticky="nsew", padx=6, pady=6)
        info_frame = ttk.LabelFrame(viewer, text="Assembly information / legend", padding=8)
        info_frame.grid(row=1, column=1, sticky="ns", padx=(0, 6), pady=6)
        info_frame.columnconfigure(0, weight=1)
        info_frame.rowconfigure(0, weight=1)
        self.global_viewer_info = tk.StringVar(value="Preparing global assembly view...")
        try:
            info_text = tk.Text(info_frame, width=38, height=32, wrap="word", relief="flat", borderwidth=0,
                                background=_ttk_background(info_frame), takefocus=True)
            info_text.grid(row=0, column=0, sticky="nsew")
            info_scroll = ttk.Scrollbar(info_frame, orient="vertical", command=info_text.yview)
            info_scroll.grid(row=0, column=1, sticky="ns")
            info_text.configure(yscrollcommand=info_scroll.set, state="disabled")
            self.global_viewer_info_text = info_text
        except Exception as exc:
            # A platform-specific Text/style failure must be visible in the
            # canvas and must not prevent the controls or geometry from being
            # created.
            self.global_viewer_info_text = None
            self._global_viewer_render_error(canvas, exc)
        controls = ttk.Frame(viewer)
        controls.grid(row=2, column=0, columnspan=2, sticky="ew", padx=6, pady=(0, 6))
        self.global_viewer_buttons = {}
        for label, callback, padding in (
            ("Reset view", self._reset_global_assembly_view, (0, 0)),
            ("Isometric", self._global_view_isometric, (8, 0)),
            ("View along X", self._global_view_along_x, (4, 0)),
            ("View along Y", self._global_view_along_y, (4, 0)),
            ("View along Z", self._global_view_along_z, (4, 0)),
            ("View XY", self._global_view_xy, (4, 0)),
            ("View XZ", self._global_view_xz, (4, 0)),
            ("View YZ", self._global_view_yz, (4, 0)),
            ("Particle path through detectors", self._global_view_particle_path, (4, 0)),
            ("Fit assembly", self._fit_global_assembly, (4, 0)),
        ):
            button = ttk.Button(controls, text=label, command=callback)
            button.pack(side="left", padx=padding)
            self.global_viewer_buttons[label] = button
        ttk.Checkbutton(controls, text="Show compact callouts", variable=self.show_3d_labels,
                        command=self._draw_global_assembly_3d_view).pack(side="left", padx=(12, 0))
        ttk.Label(controls, text="Drag to rotate • wheel to zoom", foreground="#555555").pack(side="left", padx=12)
        ttk.Button(controls, text="Close", command=viewer.destroy).pack(side="right")
        canvas.bind("<Configure>", lambda event, v=viewer, c=canvas: self._schedule_global_assembly_3d_draw(v, c, event))
        canvas.bind("<ButtonPress-1>", self._global_assembly_viewer_press)
        canvas.bind("<B1-Motion>", self._global_assembly_viewer_drag)
        canvas.bind("<ButtonRelease-1>", lambda _event: setattr(self, "global_viewer_drag", None))
        canvas.bind("<MouseWheel>", self._global_assembly_viewer_zoom)
        canvas.bind("<Button-4>", self._global_assembly_viewer_zoom)
        canvas.bind("<Button-5>", self._global_assembly_viewer_zoom)
        viewer.protocol("WM_DELETE_WINDOW", lambda v=viewer: self._close_global_assembly_3d_viewer(v))
        viewer.update_idletasks()
        self._schedule_global_assembly_3d_draw(viewer, canvas)

    def _build_run_tab(self, parent: ttk.Frame) -> None:
        # Keep the two top panels as independent, bounded cells.  In v44 the
        # button/status row's requested width could make a long worker phase
        # resize this grid and squeeze the Simulation parameters column.
        parent.columnconfigure(0, weight=3, uniform="run_top_panels", minsize=460)
        parent.columnconfigure(1, weight=2, uniform="run_top_panels", minsize=360)
        parent.rowconfigure(0, weight=0, minsize=190)
        parent.rowconfigure(2, weight=3)
        parent.rowconfigure(4, weight=1)

        settings = ttk.LabelFrame(parent, text="Simulation parameters", padding=5)
        settings.grid(row=0, column=0, sticky="nsew", padx=(0, 5))
        settings.columnconfigure(1, weight=1)
        ttk.Label(settings, text="Output directory").grid(row=0, column=0, sticky="w", padx=(0, 8), pady=3)
        ttk.Entry(settings, textvariable=self.values["output_directory"]).grid(row=0, column=1, sticky="ew", pady=3)
        ttk.Button(settings, text="Browse", command=self._browse_output).grid(row=0, column=2, padx=(5, 0), pady=3)
        self._build_grid(settings, [
            ("Events", "events"), ("Photons per event", "photons_per_event"),
            ("Photon sampling mode", "sampling_mode"), ("Seed", "seed"),
            ("Track margin (mm)", "track_margin_mm"), ("Display photons", "display_photons"),
            ("Detectors to simulate", "detector_run_mode"),
            ("Inclination from +Z / track tilt (deg)", "inclination_deg"),
            ("Maximum photon bounces (selected detector(s))", "maximum_photon_bounces"),
        ], start_row=1)
        gamma_frame = ttk.LabelFrame(settings, text="Gamma transport (only when Gamma transport is selected)", padding=4)
        gamma_frame.grid(row=11, column=0, columnspan=3, sticky="ew", pady=(5, 0))
        gamma_frame.columnconfigure(1, weight=1)
        ttk.Label(gamma_frame, text="Approved source profile").grid(row=0, column=0, sticky="w", padx=(0, 6), pady=2)
        ttk.Entry(gamma_frame, textvariable=self.values["gamma_source_profile_path"]).grid(row=0, column=1, sticky="ew", pady=2)
        ttk.Button(gamma_frame, text="Browse", command=self._browse_gamma_profile).grid(row=0, column=2, padx=(5, 0), pady=2)
        line_fields = (
            ("511", "gamma_line_weight_511_keV"), ("835", "gamma_line_weight_835_keV"),
            ("1173", "gamma_line_weight_1173_keV"), ("1275", "gamma_line_weight_1275_keV"),
            ("1332", "gamma_line_weight_1332_keV"), ("1461", "gamma_line_weight_1461_keV"),
        )
        for index, (energy, key) in enumerate(line_fields):
            row, column = 1 + index // 3, (index % 3) * 2
            ttk.Label(gamma_frame, text=f"Line {energy} keV weight").grid(row=row, column=column, sticky="w", padx=(0, 4), pady=2)
            ttk.Entry(gamma_frame, textvariable=self.values[key], width=10).grid(row=row, column=column + 1, sticky="ew", padx=(0, 8), pady=2)
        ttk.Label(gamma_frame, text="All line weights default to 0/unset. Continuum follows the approved relative profile. Weights are not calibrated intensities.", foreground="#8a5200", wraplength=470, justify="left").grid(row=3, column=0, columnspan=3, sticky="w", pady=(3, 2))
        ttk.Label(gamma_frame, text="Maximum interactions / light samples per run").grid(row=4, column=0, sticky="w", padx=(0, 6), pady=2)
        ttk.Entry(gamma_frame, textvariable=self.values["gamma_max_interactions"], width=10).grid(row=4, column=1, sticky="w", pady=2)
        ttk.Entry(gamma_frame, textvariable=self.values["gamma_max_light_samples"], width=10).grid(row=4, column=2, sticky="e", pady=2)
        ttk.Label(
            settings,
            text="All combines Cherenkov + Scintillator optical transport. Gamma transport is a separate explicit run choice; Gamma Results keep the last completed Gamma run across later optical runs and startup. Gamma mode uses the approved collimated +Z source and current assembly snapshot. Its expected local light is not transported or detected signal.",
            foreground="#245a35", wraplength=470, justify="left",
        ).grid(row=5, column=0, columnspan=3, sticky="w", pady=(7, 2))

        action_frame = ttk.LabelFrame(parent, text="Run controls and status", padding=5)
        action_frame.grid(row=0, column=1, sticky="nsew")
        action_frame.columnconfigure(0, weight=1)
        for row, minimum in ((0, 32), (1, 34), (2, 22), (3, 34), (4, 30), (5, 28)):
            action_frame.rowconfigure(row, weight=0, minsize=minimum)
        action_frame.grid_propagate(False)
        action_frame.configure(width=430, height=190)
        button_row = ttk.Frame(action_frame)
        button_row.grid(row=0, column=0, sticky="ew")
        button_row.grid_propagate(False)
        button_row.configure(width=418)
        self.run_button = ttk.Button(button_row, text="Run Simulation", command=self._run)
        self.run_button.pack(side="left")
        self.stop_button = ttk.Button(button_row, text="Stop Simulation", command=self._stop, state="disabled")
        self.stop_button.pack(side="left", padx=5)
        ttk.Button(button_row, text="Open Output Folder", command=self._open_output_folder).pack(side="left", padx=2)
        ttk.Button(button_row, text="Open Latest Run Folder", command=self._open_latest_run_folder).pack(side="left", padx=2)
        # A fixed-height clipped status cell prevents worker phase text from
        # changing the requested height of the action panel.  The progress
        # bar is always the next fixed row, so its y position is invariant.
        status_cell = ttk.Frame(action_frame, width=418, height=34)
        status_cell.grid(row=1, column=0, sticky="ew", pady=(8, 2))
        status_cell.grid_propagate(False)
        self.status_label = ttk.Label(
            status_cell, textvariable=self.status, anchor="w", justify="left",
            wraplength=398, width=57,
        )
        self.status_label.grid(row=0, column=0, sticky="ew")
        self.progress = ttk.Progressbar(action_frame, mode="determinate", maximum=100.0)
        self.progress.grid(row=2, column=0, sticky="ew", pady=(2, 6))
        self.progress_text = tk.StringVar(value="Ready")
        self.progress_label = ttk.Label(
            action_frame, textvariable=self.progress_text, foreground="#245a35",
            anchor="w", justify="left", wraplength=404,
        )
        self.progress_label.grid(row=3, column=0, sticky="ew", pady=(0, 3))
        event_row = ttk.Frame(action_frame)
        event_row.grid(row=4, column=0, sticky="w")
        ttk.Label(event_row, text="Event index").pack(side="left")
        self.bottom_event_entry = ttk.Entry(event_row, textvariable=self.selected_event, width=6)
        self.bottom_event_entry.pack(side="left", padx=5)
        self.bottom_event_entry.bind("<Return>", lambda _event: self._show_selected_event())
        self.bottom_show_event_button = ttk.Button(event_row, text="Show Selected Event", command=self._show_selected_event)
        self.bottom_show_event_button.pack(side="left")
        self.top_event_3d_button = ttk.Button(event_row, text="Interactive 3D Viewer", command=self._open_selected_event_3d_viewer)
        self.top_event_3d_button.pack(side="left", padx=(5, 0))
        ttk.Button(action_frame, text="Quit", command=self.destroy).grid(row=5, column=0, sticky="e", pady=(7, 0))

        ttk.Label(parent, text="Event/run results", font=("TkDefaultFont", 10, "bold")).grid(row=1, column=0, columnspan=2, sticky="w", pady=(7, 2))
        results_frame = ttk.Frame(parent)
        results_frame.grid(row=2, column=0, columnspan=2, sticky="nsew")
        results_frame.columnconfigure(0, weight=3)
        results_frame.columnconfigure(1, weight=2)
        results_frame.rowconfigure(0, weight=1)
        self._build_event_results(results_frame)

        viewer_frame = ttk.LabelFrame(parent, text="Selected event viewers", padding=5)
        viewer_frame.grid(row=3, column=0, columnspan=2, sticky="ew", pady=(5, 0))
        viewer_frame.columnconfigure(0, weight=1)
        button_row = ttk.Frame(viewer_frame)
        button_row.grid(row=0, column=0, sticky="ew")
        ttk.Button(button_row, text="Show static event image", command=self._show_selected_event).pack(side="left")
        ttk.Button(button_row, text="Open interactive 3D event viewer", command=self._open_selected_event_3d_viewer).pack(side="left", padx=(6, 0))
        filter_frame = ttk.LabelFrame(viewer_frame, text="Photon categories shown in both viewers", padding=(4, 2))
        filter_frame.grid(row=1, column=0, sticky="ew", pady=(5, 0))
        for column, (key, label) in enumerate(EVENT_VIEW_CATEGORY_LABELS.items()):
            ttk.Checkbutton(
                filter_frame,
                text=label,
                variable=self.event_view_category_vars[key],
                command=self._event_view_filters_changed,
            ).grid(row=column // 5, column=column % 5, sticky="w", padx=(0, 10), pady=1)
        ttk.Label(
            filter_frame,
            text="All photon paths overrides the outcome filters; turn it off to use the independent categories.",
            foreground="#555555",
        ).grid(row=2, column=0, columnspan=5, sticky="w", pady=(2, 0))
        self.event_view_status = tk.StringVar(value="Select an event row, then open its static image or interactive 3D viewer.")
        ttk.Label(viewer_frame, textvariable=self.event_view_status, foreground="#555555").grid(row=2, column=0, sticky="w", pady=(3, 0))

        log_frame = ttk.LabelFrame(parent, text="Run log and summary", padding=3)
        log_frame.grid(row=4, column=0, columnspan=2, sticky="nsew", pady=(5, 0))
        log_frame.columnconfigure(0, weight=1)
        log_frame.rowconfigure(0, weight=1)
        self.log_text = tk.Text(log_frame, height=4, wrap="word", state="disabled")
        self.log_text.grid(row=0, column=0, sticky="nsew")
        scrollbar = ttk.Scrollbar(log_frame, orient="vertical", command=self.log_text.yview)
        scrollbar.grid(row=0, column=1, sticky="ns")
        self.log_text.configure(yscrollcommand=scrollbar.set)

    def _build_event_results(self, parent: ttk.Frame) -> None:
        table_frame = ttk.LabelFrame(parent, text="Per-event results", padding=3)
        table_frame.grid(row=0, column=0, sticky="nsew", padx=(0, 5))
        table_frame.columnconfigure(0, weight=1)
        table_frame.rowconfigure(0, weight=1)
        columns = ("event_id", "mode", "expected", "generated", "transported", "tir_plus", "crossings", "valid_hits", "detected", "plus_x", "minus_x", "mean_time", "requested", "attempts", "tir_minus", "losses", "scint_generated", "scint_readout", "scint_detected", "mylar_reflections", "mylar_losses")
        self.event_tree = ttk.Treeview(table_frame, columns=columns, show="headings")
        headings = {
            "event_id": "Event", "mode": "Mode", "expected": "FT yield", "requested": "Requested",
            "generated": "Generated", "transported": "Transported", "attempts": "Detector attempts",
            "tir_plus": "TIR +x", "tir_minus": "TIR -x", "crossings": "Face crossings",
            "valid_hits": "Valid window hits",
            "losses": "Detector losses", "detected": "Detected", "plus_x": "Detected +x",
            "minus_x": "Detected -x", "mean_time": "Mean arrival (ns)",
            "scint_generated": "Scint generated", "scint_readout": "Scint readout", "scint_detected": "Scint detected", "mylar_reflections": "Mylar refl.", "mylar_losses": "Mylar loss",
        }
        widths = {"event_id": 45, "mode": 60, "expected": 65, "requested": 70, "generated": 65, "transported": 75, "attempts": 92, "tir_plus": 55, "tir_minus": 55, "crossings": 78, "valid_hits": 92, "losses": 84, "detected": 65, "plus_x": 75, "minus_x": 75, "mean_time": 100, "scint_generated": 90, "scint_readout": 82, "scint_detected": 85, "mylar_reflections": 75, "mylar_losses": 70}
        for column in columns:
            self.event_tree.heading(column, text=headings[column])
            self.event_tree.column(column, width=widths[column], anchor="center", stretch=False)
        self.event_tree.grid(row=0, column=0, sticky="nsew")
        table_scroll = ttk.Scrollbar(table_frame, orient="vertical", command=self.event_tree.yview)
        table_scroll.grid(row=0, column=1, sticky="ns")
        table_hscroll = ttk.Scrollbar(table_frame, orient="horizontal", command=self.event_tree.xview)
        table_hscroll.grid(row=1, column=0, sticky="ew")
        self.event_tree.configure(yscrollcommand=table_scroll.set, xscrollcommand=table_hscroll.set)
        self.event_tree.bind("<<TreeviewSelect>>", self._table_selection_changed)

        detail_frame = ttk.LabelFrame(parent, text="Selected event details", padding=5)
        detail_frame.grid(row=0, column=1, sticky="nsew")
        detail_frame.columnconfigure(0, weight=1)
        detail_frame.rowconfigure(0, weight=1)
        self.details_text = tk.Text(detail_frame, height=10, width=42, wrap="word", state="disabled")
        self.details_text.grid(row=0, column=0, sticky="nsew")
        details_scroll = ttk.Scrollbar(detail_frame, orient="vertical", command=self.details_text.yview)
        details_scroll.grid(row=0, column=1, sticky="ns")
        self.details_text.configure(yscrollcommand=details_scroll.set)
        ttk.Label(parent, textvariable=self.total_summary_text, foreground="#333333").grid(row=1, column=0, columnspan=2, sticky="w", pady=(4, 0))

    def _build_results_tab(self, parent: ttk.Frame) -> None:
        """Build the detector-owned readout arrival analysis workspace."""

        parent.columnconfigure(0, weight=0, minsize=265)
        parent.columnconfigure(1, weight=1)
        parent.rowconfigure(1, weight=1)
        parent.rowconfigure(2, weight=1)
        header = ttk.Frame(parent, padding=(4, 4, 4, 2))
        header.grid(row=0, column=0, columnspan=2, sticky="ew")
        header.columnconfigure(3, weight=1)
        ttk.Label(header, text="Results: readout arrivals and timing", font=("TkDefaultFont", 11, "bold")).grid(row=0, column=0, columnspan=2, sticky="w")
        ttk.Label(header, text="Event selector (event_id)").grid(row=1, column=0, sticky="w", pady=(5, 0))
        self.results_event_selector = ttk.Combobox(header, textvariable=self.selected_event, state="readonly", width=12, values=())
        self.results_event_selector.grid(row=1, column=1, sticky="w", padx=(6, 12), pady=(5, 0))
        self.results_event_selector.bind("<<ComboboxSelected>>", self._results_event_selected)
        self.results_event_selector.bind("<Return>", self._results_event_selected)
        ttk.Label(header, textvariable=tk.StringVar(value="Arrival = physical entry into the configured finite sensor window before angular acceptance, collection, or QE; units ns."), foreground="#245a35", wraplength=780, justify="left").grid(row=1, column=2, columnspan=2, sticky="ew", pady=(5, 0))
        self.results_status = tk.StringVar(value="No completed run results")
        ttk.Label(header, textvariable=self.results_status, foreground="#8a1c1c", wraplength=1040, justify="left").grid(row=2, column=0, columnspan=4, sticky="w", pady=(3, 0))

        selected_frame = ttk.LabelFrame(parent, text="Selected event", padding=5)
        selected_frame.grid(row=1, column=0, columnspan=2, sticky="nsew", padx=4, pady=(2, 4))
        selected_frame.columnconfigure(0, weight=0, minsize=250)
        selected_frame.columnconfigure(1, weight=1)
        selected_frame.rowconfigure(0, weight=1)
        selected_info = ttk.Frame(selected_frame)
        selected_info.grid(row=0, column=0, sticky="nsew", padx=(0, 8))
        selected_info.columnconfigure(0, weight=1)
        self.results_selected_card_vars = {}
        for row, detector in enumerate(("cherenkov", "scintillator")):
            card = ttk.LabelFrame(selected_info, text=detector.title(), padding=5)
            card.grid(row=row, column=0, sticky="ew", pady=(0, 5))
            card.columnconfigure(0, weight=1)
            variable = tk.StringVar(value="No completed run results")
            self.results_selected_card_vars[detector] = variable
            ttk.Label(card, textvariable=variable, justify="left", anchor="nw", wraplength=235).grid(row=0, column=0, sticky="ew")
        self.results_selected_figure, self.results_selected_canvas = self._make_results_figure(selected_frame, row=0, column=1)

        aggregate_frame = ttk.LabelFrame(parent, text="All events", padding=5)
        aggregate_frame.grid(row=2, column=0, columnspan=2, sticky="nsew", padx=4, pady=(0, 4))
        aggregate_frame.columnconfigure(0, weight=0, minsize=250)
        aggregate_frame.columnconfigure(1, weight=1)
        aggregate_frame.rowconfigure(0, weight=1)
        aggregate_info = ttk.Frame(aggregate_frame)
        aggregate_info.grid(row=0, column=0, sticky="nsew", padx=(0, 8))
        aggregate_info.columnconfigure(0, weight=1)
        self.results_aggregate_card_vars = {}
        for row, detector in enumerate(("cherenkov", "scintillator")):
            card = ttk.LabelFrame(aggregate_info, text=detector.title(), padding=5)
            card.grid(row=row, column=0, sticky="ew", pady=(0, 5))
            card.columnconfigure(0, weight=1)
            variable = tk.StringVar(value="No completed run results")
            self.results_aggregate_card_vars[detector] = variable
            ttk.Label(card, textvariable=variable, justify="left", anchor="nw", wraplength=235).grid(row=0, column=0, sticky="ew")
        self.results_aggregate_figure, self.results_aggregate_canvas = self._make_results_figure(aggregate_frame, row=0, column=1)

    def _build_gamma_results_tab(self, parent: ttk.Frame) -> None:
        """Build the persisted Gamma results and event-viewer workspace."""

        parent.columnconfigure(0, weight=1)
        parent.rowconfigure(3, weight=1)
        header = ttk.LabelFrame(parent, text="Gamma transport results", padding=7)
        header.grid(row=0, column=0, sticky="ew", padx=5, pady=5)
        header.columnconfigure(1, weight=1)
        header.columnconfigure(2, weight=0)
        ttk.Label(header, text="IMPORTANT: local expected light only — no optical-photon transport, no fiber transport, no sensor arrivals/detections.", foreground="#a32020", wraplength=980, justify="left").grid(row=0, column=0, columnspan=3, sticky="w", pady=(0, 5))
        ttk.Label(header, text="The All run is the Cherenkov + Scintillator optical combination. Gamma transport is a separate run mode; this tab retains its last completed result independently.", foreground="#245a35", wraplength=980, justify="left").grid(row=1, column=0, columnspan=3, sticky="w")
        self.gamma_results_status = tk.StringVar(value="No completed gamma run")
        ttk.Label(header, textvariable=self.gamma_results_status, foreground="#245a35", wraplength=980, justify="left").grid(row=2, column=0, columnspan=3, sticky="w")
        ttk.Label(header, text="Event").grid(row=3, column=0, sticky="w", pady=(5, 0))
        self.gamma_results_event_selector = ttk.Combobox(header, textvariable=self.gamma_selected_event, state="readonly", width=10, values=())
        self.gamma_results_event_selector.grid(row=3, column=1, sticky="w", padx=(6, 0), pady=(5, 0))
        self.gamma_results_event_selector.bind("<<ComboboxSelected>>", self._gamma_results_event_selected)
        self.gamma_global_3d_button = ttk.Button(
            header, text="Open interactive global 3D Gamma view",
            command=self._open_gamma_global_3d_viewer,
        )
        self.gamma_global_3d_button.grid(row=3, column=2, sticky="e", padx=(12, 0), pady=(5, 0))
        self.gamma_results_summary_text = tk.StringVar(value="No completed gamma run")
        ttk.Label(header, textvariable=self.gamma_results_summary_text, justify="left", wraplength=980).grid(row=4, column=0, columnspan=3, sticky="w", pady=(5, 0))
        aggregate = ttk.LabelFrame(parent, text="Aggregate by assembly / volume", padding=5)
        aggregate.grid(row=1, column=0, sticky="ew", padx=5, pady=(0, 5))
        aggregate.columnconfigure(0, weight=1)
        aggregate.columnconfigure(1, weight=0)
        self.gamma_results_aggregate_text = tk.Text(aggregate, height=6, wrap="word", state="disabled")
        self.gamma_results_aggregate_text.grid(row=0, column=0, sticky="ew")
        aggregate_scrollbar = ttk.Scrollbar(aggregate, orient="vertical", command=self.gamma_results_aggregate_text.yview)
        aggregate_scrollbar.grid(row=0, column=1, sticky="ns")
        self.gamma_results_aggregate_text.configure(yscrollcommand=aggregate_scrollbar.set)
        event_table_frame = ttk.LabelFrame(parent, text="Per-event results", padding=5)
        event_table_frame.grid(row=2, column=0, sticky="ew", padx=5, pady=(0, 5))
        event_table_frame.columnconfigure(0, weight=1)
        event_columns = (
            "event_id", "status", "source_energy", "source_component", "interaction_count",
            "processes", "deposits", "charged_tracks", "local_light",
        )
        self.gamma_results_event_table = ttk.Treeview(
            event_table_frame, columns=event_columns, show="headings", height=5, selectmode="browse",
        )
        event_headings = {
            "event_id": "Event ID", "status": "Status", "source_energy": "Source energy (keV)",
            "source_component": "Source component", "interaction_count": "Interactions",
            "processes": "Processes", "deposits": "Deposits (records / keV)",
            "charged_tracks": "Charged tracks", "local_light": "Local light (expected photons)",
        }
        event_widths = {
            "event_id": 62, "status": 100, "source_energy": 112, "source_component": 105,
            "interaction_count": 88, "processes": 170, "deposits": 150,
            "charged_tracks": 92, "local_light": 215,
        }
        for column in event_columns:
            self.gamma_results_event_table.heading(column, text=event_headings[column])
            self.gamma_results_event_table.column(column, width=event_widths[column], minwidth=55, anchor="w")
        self.gamma_results_event_table.grid(row=0, column=0, sticky="ew")
        table_xscroll = ttk.Scrollbar(event_table_frame, orient="horizontal", command=self.gamma_results_event_table.xview)
        table_xscroll.grid(row=1, column=0, sticky="ew")
        table_yscroll = ttk.Scrollbar(event_table_frame, orient="vertical", command=self.gamma_results_event_table.yview)
        table_yscroll.grid(row=0, column=1, sticky="ns")
        self.gamma_results_event_table.configure(xscrollcommand=table_xscroll.set, yscrollcommand=table_yscroll.set)
        self.gamma_results_event_table.bind("<<TreeviewSelect>>", self._gamma_results_table_selected)
        event_frame = ttk.LabelFrame(parent, text="Selected event: source, interactions, deposits, secondaries, local expected light", padding=5)
        event_frame.grid(row=3, column=0, sticky="nsew", padx=5, pady=(0, 5))
        event_frame.columnconfigure(0, weight=3)
        event_frame.columnconfigure(1, weight=2)
        event_frame.rowconfigure(0, weight=1)
        viewer_frame = ttk.LabelFrame(event_frame, text="Gamma event viewer · global X/Y/Z in mm", padding=3)
        viewer_frame.grid(row=0, column=0, sticky="nsew", padx=(0, 5))
        viewer_frame.columnconfigure(0, weight=1)
        viewer_frame.rowconfigure(0, weight=1)
        self.gamma_event_figure, self.gamma_event_canvas = self._make_results_figure(viewer_frame, row=0, column=0)
        details_frame = ttk.LabelFrame(event_frame, text="Event data", padding=3)
        details_frame.grid(row=0, column=1, sticky="nsew")
        details_frame.columnconfigure(0, weight=1)
        details_frame.rowconfigure(0, weight=1)
        self.gamma_results_event_text = tk.Text(details_frame, wrap="word", state="disabled")
        self.gamma_results_event_text.grid(row=0, column=0, sticky="nsew")
        scrollbar = ttk.Scrollbar(details_frame, orient="vertical", command=self.gamma_results_event_text.yview)
        scrollbar.grid(row=0, column=1, sticky="ns")
        self.gamma_results_event_text.configure(yscrollcommand=scrollbar.set)
        self._scrollable_mousewheel_passthrough.extend((
            self.gamma_results_aggregate_text,
            self.gamma_results_event_table,
            self.gamma_results_event_text,
        ))

    @staticmethod
    def _set_readonly_text(widget, text: str) -> None:
        widget.configure(state="normal")
        widget.delete("1.0", "end")
        widget.insert("1.0", text)
        widget.configure(state="disabled")

    def _gamma_results_event_selected(self, _event=None) -> None:
        self._refresh_gamma_results_view()

    def _gamma_results_table_selected(self, _event=None) -> None:
        table = _get_instance_attr(self, "gamma_results_event_table", None)
        if table is None:
            return
        selected = table.selection()
        if not selected:
            return
        try:
            event_id = int(selected[0])
        except (TypeError, ValueError):
            return
        self.gamma_selected_event.set(str(event_id))
        self._refresh_gamma_results_view()

    def _refresh_gamma_results_view(self) -> None:
        if not _has_instance_attr(self, "gamma_results_status"):
            return
        payload = _get_instance_attr(self, "latest_gamma_payload", None)
        if not isinstance(payload, Mapping):
            self.gamma_results_status.set("No completed gamma run")
            self.gamma_results_summary_text.set("No completed gamma run")
            self._set_readonly_text(self.gamma_results_aggregate_text, "No completed gamma run.")
            self._set_readonly_text(self.gamma_results_event_text, "No completed gamma run.")
            self.gamma_results_event_selector.configure(values=())
            table = _get_instance_attr(self, "gamma_results_event_table", None)
            if table is not None:
                table.delete(*table.get_children())
                self._gamma_event_table_payload = None
            draw_gamma_event_view(self.gamma_event_figure, None)
            self.gamma_event_canvas.draw_idle()
            QuartzCherenkovGUI._refresh_gamma_global_3d_view(self)
            return
        summary = payload.get("summary", {})
        events = list(payload.get("events", ()))
        table_rows = gamma_event_table_rows(events)
        event_ids = [int(item["event_id"]) for item in table_rows]
        self.gamma_results_event_selector.configure(values=tuple(str(value) for value in event_ids))
        table = _get_instance_attr(self, "gamma_results_event_table", None)
        if table is not None:
            row_ids = [str(row["event_id"]) for row in table_rows]
            rebuild_table = (
                _get_instance_attr(self, "_gamma_event_table_payload", None) is not payload
                or list(table.get_children()) != row_ids
            )
            if rebuild_table:
                table.delete(*table.get_children())
                for row in table_rows:
                    table.insert("", "end", iid=str(row["event_id"]), values=(
                        row["event_id"], row["status"], f"{row['source_energy_keV']:.6g}",
                        row["source_component"], row["interaction_count"],
                        row["interaction_process_summary"],
                        f"{row['deposit_record_count']} / {row['deposit_total_keV']:.6g}",
                        row["charged_track_count"], row["local_light_summary"],
                    ))
                self._gamma_event_table_payload = payload
        try:
            event_id = int(self.gamma_selected_event.get())
        except (TypeError, ValueError):
            event_id = event_ids[0] if event_ids else 0
        if event_ids and event_id not in event_ids:
            event_id = event_ids[0]
        self.gamma_selected_event.set(str(event_id))
        try:
            self.gamma_results_event_selector.set(str(event_id))
        except Exception:
            pass
        if table is not None and event_ids and list(table.selection()) != [str(event_id)]:
            table.selection_set(str(event_id))
        if table is not None and event_ids:
            table.see(str(event_id))
        direction = summary.get("source_direction_global", [0, 0, 1])
        active_lines = summary.get("active_line_energies_keV", [])
        active_weights = summary.get("active_line_weights", [])
        geometry_snapshot = payload.get("geometry_snapshot")
        geometry_message = (
            "geometry snapshot: complete Global Assembly scene"
            if gamma_geometry_view_data(geometry_snapshot)["components"]
            else "geometry snapshot unavailable for this saved result"
        )
        self.gamma_results_status.set(
            f"Model status: {payload.get('model_status', summary.get('model_status', 'unknown'))} · "
            f"events={summary.get('events', len(events))} · source +Z={direction} · "
            f"{geometry_message} · all light values below are expected local production, not transported/detected signal."
        )
        self.gamma_results_summary_text.set(
            f"Profile: {summary.get('source_profile', 'n/a')}\n"
            f"Active discrete lines (user weights only): {list(zip(active_lines, active_weights)) or 'none; continuum-only default'}\n"
            f"Event accounting: requested/completed/failed/excluded = {summary.get('event_accounting', {}).get('requested', summary.get('events', 0))} / {summary.get('event_accounting', {}).get('completed', summary.get('events', 0))} / {summary.get('event_accounting', {}).get('failed', 0)} / {summary.get('event_accounting', {}).get('excluded', 0)}; exclusion is not applicable to this no-sensor stage.\n"
            f"Total deposited: {float(summary.get('total_deposited_energy_keV', 0.0)):.6g} keV · "
            f"charged secondary tracks: {int(summary.get('total_charged_secondary_tracks', 0))} · "
            f"expected local Cherenkov: {float(summary.get('total_expected_cherenkov_photons', 0.0)):.6g} · "
            f"expected local scintillation: {float(summary.get('total_expected_scintillation_photons', 0.0)):.6g}"
        )
        aggregate_lines = ["Units: deposited energy keV/MeV; light counts are EXPECTED LOCAL photons only."]
        for assembly, item in sorted((summary.get("by_assembly", {}) or {}).items()):
            aggregate_lines.append(
                f"{assembly}: deposited={float(item.get('deposited_energy_keV', 0.0)):.6g} keV; "
                f"interactions={item.get('interactions', {})}; "
                f"local Cerenkov expected={float(item.get('expected_cherenkov_photons', 0.0)):.6g}; "
                f"local scintillation expected={float(item.get('expected_scintillation_photons', 0.0)):.6g}"
            )
        for volume, item in sorted((summary.get("by_volume", {}) or {}).items()):
            aggregate_lines.append(
                f"  {volume} [{item.get('material', 'unknown')}]: {float(item.get('deposited_energy_keV', 0.0)):.6g} keV; {item.get('interactions', {})}"
            )
        self._set_readonly_text(self.gamma_results_aggregate_text, "\n".join(aggregate_lines))
        selected = next((item for item in events if isinstance(item, Mapping) and int(_gamma_number(item.get("event_id", -1), -1)) == event_id), None)
        if selected is None:
            text = "No selected event."
            if not gamma_geometry_view_data(geometry_snapshot)["components"]:
                text += "\nDetector geometry is unavailable because this saved result has no v72 run geometry snapshot."
            self._set_readonly_text(self.gamma_results_event_text, text)
            draw_gamma_event_view(self.gamma_event_figure, None, geometry_snapshot)
            self.gamma_event_canvas.draw_idle()
            QuartzCherenkovGUI._refresh_gamma_global_3d_view(self)
            return
        deposits = [item for item in (selected.get("deposits", ()) or ()) if isinstance(item, Mapping)]
        tracks = [item for item in (selected.get("charged_tracks", ()) or ()) if isinstance(item, Mapping)]
        light_available = isinstance(selected.get("light"), (list, tuple))
        lights = [item for item in (selected.get("light", ()) or ()) if isinstance(item, Mapping)]
        selected_row = next((row for row in table_rows if row["event_id"] == event_id), {})
        event_lines = [
            f"Event {event_id} — {selected_row.get('status', 'UNKNOWN')}: source energy={_gamma_number(selected.get('source_energy_keV', 0.0)):.6g} keV; source component={selected.get('source_component', 'unknown')}",
            f"Interactions: {selected_row.get('interaction_count', 0)}; processes: {selected_row.get('interaction_process_summary', 'none')}",
            f"Deposits: {len(deposits)} records; total={sum(float(item.get('energy_keV', 0.0)) for item in deposits):.6g} keV",
            f"Charged secondaries: {len(tracks)} tracks; escaped charged energy={sum(float(item.get('escaped_energy_keV', 0.0)) for item in tracks):.6g} keV",
            "Expected local light (NOT transported/detected):",
        ]
        for item in lights:
            event_lines.append(
                f"  {item.get('mechanism')}: {float(item.get('expected_photons', 0.0)):.6g} expected photons in {item.get('volume')} "
                f"from {float(item.get('deposited_energy_keV', 0.0)):.6g} keV over {float(item.get('path_mm', 0.0)):.6g} mm"
            )
        if not light_available:
            event_lines.append("  Local light data is unavailable in this saved event.")
        if not gamma_geometry_view_data(geometry_snapshot)["components"]:
            event_lines.append("Detector geometry is unavailable because this saved result has no v72 run geometry snapshot.")
        event_lines.append(f"Escaped gamma energy: {float(selected.get('escaped_gamma_energy_keV', 0.0)):.6g} keV")
        self._set_readonly_text(self.gamma_results_event_text, "\n".join(event_lines))
        draw_gamma_event_view(self.gamma_event_figure, selected, geometry_snapshot)
        self.gamma_event_canvas.draw_idle()
        QuartzCherenkovGUI._refresh_gamma_global_3d_view(self)

    def _open_gamma_global_3d_viewer(self) -> None:
        """Open the independent, interactive global Gamma event window."""

        payload = _get_instance_attr(self, "latest_gamma_payload", None)
        if not isinstance(payload, Mapping) or not isinstance(payload.get("events"), (list, tuple)):
            messagebox.showinfo(
                "Gamma global 3D view", "No completed Gamma result is available in Gamma Results yet.",
                parent=self,
            )
            return
        viewer = _get_instance_attr(self, "gamma_global_3d_viewer", None)
        try:
            if viewer is not None and viewer.winfo_exists():
                viewer.lift()
                viewer.focus_force()
                self._refresh_gamma_global_3d_view()
                return
        except tk.TclError:
            pass

        viewer = tk.Toplevel(self)
        viewer.title("Global Gamma 3D Viewer")
        viewer.geometry("1120x760")
        viewer.minsize(720, 500)
        viewer.columnconfigure(0, weight=1)
        viewer.columnconfigure(1, weight=0)
        viewer.rowconfigure(0, weight=1)
        canvas = tk.Canvas(viewer, background="white", highlightthickness=1, highlightbackground="#909090")
        canvas.grid(row=0, column=0, sticky="nsew", padx=7, pady=7)
        info = ttk.LabelFrame(viewer, text="Global Gamma view", padding=8)
        info.grid(row=0, column=1, sticky="ns", padx=(0, 7), pady=7)
        info.columnconfigure(0, weight=1)
        info.rowconfigure(0, weight=1)
        info_text = tk.Text(
            info, width=34, height=28, wrap="word", relief="flat", borderwidth=0,
            background=_ttk_background(info), takefocus=True,
        )
        info_text.grid(row=0, column=0, sticky="nsew")
        info_scroll = ttk.Scrollbar(info, orient="vertical", command=info_text.yview)
        info_scroll.grid(row=0, column=1, sticky="ns")
        info_text.configure(yscrollcommand=info_scroll.set, state="disabled")
        self.gamma_global_3d_info_text = info_text

        controls = ttk.Frame(viewer)
        controls.grid(row=1, column=0, columnspan=2, sticky="ew", padx=7, pady=(0, 7))
        ttk.Label(controls, text="Event").pack(side="left")
        self.gamma_global_3d_event_selector = ttk.Combobox(
            controls, textvariable=self.gamma_global_3d_event_choice,
            state="readonly", width=9,
        )
        self.gamma_global_3d_event_selector.pack(side="left", padx=(5, 10))
        self.gamma_global_3d_event_selector.bind("<<ComboboxSelected>>", self._gamma_global_3d_event_selected)
        ttk.Button(controls, text="Reset view", command=self._reset_gamma_global_3d_view).pack(side="left")
        ttk.Checkbutton(
            controls, text="Show geometry", variable=self.gamma_global_3d_show_geometry,
            command=self._draw_gamma_global_3d_view,
        ).pack(side="left", padx=(10, 0))
        ttk.Checkbutton(
            controls, text="Show source", variable=self.gamma_global_3d_show_source,
            command=self._draw_gamma_global_3d_view,
        ).pack(side="left", padx=(7, 0))
        ttk.Checkbutton(
            controls, text="Show tracks / deposits", variable=self.gamma_global_3d_show_tracks,
            command=self._draw_gamma_global_3d_view,
        ).pack(side="left", padx=(7, 0))
        ttk.Label(
            controls, text="Drag to rotate · mouse wheel to zoom · global X/Y/Z in mm",
            foreground="#555555",
        ).pack(side="left", padx=(12, 0))
        ttk.Button(controls, text="Close", command=self._close_gamma_global_3d_viewer).pack(side="right")

        self.gamma_global_3d_viewer = viewer
        self.gamma_global_3d_canvas = canvas
        canvas.bind("<Configure>", lambda _event: self._draw_gamma_global_3d_view())
        canvas.bind("<ButtonPress-1>", self._gamma_global_3d_press)
        canvas.bind("<B1-Motion>", self._gamma_global_3d_drag_motion)
        canvas.bind("<MouseWheel>", self._gamma_global_3d_zoom_event)
        canvas.bind("<Button-4>", self._gamma_global_3d_zoom_event)
        canvas.bind("<Button-5>", self._gamma_global_3d_zoom_event)
        viewer.protocol("WM_DELETE_WINDOW", self._close_gamma_global_3d_viewer)
        self._refresh_gamma_global_3d_view()

    def _close_gamma_global_3d_viewer(self) -> None:
        viewer = _get_instance_attr(self, "gamma_global_3d_viewer", None)
        self.gamma_global_3d_viewer = None
        self.gamma_global_3d_canvas = None
        try:
            if viewer is not None and viewer.winfo_exists():
                viewer.destroy()
        except tk.TclError:
            pass

    def _gamma_global_3d_event_selected(self, _event=None) -> None:
        choice = _get_instance_attr(self, "gamma_global_3d_event_choice", None)
        if choice is None:
            return
        try:
            event_id = int(choice.get())
        except (TypeError, ValueError, tk.TclError):
            return
        payload = _get_instance_attr(self, "latest_gamma_payload", None)
        event_ids = [row["event_id"] for row in gamma_event_table_rows(payload.get("events", ()))] if isinstance(payload, Mapping) else []
        if event_id not in event_ids:
            return
        self.gamma_selected_event.set(str(event_id))
        table = _get_instance_attr(self, "gamma_results_event_table", None)
        if table is not None:
            try:
                table.selection_set(str(event_id))
                table.see(str(event_id))
            except tk.TclError:
                pass
        selector = _get_instance_attr(self, "gamma_results_event_selector", None)
        if selector is not None:
            try:
                selector.set(str(event_id))
            except (AttributeError, tk.TclError):
                pass
        self._refresh_gamma_results_view()

    def _refresh_gamma_global_3d_view(self) -> None:
        """Refresh the detached scene while preserving its camera and filters."""

        canvas = _get_instance_attr(self, "gamma_global_3d_canvas", None)
        if canvas is None:
            return
        try:
            if not canvas.winfo_exists():
                return
        except tk.TclError:
            return
        payload = _get_instance_attr(self, "latest_gamma_payload", None)
        event_id_text = _get_instance_attr(self, "gamma_selected_event", None)
        events = list(payload.get("events", ())) if isinstance(payload, Mapping) else []
        rows = gamma_event_table_rows(events)
        event_ids = [int(row["event_id"]) for row in rows]
        selector = _get_instance_attr(self, "gamma_global_3d_event_selector", None)
        if selector is not None:
            selector.configure(values=tuple(str(value) for value in event_ids))
        choice = _get_instance_attr(self, "gamma_global_3d_event_choice", None)
        try:
            selected_id = int(event_id_text.get()) if event_id_text is not None else (event_ids[0] if event_ids else 0)
        except (TypeError, ValueError, tk.TclError):
            selected_id = event_ids[0] if event_ids else 0
        if event_ids and selected_id not in event_ids:
            selected_id = event_ids[0]
        if choice is not None:
            choice.set(str(selected_id))
        event = next((item for item in events if isinstance(item, Mapping) and int(_gamma_number(item.get("event_id", -1), -1)) == selected_id), None)
        summary = payload.get("summary", {}) if isinstance(payload, Mapping) else {}
        self.gamma_global_3d_scene = build_gamma_global_view_scene(
            event, payload.get("geometry_snapshot") if isinstance(payload, Mapping) else None,
            summary.get("source_direction_global", (0.0, 0.0, 1.0)) if isinstance(summary, Mapping) else (0.0, 0.0, 1.0),
        )
        self._draw_gamma_global_3d_view()

    def _reset_gamma_global_3d_view(self) -> None:
        self.gamma_global_3d_yaw = -0.96
        self.gamma_global_3d_pitch = 0.38
        self.gamma_global_3d_zoom = 1.0
        self._draw_gamma_global_3d_view()

    def _gamma_global_3d_press(self, event) -> None:
        self.gamma_global_3d_drag = (event.x, event.y, self.gamma_global_3d_yaw, self.gamma_global_3d_pitch)

    def _gamma_global_3d_drag_motion(self, event) -> None:
        drag = _get_instance_attr(self, "gamma_global_3d_drag", None)
        if drag is None:
            return
        start_x, start_y, start_yaw, start_pitch = drag
        self.gamma_global_3d_yaw = start_yaw + (event.x - start_x) * 0.012
        self.gamma_global_3d_pitch = max(-1.45, min(1.45, start_pitch + (event.y - start_y) * 0.012))
        self._draw_gamma_global_3d_view()

    def _gamma_global_3d_zoom_event(self, event) -> None:
        delta = getattr(event, "delta", 0)
        if delta == 0:
            delta = 120 if getattr(event, "num", 4) == 4 else -120
        self.gamma_global_3d_zoom *= 1.12 if delta > 0 else 1.0 / 1.12
        self.gamma_global_3d_zoom = max(0.2, min(8.0, self.gamma_global_3d_zoom))
        self._draw_gamma_global_3d_view()

    def _draw_gamma_global_3d_view(self) -> None:
        canvas = _get_instance_attr(self, "gamma_global_3d_canvas", None)
        scene = _get_instance_attr(self, "gamma_global_3d_scene", None)
        if canvas is None or scene is None:
            return
        try:
            if not canvas.winfo_exists():
                return
        except tk.TclError:
            return
        canvas.delete("all")
        width, height = max(1, int(canvas.winfo_width())), max(1, int(canvas.winfo_height()))
        selected = _get_instance_attr(self, "gamma_selected_event", None)
        event_text = selected.get() if selected is not None else "?"
        canvas.create_text(12, 10, anchor="nw", text=f"GLOBAL GAMMA 3D · event {event_text} · {scene['status']}", font=("TkDefaultFont", 12, "bold"), fill="#1b2730")
        geometry_visible = bool(self.gamma_global_3d_show_geometry.get())
        tracks_visible = bool(self.gamma_global_3d_show_tracks.get())
        source_visible = bool(self.gamma_global_3d_show_source.get())
        geometry = scene["geometry_components"] if geometry_visible else []
        gamma_segments = scene["gamma_segments"] if tracks_visible else []
        charged_segments = scene["charged_segments"] if tracks_visible else []
        deposits = scene["deposits"] if tracks_visible else []
        source_point = scene["source_position_mm"] if source_visible else None
        source_end = scene["source_ray_end_mm"] if source_visible else None

        visible_points = []
        for component in geometry:
            visible_points.extend(point for face in component["faces_mm"] for point in face)
        for segment in gamma_segments + charged_segments:
            visible_points.extend((segment["start_mm"], segment["end_mm"]))
        visible_points.extend(item["position_mm"] for item in deposits)
        if source_point is not None:
            visible_points.append(source_point)
        if source_end is not None:
            visible_points.append(source_end)
        if not visible_points:
            visible_points = [[0.0, 0.0, 0.0]]
        center = [sum(point[axis] for point in visible_points) / len(visible_points) for axis in range(3)]
        extent = max(max(max(point[axis] for point in visible_points) - min(point[axis] for point in visible_points) for axis in range(3)), 1.0)
        plot_top, plot_bottom = 36.0, max(40.0, height - 30.0)
        plot_height = max(plot_bottom - plot_top, 1.0)
        scale = min(max(width - 40.0, 1.0), plot_height) * 0.72 / extent * self.gamma_global_3d_zoom
        yaw, pitch = self.gamma_global_3d_yaw, self.gamma_global_3d_pitch
        cos_yaw, sin_yaw = math.cos(yaw), math.sin(yaw)
        cos_pitch, sin_pitch = math.cos(pitch), math.sin(pitch)
        center_x, center_y = width / 2.0, (plot_top + plot_bottom) / 2.0

        def project(point):
            x, y, z = (float(point[i]) - center[i] for i in range(3))
            rx, ry = cos_yaw * x - sin_yaw * y, sin_yaw * x + cos_yaw * y
            sy, depth = cos_pitch * ry - sin_pitch * z, sin_pitch * ry + cos_pitch * z
            return center_x + rx * scale, center_y - sy * scale, depth

        polygons = []
        if geometry_visible:
            for component in geometry:
                for face in component["faces_mm"]:
                    projected = [project(point) for point in face]
                    polygons.append((sum(point[2] for point in projected) / len(projected), component, projected))
        for _depth, component, projected in sorted(polygons, key=lambda row: row[0]):
            coords = [coordinate for point in projected for coordinate in point[:2]]
            canvas.create_polygon(
                *coords, fill=component["color"], stipple="gray25",
                outline="#45515a", width=1,
            )
        if tracks_visible:
            for segment in gamma_segments:
                start, end = project(segment["start_mm"]), project(segment["end_mm"])
                transmitted = segment["interaction"] == "transmitted"
                canvas.create_line(*start[:2], *end[:2], fill="#414141" if transmitted else "#111111", width=2, dash=(5, 3) if transmitted else ())
            kind_colors = {"electron": "#147d64", "positron": "#8246af"}
            for segment in charged_segments:
                start, end = project(segment["start_mm"]), project(segment["end_mm"])
                canvas.create_line(*start[:2], *end[:2], fill=kind_colors.get(segment["kind"], "#147d64"), width=2)
            deposit_colors = {"compton": "#e69f00", "photoelectric": "#d62728", "pair": "#8246af"}
            for item in deposits:
                x, y, _depth = project(item["position_mm"])
                radius = 5
                canvas.create_oval(x-radius, y-radius, x+radius, y+radius, fill=deposit_colors.get(item["interaction"], "#e69f00"), outline="white", width=1)
        if source_point is not None:
            x, y, _depth = project(source_point)
            if source_end is not None:
                ex, ey, _ = project(source_end)
                canvas.create_line(x, y, ex, ey, fill="#0072b2", width=2, arrow="last")
            canvas.create_oval(x-6, y-6, x+6, y+6, fill="#0072b2", outline="white", width=2)
            label = f"Source · {scene['source_energy_keV']:.5g} keV · {scene['source_component']}"
            canvas.create_text(x+9, y-8, anchor="sw", text=label, fill="#005482", font=("TkDefaultFont", 9, "bold"))

        # Keep a small global-axis triad visible regardless of zoom or camera
        # orientation so the detached canvas remains legible as an X/Y/Z view.
        axis_origin = (48.0, max(48.0, height - 48.0))
        axis_length = 27.0
        for vector, color, label in (((1.0, 0.0, 0.0), "#d62728", "X"),
                                     ((0.0, 1.0, 0.0), "#2ca02c", "Y"),
                                     ((0.0, 0.0, 1.0), "#1f77b4", "Z")):
            vx = cos_yaw * vector[0] - sin_yaw * vector[1]
            vy = cos_pitch * (sin_yaw * vector[0] + cos_yaw * vector[1]) - sin_pitch * vector[2]
            endpoint = (axis_origin[0] + vx * axis_length, axis_origin[1] - vy * axis_length)
            canvas.create_line(*axis_origin, *endpoint, fill=color, width=2, arrow="last")
            canvas.create_text(endpoint[0] + 4, endpoint[1] - 4, anchor="sw", text=label, fill=color, font=("TkDefaultFont", 9, "bold"))

        if not geometry and geometry_visible:
            canvas.create_text(12, height-12, anchor="sw", text="Saved detector geometry snapshot unavailable for this Gamma result", fill="#8a3a20")
        if not gamma_segments and not charged_segments and not deposits and tracks_visible:
            canvas.create_text(center_x, center_y, text="No Gamma or charged tracks / deposits in this event", fill="#5a6268")
        legend = [
            "Blue arrow: source direction",
            "Black / dashed gray: interacting / transmitted gamma",
            "Green / purple: electron / positron track",
            "Orange / red / purple markers: Compton / photoelectric / pair deposit",
        ]
        if scene["source_position_inferred"] and source_visible:
            legend.append("Source position inferred from saved geometry bounds")
        legend.extend([
            f"Saved geometry components: {len(scene['geometry_components'])}",
            "Camera: drag to rotate · wheel to zoom · global coordinates (mm)",
        ])
        self._set_readonly_text(self.gamma_global_3d_info_text, "\n".join(legend))

    @staticmethod
    def _make_results_figure(parent, row: int, column: int):
        from matplotlib.backends.backend_tkagg import FigureCanvasTkAgg
        from matplotlib.figure import Figure
        figure = Figure(figsize=(6.4, 2.55), dpi=100)
        canvas = FigureCanvasTkAgg(figure, master=parent)
        canvas.get_tk_widget().grid(row=row, column=column, sticky="nsew")
        return figure, canvas

    def _available_result_event_ids(self) -> list[int]:
        if self.latest_results:
            return [int(result.event_id) for result in self.latest_results]
        summary = self.latest_summary or {}
        return [int(item.get("event_id", index)) for index, item in enumerate(summary.get("event_summaries", ())) if isinstance(item, Mapping)]

    def _results_event_selected(self, _event=None) -> None:
        value = self.selected_event.get().strip()
        try:
            event_id = int(value)
        except (TypeError, ValueError):
            return
        self._select_event_id(event_id)

    def _select_event_id(self, event_id: int) -> None:
        event_ids = self._available_result_event_ids()
        if event_ids and int(event_id) not in event_ids:
            event_id = event_ids[0]
        self.selected_event.set(str(event_id))
        if _has_instance_attr(self, "event_tree"):
            for row_index, row in enumerate(_get_instance_attr(self, "_visible_event_rows", ())):
                if int(row.get("event_id", row_index)) == int(event_id):
                    iid = f"event-{row_index}"
                    self.event_tree.selection_set(iid)
                    self.event_tree.focus(iid)
                    self._update_event_details(row_index)
                    break
        self._refresh_results_tab()

    @staticmethod
    def _results_card_text(detector: str, item: Mapping[str, object], *, aggregate: bool = False) -> str:
        if not bool(item.get("included", False)):
            return "Not included in run"
        mean = item.get("mean_arrival_ns")
        median = item.get("median_arrival_ns")
        scope = "all events" if aggregate else "selected event"
        return (
            f"{scope}\n"
            f"Generated: {int(item.get('generated', 0))}\n"
            f"Arrived at readout: {int(item.get('arrivals', 0))}\n"
            f"Detected by sensor: {int(item.get('detected', 0))}\n"
            f"Mean / median arrival: {('n/a' if mean is None else f'{float(mean):.4f}')} / "
            f"{('n/a' if median is None else f'{float(median):.4f}')} ns"
        )

    def _draw_results_figure(self, figure, canvas, analysis: Mapping[str, object], *, aggregate: bool) -> None:
        figure.clear()
        axis = figure.add_subplot(111)
        section = analysis.get("aggregate" if aggregate else "selected", analysis)
        histograms = section.get("histograms", {})
        edges = np.asarray(histograms.get("bin_edges_ns", ()), dtype=float)
        colors = {"cherenkov": "#d95f02", "scintillator": "#1b9e77"}
        labels = {"cherenkov": "Cherenkov", "scintillator": "Scintillator"}
        plotted = 0
        for detector in ("cherenkov", "scintillator"):
            item = section.get("detectors", {}).get(detector, {})
            values = item.get("arrival_times_ns", ())
            if not bool(item.get("included", False)) or not values or len(edges) < 2:
                continue
            histogram = histograms.get("histograms", {}).get(detector, {})
            y_values = histogram.get("normalized_density_per_ns" if aggregate else "counts", ())
            if not y_values:
                continue
            centers = (edges[:-1] + edges[1:]) / 2.0
            axis.step(centers, y_values, where="mid", color=colors[detector], linewidth=2.0, label=labels[detector])
            plotted += 1
        axis.set_xlabel("Arrival time at readout (ns)")
        axis.set_ylabel("Probability density (1/ns)" if aggregate else "Photons / bin")
        axis.set_title("All-event normalized arrival distribution" if aggregate else "Selected-event arrival-time comparison")
        if plotted:
            axis.legend(loc="best", frameon=False)
            axis.grid(True, alpha=0.22)
        else:
            axis.text(0.5, 0.5, "No readout arrival records", transform=axis.transAxes, ha="center", va="center", color="#555555")
            axis.grid(False)
        figure.tight_layout()
        canvas.draw_idle()

    def _refresh_results_tab(self) -> None:
        if not _has_instance_attr(self, "results_status"):
            return
        event_ids = self._available_result_event_ids()
        if self.results_event_selector is not None:
            self.results_event_selector.configure(values=tuple(str(value) for value in event_ids))
        if not event_ids:
            self.results_status.set("No completed run results")
            for variables in (self.results_selected_card_vars, self.results_aggregate_card_vars):
                for variable in variables.values():
                    variable.set("No completed run results")
            empty = {"selected": {"detectors": {}, "histograms": {}}, "aggregate": {"detectors": {}, "histograms": {}}}
            self._draw_results_figure(self.results_selected_figure, self.results_selected_canvas, empty, aggregate=False)
            self._draw_results_figure(self.results_aggregate_figure, self.results_aggregate_canvas, empty, aggregate=True)
            return
        try:
            event_id = int(self.selected_event.get())
        except (TypeError, ValueError):
            event_id = event_ids[0]
        if event_id not in event_ids:
            event_id = event_ids[0]
        self.selected_event.set(str(event_id))
        if self.latest_results:
            analysis = simulator.readout_arrival_analysis(self.latest_results, event_id=event_id)
        elif self.latest_readout_payload is not None:
            analysis = simulator.readout_arrival_analysis_from_rows(
                self.latest_readout_payload.get("rows", ()),
                (self.latest_summary or {}).get("event_summaries", ()),
                self.latest_readout_payload.get("detector_mode", (self.latest_summary or {}).get("detector_mode", "cherenkov")),
                event_id=event_id,
            )
        else:
            self.results_status.set("No readout arrival artifact is available for this saved run")
            return
        self.latest_arrival_analysis = analysis
        self.results_status.set(
            f"Event {event_id}; arrival counts are finite active-window crossings before QE/collection response. "
            "Aggregate curves are normalized to unit area (probability density per ns)."
        )
        for detector, variable in self.results_selected_card_vars.items():
            variable.set(self._results_card_text(detector, analysis["selected"]["detectors"].get(detector, {})))
        for detector, variable in self.results_aggregate_card_vars.items():
            variable.set(self._results_card_text(detector, analysis["aggregate"]["detectors"].get(detector, {}), aggregate=True))
        self._draw_results_figure(self.results_selected_figure, self.results_selected_canvas, analysis, aggregate=False)
        self._draw_results_figure(self.results_aggregate_figure, self.results_aggregate_canvas, analysis, aggregate=True)

    def _build_particle_tab(self, parent: ttk.Frame) -> None:
        parent.columnconfigure(1, weight=1)
        self._build_grid(parent, [
            ("Particle type", "particle_type"), ("Mass (MeV/c2)", "mass_mev_c2"),
            ("Kinetic energy (MeV)", "kinetic_mev"), ("Entry x (mm)", "entry_x_mm"),
            ("Entry y (mm)", "entry_y_mm"), ("Inclination from +Z / track tilt (deg)", "inclination_deg"),
            ("Azimuth around +Z (deg, compatibility)", "azimuth_deg"),
        ])
        ttk.Label(parent, text="Run angle convention: zero means global +Z; positive tilts toward +X in the X-Z plane. Entry X/Y remain available for the Cherenkov assembly; the canonical All path derives the second assembly impact from the shared global line and reports a clear preflight miss when needed.", foreground="#245a35", wraplength=700, justify="left").grid(row=8, column=0, columnspan=2, sticky="w", pady=(12, 0))

    def _build_detector_tab(self, parent: ttk.Frame) -> None:
        # Sensor/readout characteristics are owned by the detector-specific
        # Cherenkov and Scintillator Sensor tabs.  Run mirrors two separate
        # read-only reports and never creates a generic mixed detector panel.
        parent.columnconfigure(0, weight=1)
        ttk.Label(parent, text="These Run views are read-only mirrors. Edit geometry/readout only in the matching Cherenkov or Scintillator tabs; the two physical generation models are never combined.", wraplength=980, justify="left", foreground="#245a35").grid(row=0, column=0, sticky="w", pady=(0, 7))
        self.detector_mode_text = tk.StringVar(value="Detector-specific transport diagnostics ready")
        self.detector_overview_text = tk.StringVar(value="")
        self.detector_mirror_text = tk.StringVar(value="")
        for row, (assembly, title) in enumerate((("cherenkov", "Cherenkov Transport Diagnostic"), ("scintillator", "Scintillator Transport Diagnostic")), start=1):
            frame = ttk.LabelFrame(parent, text=title, padding=7)
            frame.grid(row=row, column=0, sticky="ew", pady=(0, 7))
            frame.columnconfigure(0, weight=1)
            text = tk.StringVar(value="Preparing…")
            setattr(self, f"run_{assembly}_transport_diagnostic_text", text)
            ttk.Label(frame, textvariable=text, justify="left", wraplength=980).grid(row=0, column=0, sticky="ew")
            ttk.Button(frame, text=f"Refresh {assembly.title()} Diagnostic", command=lambda a=assembly: self._refresh_transport_diagnostic(a)).grid(row=1, column=0, sticky="w", pady=(5, 0))
        self._refresh_run_transport_diagnostics()

    def _refresh_run_transport_diagnostics(self) -> None:
        for assembly in ("cherenkov", "scintillator"):
            self._refresh_transport_diagnostic(assembly)
            source = _get_instance_attr(self, f"{assembly}_transport_diagnostic_text", None)
            target = _get_instance_attr(self, f"run_{assembly}_transport_diagnostic_text", None)
            if source is not None and target is not None:
                target.set(source.get())

    def _build_active_area_controls(self, parent: ttk.Frame, prefix: str, start_row: int) -> None:
        """Build shape-specific, editable active-area fields for one detector."""

        if not _has_instance_attr(self, "_detector_area_widgets"):
            self._detector_area_widgets = {}
        fields = (
            ("Shape", "shape", "shape"),
            ("Active width Y (mm)", "active_width_mm", "rectangular"),
            ("Active height Z (mm)", "active_height_mm", "rectangular"),
            ("Active radius (mm)", "active_radius_mm", "circular"),
            ("Semi-axis A / Y (mm)", "active_semi_axis_a_mm", "elliptical"),
            ("Semi-axis B / Z (mm)", "active_semi_axis_b_mm", "elliptical"),
            ("Offset Y from portal centre (mm)", "active_offset_y_mm", "always"),
            ("Offset Z from portal centre (mm)", "active_offset_z_mm", "always"),
            ("Rotation in Y/Z plane (deg)", "active_rotation_deg", "always"),
        )
        for offset, (label, suffix, kind) in enumerate(fields):
            key = f"{prefix}_{suffix}"
            row = start_row + offset
            ttk.Label(parent, text=label).grid(row=row, column=0, sticky="w", padx=(0, 8), pady=2)
            if kind == "shape":
                widget = ttk.Combobox(parent, textvariable=self.values[key], values=DETECTOR_ACTIVE_AREA_SHAPES, state="readonly")
            else:
                widget = ttk.Entry(parent, textvariable=self.values[key])
            widget.grid(row=row, column=1, sticky="ew", pady=2)
            self._detector_area_widgets[key] = widget
            self.values[key].trace_add("write", self._detector_area_changed)
        ttk.Label(parent, text="The active area is checked against the actual terminal Y/Z portal before transport; offsets and rotation are relative to its centre.", foreground="#245a35", wraplength=600, justify="left").grid(row=start_row + len(fields), column=0, columnspan=2, sticky="w", pady=(4, 2))

    def _update_detector_area_field_states(self) -> None:
        if not _has_instance_attr(self, "_detector_area_widgets"):
            return
        mode = _value_text(self.values, "opposite_end_mode", "open").lower()
        for prefix in ("detector", "opposite_detector"):
            shape = _value_text(self.values, f"{prefix}_shape", "rectangular").lower()
            for suffix, enabled_shapes in {
                "active_width_mm": {"rectangular"},
                "active_height_mm": {"rectangular"},
                "active_radius_mm": {"circular"},
                "active_semi_axis_a_mm": {"elliptical"},
                "active_semi_axis_b_mm": {"elliptical"},
            }.items():
                widget = self._detector_area_widgets.get(f"{prefix}_{suffix}")
                if widget is not None:
                    enabled = prefix == "detector" or mode == "detector"
                    widget.configure(state="normal" if enabled and shape in enabled_shapes else "disabled")
            shape_widget = self._detector_area_widgets.get(f"{prefix}_shape")
            if shape_widget is not None:
                shape_widget.configure(state="readonly" if prefix == "detector" or mode == "detector" else "disabled")
            for suffix in ("active_offset_y_mm", "active_offset_z_mm", "active_rotation_deg"):
                widget = self._detector_area_widgets.get(f"{prefix}_{suffix}")
                if widget is not None:
                    widget.configure(state="normal" if prefix == "detector" or mode == "detector" else "disabled")

    def _detector_area_changed(self, *_args) -> None:
        self._update_detector_area_field_states()
        self._update_detector_active_area()

    def _build_geometry_tab(self, parent: ttk.Frame) -> None:
        parent.columnconfigure(0, weight=0)
        parent.columnconfigure(1, weight=1)
        parent.rowconfigure(0, weight=1)
        sections_frame = ttk.LabelFrame(parent, text="Geometry sections", padding=4, width=250)
        sections_frame.grid(row=0, column=0, sticky="ns", padx=(0, 5))
        sections_frame.rowconfigure(0, weight=1)
        sections_frame.columnconfigure(0, weight=1)
        self.geometry_tree = ttk.Treeview(sections_frame, columns=("name", "shape", "material", "dimensions"), show="headings", selectmode="browse", height=12)
        for column, heading, width in (("name", "Name", 100), ("shape", "Shape", 65), ("material", "Material", 90), ("dimensions", "X × Y × Z (mm)", 115)):
            self.geometry_tree.heading(column, text=heading)
            self.geometry_tree.column(column, width=width, anchor="w", stretch=False)
        self.geometry_tree.grid(row=0, column=0, sticky="nsew")
        section_scroll = ttk.Scrollbar(sections_frame, orient="vertical", command=self.geometry_tree.yview)
        section_scroll.grid(row=0, column=1, sticky="ns")
        self.geometry_tree.configure(yscrollcommand=section_scroll.set)
        self.geometry_tree.bind("<<TreeviewSelect>>", self._geometry_selection_changed)
        section_buttons = ttk.Frame(sections_frame)
        section_buttons.grid(row=1, column=0, columnspan=2, sticky="ew", pady=(5, 0))
        ttk.Button(section_buttons, text="Add", command=self._add_geometry_section).pack(side="left")
        ttk.Button(section_buttons, text="Remove", command=self._remove_geometry_section).pack(side="left", padx=4)
        ttk.Button(section_buttons, text="Save Geometry", command=self._save_geometry).pack(side="left", padx=4)
        ttk.Button(section_buttons, text="Load Geometry", command=self._load_geometry).pack(side="left")

        workspace = ttk.Frame(parent)
        workspace.grid(row=0, column=1, sticky="nsew")
        workspace.columnconfigure(0, weight=1)
        workspace.rowconfigure(2, weight=1)
        self.geometry_support_text = tk.StringVar()
        ttk.Label(workspace, textvariable=self.geometry_support_text, foreground="#8a1c1c", wraplength=760, justify="left").grid(row=0, column=0, sticky="ew", pady=(0, 4))
        properties = ttk.LabelFrame(workspace, text="Selected section properties", padding=4)
        properties.grid(row=1, column=0, sticky="ew", pady=(0, 5))
        properties.columnconfigure(0, weight=1)
        properties.rowconfigure(0, weight=1)
        self.geometry_properties_canvas = tk.Canvas(properties, height=176, highlightthickness=0)
        self.geometry_properties_canvas.grid(row=0, column=0, sticky="nsew")
        property_scroll = ttk.Scrollbar(properties, orient="vertical", command=self.geometry_properties_canvas.yview)
        property_scroll.grid(row=0, column=1, sticky="ns")
        self.geometry_properties_canvas.configure(yscrollcommand=property_scroll.set)
        self.geometry_properties_inner = ttk.Frame(self.geometry_properties_canvas)
        self.geometry_properties_window = self.geometry_properties_canvas.create_window((0, 0), window=self.geometry_properties_inner, anchor="nw")
        self.geometry_properties_inner.bind("<Configure>", lambda _event: self.geometry_properties_canvas.configure(scrollregion=self.geometry_properties_canvas.bbox("all")))
        self.geometry_properties_canvas.bind("<Configure>", lambda event: self.geometry_properties_canvas.itemconfigure(self.geometry_properties_window, width=event.width))
        self._build_geometry_properties()

        preview_frame = ttk.LabelFrame(workspace, text="Synchronized geometry previews", padding=4)
        preview_frame.grid(row=2, column=0, sticky="nsew")
        preview_frame.columnconfigure(0, weight=1)
        preview_frame.columnconfigure(1, weight=1)
        preview_frame.rowconfigure(0, weight=1)
        preview_frame.rowconfigure(1, weight=1)
        for index, (key, title) in enumerate((("xy", "XY  (x / y)"), ("xz", "XZ  (x / z)"), ("yz", "YZ  (y / z)"), ("iso", "3D  isometric"))):
            row, column = divmod(index, 2)
            pane = ttk.Frame(preview_frame)
            pane.grid(row=row, column=column, sticky="nsew", padx=2, pady=2)
            pane.columnconfigure(0, weight=1)
            pane.rowconfigure(1, weight=1)
            ttk.Label(pane, text=title).grid(row=0, column=0, sticky="w")
            canvas = tk.Canvas(pane, background="white", highlightthickness=1, highlightbackground="#b0b0b0", width=300, height=190)
            canvas.grid(row=1, column=0, sticky="nsew")
            canvas.bind("<Configure>", lambda _event, view=key: self._refresh_geometry_view(view))
            self.geometry_canvases[key] = canvas
        ttk.Button(workspace, text="Open 3D Viewer", command=self._open_geometry_3d_viewer).grid(row=3, column=0, sticky="e", pady=(5, 0))
        self._refresh_geometry_tree()

    def _build_grid(self, parent: ttk.Frame, fields, start_row: int = 0) -> None:
        parent.columnconfigure(1, weight=1)
        for offset, (label, key) in enumerate(fields):
            row = start_row + offset
            label_widget = ttk.Label(parent, text=label)
            label_widget.grid(row=row, column=0, sticky="w", padx=(0, 8), pady=3)
            if key == "sampling_mode":
                widget = ttk.Combobox(parent, textvariable=self.values[key], values=list(SAMPLING_MODE_LABELS.values()), state="readonly")
            elif key == "detector_run_mode":
                # v62 optical labels remain unchanged; Gamma is an explicit
                # additional GUI mode with its own worker/result branch.
                widget = ttk.Combobox(parent, textvariable=self.values[key], values=RUN_MODE_COMBO_VALUES, state="readonly")
            elif key in {"particle_type", "opposite_end_mode"}:
                values = list(PARTICLE_MASSES_MEV_C2) + ["custom"] if key == "particle_type" else ["open", "mirror", "detector"]
                widget = ttk.Combobox(parent, textvariable=self.values[key], values=values, state="readonly")
            else:
                widget = ttk.Entry(parent, textvariable=self.values[key])
            widget.grid(row=row, column=1, sticky="ew", pady=3)
            self._field_widgets[key] = widget
            self._field_labels[key] = label_widget
            if key == "photons_per_event":
                self.photon_count_label = label_widget
                self.photon_count_entry = widget

    def _particle_type_changed(self, *_args) -> None:
        particle_type = self.values["particle_type"].get()
        if particle_type in PARTICLE_MASSES_MEV_C2:
            self.values["mass_mev_c2"].set(str(PARTICLE_MASSES_MEV_C2[particle_type]))

    def _mode_changed(self, *_args) -> None:
        mode = self.values["opposite_end_mode"].get().strip().lower()
        self._update_detector_area_field_states()
        mirror_state = "normal" if mode == "mirror" else "disabled"
        detector_state = "normal" if mode == "detector" else "disabled"
        for key in ("opposite_detector_index", "opposite_detector_acceptance_deg", "opposite_detector_qe", "opposite_collection_efficiency"):
            widget = self._field_widgets.get(key)
            if widget is not None:
                widget.configure(state=detector_state)
        mirror_widget = self._field_widgets.get("mirror_reflectivity")
        if mirror_widget is not None:
            mirror_widget.configure(state=mirror_state)
        mirror_frame = _get_instance_attr(self, "detector_mirror_frame", None)
        second_frame = _get_instance_attr(self, "detector_second_frame", None)
        if mirror_frame is not None and second_frame is not None:
            mirror_frame.grid_remove()
            second_frame.grid_remove()
            if mode == "mirror":
                mirror_frame.grid(row=0, column=0, sticky="ew")
            elif mode == "detector":
                second_frame.grid(row=0, column=0, sticky="ew")
        if mode == "mirror":
            self.detector_mode_text.set("-x mirror is active; mirror reflectivity controls the -x portal.")
            self.detector_overview_text.set("Active mode: MIRROR\nThe -x portal is reflective. The mirror configuration below is used by the simulator.")
            self.detector_mirror_text.set("Specular reflection at the -x portal; photons can return into the active solid when reflected.")
        elif mode == "detector":
            self.detector_mode_text.set("-x second detector is active; blank fields inherit the corresponding +x values.")
            self.detector_overview_text.set("Active mode: SECOND DETECTOR (-x)\nThe second detector panel below supplies optional overrides for the -x portal.")
        else:
            self.detector_mode_text.set("-x is open; no mirror or second-detector controls are active.")
            self.detector_overview_text.set("Active mode: OPEN\nPhotons leaving the -x portal are handled as an external open interface.")
        if _has_instance_attr(self, "detector_mirror_canvas"):
            self._draw_mirror_portal()
        if _has_instance_attr(self, "detector_minus_face_canvas"):
            self._update_detector_active_area()
        if mode == "mirror":
            self.status.set("-x mirror mode selected")
        elif mode == "detector":
            self.status.set("-x second detector mode selected")
        elif self.worker is None or not self.worker.is_alive():
            self.status.set("Ready")

    def _sampling_mode_changed(self, *_args) -> None:
        value = self.values["sampling_mode"].get()
        mode = SAMPLING_MODE_BY_LABEL.get(value, value.lower())
        if _has_instance_attr(self, "photon_count_entry"):
            self.photon_count_entry.configure(state="disabled" if mode == "poisson" else "normal")
            self.photon_count_label.configure(text="Fixed test sample photons (ignored in Poisson mode)" if mode == "poisson" else "Photons per event")
        if mode == "poisson":
            self.status.set("Poisson mode: fixed sample field ignored")
        elif self.worker is None or not self.worker.is_alive():
            self.status.set("Ready")

    def _apply_coupled_preset(self) -> None:
        values = {
            "events": "100", "photons_per_event": "1000", "sampling_mode": SAMPLING_MODE_LABELS["fixed"],
            "detector_index": "1.40", "detector_acceptance_deg": "90", "detector_qe": "1.0", "collection_efficiency": "1.0",
            "sensor_1_index": "1.40", "sensor_1_acceptance_deg": "90", "sensor_1_qe": "1.0", "sensor_1_collection": "1.0",
        }
        for key, value in values.items():
            self.values[key].set(value)
        self.status.set("Coupled detector quick test loaded; -x mode unchanged")
        self._append_log("Loaded coupled detector quick test: index 1.40, acceptance 90 deg, QE 1, collection 1")

    def _save_detector_profile(self) -> None:
        path = filedialog.asksaveasfilename(
            title="Save detector profile",
            defaultextension=".json",
            filetypes=(("Detector JSON", "*.json"), ("All files", "*.*")),
        )
        if not path:
            return
        try:
            save_detector_profile(path, self.values)
            self.detector_profile_text.set(f"Current detector profile: {Path(path).name}")
            self.status.set(f"Detector profile saved: {Path(path).name}")
            self._append_log(f"Saved detector profile: {Path(path).resolve()}")
        except (OSError, TypeError, ValueError) as exc:
            messagebox.showerror("Save detector profile", str(exc), parent=self)

    def _load_detector_profile(self) -> None:
        path = filedialog.askopenfilename(
            title="Load detector profile",
            filetypes=(("Detector JSON", "*.json"), ("All files", "*.*")),
        )
        if not path:
            return
        try:
            load_detector_profile(path, self.values)
            self._mode_changed()
            self.detector_profile_text.set(f"Current detector profile: {Path(path).name}")
            self.status.set(f"Detector profile loaded: {Path(path).name}")
            self._append_log(f"Loaded detector profile: {Path(path).resolve()}")
        except (OSError, json.JSONDecodeError, TypeError, ValueError) as exc:
            # load_detector_profile validates before applying, so this branch
            # leaves all current detector StringVars untouched.
            messagebox.showerror("Load detector profile", str(exc), parent=self)

    def _save_sensor_profile(self) -> None:
        path = filedialog.asksaveasfilename(title="Save beam-monitor sensors", defaultextension=".json", filetypes=(("Sensors JSON", "*.json"), ("All files", "*.*")))
        if not path:
            return
        try:
            save_component_profile(path, sensor_profile_from_values(self.values, self.assembly_models))
            self.status.set(f"Sensors saved: {Path(path).name}")
        except (OSError, TypeError, ValueError) as exc:
            messagebox.showerror("Save sensors", str(exc), parent=self)

    def _load_sensor_profile(self) -> None:
        path = filedialog.askopenfilename(title="Load beam-monitor sensors", filetypes=(("Sensors JSON", "*.json"), ("All files", "*.*")))
        if not path:
            return
        try:
            profile = json.loads(Path(path).expanduser().resolve().read_text(encoding="utf-8"))
            apply_sensor_profile(self.values, profile)
            self.status.set(f"Sensors loaded: {Path(path).name}")
        except (OSError, json.JSONDecodeError, TypeError, ValueError) as exc:
            messagebox.showerror("Load sensors", str(exc), parent=self)

    def _browse_output(self) -> None:
        selected = filedialog.askdirectory(initialdir=self.values["output_directory"].get() or str(Path.cwd()))
        if selected:
            self.values["output_directory"].set(selected)

    def _browse_gamma_profile(self) -> None:
        selected = filedialog.askopenfilename(
            title="Select approved gamma source profile",
            initialdir=str(Path(self.values["gamma_source_profile_path"].get()).expanduser().parent),
            filetypes=(("Gamma source profile", "*.json"), ("All files", "*.*")),
        )
        if selected:
            self.values["gamma_source_profile_path"].set(selected)

    def _append_log(self, message: str) -> None:
        self.log_text.configure(state="normal")
        self.log_text.insert("end", message.rstrip() + "\n")
        self.log_text.see("end")
        self.log_text.configure(state="disabled")

    def _show_run_phase(self, phase: str) -> None:
        """Make a pre-worker phase visible before continuing synchronous setup."""

        self.run_state = "validating" if phase == "Validating" else "running"
        self.status.set(phase)
        self.progress_text.set(f"{phase} · 0.0% complete, elapsed 0.0 s")
        self.update_idletasks()

    def _return_to_idle_after_run_rejection(self) -> None:
        self.run_state = "idle"
        self.progress.configure(value=0.0)
        self.status.set("Ready")
        self.progress_text.set("Ready")

    def _run(self) -> None:
        if self.worker is not None and self.worker.is_alive():
            return
        self.progress_started_at = time.monotonic()
        self.progress.configure(value=0.0)
        self._show_run_phase("Preparing run")
        geometry_error = self._validate_geometry_for_run()
        if geometry_error:
            messagebox.showerror("Geometry not supported by transport", geometry_error, parent=self)
            self._append_log("ERROR: " + geometry_error)
            self._return_to_idle_after_run_rejection()
            return
        try:
            simulation_geometry = self._geometry_for_simulation()
        except (TypeError, ValueError) as exc:
            messagebox.showerror("Invalid simulation geometry", str(exc), parent=self)
            self._append_log(f"ERROR: {exc}")
            self._return_to_idle_after_run_rejection()
            return
        raw = {key: variable.get() for key, variable in self.values.items()}
        try:
            settings = validate_settings(raw)
            # Current runs use the complete canonical assembly validator.  The
            # simulator's legacy active-area validator remains available for
            # direct legacy BarGeometry/SolidGeometry callers, but must not be
            # a second authority for the editor's separate photosensor model.
            validate_canonical_run_assemblies(self.values, self.assembly_models, settings.detector_run_mode)
            run_directory = create_run_directory(settings.output_directory)
            # Serialize/deep-copy every GUI-owned input before creating the
            # worker.  The worker will reconstruct private models from these
            # data-only payloads and never read StringVars, widgets, or the
            # live assembly_models mapping.
            snapshot = _make_run_snapshot(settings, simulation_geometry, self.assembly_models)
        except Exception as exc:
            messagebox.showerror("Invalid settings", str(exc), parent=self)
            self._append_log(f"ERROR: {exc}")
            self._return_to_idle_after_run_rejection()
            return
        self.last_run_directory = run_directory
        mode_label = detector_run_mode_label(settings.detector_run_mode)
        self.stop_requested.clear()
        self.progress_started_at = time.monotonic()
        self.progress_total_events = settings.events
        self.progress_event_index = 0
        self.progress_event_fraction = 0.0
        self.progress.configure(value=0.0)
        self.progress_text.set(f"Starting {mode_label}: event 1/{settings.events} (0.0%), elapsed 0.0 s, ETA n/a")
        self.run_button.configure(state="disabled")
        self.stop_button.configure(state="normal")
        self.run_state = "running"
        self._show_run_phase(f"Starting worker · {mode_label}")
        self.progress_text.set(f"Starting worker · {mode_label} · 0.0% complete, elapsed 0.0 s")
        self._append_log(f"Output directory: {run_directory}")
        self._append_log(f"Detector mode: {mode_label} ({settings.detector_run_mode})")
        self.active_detector_mode_label = mode_label
        work_queue = self.queue
        stop_event = self.stop_requested
        self._worker_terminal_seen = False
        self.worker = threading.Thread(
            target=self._worker_run,
            args=(snapshot, run_directory, work_queue, stop_event),
            daemon=True,
            name="quartz-simulation-worker",
        )
        try:
            self.worker.start()
        except Exception as exc:
            self._finish_error(f"Could not start simulation worker: {type(exc).__name__}: {exc}")

    def _validate_geometry_for_run(self) -> Optional[str]:
        mode = detector_run_mode_from_values(self.values)
        selected = {"cherenkov", "scintillator"} if mode in {"combined", GAMMA_RUN_MODE} else {mode}
        if "cherenkov" in selected:
            model = self.assembly_models.get("cherenkov")
            if not isinstance(model, AssemblyGeometryModel) or not model.sections:
                return "Run requires at least one Cherenkov geometry section."
            try:
                model.guide_enabled = _bool_text(self.values, "cherenkov_guide_enabled")
                model.mylar_enabled = _bool_text(self.values, "cherenkov_mylar_enabled")
                model.mylar_reflectivity = float(_value_text(self.values, "cherenkov_mylar_reflectivity_percent", "95")) / 100.0
                model.to_simulation_geometry().validate()
            except (TypeError, ValueError) as exc:
                return f"Run blocked: Cherenkov assembly geometry is invalid: {exc}"
        if "scintillator" in selected:
            model = self.assembly_models.get("scintillator")
            if not isinstance(model, AssemblyGeometryModel) or not model.sections:
                return "Run requires at least one Scintillator geometry section."
            try:
                sensor = complete_sensor_readout_from_values(self.values, "scintillator", model)
                model.validate_complete_assembly(
                    sensor_shape=str(sensor["shape"]),
                    sensor_width_mm=float(sensor["dimensions_mm"]["width"]),
                    sensor_height_mm=float(sensor["dimensions_mm"]["height"]),
                    sensor_depth_mm=float(sensor["dimensions_mm"]["depth"]),
                    sensor_radius_mm=float(sensor["readout_radius_mm"]),
                    sensor_semi_axis_a_mm=float(sensor["readout_semiaxes_mm"]["a"]),
                    sensor_semi_axis_b_mm=float(sensor["readout_semiaxes_mm"]["b"]),
                )
            except (TypeError, ValueError) as exc:
                return f"Run blocked: Scintillator assembly geometry is invalid: {exc}"
        return None

    def _geometry_for_simulation(self):
        if detector_run_mode_from_values(self.values) == "scintillator":
            # The simulator's scintillator-only branch owns its optical
            # transport geometry.  A neutral valid geometry keeps the shared
            # snapshot/result contract without validating or simulating a
            # hidden Cherenkov detector.
            return simulator.BarGeometry()
        model = self._assembly_model("cherenkov")
        model.guide_enabled = _bool_text(self.values, "cherenkov_guide_enabled")
        model.mylar_enabled = _bool_text(self.values, "cherenkov_mylar_enabled")
        model.mylar_reflectivity = float(_value_text(self.values, "cherenkov_mylar_reflectivity_percent", "95")) / 100.0
        model.terminal_plus_surface = _value_text(self.values, "terminal_plus_surface", model.terminal_plus_surface)
        model.terminal_minus_surface = _value_text(self.values, "terminal_minus_surface", model.terminal_minus_surface)
        sensor = complete_sensor_readout_from_values(self.values, "cherenkov", model)
        geometry = model.to_simulation_geometry(
            sensor_shape=str(sensor["shape"]),
            sensor_width_mm=float(sensor["dimensions_mm"]["width"]),
            sensor_height_mm=float(sensor["dimensions_mm"]["height"]),
            sensor_depth_mm=float(sensor["dimensions_mm"]["depth"]),
            sensor_radius_mm=float(sensor["readout_radius_mm"]),
            sensor_semi_axis_a_mm=float(sensor["readout_semiaxes_mm"]["a"]),
            sensor_semi_axis_b_mm=float(sensor["readout_semiaxes_mm"]["b"]),
        )
        geometry.validate()
        return geometry

    def _stop(self) -> None:
        if self.worker is None or not self.worker.is_alive():
            return
        self.stop_requested.set()
        self.stop_button.configure(state="disabled")
        self.status.set("Stop requested; finishing the current photon/event safely")
        self.progress_text.set("Stop requested; finishing the current photon/event safely")
        self._append_log("Stop requested. The simulator will return completed/partial results cooperatively.")

    def _worker_run(
        self,
        snapshot: RunSnapshot | RunSettings,
        run_directory: Path,
        work_queue: Optional[queue.Queue] = None,
        stop_event: Optional[threading.Event] = None,
    ) -> None:
        """Run from detached data and publish only queue messages.

        This method intentionally has no access path to Tk variables/widgets
        or the live assembly model mapping.  Even the queue and stop event are
        passed as thread-safe objects rather than looked up through ``self``.
        """

        # Keep direct/headless callers of the pre-v41 private method working;
        # the GUI path above always supplies a RunSnapshot, queue, and event.
        if not isinstance(snapshot, RunSnapshot):
            legacy_geometry = work_queue
            work_queue = _get_instance_attr(self, "queue")
            stop_event = stop_event or _get_instance_attr(self, "stop_requested", threading.Event())
            snapshot = _make_run_snapshot(
                snapshot,
                legacy_geometry,
                _get_instance_attr(self, "assembly_models", {}),
            )
        if work_queue is None:
            raise TypeError("worker queue is required")
        if stop_event is None:
            stop_event = threading.Event()

        terminal_lock = threading.Lock()
        terminal_sent = False
        heartbeat_stop = threading.Event()
        phase_lock = threading.Lock()
        phase_state = {"name": "Worker started", "started_at": time.monotonic()}

        def enqueue(kind: str, payload) -> bool:
            try:
                work_queue.put((kind, payload))
                return True
            except Exception as exc:
                _append_persistent_run_log(
                    run_directory / "run.log",
                    f"GUI queue {kind} callback failed: {type(exc).__name__}: {exc}",
                )
                return False

        def log(message: str) -> None:
            enqueue("log", str(message))

        def status(message: str) -> None:
            nonlocal last_status_emit, last_status_message
            message = str(message)
            with phase_lock:
                phase_state["name"] = message
                phase_state["started_at"] = time.monotonic()
            now = time.monotonic()
            if message == last_status_message and now - last_status_emit < 0.08:
                return
            last_status_emit = now
            last_status_message = message
            enqueue("status", message)

        def progress(event_index: int, total_events: int, event_fraction: float) -> None:
            nonlocal last_progress_emit
            now = time.monotonic()
            force = event_fraction >= 1.0
            if not force and now - last_progress_emit < 0.08:
                return
            last_progress_emit = now
            enqueue("progress", (event_index, total_events, event_fraction))

        def heartbeat() -> None:
            while not heartbeat_stop.wait(3.0):
                with phase_lock:
                    name = str(phase_state["name"])
                    elapsed = max(0.0, time.monotonic() - float(phase_state["started_at"]))
                enqueue("status", f"{name} (active {elapsed:.1f} s)")

        def send_terminal(kind: str, payload) -> bool:
            nonlocal terminal_sent
            with terminal_lock:
                if terminal_sent:
                    return False
                terminal_sent = True
            heartbeat_stop.set()
            return enqueue(kind, payload)

        last_progress_emit = 0.0
        last_status_emit = 0.0
        last_status_message = ""
        heartbeat_thread = None
        try:
            # This is the first worker action, before snapshot reconstruction
            # or any geometry preparation.  It proves that Thread.start()
            # succeeded and gives the Tk queue a distinct phase to display.
            _append_persistent_run_log(run_directory / "run.log", "PHASE: Worker started")
            enqueue("status", "Worker started")
            heartbeat_thread = threading.Thread(
                target=heartbeat,
                name="quartz-simulation-heartbeat",
                daemon=True,
            )
            heartbeat_thread.start()
            settings = _settings_from_snapshot(snapshot)
            assembly_models = _assembly_models_from_snapshot(snapshot)
            geometry = _transport_geometry_from_snapshot(snapshot.transport_geometry)
            payload = run_simulation_to_directory(
                settings,
                run_directory,
                log=log,
                geometry=geometry,
                progress_callback=progress,
                status_callback=status,
                stop_event=stop_event,
                assembly_models=assembly_models,
            )
            send_terminal("done", payload)
        except Exception as exc:
            message = f"{type(exc).__name__}: {exc}"
            _append_persistent_run_log(run_directory / "run.log", f"ERROR: {message}")
            send_terminal("error", message)
        finally:
            heartbeat_stop.set()
            # This is a last-resort guard for failures while reconstructing
            # inputs or publishing the normal terminal payload.  The lock in
            # send_terminal guarantees exactly one terminal queue message.
            with terminal_lock:
                already_sent = terminal_sent
            if not already_sent:
                message = "Worker terminated before reporting a terminal state"
                _append_persistent_run_log(run_directory / "run.log", f"ERROR: {message}")
                send_terminal("error", message)

    def _update_progress(self, payload) -> None:
        try:
            event_index, total_events, event_fraction = payload
            event_index = int(event_index)
            total_events = max(1, int(total_events))
            event_fraction = max(0.0, min(1.0, float(event_fraction)))
        except (TypeError, ValueError):
            return
        self.progress_event_index = event_index
        self.progress_total_events = total_events
        self.progress_event_fraction = event_fraction
        overall_fraction = max(0.0, min(1.0, (event_index + event_fraction) / total_events))
        elapsed = max(0.0, time.monotonic() - self.progress_started_at) if self.progress_started_at else 0.0
        eta = elapsed * (1.0 - overall_fraction) / overall_fraction if overall_fraction > 1e-9 else None
        eta_text = f"{eta:.1f} s" if eta is not None else "n/a"
        mode_label = _get_instance_attr(self, "active_detector_mode_label", "selected detector mode")
        text = f"Running {mode_label}: event {min(event_index + 1, total_events)}/{total_events} ({overall_fraction * 100.0:.1f}%), elapsed {elapsed:.1f} s, ETA {eta_text}"
        self.progress.configure(value=overall_fraction * 100.0)
        self.progress_text.set(text)
        self.status.set(text)

    def _update_status(self, phase: str) -> None:
        """Display worker phase/heartbeat while preserving honest percentage."""

        phase = str(phase).strip() or "Running"
        elapsed = max(0.0, time.monotonic() - self.progress_started_at) if self.progress_started_at else 0.0
        overall_fraction = max(0.0, min(1.0, (self.progress_event_index + self.progress_event_fraction) / max(self.progress_total_events, 1)))
        mode_label = _get_instance_attr(self, "active_detector_mode_label", "selected detector mode")
        text = f"{phase} · {mode_label} · {overall_fraction * 100.0:.1f}% complete, elapsed {elapsed:.1f} s"
        self.progress_text.set(text)
        self.status.set(text)

    def _poll_queue(self) -> None:
        self._queue_poll_scheduled = False
        try:
            while True:
                kind, payload = self.queue.get_nowait()
                try:
                    if kind == "log":
                        self._append_log(str(payload))
                    elif kind == "progress":
                        self._update_progress(payload)
                    elif kind == "status":
                        self._update_status(str(payload))
                    elif kind == "done":
                        self._worker_terminal_seen = True
                        try:
                            self._finish_run(payload)
                        except Exception as exc:
                            # A malformed result or a GUI-only finalization
                            # problem must still return the Run button to idle.
                            self._finish_error(f"GUI finalization failed: {type(exc).__name__}: {exc}")
                    elif kind == "error":
                        self._worker_terminal_seen = True
                        self._finish_error(str(payload))
                except Exception as exc:
                    # One malformed/logging message must not strand all later
                    # messages or prevent the next after() poll from being
                    # scheduled.
                    try:
                        self._append_log(f"ERROR: queue message {kind!r} failed: {type(exc).__name__}: {exc}")
                    except Exception:
                        pass
        except queue.Empty:
            pass
        finally:
            # A worker that died without a terminal message is a real startup
            # failure, not a slow run.  Normal workers always enqueue done or
            # error from their finally path before this check can fire.
            worker = _get_instance_attr(self, "worker", None)
            if (
                self.run_state == "running"
                and worker is not None
                and not worker.is_alive()
                and not _get_instance_attr(self, "_worker_terminal_seen", False)
                and self.queue.empty()
            ):
                self._worker_terminal_seen = True
                try:
                    self._finish_error("Simulation worker terminated without a result or error message")
                except Exception:
                    pass
            self._schedule_queue_poll()

    def _set_details_text(self, text: str) -> None:
        self.details_text.configure(state="normal")
        self.details_text.delete("1.0", "end")
        self.details_text.insert("1.0", text)
        self.details_text.configure(state="disabled")

    def _populate_event_table(self, results=None, summary: Optional[Mapping[str, object]] = None) -> None:
        all_rows = event_table_rows(results) if results is not None else event_table_rows_from_summary(summary or {})
        # The complete event summaries remain in summary.json and the full
        # SimulationResult list remains in memory for event selection.  Tk
        # only needs a bounded, responsive view; inserting millions of rows on
        # the main thread was a second post-worker freeze in large runs.
        rows = all_rows[:MAX_DISPLAYED_EVENT_ROWS]
        self._visible_event_rows = rows
        for item in self.event_tree.get_children():
            self.event_tree.delete(item)
        for row_index, row in enumerate(rows):
            mean_time = row.get("mean_detected_arrival_ns")
            values = (
                row.get("event_id", row_index), detector_run_mode_label(row.get("detector_mode", "cherenkov")), f"{float(row.get('physical_expected_photons', 0.0)):.2f}",
                row.get("generated_photons", 0), row.get("transported_photons", 0), row.get("tir_plus_x", 0),
                row.get("detector_face_crossings", 0), row.get("valid_sensor_hits", 0), row.get("detected_photons", 0),
                row.get("detected_plus_x", 0), row.get("detected_minus_x", 0), "n/a" if mean_time is None else f"{float(mean_time):.4f}",
                row.get("requested_photons", 0), row.get("detector_attempts", 0), row.get("tir_minus_x", 0), row.get("detector_losses", 0),
                row.get("scintillator_generated_photons", 0), row.get("scintillator_photons_reaching_readout", 0), row.get("scintillator_detected_photons", 0), row.get("scintillator_mylar_reflections", 0), row.get("scintillator_mylar_losses", 0),
            )
            self.event_tree.insert("", "end", iid=f"event-{row_index}", values=values)
        if rows:
            first_event_id = int(rows[0].get("event_id", 0))
            self.selected_event.set(str(first_event_id))
            self.event_tree.selection_set("event-0")
            self.event_tree.focus("event-0")
            self._update_event_details(0)
        else:
            self._set_details_text("No event rows are available yet.")
        if len(all_rows) > len(rows):
            self._append_log(
                f"Event table shows first {len(rows)} of {len(all_rows)} events; full per-event summaries remain in summary.json."
            )

    def _update_event_details(self, index: int) -> None:
        if self.latest_results and 0 <= index < len(self.latest_results):
            display = self.latest_settings.display_photons if self.latest_settings else None
            self._set_details_text(format_event_details(self.latest_results[index], display))
            return
        rows = _get_instance_attr(self, "_visible_event_rows", [])
        if 0 <= index < len(rows):
            row = rows[index]
            self._set_details_text(
                f"Event {row.get('event_id', index)}\nSampling mode: {row.get('sampling_mode', 'fixed')}\n"
                f"Physical expected yield: {row.get('physical_expected_photons', 0.0)} photons/event\n"
                f"Requested/generated/transported: {row.get('requested_photons', 0)} / {row.get('generated_photons', 0)} / {row.get('transported_photons', 0)}\n"
                f"Detector attempts: {row.get('detector_attempts', 0)}; face crossings: {row.get('detector_face_crossings', 0)}\n"
                f"TIR +x / -x: {row.get('tir_plus_x', 0)} / {row.get('tir_minus_x', 0)}\n"
                f"Detector losses: {row.get('detector_losses', 0)}\n"
                f"Detected +x / -x: {row.get('detected_plus_x', 0)} / {row.get('detected_minus_x', 0)}\n"
                f"Mean detected arrival: {row.get('mean_detected_arrival_ns') or 'n/a'} ns\n\n"
                "Full particle and photon path details are available for runs performed in this GUI process."
            )

    def _table_selection_changed(self, _event=None) -> None:
        selection = self.event_tree.selection()
        if not selection:
            return
        try:
            index = int(str(selection[0]).split("-", 1)[1])
        except (IndexError, ValueError):
            return
        event_id = int(self._visible_event_rows[index].get("event_id", index)) if 0 <= index < len(self._visible_event_rows) else index
        self.selected_event.set(str(event_id))
        self._update_event_details(index)
        self._refresh_open_event_viewer(index)
        self._refresh_results_tab()

    def _refresh_open_event_viewer(self, index: int) -> None:
        """Refresh event data while retaining the user's camera controls."""

        viewer = self.interactive_event_viewer
        if viewer is None or not self.latest_results or not (0 <= index < len(self.latest_results)):
            return
        try:
            if not viewer.winfo_exists():
                return
        except tk.TclError:
            return
        result = self.latest_results[index]
        mode = getattr(result, "detector_mode", "auto")
        if mode != _get_instance_attr(self, "interactive_event_mode", mode):
            self.interactive_event_yaw = 0.75
            self.interactive_event_pitch = 0.45
            self.interactive_event_zoom = 1.0
            self.interactive_event_mode = mode
        self.interactive_event_result = result
        self.interactive_event_scene = self._interactive_scene(result, self._event_view_categories())
        self.interactive_event_filter_message = (
            "No photon categories selected; solid and particle track remain visible."
            if not self._event_view_categories() else
            "No photon paths match the selected categories."
            if not self.interactive_event_scene["photon_lines"] else ""
        )
        self._draw_interactive_event_view()
        self.event_view_status.set(f"Interactive viewer updated: event {index}")

    def _open_image(self, path: Path) -> None:
        path = Path(path).expanduser().resolve()
        self._append_log(f"Event plot: {path}")
        try:
            if sys.platform.startswith("linux"):
                subprocess.Popen(["xdg-open", str(path)])
            elif sys.platform == "darwin":
                subprocess.Popen(["open", str(path)])
            else:
                os.startfile(str(path))
            self.status.set(f"Opened event plot: {path.name}")
        except Exception as exc:
            self.status.set("Event plot saved")
            self._append_log(f"Open viewer unavailable: {exc}")
            messagebox.showinfo("Event plot saved", f"The selected event plot was saved at:\n{path}\n\nOpen it with an image viewer.", parent=self)

    def _event_view_categories(self) -> set[str]:
        """Return the current persistent viewer checkbox selection."""

        return {
            key for key, variable in self.event_view_category_vars.items()
            if bool(variable.get())
        }

    def _event_view_filters_changed(self) -> None:
        """Redraw an already-open interactive viewer after a checkbox change."""

        viewer = self.interactive_event_viewer
        result = self.interactive_event_result
        if viewer is None or result is None:
            return
        try:
            if not viewer.winfo_exists():
                return
        except tk.TclError:
            return
        categories = self._event_view_categories()
        self.interactive_event_scene = self._interactive_scene(result, categories)
        self.interactive_event_filter_message = (
            "No photon categories selected; solid and particle track remain visible."
            if not categories else
            "No photon paths match the selected categories."
            if not self.interactive_event_scene["photon_lines"] else ""
        )
        self._draw_interactive_event_view()
        self.event_view_status.set("Viewer filters updated")

    def _result_index_for_event_id(self, event_id: int) -> Optional[int]:
        for index, result in enumerate(self.latest_results):
            if int(getattr(result, "event_id", index)) == int(event_id):
                return index
        # Compatibility with the old Event index entry when event IDs are
        # contiguous but a caller supplied the row index.
        if 0 <= int(event_id) < len(self.latest_results):
            return int(event_id)
        return None

    def _show_selected_event(self) -> None:
        try:
            index = int(self.selected_event.get().strip())
        except ValueError:
            messagebox.showerror("Event selection", "Event index must be a non-negative integer", parent=self)
            return
        if index < 0:
            messagebox.showerror("Event selection", "Event index must be non-negative", parent=self)
            return
        if self.latest_results:
            result_index = self._result_index_for_event_id(index)
            if result_index is None:
                messagebox.showerror("Event selection", f"Choose an event_id from {self._available_result_event_ids()[0]} to {self._available_result_event_ids()[-1]}", parent=self)
                return
            self._update_event_details(result_index)
            self._refresh_open_event_viewer(result_index)
            display_count = self.latest_settings.display_photons if self.latest_settings else int(default_values()["display_photons"])
            track_margin_mm = self.latest_settings.track_margin_mm if self.latest_settings else 20.0
            try:
                run_directory = self.last_run_directory or Path(self.values["output_directory"].get()).expanduser().resolve()
                output_path = create_selected_event_plot(
                    self.latest_results,
                    result_index,
                    run_directory,
                    display_count,
                    track_margin_mm,
                    self._event_view_categories(),
                    canonical_scene=self.latest_event_global_scene,
                )
            except Exception as exc:
                self._append_log(f"ERROR: could not create selected event plot: {exc}")
                messagebox.showerror("Event plot", str(exc), parent=self)
                return
            self._open_image(output_path)
            return
        self._update_event_details(index)
        base = Path(self.values["output_directory"].get()).expanduser().resolve()
        latest = latest_run_directory(base)
        if latest is None:
            messagebox.showinfo("Event selection", "Run a simulation first. No prior quartz_run_* directory was found.", parent=self)
            return
        self.last_run_directory = latest
        self.latest_run_text.set(f"Latest run: {latest}")
        if self._event_view_categories() != {"photon_paths"}:
            messagebox.showinfo(
                "Event filters",
                "Saved event images contain the original paths. Run the simulation in this GUI process to apply the selected photon-category filters.",
                parent=self,
            )
            return
        selected_path = latest / f"event_display_event_{index}.png"
        if not selected_path.exists():
            selected_path = latest / f"event_display_event_{index}_selected.png"
        if selected_path.exists():
            self._open_image(selected_path)
        else:
            messagebox.showinfo("Event plot", f"No saved plot for event {index} was found in:\n{latest}\n\nRun the simulation in this GUI to retain full in-memory event paths.", parent=self)

    @staticmethod
    def _event_geometry_rings(geometry):
        """Return geometry rings in mm for the interactive event canvas."""

        if isinstance(geometry, (simulator.SolidGeometry, simulator.CompoundGeometry)):
            return [[tuple(float(value) * 1e3 for value in point) for point in ring] for ring in geometry.surface_rings()]
        hx, hy, hz = geometry.half_x * 1e3, geometry.half_y * 1e3, geometry.half_z * 1e3
        return [
            [(-hx, -hy, -hz), (hx, -hy, -hz), (hx, hy, -hz), (-hx, hy, -hz)],
            [(-hx, -hy, hz), (hx, -hy, hz), (hx, hy, hz), (-hx, hy, hz)],
        ]

    def _interactive_scene(self, result, selected_categories=None):
        canonical = _get_instance_attr(self, "interactive_event_global_scene", None)
        if canonical is None:
            canonical = global_assembly_scene(self.values, self._geometry_for_simulation(), self.assembly_models)
        mode = _get_instance_attr(self, "interactive_event_mode", None)
        display_count = self.latest_settings.display_photons if self.latest_settings else 40
        margin = self.latest_settings.track_margin_mm if self.latest_settings else 20.0
        return build_interactive_event_scene(
            result, canonical, detector_mode=mode, selected_categories=selected_categories,
            display_photons=display_count, track_margin_mm=margin,
        )

    def _open_selected_event_3d_viewer(self) -> None:
        try:
            index = int(self.selected_event.get().strip())
        except (AttributeError, TypeError, ValueError):
            messagebox.showerror("Event selection", "Event index must be a non-negative integer", parent=self)
            return
        result_index = self._result_index_for_event_id(index) if self.latest_results else None
        if index < 0 or not self.latest_results or result_index is None:
            messagebox.showinfo("Interactive event viewer", "Run a simulation in this GUI first; the interactive viewer requires in-memory photon paths.", parent=self)
            return
        if self.interactive_event_viewer is not None and self.interactive_event_viewer.winfo_exists():
            self.interactive_event_viewer.destroy()
        viewer = tk.Toplevel(self)
        viewer.title(f"Interactive 3D event viewer — event {index}")
        viewer.geometry("860x650")
        viewer.minsize(560, 420)
        viewer.columnconfigure(0, weight=1)
        viewer.columnconfigure(1, weight=0)
        viewer.rowconfigure(0, weight=1)
        canvas = tk.Canvas(viewer, background="white", highlightthickness=1, highlightbackground="#909090")
        canvas.grid(row=0, column=0, sticky="nsew", padx=6, pady=6)
        info_frame = ttk.LabelFrame(viewer, text="Event legend / details", padding=8)
        info_frame.grid(row=0, column=1, sticky="ns", padx=(0, 6), pady=6)
        info_frame.columnconfigure(0, weight=1)
        info_frame.rowconfigure(0, weight=1)
        self.interactive_event_info = tk.StringVar(value="Preparing event view...")
        event_info_text = tk.Text(info_frame, width=36, height=28, wrap="word", relief="flat", borderwidth=0,
                                  background=_ttk_background(info_frame), takefocus=True)
        event_info_text.grid(row=0, column=0, sticky="nsew")
        event_info_scroll = ttk.Scrollbar(info_frame, orient="vertical", command=event_info_text.yview)
        event_info_scroll.grid(row=0, column=1, sticky="ns")
        event_info_text.configure(yscrollcommand=event_info_scroll.set, state="disabled")
        self.interactive_event_info_text = event_info_text
        controls = ttk.Frame(viewer)
        controls.grid(row=1, column=0, columnspan=2, sticky="ew", padx=6, pady=(0, 6))
        ttk.Button(controls, text="Reset view", command=self._reset_interactive_event_view).pack(side="left")
        ttk.Checkbutton(controls, text="Show compact callouts", variable=self.show_3d_labels,
                        command=self._draw_interactive_event_view).pack(side="left", padx=(12, 0))
        ttk.Label(controls, text="Drag to rotate • mouse wheel to zoom • colors show photon outcomes", foreground="#555555").pack(side="left", padx=(10, 0))
        ttk.Button(controls, text="Close", command=viewer.destroy).pack(side="right")
        self.interactive_event_viewer = viewer
        self.interactive_event_canvas = canvas
        self.interactive_event_result = self.latest_results[result_index]
        self.interactive_event_mode = getattr(self.interactive_event_result, "detector_mode", "auto")
        self.interactive_event_global_scene = _get_instance_attr(self, "latest_event_global_scene", None)
        if self.interactive_event_global_scene is None:
            self.interactive_event_global_scene = global_assembly_scene(self.values, self._geometry_for_simulation(), self.assembly_models)
        self.interactive_event_yaw = 0.75
        self.interactive_event_pitch = 0.45
        self.interactive_event_zoom = 1.0
        categories = self._event_view_categories()
        self.interactive_event_scene = self._interactive_scene(self.interactive_event_result, categories)
        self.interactive_event_filter_message = (
            "No photon categories selected; solid and particle track remain visible."
            if not categories else
            "No photon paths match the selected categories."
            if not self.interactive_event_scene["photon_lines"] else ""
        )
        self.interactive_event_drag = None
        canvas.bind("<Configure>", lambda _event: self._draw_interactive_event_view())
        canvas.bind("<ButtonPress-1>", self._interactive_event_press)
        canvas.bind("<B1-Motion>", self._interactive_event_drag_motion)
        canvas.bind("<MouseWheel>", self._interactive_event_zoom_event)
        canvas.bind("<Button-4>", self._interactive_event_zoom_event)
        canvas.bind("<Button-5>", self._interactive_event_zoom_event)
        viewer.protocol("WM_DELETE_WINDOW", viewer.destroy)
        self._draw_interactive_event_view()
        self.event_view_status.set(f"Interactive viewer open: event {index}")

    def _reset_interactive_event_view(self) -> None:
        self.interactive_event_yaw = 0.75
        self.interactive_event_pitch = 0.45
        self.interactive_event_zoom = 1.0
        self._draw_interactive_event_view()

    def _interactive_event_press(self, event) -> None:
        self.interactive_event_drag = (event.x, event.y, self.interactive_event_yaw, self.interactive_event_pitch)

    def _interactive_event_drag_motion(self, event) -> None:
        if self.interactive_event_drag is None:
            return
        start_x, start_y, start_yaw, start_pitch = self.interactive_event_drag
        self.interactive_event_yaw = start_yaw + (event.x - start_x) * 0.012
        self.interactive_event_pitch = max(-1.45, min(1.45, start_pitch + (event.y - start_y) * 0.012))
        self._draw_interactive_event_view()

    def _interactive_event_zoom_event(self, event) -> None:
        delta = getattr(event, "delta", 0)
        if delta == 0:
            delta = 120 if getattr(event, "num", 4) == 4 else -120
        self.interactive_event_zoom *= 1.12 if delta > 0 else 1.0 / 1.12
        self.interactive_event_zoom = max(0.2, min(8.0, self.interactive_event_zoom))
        self._draw_interactive_event_view()

    def _draw_interactive_event_view(self) -> None:
        canvas = self.interactive_event_canvas
        scene = _get_instance_attr(self, "interactive_event_scene", None)
        if canvas is None or scene is None or not canvas.winfo_exists():
            return
        canvas.delete("all")
        center = tuple(scene["center_mm"])
        extent = float(scene["extent_mm"])
        # Fit against the actual canvas dimensions after reserving the
        # responsive title/status and legend bands.
        width = max(1, int(canvas.winfo_width()))
        height = max(1, int(canvas.winfo_height()))
        title_font = tkfont.Font(root=canvas, family="TkDefaultFont", size=12, weight="bold") if tkfont else ("TkDefaultFont", 12, "bold")
        status_font = tkfont.Font(root=canvas, family="TkDefaultFont", size=9) if tkfont else ("TkDefaultFont", 9)
        title_measure = title_font.measure if tkfont else (lambda text: _event_viewer_default_measure(text, 12.0))
        status_measure = status_font.measure if tkfont else (lambda text: _event_viewer_default_measure(text, 9.0))
        selected_event = _get_instance_attr(self, "selected_event", None)
        selected_event_text = selected_event.get() if selected_event is not None else "?"
        filter_message = _get_instance_attr(self, "interactive_event_filter_message", "")
        layout = interactive_event_view_layout(
            width, height,
            f"INTERACTIVE EVENT 3D • event {selected_event_text}",
            filter_message or "C/S assemblies • selected photon categories • drag to rotate",
            title_measure=title_measure, status_measure=status_measure,
            title_line_height=title_font.metrics("linespace") if tkfont else 18.0,
            status_line_height=status_font.metrics("linespace") if tkfont else 13.0,
        )
        scene_left, scene_top, scene_right, scene_bottom = layout["scene_rect"]
        scene_width = max(scene_right - scene_left, 1.0)
        scene_height = max(scene_bottom - scene_top, 1.0)
        scene_center_x = (scene_left + scene_right) / 2.0
        scene_center_y = (scene_top + scene_bottom) / 2.0
        scale = min(scene_width, scene_height) * 0.72 / max(extent, 1e-9) * self.interactive_event_zoom
        label_font = tkfont.Font(root=canvas, family="TkDefaultFont", size=9, weight="bold") if tkfont else ("TkDefaultFont", 9, "bold")
        label_measure = (lambda text: (label_font.measure(text) + 8, label_font.metrics("linespace") + 2)) if tkfont else None
        cos_yaw, sin_yaw = math.cos(self.interactive_event_yaw), math.sin(self.interactive_event_yaw)
        cos_pitch, sin_pitch = math.cos(self.interactive_event_pitch), math.sin(self.interactive_event_pitch)

        def project(point):
            x = point[0] - center[0]
            y = point[1] - center[1]
            z = point[2] - center[2]
            rotated_x = cos_yaw * x - sin_yaw * y
            rotated_y = sin_yaw * x + cos_yaw * y
            screen_y = cos_pitch * rotated_y - sin_pitch * z
            depth = sin_pitch * rotated_y + cos_pitch * z
            return scene_center_x + rotated_x * scale, scene_center_y - screen_y * scale, depth

        callout_items = []
        for component in scene["geometry_components"]:
            projected_component_points = [
                project(point)
                for face in component.get("faces_mm", ())
                for point in face
            ]
            if not projected_component_points:
                continue
            target = (
                sum(point[0] for point in projected_component_points) / len(projected_component_points),
                sum(point[1] for point in projected_component_points) / len(projected_component_points),
            )
            callout_items.append({
                "component_id": component.get("component_id", component.get("label", "component")),
                "assembly_type": component.get("assembly_type"),
                "kind": component.get("kind", "geometry"),
                "label": detector_callout_label(component),
                "detail": component.get("detail", component.get("label", "component")),
                "target": target,
            })

        # Faces are drawn back-to-front with detector-specific assembly colors.
        projected_faces = []
        for face, color, assembly_type, kind in scene["geometry_faces"]:
            points = [project(point) for point in face]
            projected_faces.append((sum(point[2] for point in points) / len(points), points, color, assembly_type, kind))
        for _depth, points, color, _assembly_type, kind in sorted(projected_faces, key=lambda item: item[0]):
            fill = self._geometry_viewer_color(color, 0.86 if kind in {"geometry", "guide"} else 1.05)
            canvas.create_polygon([coordinate for point in points for coordinate in point[:2]], fill=fill, outline="", stipple="gray50")
        for line, assembly_type, kind in scene["geometry_lines"]:
            points = [project(point) for point in line]
            edge = "#263238" if kind not in {"mylar_overlay"} else "#8c6d1f"
            canvas.create_line([coordinate for point in points for coordinate in point[:2]], fill=edge, width=2)
        track = scene["track_mm"]
        if len(track) >= 2:
            track_points = [project(point) for point in track]
            canvas.create_line([coordinate for point in track_points for coordinate in point[:2]], fill="#111111", width=3)
            canvas.create_oval(track_points[0][0] - 4, track_points[0][1] - 4, track_points[0][0] + 4, track_points[0][1] + 4, fill="#00a000", outline="")
            canvas.create_rectangle(track_points[-1][0] - 4, track_points[-1][1] - 4, track_points[-1][0] + 4, track_points[-1][1] + 4, fill="#a00000", outline="")
        for marker in scene["impact_markers"]:
            point = project(marker["point_mm"])
            color = "#1f77b4" if marker["assembly_type"] == "cherenkov" else "#17a673"
            canvas.create_oval(point[0] - 4, point[1] - 4, point[0] + 4, point[1] + 4, fill=color, outline="#17324d")
        for line in scene["photon_lines"]:
            points = [project(point) for point in line["points_mm"]]
            canvas.create_line([coordinate for point in points for coordinate in point[:2]], fill=line["color"], width=1)
            if line.get("detected_marker") and points:
                point = points[-1]
                canvas.create_polygon(
                    point[0], point[1] - 4, point[0] + 4, point[1],
                    point[0], point[1] + 4, point[0] - 4, point[1],
                    fill=line["color"], outline="#111827",
                )
        label_layout = collision_aware_callout_layout(
            callout_items, width, height,
            occupied_rectangles=layout["occupied_rectangles"],
            show_labels=self._three_d_labels_enabled(), text_measure=label_measure,
        )
        self.interactive_event_label_layout = label_layout
        for callout in label_layout["callouts"]:
            start_x, start_y = callout["leader_start"]
            target_x, target_y = callout["target"]
            canvas.create_line(start_x, start_y, target_x, target_y,
                               fill="#455a64", dash=(3, 2), width=1)
            canvas.create_text(callout["position"][0], callout["position"][1],
                               anchor="nw", text=callout["label"], fill="#17324d",
                               font=label_font)
        legend_lines = compact_scene_legend_lines(callout_items, label_layout["omitted"], scene.get("mode", "combined"))
        legend_x1, legend_y1, legend_x2, legend_y2 = layout["legend_rect"]
        canvas.create_rectangle(legend_x1, legend_y1, legend_x2, legend_y2,
                                fill="#ffffff", outline="#b0b0b0")
        canvas_legend_lines = (
            f"MODE: {scene.get('mode', 'combined').upper()}",
            "C/S geometry + sensor readouts",
            "black = particle track",
            "colored = photon outcomes",
            "compact labels prefer sensors/readouts",
        )
        for row, label in enumerate(canvas_legend_lines):
            canvas.create_text(legend_x1 + 7, legend_y1 + 7 + 14 * row,
                               anchor="nw", text=label, fill="#333333",
                               font=("TkDefaultFont", 8, "bold" if row == 0 else "normal"))
        header_x1, header_y1, header_x2, header_y2 = layout["header_rect"]
        canvas.create_rectangle(header_x1, header_y1, header_x2, header_y2,
                                fill="#f5f8fb", outline="")
        title_x = (layout["title_rect"][0] + layout["title_rect"][2]) / 2.0
        canvas.create_text(title_x, layout["title_rect"][1], anchor="n", text=layout["title"],
                           fill="#17324d", font=title_font)
        status_color = "#8a1c1c" if _get_instance_attr(self, "interactive_event_filter_message", "") else "#8a5200"
        canvas.create_text(layout["status_rect"][0], layout["status_rect"][1], anchor="nw",
                           text=layout["status"], fill=status_color, font=status_font,
                           width=max(24.0, layout["status_rect"][2] - layout["status_rect"][0]),
                           justify="left")
        if _has_instance_attr(self, "interactive_event_info"):
            status_detail = layout["full_status"] or "Photon paths are shown for the selected categories."
            self._set_viewer_details_text(
                _get_instance_attr(self, "interactive_event_info_text", None),
                f"MODE: {scene['mode']}\n\nSTATUS\n{status_detail}\n\n" +
                "\n".join(legend_lines) + "\n\n" +
                f"Callouts shown: {len(label_layout['callouts'])}\n"
                f"In legend: {len(label_layout['omitted'])}\n"
                f"Photon paths shown: {len(scene['photon_lines'])}"
            )

    def _load_latest_summary(self) -> None:
        output_root = Path(self.values["output_directory"].get()).expanduser().resolve()
        gamma_payload, _gamma_run, _gamma_artifact = load_latest_gamma_result(output_root)
        self.latest_gamma_payload = gamma_payload
        if isinstance(gamma_payload, Mapping):
            first_event_id = next((item.get("event_id") for item in gamma_payload.get("events", ()) if isinstance(item, Mapping)), 0)
            self.gamma_selected_event.set(str(first_event_id))
        latest = latest_run_directory(output_root)
        if latest is None or not (latest / "summary.json").exists():
            self._refresh_gamma_results_view()
            return
        try:
            summary = json.loads((latest / "summary.json").read_text(encoding="utf-8"))
        except (OSError, json.JSONDecodeError) as exc:
            self._append_log(f"Could not load latest summary: {exc}")
            return
        self.last_run_directory = latest
        self.latest_run_text.set(f"Latest run: {latest}")
        self.latest_summary = summary.get("summary", {}) if isinstance(summary.get("summary", {}), Mapping) else {}
        self.latest_readout_payload = None
        readout_path = summary.get("readout_arrivals_path") or str(latest / "readout_arrivals.json")
        try:
            candidate = Path(str(readout_path))
            if not candidate.is_absolute():
                candidate = latest / candidate
            if candidate.exists():
                payload = json.loads(candidate.read_text(encoding="utf-8"))
                if payload.get("schema") == "readout-arrivals-v1":
                    self.latest_readout_payload = payload
        except (OSError, json.JSONDecodeError, TypeError, ValueError) as exc:
            self._append_log(f"Could not load readout arrival artifact: {exc}")
        self.latest_event_global_scene = summary.get("event_view_scene")
        self._populate_event_table(summary=self.latest_summary)
        aggregate = summary_with_diagnostic_defaults(summary.get("summary", {}))
        if isinstance(aggregate, Mapping):
            self.total_summary_text.set(f"Latest saved run: {aggregate.get('events', 0)} events, detector mode {detector_run_mode_label(aggregate.get('detector_mode', 'cherenkov'))}, generated {aggregate.get('generated_photons', 0)}, scintillator generated {aggregate.get('scintillator_generated_photons', 0)}, detected {aggregate.get('detected_photons', 0)}, scintillator detected {aggregate.get('scintillator_detected_photons', 0)}")
        self._refresh_results_tab()
        self._refresh_gamma_results_view()

    def _open_latest_run_folder(self) -> None:
        directory = self.last_run_directory or latest_run_directory(Path(self.values["output_directory"].get()))
        if directory is None:
            messagebox.showerror("Latest run", "No quartz_run_* output directory was found yet.", parent=self)
            return
        self.last_run_directory = directory
        self.latest_run_text.set(f"Latest run: {directory}")
        self._open_output_folder()

    def _finish_run(self, payload: Mapping[str, object]) -> None:
        self._worker_terminal_seen = True
        self.run_state = "stopped" if bool(payload.get("stopped")) or self.stop_requested.is_set() else "completed"
        self.run_button.configure(state="normal")
        self.stop_button.configure(state="disabled")
        self.progress.stop()
        summary = payload["summary"]
        display_summary = summary_with_diagnostic_defaults(summary) if isinstance(summary, Mapping) else summary
        mode = detector_run_mode(display_summary.get("detector_mode", "cherenkov")) if isinstance(display_summary, Mapping) else "cherenkov"
        mode_label = detector_run_mode_label(mode)
        self.last_run_directory = Path(str(payload["run_directory"])).expanduser().resolve()
        self.latest_run_text.set(f"Latest run: {self.last_run_directory}")
        live_results = payload.get("_results", [])
        self.latest_results = list(live_results) if isinstance(live_results, list) else []
        self.latest_summary = summary if isinstance(summary, Mapping) else {}
        gamma_run_stopped = bool(payload.get("stopped")) or self.stop_requested.is_set()
        previous_gamma_payload = self.latest_gamma_payload
        self.latest_gamma_payload = retain_latest_completed_gamma_payload(previous_gamma_payload, payload)
        if self.latest_gamma_payload is not previous_gamma_payload:
            gamma_events = self.latest_gamma_payload.get("events", ())
            first_event_id = next((item.get("event_id") for item in gamma_events if isinstance(item, Mapping)), 0)
            self.gamma_selected_event.set(str(first_event_id))
        self.latest_readout_payload = None
        readout_path = payload.get("readout_arrivals_path")
        if readout_path:
            try:
                self.latest_readout_payload = json.loads(Path(str(readout_path)).read_text(encoding="utf-8"))
            except (OSError, json.JSONDecodeError, TypeError, ValueError) as exc:
                self._append_log(f"Could not load in-memory readout arrival artifact: {exc}")
        settings_object = payload.get("_settings_object")
        self.latest_settings = settings_object if isinstance(settings_object, RunSettings) else None
        self.latest_event_global_scene = payload.get("_event_view_scene") or payload.get("event_view_scene")
        if self.interactive_event_viewer is not None:
            try:
                if self.interactive_event_viewer.winfo_exists():
                    self.interactive_event_viewer.destroy()
            except tk.TclError:
                pass
            self.interactive_event_viewer = None
            self.interactive_event_canvas = None
        # Use compact aggregate event summaries for the table.  The old path
        # rescanned every raw photon record in the Tk thread; raw paths are
        # retained in latest_results solely for the selected-event viewer.
        if isinstance(summary, Mapping):
            self._populate_event_table(summary=summary)
        self._refresh_results_tab()
        self._refresh_gamma_results_view()
        if isinstance(display_summary, Mapping):
            if mode == GAMMA_RUN_MODE:
                if gamma_run_stopped:
                    self.total_summary_text.set("Gamma run stopped; the last completed Gamma result is retained in Gamma Results.")
                else:
                    self.total_summary_text.set(
                        f"Gamma total: {display_summary.get('events', 0)} events; deposited {float(display_summary.get('total_deposited_energy_keV', 0.0)):.6g} keV; "
                        f"charged secondaries {display_summary.get('total_charged_secondary_tracks', 0)}; local expected light only. See Gamma Results."
                    )
            else:
                self.total_summary_text.set(f"Total: {display_summary.get('events', 0)} events, detector mode {mode_label}, requested {display_summary.get('requested_photons', 0)}, generated {display_summary.get('generated_photons', 0)}, scintillator generated {display_summary.get('scintillator_generated_photons', 0)}, detected {display_summary.get('detected_photons', 0)}, scintillator detected {display_summary.get('scintillator_detected_photons', 0)}")
        was_stopped = gamma_run_stopped
        if was_stopped:
            self.status.set(f"Stopped cooperatively · {mode_label}")
            self.progress_text.set(f"Stopped ({mode_label}): retained {len(self.latest_results)} completed/partial event(s)")
        else:
            self.progress.configure(value=100.0)
            self.status.set(f"Completed · {mode_label}")
            self.progress_text.set(f"Completed ({mode_label}): {display_summary.get('events', 0) if isinstance(display_summary, Mapping) else 0} event(s)")
        self._append_log(f"Completed: {payload['run_directory']}")
        event_zero_path = self.last_run_directory / "event_display_event_0.png"
        if event_zero_path.exists():
            self._append_log(f"Event 0 plot ready: {event_zero_path}")
        if isinstance(display_summary, Mapping):
            if mode == GAMMA_RUN_MODE:
                self._append_log(f"Gamma summary: events={display_summary.get('events', 0)}, deposited_keV={display_summary.get('total_deposited_energy_keV', 0.0)}, charged_secondaries={display_summary.get('total_charged_secondary_tracks', 0)}, expected_local_light_only=true")
            else:
                self._append_log(f"Summary: detector_mode={mode_label}, events={display_summary.get('events', 0)}, generated={display_summary.get('generated_photons', 0)}, scint_generated={display_summary.get('scintillator_generated_photons', 0)}, detected={display_summary.get('detected_photons', 0)}, scint_detected={display_summary.get('scintillator_detected_photons', 0)}")
            if display_summary.get("tir_plus_x", 0) and not display_summary.get("detector_face_crossings", 0):
                self._append_log("DIAGNOSTIC: photons are reaching +x but reflecting by total internal reflection. Increase the +x coupling index; QE does not overcome TIR.")
        root_error = payload.get("root_error")
        if root_error:
            # ROOT is an optional serialization product.  It must not turn a
            # completed in-memory simulation into a GUI error or hide the
            # preserved JSON/PNG outputs behind a misleading physics dialog.
            self.status.set("Completed (ROOT optional output unavailable)")
            self._append_log("ROOT OUTPUT ONLY: " + str(root_error))
            self.progress_text.set(f"Completed ({mode_label}); JSON/PNG outputs saved; ROOT unavailable")
            messagebox.showinfo(
                "Simulation complete",
                f"Outputs saved in:\n{payload['run_directory']}\n\nOptional ROOT output was not saved; see run.log for details.",
                parent=self,
            )
        else:
            messagebox.showinfo("Simulation stopped" if was_stopped else "Simulation complete", f"Outputs saved in:\n{payload['run_directory']}", parent=self)

    def _finish_error(self, message: str) -> None:
        self._worker_terminal_seen = True
        self.run_state = "error"
        self.run_button.configure(state="normal")
        self.stop_button.configure(state="disabled")
        self.progress.stop()
        self.progress.configure(value=0.0)
        self.status.set("Error")
        self.progress_text.set("Error: simulation did not complete")
        self._append_log("ERROR: " + message)
        messagebox.showerror("Simulation error", message, parent=self)

    def _open_output_folder(self) -> None:
        directory = self.last_run_directory or Path(self.values["output_directory"].get()).expanduser().resolve()
        if not directory.exists():
            messagebox.showerror("Output folder", f"Folder does not exist yet:\n{directory}", parent=self)
            return
        try:
            if sys.platform.startswith("linux"):
                subprocess.Popen(["xdg-open", str(directory)])
            elif sys.platform == "darwin":
                subprocess.Popen(["open", str(directory)])
            else:
                os.startfile(str(directory))
        except Exception as exc:
            messagebox.showerror("Output folder", f"Could not open {directory}: {exc}", parent=self)

    # Geometry editor -------------------------------------------------

    @staticmethod
    def _geometry_shape_key(value: str) -> str:
        value = str(value).strip()
        key = value.lower()
        return GEOMETRY_SHAPE_BY_LABEL.get(value, GEOMETRY_SHAPE_ALIASES.get(key, key))

    @staticmethod
    def _geometry_shape_label(shape: str) -> str:
        return GEOMETRY_SHAPE_LABELS.get(shape, shape)

    @staticmethod
    def _geometry_dimension_summary(section: GeometrySection) -> str:
        if section.shape == "box":
            return f"{section.length_x_mm:g} × {section.width_y_mm:g} × {section.thickness_z_mm:g}"
        if section.shape == "tapered_rectangular_prism":
            return f"X {section.length_x_mm:g}; Y {section.width_y_start_mm:g}→{section.width_y_end_mm:g}; Z {section.thickness_z_start_mm:g}→{section.thickness_z_end_mm:g}"
        if section.shape == "cylinder":
            return f"R {section.radius_mm:g}; L {section.length_mm:g}; {section.axis}"
        if section.shape == "elliptical_cylinder":
            return f"A {section.semi_axis_a_mm:g}; B {section.semi_axis_b_mm:g}; L {section.length_mm:g}; {section.axis}"
        return f"R {section.radius_start_mm:g}→{section.radius_end_mm:g}; L {section.length_mm:g}; {section.axis}"

    def _update_geometry_support_status(self) -> None:
        if not _has_instance_attr(self, "geometry_support_text"):
            return
        try:
            geometry = self._geometry_for_simulation()
            geometry.validate()
        except (TypeError, ValueError) as exc:
            names = ", ".join(section.name for section in self.geometry_sections) or "no sections"
            self.geometry_support_text.set(f"Preview only / invalid for Run: {names}. {exc}")
            return
        if isinstance(geometry, simulator.CompoundGeometry):
            description = f"{len(geometry.sections)} adjacent X sections"
            materials = ", ".join(section.normalized_material() for section in geometry.sections)
            self.geometry_support_text.set(
                f"Simulation supported: CompoundGeometry ({description}; materials {materials}). "
                "Internal interfaces use Fresnel/Snell transport; terminal +/-X portals are available."
            )
        else:
            portal_note = " +/-X planar detector portals available." if geometry.detector_portals_supported else " No +/-X planar detector portal for this axis; end-detector logic is inactive."
            self.geometry_support_text.set(f"Simulation supported: {self._geometry_shape_label(geometry.shape)}, material {geometry.material}.{portal_note}")

    def _refresh_geometry_tree(self) -> None:
        self._suppress_geometry_selection_event = True
        try:
            for item in self.geometry_tree.get_children():
                self.geometry_tree.delete(item)
            for index, section in enumerate(self.geometry_sections):
                self.geometry_tree.insert("", "end", iid=str(index), values=(section.name, self._geometry_shape_label(section.shape), section.material, self._geometry_dimension_summary(section)))
            if self.geometry_sections:
                self.selected_geometry_index = min(self.selected_geometry_index, len(self.geometry_sections) - 1)
                iid = str(self.selected_geometry_index)
                self.geometry_tree.selection_set(iid)
                self.geometry_tree.focus(iid)
        finally:
            self._suppress_geometry_selection_event = False
        self._update_geometry_support_status()

    def _geometry_selection_changed(self, _event=None) -> None:
        if self._suppress_geometry_selection_event:
            return
        selection = self.geometry_tree.selection()
        if not selection:
            return
        try:
            self.selected_geometry_index = int(selection[0])
        except (TypeError, ValueError):
            return
        self._build_geometry_properties()
        self._refresh_geometry_views()

    def _add_geometry_section(self) -> None:
        number = len(self.geometry_sections) + 1
        # Keep the new default section contiguous with the 100 x 20 x 10 mm
        # bar so Add Section immediately produces a valid compound preview.
        self.geometry_sections.append(GeometrySection(f"Section {number}", length_x_mm=20.0, width_y_mm=20.0, thickness_z_mm=10.0, center_x_mm=60.0))
        self.selected_geometry_index = len(self.geometry_sections) - 1
        self._refresh_geometry_tree()
        self._build_geometry_properties()
        self._refresh_geometry_views()

    def _remove_geometry_section(self) -> None:
        if len(self.geometry_sections) <= 1:
            messagebox.showinfo("Remove section", "Keep at least one rectangular section in the geometry editor.", parent=self)
            return
        del self.geometry_sections[self.selected_geometry_index]
        self.selected_geometry_index = max(0, min(self.selected_geometry_index, len(self.geometry_sections) - 1))
        self._refresh_geometry_tree()
        self._build_geometry_properties()
        self._refresh_geometry_views()

    def _build_geometry_properties(self) -> None:
        for child in self.geometry_properties_inner.winfo_children():
            child.destroy()
        if not self.geometry_sections:
            return
        section = self.geometry_sections[self.selected_geometry_index]
        self.geometry_property_vars = {}
        self.geometry_property_widgets = {}
        fields = [("Name", "name", "entry"), ("Shape", "shape", "shape"), ("Material", "material", "material")]
        if section.shape == "box":
            fields += [("Length X (mm)", "length_x_mm", "entry"), ("Width Y (mm)", "width_y_mm", "entry"), ("Thickness Z (mm)", "thickness_z_mm", "entry")]
        elif section.shape == "tapered_rectangular_prism":
            fields += [
                ("X length (mm)", "length_x_mm", "entry"),
                ("Y start width (mm)", "width_y_start_mm", "entry"), ("Y end width (mm)", "width_y_end_mm", "entry"),
                ("Z start thickness (mm)", "thickness_z_start_mm", "entry"), ("Z end thickness (mm)", "thickness_z_end_mm", "entry"),
            ]
        elif section.shape == "cylinder":
            fields += [("Axis", "axis", "axis"), ("Radius (mm)", "radius_mm", "entry"), ("Length (mm)", "length_mm", "entry")]
        elif section.shape == "elliptical_cylinder":
            fields += [("Axis", "axis", "axis"), ("Semi-axis A (mm)", "semi_axis_a_mm", "entry"), ("Semi-axis B (mm)", "semi_axis_b_mm", "entry"), ("Length (mm)", "length_mm", "entry")]
        elif section.shape == "truncated_cone":
            fields += [("Axis", "axis", "axis"), ("Start radius (mm)", "radius_start_mm", "entry"), ("End radius (mm)", "radius_end_mm", "entry"), ("Length (mm)", "length_mm", "entry")]
        fields += [("Center X (mm)", "center_x_mm", "entry"), ("Center Y (mm)", "center_y_mm", "entry"), ("Center Z (mm)", "center_z_mm", "entry")]
        for row, (label, key, field_kind) in enumerate(fields):
            ttk.Label(self.geometry_properties_inner, text=label).grid(row=row, column=0, sticky="w", padx=(2, 10), pady=2)
            value = self._geometry_shape_label(section.shape) if key == "shape" else getattr(section, key)
            variable = tk.StringVar(value=str(value))
            self.geometry_property_vars[key] = variable
            if field_kind == "shape":
                widget = ttk.Combobox(self.geometry_properties_inner, textvariable=variable, values=tuple(GEOMETRY_SHAPE_LABELS.values()), state="readonly")
            elif field_kind == "material":
                widget = ttk.Combobox(self.geometry_properties_inner, textvariable=variable, values=("fused silica", "quartz", "air", "detector glass", "custom"), state="readonly")
            elif field_kind == "axis":
                widget = ttk.Combobox(self.geometry_properties_inner, textvariable=variable, values=GEOMETRY_AXES, state="readonly")
            else:
                # Explicitly normal: dimension fields must remain editable.
                widget = ttk.Entry(self.geometry_properties_inner, textvariable=variable, state="normal")
            widget.grid(row=row, column=1, sticky="ew", padx=(0, 2), pady=2)
            self.geometry_property_widgets[key] = widget
            if field_kind == "entry":
                widget.bind("<Return>", self._geometry_property_changed)
                widget.bind("<FocusOut>", self._geometry_property_changed)
            else:
                widget.bind("<<ComboboxSelected>>", self._geometry_property_changed)
        self.geometry_properties_inner.columnconfigure(1, weight=1)
        note = "Box uses X length, Y width, Z thickness. Tapered prisms use start/end Y and Z dimensions. Cylinders and cones use an axis plus the shown radii/axes and length. All values are mm; centers are offsets from the preview origin."
        ttk.Label(self.geometry_properties_inner, text=note, foreground="#245a35", wraplength=700, justify="left").grid(row=len(fields), column=0, columnspan=2, sticky="w", pady=(7, 2))
        ttk.Button(self.geometry_properties_inner, text="Apply / Update", command=self._geometry_property_changed).grid(row=len(fields) + 1, column=1, sticky="e", pady=(4, 2))
        self._update_geometry_support_status()

    def _geometry_property_changed(self, _event=None) -> Optional[str]:
        if self._loading_geometry_properties or not self.geometry_sections:
            return
        section = self.geometry_sections[self.selected_geometry_index]
        numeric = {
            "length_x_mm", "width_y_mm", "thickness_z_mm", "width_y_start_mm", "width_y_end_mm",
            "thickness_z_start_mm", "thickness_z_end_mm", "radius_mm", "length_mm", "semi_axis_a_mm",
            "semi_axis_b_mm", "radius_start_mm", "radius_end_mm", "center_x_mm", "center_y_mm", "center_z_mm",
        }
        original = section.to_dict()
        shape_changed = False
        try:
            for key, variable in self.geometry_property_vars.items():
                if key == "shape":
                    new_shape = self._geometry_shape_key(variable.get())
                    shape_changed = new_shape != section.shape
                    setattr(section, key, new_shape)
                elif key == "axis":
                    setattr(section, key, variable.get().upper())
                elif key in numeric:
                    value = float(variable.get())
                    setattr(section, key, value)
                else:
                    setattr(section, key, variable.get())
            section.validate()
        except (TypeError, ValueError) as exc:
            for key, value in original.items():
                setattr(section, key, value)
            for key, variable in self.geometry_property_vars.items():
                if key in original:
                    variable.set(self._geometry_shape_label(original[key]) if key == "shape" else str(original[key]))
            self.geometry_support_text.set(f"Invalid geometry value; last valid value restored: {exc}")
            return "break" if getattr(_event, "keysym", "") == "Return" else None
        self._refresh_geometry_tree()
        if shape_changed:
            self._build_geometry_properties()
        self._refresh_geometry_views()
        return "break" if getattr(_event, "keysym", "") == "Return" else None

    def _save_geometry(self) -> None:
        path = filedialog.asksaveasfilename(title="Save geometry", defaultextension=".json", filetypes=(("Geometry JSON", "*.json"), ("All files", "*.*")))
        if not path:
            return
        payload = {"format": "quartz-cherenkov-gui-geometry", "version": 1, "sections": [section.to_dict() for section in self.geometry_sections]}
        try:
            Path(path).write_text(json.dumps(payload, indent=2) + "\n", encoding="utf-8")
            self.status.set(f"Geometry saved: {Path(path).name}")
            self._append_log(f"Saved geometry: {Path(path).resolve()}")
        except OSError as exc:
            messagebox.showerror("Save geometry", str(exc), parent=self)

    def _load_geometry(self) -> None:
        path = filedialog.askopenfilename(title="Load geometry", filetypes=(("Geometry JSON", "*.json"), ("All files", "*.*")))
        if not path:
            return
        try:
            payload = json.loads(Path(path).read_text(encoding="utf-8"))
            raw_sections = payload.get("sections", payload) if isinstance(payload, Mapping) else payload
            if not isinstance(raw_sections, list) or not raw_sections:
                raise ValueError("The file does not contain a non-empty sections list")
            sections = [GeometrySection.from_dict(item, index) for index, item in enumerate(raw_sections) if isinstance(item, Mapping)]
            if not sections:
                raise ValueError("The file does not contain valid geometry sections")
            self.geometry_sections = sections
            self.selected_geometry_index = 0
            self._refresh_geometry_tree()
            self._build_geometry_properties()
            self._refresh_geometry_views()
            self.status.set(f"Geometry loaded: {Path(path).name}")
            self._append_log(f"Loaded geometry: {Path(path).resolve()}")
        except (OSError, json.JSONDecodeError, ValueError, TypeError) as exc:
            messagebox.showerror("Load geometry", str(exc), parent=self)

    @staticmethod
    def _detector_portal_profile(section: GeometrySection, face: str = "plus_x"):
        """Return (available, width_y_mm, height_z_mm, kind, description)."""

        positive_x = face == "plus_x"
        if section.shape == "box":
            width, height, kind = section.width_y_mm, section.thickness_z_mm, "rectangle"
        elif section.shape == "tapered_rectangular_prism":
            width = section.width_y_end_mm if positive_x else section.width_y_start_mm
            height = section.thickness_z_end_mm if positive_x else section.thickness_z_start_mm
            kind = "rectangle"
        elif section.axis == "X" and section.shape == "cylinder":
            width = height = 2.0 * section.radius_mm
            kind = "ellipse"
        elif section.axis == "X" and section.shape == "elliptical_cylinder":
            width, height, kind = 2.0 * section.semi_axis_a_mm, 2.0 * section.semi_axis_b_mm, "ellipse"
        elif section.axis == "X" and section.shape == "truncated_cone":
            radius = section.radius_end_mm if positive_x else section.radius_start_mm
            width = height = 2.0 * radius
            kind = "ellipse"
        else:
            return False, 0.0, 0.0, "none", f"no planar {face.replace('_', ' ')} portal for {section.shape} axis {section.axis}"
        side = "+x" if positive_x else "-x"
        return True, width, height, kind, f"{kind} {side} portal: Y {width:g} mm × Z {height:g} mm"

    def _terminal_editor_section(self, face: str = "plus_x") -> GeometrySection:
        """Return the section owning the external +/-X portal."""

        if not self.geometry_sections:
            return GeometrySection("Quartz bar")
        return self.geometry_sections[-1] if face == "plus_x" else self.geometry_sections[0]

    def _update_detector_active_area(self) -> None:
        if not _has_instance_attr(self, "detector_active_area_text"):
            return
        plus_section = self._terminal_editor_section("plus_x")
        _, _, _, _, plus_description = self._detector_portal_profile(plus_section, "plus_x")
        try:
            plus_spec = detector_active_area_from_values(self.values, "detector", "First detector")
            plus_area_text = self._active_area_description(plus_spec)
        except ValueError as exc:
            plus_area_text = f"Invalid active area: {exc}"
        self.detector_active_area_text.set(f"Active area: {plus_area_text}\nPortal: {plus_description}\nTerminal +x section: {plus_section.name}")
        if _has_instance_attr(self, "detector_face_canvas"):
            self._draw_detector_face(self.detector_face_canvas, "plus_x")
        minus_section = self._terminal_editor_section("minus_x")
        minus_available, _, _, _, minus_description = self._detector_portal_profile(minus_section, "minus_x")
        if _has_instance_attr(self, "detector_minus_area_text"):
            try:
                minus_spec = detector_active_area_from_values(self.values, "opposite_detector", "Second detector")
                minus_area_text = self._active_area_description(minus_spec)
            except ValueError as exc:
                minus_area_text = f"Invalid active area: {exc}"
            self.detector_minus_area_text.set(f"Active area: {minus_area_text}\nPortal: {minus_description}\nTerminal -x section: {minus_section.name}")
        if _has_instance_attr(self, "detector_minus_face_canvas"):
            self._draw_detector_face(self.detector_minus_face_canvas, "minus_x")
        if _has_instance_attr(self, "detector_mirror_text"):
            self.detector_mirror_text.set(
                "Specular reflection at the -x portal; photons can return into the active solid."
                if minus_available else f"Mirror configuration unavailable: {minus_description}."
            )
        if _has_instance_attr(self, "detector_mirror_canvas"):
            self._draw_mirror_portal()

    @staticmethod
    def _active_area_description(spec: Mapping[str, object]) -> str:
        shape = str(spec["shape"])
        if shape == "rectangular":
            dimensions = f"Y {float(spec['active_width_mm']):g} × Z {float(spec['active_height_mm']):g} mm"
        elif shape == "circular":
            dimensions = f"radius {float(spec['active_radius_mm']):g} mm"
        else:
            dimensions = f"semi-axes {float(spec['active_semi_axis_a_mm']):g} × {float(spec['active_semi_axis_b_mm']):g} mm"
        return f"{shape}, {dimensions}, offset ({float(spec['active_offset_y_mm']):g}, {float(spec['active_offset_z_mm']):g}) mm, rotation {float(spec['active_rotation_deg']):g}°"

    def _draw_detector_face(self, canvas: tk.Canvas, face: str = "plus_x") -> None:
        canvas.delete("all")
        section = self._terminal_editor_section(face)
        portal, width, height, profile, description = self._detector_portal_profile(section, face)
        if not portal:
            canvas.create_text(125, 55, text="No planar portal", fill="#8a1c1c")
            canvas.create_text(125, 78, text=description, fill="#555555")
            return
        width, height = max(width, 1.0), max(height, 1.0)
        scale = min(190 / width, 85 / height)
        x0, y0 = 125 - width * scale / 2, 65 - height * scale / 2
        x1, y1 = 125 + width * scale / 2, 65 + height * scale / 2
        if profile == "ellipse":
            canvas.create_oval(x0, y0, x1, y1, fill="#f4b183", outline="#7f5539", width=2)
        else:
            canvas.create_rectangle(x0, y0, x1, y1, fill="#f4b183", outline="#7f5539", width=2)
        canvas.create_text(125, 15, text=f"{face.replace('_', ' ')} detector portal", fill="#333333")
        prefix = "detector" if face == "plus_x" else "opposite_detector"
        try:
            points = detector_area_preview_geometry(self.values, prefix)
            screen_points = [(125 + y * scale, 65 - z * scale) for y, z in points]
            flat = [coordinate for point in screen_points for coordinate in point]
            if screen_points[0] != screen_points[-1]:
                flat += [screen_points[0][0], screen_points[0][1]]
            canvas.create_polygon(flat, fill="#5b9bd5", outline="#1f4e79", width=2, stipple="gray25")
            area_spec = detector_active_area_from_values(self.values, prefix)
            canvas.create_text(125, y1 + 31, text=self._active_area_description(area_spec), fill="#1f4e79", width=230)
        except ValueError as exc:
            canvas.create_text(125, y1 + 28, text=f"Invalid active area: {exc}", fill="#8a1c1c", width=230)
        canvas.create_text(125, y1 + 16, text=f"Portal Y = {width:g} mm", fill="#333333")
        canvas.create_text(x1 + 35, 65, text=f"Portal Z = {height:g} mm", fill="#333333", angle=90)

    def _draw_mirror_portal(self) -> None:
        canvas = _get_instance_attr(self, "detector_mirror_canvas", None)
        if canvas is None:
            return
        canvas.delete("all")
        section = self.geometry_sections[0] if self.geometry_sections else GeometrySection("Quartz bar")
        portal, width, height, profile, description = self._detector_portal_profile(section, "minus_x")
        if not portal:
            canvas.create_text(125, 45, text="No planar -x mirror portal", fill="#8a1c1c")
            return
        width, height = max(width, 1.0), max(height, 1.0)
        scale = min(190 / width, 55 / height)
        x0, y0 = 125 - width * scale / 2, 48 - height * scale / 2
        x1, y1 = 125 + width * scale / 2, 48 + height * scale / 2
        if profile == "ellipse":
            canvas.create_oval(x0, y0, x1, y1, fill="#d9d9d9", outline="#444444", width=2)
        else:
            canvas.create_rectangle(x0, y0, x1, y1, fill="#d9d9d9", outline="#444444", width=2)
        for offset in range(-80, 90, 14):
            canvas.create_line(125 + offset, y1, 125 + offset + 20, y0, fill="#777777")
        canvas.create_text(125, 12, text="-x mirror portal", fill="#333333")

    def _refresh_geometry_views(self) -> None:
        if _has_instance_attr(self, "geometry_canvases"):
            for key in self.geometry_canvases:
                self._refresh_geometry_view(key)
        self._update_detector_active_area()
        self._refresh_cherenkov_viewer()

    def _refresh_geometry_view(self, view: str) -> None:
        canvas = self.geometry_canvases.get(view)
        if canvas is not None and canvas.winfo_exists():
            self._draw_geometry_canvas(canvas, view)

    def _geometry_color(self, section: GeometrySection) -> str:
        return self.GEOMETRY_MATERIAL_COLORS.get(section.material.lower(), "#c5b0d5")

    @staticmethod
    def _convex_hull(points):
        """Return the monotonic-chain hull for a small projected point set."""

        unique = sorted(set((round(float(x), 8), round(float(y), 8)) for x, y in points))
        if len(unique) <= 2:
            return unique

        def cross(origin, first, second):
            return (first[0] - origin[0]) * (second[1] - origin[1]) - (first[1] - origin[1]) * (second[0] - origin[0])

        lower = []
        for point in unique:
            while len(lower) >= 2 and cross(lower[-2], lower[-1], point) <= 0:
                lower.pop()
            lower.append(point)
        upper = []
        for point in reversed(unique):
            while len(upper) >= 2 and cross(upper[-2], upper[-1], point) <= 0:
                upper.pop()
            upper.append(point)
        return lower[:-1] + upper[:-1]

    @staticmethod
    def _section_rings(section: GeometrySection, samples: int = 32):
        """Return sampled surface rings in mm for all editor primitives."""

        cx, cy, cz = section.center_x_mm, section.center_y_mm, section.center_z_mm
        if section.shape == "box":
            hx, hy, hz = section.length_x_mm / 2.0, section.width_y_mm / 2.0, section.thickness_z_mm / 2.0
            return [
                [(cx - hx, cy - hy, cz - hz), (cx + hx, cy - hy, cz - hz), (cx + hx, cy + hy, cz - hz), (cx - hx, cy + hy, cz - hz)],
                [(cx - hx, cy - hy, cz + hz), (cx + hx, cy - hy, cz + hz), (cx + hx, cy + hy, cz + hz), (cx - hx, cy + hy, cz + hz)],
            ]
        if section.shape == "tapered_rectangular_prism":
            x0, x1 = cx - section.length_x_mm / 2.0, cx + section.length_x_mm / 2.0
            def rectangular_ring(x, width, thickness):
                hy, hz = width / 2.0, thickness / 2.0
                return [(x, cy - hy, cz - hz), (x, cy + hy, cz - hz), (x, cy + hy, cz + hz), (x, cy - hy, cz + hz)]
            return [rectangular_ring(x0, section.width_y_start_mm, section.thickness_z_start_mm), rectangular_ring(x1, section.width_y_end_mm, section.thickness_z_end_mm)]

        if section.shape == "cylinder":
            radius_a = radius_b = section.radius_mm
        elif section.shape == "elliptical_cylinder":
            radius_a, radius_b = section.semi_axis_a_mm, section.semi_axis_b_mm
        else:
            radius_a = section.radius_start_mm
            radius_b = section.radius_end_mm
        half_length = section.length_mm / 2.0

        def ring_at(axis_position, first_radius, second_radius):
            ring = []
            for index in range(samples):
                theta = 2.0 * math.pi * index / samples
                u, v = first_radius * math.cos(theta), second_radius * math.sin(theta)
                if section.axis == "X":
                    ring.append((cx + axis_position, cy + u, cz + v))
                elif section.axis == "Y":
                    ring.append((cx + u, cy + axis_position, cz + v))
                else:
                    ring.append((cx + u, cy + v, cz + axis_position))
            return ring

        if section.shape == "truncated_cone":
            return [ring_at(-half_length, radius_a, radius_a), ring_at(half_length, radius_b, radius_b)]
        return [ring_at(-half_length, radius_a, radius_b), ring_at(half_length, radius_a, radius_b)]

    @staticmethod
    def _project_geometry_point(point, view: str):
        x, y, z = point
        if view == "xy":
            return (x, y)
        if view == "xz":
            return (x, z)
        if view == "yz":
            return (y, z)
        return (x - 0.65 * y, -0.55 * z + 0.35 * y)

    def _draw_geometry_canvas(self, canvas: tk.Canvas, view: str) -> None:
        canvas.delete("all")
        width = max(int(canvas.winfo_width()), 180)
        height = max(int(canvas.winfo_height()), 130)
        if not self.geometry_sections:
            canvas.create_text(width / 2, height / 2, text="No geometry sections")
            return
        if view == "iso":
            self._draw_isometric_geometry(canvas, width, height)
        else:
            self._draw_orthographic_geometry(canvas, width, height, view)
        axis_label = {"xy": "x / y (mm)", "xz": "x / z (mm)", "yz": "y / z (mm)", "iso": "x / y / z (mm)"}[view]
        canvas.create_text(6, height - 8, anchor="w", text=axis_label, fill="#555555")

    def _draw_orthographic_geometry(self, canvas: tk.Canvas, width: int, height: int, view: str) -> None:
        hulls = []
        all_points = []
        for section_index, section in enumerate(self.geometry_sections):
            projected = [self._project_geometry_point(point, view) for ring in self._section_rings(section) for point in ring]
            hull = self._convex_hull(projected)
            if hull:
                hulls.append((section_index, section, hull))
                all_points.extend(hull)
        if not all_points:
            return
        min_a, max_a = min(point[0] for point in all_points), max(point[0] for point in all_points)
        min_b, max_b = min(point[1] for point in all_points), max(point[1] for point in all_points)
        scale = min((width - 42) / max(max_a - min_a, 1.0), (height - 34) / max(max_b - min_b, 1.0))
        def map_point(point):
            return (21 + (point[0] - min_a) * scale, height - 22 - (point[1] - min_b) * scale)
        for section_index, section, hull in hulls:
            mapped = [map_point(point) for point in hull]
            selected = section_index == self.selected_geometry_index
            canvas.create_polygon(mapped, fill=self._geometry_color(section), outline="#1f4e79" if selected else "#6b7280", width=3 if selected else 1)
            center = map_point((sum(point[0] for point in hull) / len(hull), sum(point[1] for point in hull) / len(hull)))
            canvas.create_text(center[0], center[1], text=section.name, fill="#17324d")

    def _draw_isometric_geometry(self, canvas: tk.Canvas, width: int, height: int) -> None:
        rings_by_section = [(section, self._section_rings(section)) for section in self.geometry_sections]
        all_points = [self._project_geometry_point(point, "iso") for _, rings in rings_by_section for ring in rings for point in ring]
        if not all_points:
            return
        min_a, max_a = min(point[0] for point in all_points), max(point[0] for point in all_points)
        min_b, max_b = min(point[1] for point in all_points), max(point[1] for point in all_points)
        scale = min((width - 45) / max(max_a - min_a, 1.0), (height - 38) / max(max_b - min_b, 1.0))
        def map_point(point):
            return (22 + (point[0] - min_a) * scale, height - 25 - (point[1] - min_b) * scale)
        for item_index, (section, rings) in enumerate(rings_by_section):
            mapped_rings = [[map_point(self._project_geometry_point(point, "iso")) for point in ring] for ring in rings]
            selected = item_index == self.selected_geometry_index
            line_color = "#1f4e79" if selected else "#6b7280"
            line_width = 2 if selected else 1
            for ring in mapped_rings:
                coordinates = [coordinate for point in ring + [ring[0]] for coordinate in point]
                canvas.create_line(*coordinates, fill=line_color, width=line_width)
            if len(mapped_rings) == 2 and len(mapped_rings[0]) == len(mapped_rings[1]):
                stride = max(1, len(mapped_rings[0]) // 8)
                for index in range(0, len(mapped_rings[0]), stride):
                    canvas.create_line(*mapped_rings[0][index], *mapped_rings[1][index], fill=line_color, width=line_width)
            if mapped_rings:
                canvas.create_polygon(mapped_rings[0], fill=self._geometry_color(section), outline="", stipple="gray50")
            flat = [point for ring in mapped_rings for point in ring]
            center = (sum(point[0] for point in flat) / len(flat), sum(point[1] for point in flat) / len(flat))
            canvas.create_text(center[0], center[1], text=section.name, fill="#17324d")

    # Interactive geometry viewer -----------------------------------
    # This deliberately remains a Tk Canvas software renderer.  It does not
    # depend on an interactive Matplotlib backend (the application selects
    # the non-interactive Agg backend for its saved plots), and the bindings live on the Canvas
    # that actually receives the mouse events.

    @staticmethod
    def _geometry_viewer_color(color: str, factor: float) -> str:
        """Return a shaded Tk-compatible hex color."""

        value = str(color).lstrip("#")
        if len(value) != 6:
            value = "9bbbd4"
        try:
            channels = [int(value[offset:offset + 2], 16) for offset in (0, 2, 4)]
        except ValueError:
            channels = [155, 187, 212]
        channels = [max(0, min(255, int(channel * factor))) for channel in channels]
        return "#" + "".join(f"{channel:02x}" for channel in channels)

    def _canonical_local_view_data(self, assembly: str) -> dict:
        """Build renderer faces from one canonical complete assembly scene.

        ``complete_local_assembly_scene`` supplies all state-dependent
        components (active sections, enabled guide, mirror/readout interface,
        sensor, and Mylar metadata).  The sampled faces below are only a
        projection adapter for the Tk canvas; they are never an alternate
        geometry/model source.
        """

        model = self._assembly_model(assembly)
        guide_key = f"{assembly}_guide_enabled"
        model.guide_enabled = _bool_text(self.values, guide_key)
        scene = complete_local_assembly_scene(model, self.values, assembly=assembly)
        section_palette = ("#6fa8dc", "#e69138", "#93c47d", "#8e7cc3", "#e06666", "#76a5af", "#c9b458")
        component_by_name = {
            str(item.get("name")): item for item in scene.get("physical_components", scene.get("components", ()))
            if item.get("kind") in {"geometry", "guide"}
        }
        faces = []
        edges = []
        centers = []
        sensor_labels = []
        sensor_windows = []
        mirror_labels = []
        mylar_faces = []
        all_points = []

        for section_index, section in enumerate(model.active_sections()):
            component = component_by_name.get(section.name, {})
            section_faces = [list(face) for face in component.get("faces_mm", ()) if len(face) >= 3]
            if not section_faces:
                continue
            all_points.extend(point for face in section_faces for point in face)
            flat = [point for face in section_faces for point in face]
            centers.append((section_index, section, tuple(sum(point[axis] for point in flat) / len(flat) for axis in range(3))))
            base_color = component.get("color", section_palette[section_index % len(section_palette)])
            for face_index, face in enumerate(section_faces):
                shade = 0.82 if face_index == 0 else 1.08 if face_index == 1 else 0.96 + 0.04 * (face_index % 2)
                faces.append((sum(sum(point[axis] for point in face) for axis in range(3)) / (3 * len(face)), face, section_index, base_color, shade))
                edges.append((face + [face[0]], section_index))

        for overlay in scene.get("overlay_components", ()):
            for face in overlay.get("faces_mm", ()):
                face = list(face)
                mylar_faces.append((face, overlay.get("label", "Mylar")))
                all_points.extend(face)

        # Terminal mirrors are canonical scene components too.  Draw their
        # physical YZ plate so Mirror remains visible in the interactive view.
        for terminal in scene.get("terminal_interfaces", ()):
            if terminal.get("state") != "mirror":
                continue
            x = float(terminal["x_mm"])
            cy, cz = float(terminal["center_y_mm"]), float(terminal["center_z_mm"])
            hy, hz = float(terminal["width_y_mm"]) / 2.0, float(terminal["thickness_z_mm"]) / 2.0
            face = [(x, cy - hy, cz - hz), (x, cy + hy, cz - hz), (x, cy + hy, cz + hz), (x, cy - hy, cz + hz)]
            face_index = len(centers) + 1
            faces.append((sum(sum(point[axis] for point in face) for axis in range(3)) / (3 * len(face)), face, face_index, "#d62728", 1.0))
            edges.append((face + [face[0]], face_index))
            terminal_side = "+X" if terminal.get("face") == "plus_x" else "-X"
            mirror_labels.append((f"Mirror at {terminal_side}", (x, cy, cz)))
            all_points.extend(face)

        sensor_section_index = len(model.active_sections()) + 10
        for sensor_index, sensor in enumerate(item for item in scene.get("physical_components", scene["components"]) if item.get("kind") == "sensor" and item.get("solid_faces_mm")):
            sensor_faces = sensor["solid_faces_mm"]
            all_points.extend(point for face in sensor_faces for point in face)
            sensor_labels.append((sensor.get("label", "Photosensor"), tuple(sensor["position_mm"])))
            sensor_windows.extend(sensor.get("window_faces_mm", ()))
            sensor_face_index = sensor_section_index + sensor_index
            for face in sensor_faces:
                faces.append((sum(sum(point[axis] for point in face) for axis in range(3)) / (3 * len(face)), list(face), sensor_face_index, sensor.get("color", "#9b59b6"), 1.0))
                edges.append((list(face) + [face[0]], sensor_face_index))

        if not all_points:
            all_points = [(0.0, 0.0, 0.0)]
        lows = tuple(min(point[axis] for point in all_points) for axis in range(3))
        highs = tuple(max(point[axis] for point in all_points) for axis in range(3))
        center = tuple((lows[axis] + highs[axis]) / 2.0 for axis in range(3))
        extent = max(max(highs[axis] - lows[axis] for axis in range(3)), 1.0)
        return {
            "scene": scene, "faces": faces, "edges": edges, "centers": centers,
            "sensor_labels": sensor_labels, "sensor_windows": sensor_windows,
            "mirror_labels": mirror_labels,
            "center": center, "extent": extent, "mylar_face_ids": set(),
            "mylar_faces": mylar_faces,
        }

    def _geometry_viewer_scene(self):
        """Build the live canonical Cherenkov scene for the interactive view."""

        data = self._canonical_local_view_data("cherenkov")
        self.geometry_viewer_scene = data["scene"]
        return (data["faces"], data["edges"], data["centers"],
                data["sensor_labels"], data["sensor_windows"],
                data["center"], data["extent"], data.get("mylar_faces", ()),
                data["mirror_labels"])

    def _geometry_viewer_projection(self, point, center, scale, width, height):
        """Project one millimetre-space point and return (screen x/y, depth)."""

        x = float(point[0]) - center[0]
        y = float(point[1]) - center[1]
        z = float(point[2]) - center[2]
        cos_yaw, sin_yaw = math.cos(self.geometry_viewer_yaw), math.sin(self.geometry_viewer_yaw)
        cos_pitch, sin_pitch = math.cos(self.geometry_viewer_pitch), math.sin(self.geometry_viewer_pitch)
        rotated_x = cos_yaw * x - sin_yaw * y
        rotated_y = sin_yaw * x + cos_yaw * y
        screen_y = cos_pitch * rotated_y - sin_pitch * z
        depth = sin_pitch * rotated_y + cos_pitch * z
        return width / 2.0 + rotated_x * scale, height / 2.0 - screen_y * scale, depth

    def _draw_interactive_geometry_view(self) -> None:
        canvas = self.geometry_viewer_canvas
        if canvas is None or not canvas.winfo_exists():
            return
        canvas.delete("all")
        width = max(260, canvas.winfo_width())
        height = max(200, canvas.winfo_height())
        faces, edges, centers, sensor_labels, sensor_windows, center, extent, mylar_faces, mirror_labels = self._geometry_viewer_scene()
        scale = min(width, height) * 0.38 / extent * self.geometry_viewer_zoom

        def project(point):
            return self._geometry_viewer_projection(point, center, scale, width, height)

        # Draw a world-axis triad through the scene centre before the solid.
        axis_length = extent * 0.68
        axis_colors = (("X", "#c62828", (axis_length, 0.0, 0.0)), ("Y", "#2e7d32", (0.0, axis_length, 0.0)), ("Z", "#1565c0", (0.0, 0.0, axis_length)))
        for label, color, endpoint in axis_colors:
            start = project(center)
            end = project((center[0] + endpoint[0], center[1] + endpoint[1], center[2] + endpoint[2]))
            canvas.create_line(start[0], start[1], end[0], end[1], fill=color, width=2, arrow=TK_ARROW_LAST)
            canvas.create_text(end[0] + 8, end[1], text=label, fill=color, anchor="w", font=("TkDefaultFont", 9, "bold"))
        canvas.create_oval(start[0] - 2, start[1] - 2, start[0] + 2, start[1] + 2, fill="#444444", outline="")

        # Painter's algorithm: far faces first, then crisp outlines and labels.
        projected_faces = []
        for _average, face, section_index, color, shade in faces:
            points = [project(point) for point in face]
            projected_faces.append((sum(point[2] for point in points) / len(points), points, section_index, color, shade))
        for _depth, points, _section_index, color, shade in sorted(projected_faces, key=lambda item: item[0]):
            canvas.create_polygon([coordinate for point in points for coordinate in point[:2]], fill=self._geometry_viewer_color(color, shade), outline="")
        for edge, section_index in edges:
            points = [project(point) for point in edge]
            canvas.create_line([coordinate for point in points for coordinate in point[:2]], fill="#263238", width=2, smooth=False)
        for mylar_face, mylar_label in mylar_faces:
            points = [project(point) for point in mylar_face]
            canvas.create_polygon([coordinate for point in points for coordinate in point[:2]],
                                  fill=MYLAR_DISPLAY_COLOR, outline=MYLAR_DISPLAY_EDGE,
                                  stipple="gray25")
            canvas.create_line([coordinate for point in points + [points[0]] for coordinate in point[:2]],
                               fill=MYLAR_DISPLAY_EDGE, width=2, dash=(5, 3))
            canvas.create_text(points[0][0] + 3, points[0][1] - 3, text=mylar_label,
                               fill=MYLAR_DISPLAY_EDGE, anchor="sw", font=("TkDefaultFont", 8, "bold"))

        for section_index, section, section_center in centers:
            label_point = project(section_center)
            canvas.create_text(label_point[0], label_point[1], text=section.name, fill="#17324d", font=("TkDefaultFont", 9, "bold"))
        for label, sensor_center in sensor_labels:
            label_point = project(sensor_center)
            canvas.create_text(label_point[0], label_point[1], text=label, fill="#4a235a", font=("TkDefaultFont", 9, "bold"))
        for window in sensor_windows:
            points = [project(point) for point in window]
            canvas.create_polygon([coordinate for point in points for coordinate in point[:2]],
                                  fill="#e1bee7", outline="#4a235a", width=3)

        # A compact legend keeps section names and terminal interfaces visible
        # even when several components overlap after rotation.
        legend_x, legend_y = 10, 10
        legend_height = 20 + 18 * (len(centers) + len(sensor_labels) + len(mirror_labels))
        canvas.create_rectangle(legend_x, legend_y, min(width - 8, 245), legend_y + legend_height, fill="#ffffff", outline="#b0b0b0")
        canvas.create_text(legend_x + 7, legend_y + 5, anchor="nw", text="Local geometry / interfaces", fill="#333333", font=("TkDefaultFont", 9, "bold"))
        section_palette = ("#6fa8dc", "#e69138", "#93c47d", "#8e7cc3", "#e06666", "#76a5af", "#c9b458")
        for row, (section_index, section, _section_center) in enumerate(centers):
            y = legend_y + 22 + row * 18
            color = section_palette[section_index % len(section_palette)]
            canvas.create_rectangle(legend_x + 7, y, legend_x + 18, y + 11, fill=color, outline="#4a4a4a")
            shape_label = self._geometry_shape_label(section.shape)
            canvas.create_text(legend_x + 24, y + 5, anchor="w", text=f"{section.name} — {shape_label}", fill="#333333")
        for row, (label, _sensor_center) in enumerate(sensor_labels, start=len(centers)):
            y = legend_y + 22 + row * 18
            canvas.create_rectangle(legend_x + 7, y, legend_x + 18, y + 11, fill="#9b59b6", outline="#4a235a")
            canvas.create_text(legend_x + 24, y + 5, anchor="w", text=label, fill="#4a235a")
        for row, (label, _mirror_center) in enumerate(mirror_labels, start=len(centers) + len(sensor_labels)):
            y = legend_y + 22 + row * 18
            canvas.create_rectangle(legend_x + 7, y, legend_x + 18, y + 11, fill="#d62728", outline="#7f1717")
            canvas.create_text(legend_x + 24, y + 5, anchor="w", text=label, fill="#7f1717", font=("TkDefaultFont", 9, "bold"))
        canvas.create_text(8, height - 9, anchor="sw", text="Drag: rotate   •   Wheel: zoom   •   X/Y/Z: world axes", fill="#555555")

    def _reset_geometry_view(self) -> None:
        self.geometry_viewer_yaw = 0.75
        self.geometry_viewer_pitch = 0.45
        self.geometry_viewer_zoom = 1.0
        self.geometry_viewer_drag = None
        self._draw_interactive_geometry_view()

    def _geometry_viewer_press(self, event) -> None:
        self.geometry_viewer_drag = (event.x, event.y, self.geometry_viewer_yaw, self.geometry_viewer_pitch)
        # Explicitly focus the receiving widget so a drag remains associated
        # with this viewer even when the pointer entered from another widget.
        event.widget.focus_set()

    def _geometry_viewer_drag_motion(self, event) -> None:
        if self.geometry_viewer_drag is None:
            return
        start_x, start_y, start_yaw, start_pitch = self.geometry_viewer_drag
        self.geometry_viewer_yaw = start_yaw + (event.x - start_x) * 0.012
        self.geometry_viewer_pitch = max(-1.45, min(1.45, start_pitch + (event.y - start_y) * 0.012))
        self._draw_interactive_geometry_view()

    def _geometry_viewer_release(self, _event=None) -> None:
        self.geometry_viewer_drag = None

    def _geometry_viewer_zoom_event(self, event) -> None:
        delta = getattr(event, "delta", 0)
        if not delta:
            delta = 120 if getattr(event, "num", 4) == 4 else -120
        self.geometry_viewer_zoom *= 1.12 if delta > 0 else 1.0 / 1.12
        self.geometry_viewer_zoom = max(0.2, min(8.0, self.geometry_viewer_zoom))
        self._draw_interactive_geometry_view()

    def _close_geometry_3d_viewer(self) -> None:
        canvas = self.geometry_viewer_canvas
        if canvas is not None:
            for sequence, function_id in self.geometry_viewer_bindings:
                try:
                    canvas.unbind(sequence, function_id)
                except tk.TclError:
                    pass
        self.geometry_viewer_bindings = []
        self.geometry_viewer_drag = None
        viewer = self.geometry_viewer
        self.geometry_viewer_canvas = None
        self.geometry_viewer = None
        self.geometry_viewer_scene = None
        if viewer is not None and viewer.winfo_exists():
            viewer.destroy()

    def _open_geometry_3d_viewer(self) -> None:
        if self.geometry_viewer is not None and self.geometry_viewer.winfo_exists():
            self.geometry_viewer.deiconify()
            self.geometry_viewer.lift()
            self._draw_interactive_geometry_view()
            return
        viewer = tk.Toplevel(self)
        viewer.title("Geometry 3D Viewer")
        viewer.geometry("820x620")
        viewer.minsize(560, 420)
        viewer.columnconfigure(0, weight=1)
        viewer.rowconfigure(0, weight=1)
        canvas = tk.Canvas(viewer, background="white", highlightthickness=1, highlightbackground="#909090", takefocus=True)
        canvas.grid(row=0, column=0, sticky="nsew", padx=6, pady=6)
        controls = ttk.Frame(viewer)
        controls.grid(row=1, column=0, sticky="ew", padx=6, pady=(0, 6))
        ttk.Button(controls, text="Reset view", command=self._reset_geometry_view).pack(side="left")
        ttk.Label(controls, text="Drag with left mouse button to rotate • mouse wheel to zoom", foreground="#555555").pack(side="left", padx=(10, 0))
        ttk.Button(controls, text="Close", command=self._close_geometry_3d_viewer).pack(side="right")
        self.geometry_viewer = viewer
        self.geometry_viewer_canvas = canvas
        viewer.refresh = self._refresh_cherenkov_viewer
        self.geometry_viewer_yaw = 0.75
        self.geometry_viewer_pitch = 0.45
        self.geometry_viewer_zoom = 1.0
        self.geometry_viewer_drag = None
        bindings = (
            ("<Configure>", lambda _event: self._draw_interactive_geometry_view()),
            ("<ButtonPress-1>", self._geometry_viewer_press),
            ("<B1-Motion>", self._geometry_viewer_drag_motion),
            ("<ButtonRelease-1>", self._geometry_viewer_release),
            ("<MouseWheel>", self._geometry_viewer_zoom_event),
            ("<Button-4>", self._geometry_viewer_zoom_event),
            ("<Button-5>", self._geometry_viewer_zoom_event),
        )
        self.geometry_viewer_bindings = [(sequence, canvas.bind(sequence, callback)) for sequence, callback in bindings]
        viewer.protocol("WM_DELETE_WINDOW", self._close_geometry_3d_viewer)
        self._draw_interactive_geometry_view()


def main() -> int:
    if not TK_AVAILABLE:
        print(f"GUI unavailable because tkinter could not load its system library: {TK_IMPORT_ERROR}", file=sys.stderr)
        print("Install the Linux Tk runtime for your distribution, then run: python3 quartz_cherenkov_gui.py", file=sys.stderr)
        return 2
    try:
        app = QuartzCherenkovGUI()
    except Exception as exc:
        print(f"GUI could not start: {exc}", file=sys.stderr)
        return 2
    app.mainloop()
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
