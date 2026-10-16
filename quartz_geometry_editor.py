"""Tk editor for sequential rectangle/taper sections and detector areas.

The panel is intentionally usable without a running Tk display: the
``geometry_from_rows`` and ``geometry_to_rows`` helpers are the same model
conversion used by the GUI and by headless validation.  v23 also exposes the
canonical complete-assembly scene contract used by the local XZ/XY/YZ views;
optional guides follow the explicit enable-and-configured rule.  v30 keeps
the canonical physical photosensor solid/window in the same scene contract
used by all three local views and the interactive preview; global placement
is delegated to the shared model-layer rigid transform helpers. v31 also
exposes the display-only Mylar overlays derived from those canonical faces.
"""

from __future__ import annotations

import math
from dataclasses import replace
from pathlib import Path
from typing import Any, Dict, Iterable, List, Mapping, Optional, Sequence

from quartz_geometry import (
    ActiveDetectorArea,
    AssemblyGeometryModel,
    AssemblySection,
    GeometryModel,
    GeometrySection,
    MaterialSpec,
    SECTION_SHAPES,
    derive_mylar_overlay_surfaces,
    complete_detector_assembly_geometry_from_dict,
    complete_detector_assembly_geometry_to_dict,
)


MATERIAL_PRESETS = ("fused_silica", "pmma", "custom")


def canonical_assembly_model_from_sections(
    sections: Sequence[AssemblySection],
    *,
    name: str = "Assembly",
    guide_sections: Sequence[AssemblySection] = (),
    guide_enabled: bool = False,
) -> AssemblyGeometryModel:
    """Build the same v13 section model used by both detector GUI tabs."""

    enabled = str(guide_enabled).strip().lower() in {"1", "true", "yes", "on", "enabled"}
    model = AssemblyGeometryModel(
        name=name,
        sections=list(sections),
        guide_sections=list(guide_sections),
        guide_enabled=enabled,
    )
    model.validate()
    return model


def canonical_complete_assembly_scene(model: AssemblyGeometryModel, **kwargs) -> dict[str, Any]:
    """Return the shared terminal-aware scene for non-Tk callers."""

    return model.complete_scene(**kwargs)


def canonical_sensor_components(model: AssemblyGeometryModel, **kwargs) -> list[dict[str, Any]]:
    """Return physical sensor body/window components from the shared scene."""

    scene = canonical_complete_assembly_scene(model, **kwargs)
    return [dict(item) for item in scene.get("physical_components", ()) if item.get("kind") == "sensor"]


def canonical_mylar_overlay_components(model: AssemblyGeometryModel, **kwargs) -> list[dict[str, Any]]:
    """Return display-only Mylar surfaces derived from canonical faces."""

    return [dict(item) for item in derive_mylar_overlay_surfaces(model, **kwargs)]


def save_complete_detector_assembly(model: AssemblyGeometryModel, detector_type: str,
                                    photosensor: Mapping[str, Any], **kwargs) -> dict[str, Any]:
    """Editor-facing wrapper for the complete v32 assembly serializer."""

    return complete_detector_assembly_geometry_to_dict(model, detector_type, photosensor, **kwargs)


def load_complete_detector_assembly(payload: Mapping[str, Any], detector_type: str,
                                    *, name: Optional[str] = None):
    """Editor-facing side-effect-free complete v32 parser."""

    return complete_detector_assembly_geometry_from_dict(payload, detector_type, name=name)


def reorder_canonical_sections(sections: Sequence[AssemblySection], old_index: int, new_index: int) -> list[AssemblySection]:
    """Return a reordered copy for headless editor/tests and GUI parity."""

    result = list(sections)
    if not (0 <= old_index < len(result) and 0 <= new_index < len(result)):
        raise IndexError("section reorder index out of range")
    item = result.pop(old_index)
    result.insert(new_index, item)
    return result


def _detector_profile(area: ActiveDetectorArea, width_m: float, height_m: float, points: int = 64) -> List[tuple[float, float]]:
    """Return a closed y/z profile in metres for preview and renderer parity."""

    resolved_width, resolved_height, center_y, center_z = area.resolved(width_m, height_m)
    if area.shape == "circle" and area.mode == "custom":
        radius = area.radius_m if area.radius_m is not None else resolved_width / 2.0
        return [(center_y + radius * math.cos(angle), center_z + radius * math.sin(angle)) for angle in [2.0 * math.pi * i / points for i in range(points + 1)]]
    if area.shape == "ellipse" and area.mode == "custom":
        return [(center_y + resolved_width / 2.0 * math.cos(angle), center_z + resolved_height / 2.0 * math.sin(angle)) for angle in [2.0 * math.pi * i / points for i in range(points + 1)]]
    return [
        (center_y - resolved_width / 2.0, center_z - resolved_height / 2.0),
        (center_y + resolved_width / 2.0, center_z - resolved_height / 2.0),
        (center_y + resolved_width / 2.0, center_z + resolved_height / 2.0),
        (center_y - resolved_width / 2.0, center_z + resolved_height / 2.0),
        (center_y - resolved_width / 2.0, center_z - resolved_height / 2.0),
    ]


def _draw_detector_preview(axis, geometry: GeometryModel, face: str, color: str) -> None:
    from mpl_toolkits.mplot3d.art3d import Poly3DCollection

    area = geometry.active_area(face)
    if not area.enabled:
        return
    face_width, face_height = geometry.face_dimensions(face)
    profile = _detector_profile(area, face_width, face_height)
    x_face = geometry.half_x if face == "plus_x" else -geometry.half_x
    direction = 1.0 if face == "plus_x" else -1.0
    depth = max(0.0005, min(0.003, geometry.length_x_m * 0.03))
    x0 = x_face + direction * 0.00005
    x1 = x_face + direction * depth
    front = [(x0 * 1e3, y * 1e3, z * 1e3) for y, z in profile]
    back = [(x1 * 1e3, y * 1e3, z * 1e3) for y, z in profile]
    polygons = [front, back]
    polygons.extend([[front[index], front[index + 1], back[index + 1], back[index]] for index in range(len(profile) - 1)])
    label = f"{face} detector ({area.shape}, {area.mode})"
    axis.add_collection3d(Poly3DCollection(polygons, facecolors=color, edgecolors=color, linewidths=1.6, alpha=0.34, label=label))


MATERIAL_COLORS = {"fused_silica": "#80b1d3", "pmma": "#fdb462", "custom": "#fb8072"}


def _geometry_extents(geometry: GeometryModel) -> tuple[float, float, float]:
    y_extent = max(section.dimensions_at(fraction)[0] / 2.0 for section in geometry.sections for fraction in (0.0, 1.0))
    z_extent = max(section.dimensions_at(fraction)[1] / 2.0 for section in geometry.sections for fraction in (0.0, 1.0))
    for face in ("plus_x", "minus_x"):
        area = geometry.active_area(face)
        if area.enabled:
            fw, fh = geometry.face_dimensions(face)
            aw, ah, cy, cz = area.resolved(fw, fh)
            y_extent = max(y_extent, abs(cy) + aw / 2.0)
            z_extent = max(z_extent, abs(cz) + ah / 2.0)
    return geometry.half_x, y_extent, z_extent


def _draw_orthographic_views(axes, geometry: GeometryModel) -> None:
    """Draw synchronized x-z, x-y, and y-z projections in millimetres.

    The YZ view is an end-face projection, so it is intentionally built from
    the dimensions at each real terminal (the first section at fraction 0 and
    the final section at fraction 1), rather than from nominal global fields.
    """

    from matplotlib.patches import Ellipse, Polygon, Rectangle

    x_half, y_extent, z_extent = _geometry_extents(geometry)
    axes_xz, axes_xy, axes_yz = axes
    for index, section in enumerate(geometry.sections):
        x0, x1 = geometry.x_edges[index] * 1e3, geometry.x_edges[index + 1] * 1e3
        w0, t0 = section.dimensions_at(0.0)
        w1, t1 = section.dimensions_at(1.0)
        color = MATERIAL_COLORS[section.material.preset]
        axes_xz.add_patch(Polygon([(x0, -t0 * 500), (x1, -t1 * 500), (x1, t1 * 500), (x0, t0 * 500)], closed=True, facecolor=color, edgecolor="#4b6270", alpha=0.38, linewidth=0.9))
        axes_xy.add_patch(Polygon([(x0, -w0 * 500), (x1, -w1 * 500), (x1, w1 * 500), (x0, w0 * 500)], closed=True, facecolor=color, edgecolor="#4b6270", alpha=0.38, linewidth=0.9))
        axes_xz.text((x0 + x1) / 2.0, 0.0, str(index + 1), ha="center", va="center", fontsize=7, color="#30434f")
        axes_xy.text((x0 + x1) / 2.0, 0.0, str(index + 1), ha="center", va="center", fontsize=7, color="#30434f")

    for face, color, linestyle in (("plus_x", "#1b9e77", "-"), ("minus_x", "#984ea3", "--")):
        x = geometry.half_x * 1e3 if face == "plus_x" else -geometry.half_x * 1e3
        fw, fh = geometry.face_dimensions(face)
        axes_xz.plot([x, x], [-fh * 500, fh * 500], color=color, linestyle=linestyle, linewidth=1.5)
        axes_xy.plot([x, x], [-fw * 500, fw * 500], color=color, linestyle=linestyle, linewidth=1.5)
        area = geometry.active_area(face)
        if area.enabled:
            aw, ah, cy, cz = area.resolved(fw, fh)
            axes_xz.plot([x, x], [(cz - ah / 2) * 1e3, (cz + ah / 2) * 1e3], color=color, linewidth=3.0)
            axes_xy.plot([x, x], [(cy - aw / 2) * 1e3, (cy + aw / 2) * 1e3], color=color, linewidth=3.0)

    # Always add both terminal profiles.  This is the fallback geometry even
    # when a detector area is disabled, and prevents an empty YZ axes object.
    terminal_profiles = (
        ("-x end", "minus_x", geometry.sections[0].dimensions_at(0.0), "#984ea3", "--"),
        ("+x end", "plus_x", geometry.sections[-1].dimensions_at(1.0), "#1b9e77", "-"),
    )
    yz_y_extent = max(width / 2.0 for _label, _face, (width, _height), _color, _style in terminal_profiles)
    yz_z_extent = max(height / 2.0 for _label, _face, (_width, height), _color, _style in terminal_profiles)
    for label, face, (width, thickness), color, linestyle in terminal_profiles:
        axes_yz.add_patch(Rectangle((-width * 500, -thickness * 500), width * 1e3, thickness * 1e3, fill=False, edgecolor=color, linestyle=linestyle, linewidth=1.8, label=label))
        area = geometry.active_area(face)
        if not area.enabled:
            continue
        fw, fh = geometry.face_dimensions(face)
        aw, ah, cy, cz = area.resolved(fw, fh)
        if area.shape in {"circle", "ellipse"} and area.mode == "custom":
            axes_yz.add_patch(Ellipse((cy * 1e3, cz * 1e3), aw * 1e3, ah * 1e3, fill=False, edgecolor=color, linewidth=2.5, label=f"{face} active area"))
        else:
            axes_yz.add_patch(Rectangle(((cy - aw / 2) * 1e3, (cz - ah / 2) * 1e3), aw * 1e3, ah * 1e3, fill=False, edgecolor=color, linewidth=2.5, label=f"{face} active area"))
        yz_y_extent = max(yz_y_extent, abs(cy) + aw / 2.0)
        yz_z_extent = max(yz_z_extent, abs(cz) + ah / 2.0)

    for axis in (axes_xz, axes_xy):
        axis.set_xlim(-x_half * 1.08e3, x_half * 1.08e3)
        axis.set_aspect("equal", adjustable="box")
        axis.grid(True, alpha=0.22)
    axes_xz.set_ylim(-z_extent * 1.12e3, z_extent * 1.12e3)
    axes_xy.set_ylim(-y_extent * 1.12e3, y_extent * 1.12e3)
    # Explicit non-degenerate limits keep a valid terminal rectangle visible
    # even when matplotlib has not received any detector overlay.
    axes_yz.set_xlim(-max(yz_y_extent * 1.12, 1e-9) * 1e3, max(yz_y_extent * 1.12, 1e-9) * 1e3)
    axes_yz.set_ylim(-max(yz_z_extent * 1.12, 1e-9) * 1e3, max(yz_z_extent * 1.12, 1e-9) * 1e3)
    axes_yz.set_aspect("equal", adjustable="box")
    axes_yz.grid(True, alpha=0.22)
    axes_xz.set_title("Longitudinal x-z")
    axes_xz.set_xlabel("x (mm)")
    axes_xz.set_ylabel("z (mm)")
    axes_xy.set_title("Plan x-y")
    axes_xy.set_xlabel("x (mm)")
    axes_xy.set_ylabel("y (mm)")
    axes_yz.set_title("End y-z")
    axes_yz.set_xlabel("y (mm)")
    axes_yz.set_ylabel("z (mm)")
    from matplotlib.patches import Patch
    materials = []
    for section in geometry.sections:
        if section.material.preset not in materials:
            materials.append(section.material.preset)
    axes_xz.legend(handles=[Patch(facecolor=MATERIAL_COLORS[name], alpha=0.38, label=name.replace("_", " ")) for name in materials], fontsize=6, loc="upper left")
    axes_yz.legend(fontsize=6, loc="upper left")


def _draw_3d_geometry(axis, geometry: GeometryModel) -> None:
    from mpl_toolkits.mplot3d.art3d import Poly3DCollection

    for index, section in enumerate(geometry.sections):
        x0, x1 = geometry.x_edges[index] * 1e3, geometry.x_edges[index + 1] * 1e3
        w0, t0 = section.dimensions_at(0.0)
        w1, t1 = section.dimensions_at(1.0)
        vertices = [(x0, -w0 * 500, -t0 * 500), (x0, w0 * 500, -t0 * 500), (x0, w0 * 500, t0 * 500), (x0, -w0 * 500, t0 * 500), (x1, -w1 * 500, -t1 * 500), (x1, w1 * 500, -t1 * 500), (x1, w1 * 500, t1 * 500), (x1, -w1 * 500, t1 * 500)]
        faces = ((0, 1, 2, 3), (4, 7, 6, 5), (0, 4, 5, 1), (3, 2, 6, 7), (1, 5, 6, 2), (0, 3, 7, 4))
        axis.add_collection3d(Poly3DCollection([[vertices[i] for i in face] for face in faces], facecolors=MATERIAL_COLORS[section.material.preset], alpha=.24, edgecolors="#4b6270", linewidths=.6, label=f"{section.name}: {section.shape} / {section.material.preset}"))
    for face, color in (("plus_x", "#1b9e77"), ("minus_x", "#984ea3")):
        _draw_detector_preview(axis, geometry, face, color)
    x_half, y_extent, z_extent = _geometry_extents(geometry)
    x_pad = max(0.0005, min(0.003, geometry.length_x_m * 0.03))
    axis.set_xlim((-x_half - x_pad) * 1e3, (x_half + x_pad) * 1e3)
    axis.set_ylim(-y_extent * 1.08e3, y_extent * 1.08e3)
    axis.set_zlim(-z_extent * 1.08e3, z_extent * 1.08e3)
    axis.set_box_aspect((geometry.length_x_m + 2 * x_pad, 2 * y_extent, 2 * z_extent))
    axis.set_xlabel("x (mm)")
    axis.set_ylabel("y (mm)")
    axis.set_zlabel("z (mm)")
    axis.set_title("3D geometry")
    axis.text2D(0.02, 0.02, "Sections: " + " | ".join(section.shape for section in geometry.sections), transform=axis.transAxes, fontsize=7)


def _geometry_preview_figure(geometry: GeometryModel):
    """Build the 3D view and all three synchronized orthographic views."""

    import matplotlib.pyplot as plt
    figure = plt.figure(figsize=(12, 8))
    grid = figure.add_gridspec(2, 2, wspace=0.28, hspace=0.30)
    axes_xz = figure.add_subplot(grid[0, 0])
    axes_xy = figure.add_subplot(grid[0, 1])
    axes_yz = figure.add_subplot(grid[1, 0])
    axes_3d = figure.add_subplot(grid[1, 1], projection="3d")
    _draw_orthographic_views((axes_xz, axes_xy, axes_yz), geometry)
    _draw_3d_geometry(axes_3d, geometry)
    figure.suptitle("Quartz geometry: 3D and synchronized orthographic views", fontsize=12)
    return figure


def plot_geometry_preview(geometry: GeometryModel, output_path: str | Path) -> None:
    """Save a material-colored 3D preview with visible shaped detector plates."""

    import matplotlib.pyplot as plt
    figure = _geometry_preview_figure(geometry)
    figure.savefig(output_path, dpi=160)
    plt.close(figure)


def geometry_to_rows(geometry: GeometryModel) -> List[Dict[str, Any]]:
    return [
        {
            "name": section.name,
            "shape": section.shape,
            "length_mm": section.length_m * 1e3,
            "width_start_mm": section.width_start_m * 1e3,
            "width_end_mm": (section.width_end_m or section.width_start_m) * 1e3,
            "thickness_start_mm": section.thickness_start_m * 1e3,
            "thickness_end_mm": (section.thickness_end_m or section.thickness_start_m) * 1e3,
            "material": section.material.preset,
            "custom_n": section.material.custom_n if section.material.preset == "custom" else "",
            "effective_n_450nm": section.material.effective_index(450e-9),
        }
        for section in geometry.sections
    ]


def geometry_from_rows(rows: Sequence[Mapping[str, Any]], base: Optional[GeometryModel] = None) -> GeometryModel:
    sections = []
    for index, row in enumerate(rows):
        def number(key: str, fallback: float) -> float:
            value = row.get(key, fallback)
            value = fallback if value in (None, "") else float(value)
            if not math.isfinite(value):
                raise ValueError(f"section {index + 1} {key} must be finite")
            return value
        material_name = str(row.get("material", "fused_silica")).strip().lower()
        custom = row.get("custom_n")
        material = MaterialSpec(material_name, None if custom in (None, "") else float(custom))
        width_start = number("width_start_mm", 20.0)
        width_end = number("width_end_mm", width_start)
        thickness_start = number("thickness_start_mm", 10.0)
        thickness_end = number("thickness_end_mm", thickness_start)
        raw_shape = row.get("shape", row.get("type"))
        if raw_shape in (None, ""):
            shape = "taper" if (not math.isclose(width_start, width_end) or not math.isclose(thickness_start, thickness_end)) else "rectangle"
        else:
            shape = str(raw_shape).strip().lower()
        if shape not in SECTION_SHAPES:
            raise ValueError(f"section {index + 1} shape must be one of: {', '.join(SECTION_SHAPES)}")
        if shape == "rectangle":
            width_end = width_start
            thickness_end = thickness_start
        sections.append(GeometrySection(
            length_m=number("length_mm", 100.0) * 1e-3,
            width_start_m=width_start * 1e-3,
            width_end_m=width_end * 1e-3,
            thickness_start_m=thickness_start * 1e-3,
            thickness_end_m=thickness_end * 1e-3,
            material=material,
            name=str(row.get("name", f"section_{index + 1}")),
            shape=shape,
        ))
    if not sections:
        raise ValueError("at least one geometry section is required")
    base = base or GeometryModel()
    return GeometryModel(sum(section.length_m for section in sections), base.width_y_m, base.thickness_z_m, tuple(sections), base.detector_plus_x, base.detector_minus_x)


def detector_area_from_values(values: Mapping[str, Any], prefix: str, default_face_width_mm: float, default_face_height_mm: float) -> ActiveDetectorArea:
    mode = str(values.get(f"{prefix}_mode", "match_exit_surface"))
    shape = str(values.get(f"{prefix}_shape", "full_exit_surface" if mode == "match_exit_surface" else "rectangle"))
    width = values.get(f"{prefix}_width_mm", values.get(f"{prefix}_width"))
    height = values.get(f"{prefix}_height_mm", values.get(f"{prefix}_height"))
    radius = values.get(f"{prefix}_radius_mm", values.get(f"{prefix}_radius"))
    semi_a = values.get(f"{prefix}_semi_axis_a_mm", values.get(f"{prefix}_semi_axis_a"))
    semi_b = values.get(f"{prefix}_semi_axis_b_mm", values.get(f"{prefix}_semi_axis_b"))
    raw_enabled = values.get(f"{prefix}_enabled", True)
    if isinstance(raw_enabled, str):
        enabled = raw_enabled.strip().lower() not in {"", "0", "false", "no", "off"}
    else:
        enabled = bool(raw_enabled)
    return ActiveDetectorArea(
        enabled=enabled,
        mode=mode,
        width_m=None if width in (None, "") else float(width) * 1e-3,
        height_m=None if height in (None, "") else float(height) * 1e-3,
        center_y_m=float(values.get(f"{prefix}_center_y_mm", values.get(f"{prefix}_center_y", 0.0)) or 0.0) * 1e-3,
        center_z_m=float(values.get(f"{prefix}_center_z_mm", values.get(f"{prefix}_center_z", 0.0)) or 0.0) * 1e-3,
        shape=shape,
        radius_m=None if radius in (None, "") else float(radius) * 1e-3,
        semi_axis_a_m=None if semi_a in (None, "") else float(semi_a) * 1e-3,
        semi_axis_b_m=None if semi_b in (None, "") else float(semi_b) * 1e-3,
    )


if __name__ != "__main__":
    try:
        import tkinter as tk
        from tkinter import filedialog, messagebox, ttk
    except ImportError:  # pragma: no cover - exercised only on minimal Tk installs
        tk = None
        ttk = None
        filedialog = None
        messagebox = None
else:
    tk = None
    ttk = None
    filedialog = None
    messagebox = None


if ttk is not None:
    class GeometryEditorPanel(ttk.Frame):
        """Scrollable per-section editor embedded in the main GUI."""

        def __init__(self, master: Any, geometry: Optional[GeometryModel] = None, **kwargs: Any) -> None:
            super().__init__(master, **kwargs)
            self.geometry_model = geometry or GeometryModel()
            self.rows: List[Dict[str, Any]] = geometry_to_rows(self.geometry_model)
            self.vars: List[Dict[str, Any]] = []
            self._build()
            self.set_geometry(self.geometry_model)

        def _build(self) -> None:
            self.columnconfigure(0, weight=1)
            self.rowconfigure(0, weight=0)
            self.rowconfigure(1, weight=1)
            self.rowconfigure(2, weight=0, minsize=465)

            # The toolbar and live preview stay fixed while only the editor
            # controls scroll.  This keeps shape selectors discoverable and
            # prevents long section lists from pushing the preview away.
            toolbar = ttk.Frame(self)
            toolbar.grid(row=0, column=0, sticky="ew", pady=(0, 3))
            toolbar.columnconfigure(8, weight=1)
            ttk.Button(toolbar, text="Add section", command=self.add_section).grid(row=0, column=0, padx=(0, 3))
            ttk.Button(toolbar, text="Remove last", command=self.remove_section).grid(row=0, column=1, padx=3)
            ttk.Button(toolbar, text="Load Geometry", command=self.load_geometry).grid(row=0, column=2, padx=3)
            ttk.Button(toolbar, text="Save Geometry", command=self.save_geometry).grid(row=0, column=3, padx=3)
            ttk.Button(toolbar, text="Preview PNG", command=self.preview_png).grid(row=0, column=4, padx=3)
            ttk.Button(toolbar, text="Preview 2D", command=self._refresh_preview).grid(row=0, column=5, padx=(12, 3))
            ttk.Button(toolbar, text="Open 3D Viewer", command=self.open_3d_viewer).grid(row=0, column=6, padx=3)
            ttk.Label(toolbar, text="Scroll the editor below for geometry pieces and detector settings", foreground="#4c4c4c").grid(row=0, column=8, sticky="e", padx=(12, 0))

            # The complete controls area is one scrollable canvas.  The
            # detector panels are in this same content frame, so every field
            # remains reachable after adding an arbitrary number of sections.
            scroll_frame = ttk.Frame(self)
            scroll_frame.grid(row=1, column=0, sticky="nsew")
            scroll_frame.columnconfigure(0, weight=1)
            scroll_frame.rowconfigure(0, weight=1)
            self.canvas = tk.Canvas(scroll_frame, highlightthickness=0)
            self.canvas.grid(row=0, column=0, sticky="nsew")
            scrollbar = ttk.Scrollbar(scroll_frame, orient="vertical", command=self.canvas.yview)
            scrollbar.grid(row=0, column=1, sticky="ns")
            self.canvas.configure(yscrollcommand=scrollbar.set)
            self.content = ttk.Frame(self.canvas)
            self.content.columnconfigure(0, weight=1)
            self.window_id = self.canvas.create_window((0, 0), window=self.content, anchor="nw")
            self.content.bind("<Configure>", lambda _event: self.canvas.configure(scrollregion=self.canvas.bbox("all")))
            self.canvas.bind("<Configure>", lambda event: self.canvas.itemconfigure(self.window_id, width=event.width))
            self.canvas.bind_all("<MouseWheel>", self._mousewheel, add="+")
            self.canvas.bind_all("<Button-4>", self._mousewheel, add="+")
            self.canvas.bind_all("<Button-5>", self._mousewheel, add="+")

            self.row_frame = ttk.LabelFrame(self.content, text="Geometry-piece sections (rectangle / taper / segment)", padding=4)
            self.row_frame.grid(row=0, column=0, sticky="ew", pady=4)
            self.render_rows()

            detector = ttk.LabelFrame(self.content, text="Detector active areas: independent +x and -x terminal faces")
            detector.grid(row=1, column=0, sticky="ew", pady=(3, 0))
            detector.columnconfigure(0, weight=1)
            detector.columnconfigure(1, weight=1)
            self.detector_vars: Dict[str, Any] = {}
            self.detector_widgets: Dict[str, Dict[str, Any]] = {}
            for column, (prefix, title) in enumerate((("plus_x", "+x detector active area"), ("minus_x", "-x detector active area"))):
                face = ttk.LabelFrame(detector, text=title, padding=4)
                face.grid(row=0, column=column, padx=(0 if column == 0 else 4, 0), sticky="nsew")
                face.columnconfigure(1, weight=1)
                self.detector_widgets[prefix] = {}
                values = {"mode": "match_exit_surface", "shape": "full_exit_surface", "width_mm": "", "height_mm": "", "radius_mm": "", "semi_axis_a_mm": "", "semi_axis_b_mm": "", "center_y_mm": "0", "center_z_mm": "0"}
                self.detector_vars[f"{prefix}_enabled"] = tk.BooleanVar(value=True)
                ttk.Checkbutton(face, text="Enabled", variable=self.detector_vars[f"{prefix}_enabled"]).grid(row=0, column=0, columnspan=2, sticky="w", pady=(0, 3))
                for key, value in values.items():
                    variable = tk.StringVar(value=value)
                    self.detector_vars[f"{prefix}_{key}"] = variable
                    label = {
                        "mode": "Mode",
                        "shape": "Active-area shape",
                        "width_mm": "Width (mm)",
                        "height_mm": "Height (mm)",
                        "radius_mm": "Radius (mm)",
                        "semi_axis_a_mm": "Semi-axis A (mm)",
                        "semi_axis_b_mm": "Semi-axis B (mm)",
                        "center_y_mm": "Center y offset (mm)",
                        "center_z_mm": "Center z offset (mm)",
                    }[key]
                    ttk.Label(face, text=label).grid(row=list(values).index(key) + 1, column=0, padx=(0, 5), pady=2, sticky="w")
                    if key == "mode":
                        widget = ttk.Combobox(face, textvariable=variable, values=("match_exit_surface", "custom"), state="readonly", width=18)
                    elif key == "shape":
                        widget = ttk.Combobox(face, textvariable=variable, values=("rectangle", "circle", "ellipse", "full_exit_surface"), state="readonly", width=18)
                    else:
                        widget = ttk.Entry(face, textvariable=variable, width=12)
                    widget.grid(row=list(values).index(key) + 1, column=1, padx=2, pady=2, sticky="ew")
                    self.detector_widgets[prefix][key] = widget
                    variable.trace_add("write", self._detector_value_changed)
                self.detector_vars[f"{prefix}_enabled"].trace_add("write", self._detector_value_changed)
                self.detector_widgets[prefix]["shape"].bind("<<ComboboxSelected>>", lambda _event, item=prefix: self._shape_selected(item))

            ttk.Label(
                self.content,
                text="Geometry-piece shape and detector active-area shape are separate settings. Mode match_exit_surface uses the complete terminal face; custom shapes are rectangle, circle, ellipse, or full_exit_surface. Circle uses radius; ellipse uses semiaxes A/B. Center y/z offsets are local face coordinates.",
                foreground="#245a35", wraplength=1100, justify="left",
            ).grid(row=2, column=0, sticky="w", pady=3)

            # This fixed-height panel is outside the scrolling editor and is
            # therefore visible even when many geometry sections are present.
            preview_frame = ttk.LabelFrame(self, text="Live geometry preview", padding=4)
            preview_frame.grid(row=2, column=0, sticky="nsew", pady=(4, 0))
            preview_frame.grid_propagate(False)
            preview_frame.columnconfigure(0, weight=1)
            preview_frame.rowconfigure(1, weight=1)
            ttk.Label(preview_frame, text="All views update from the same current model: longitudinal x-z, plan x-y, end y-z, and 3D. Axes are in mm.", foreground="#245a35").grid(row=0, column=0, sticky="w", pady=(0, 3))
            self.preview_host = ttk.Frame(preview_frame)
            self.preview_host.grid(row=1, column=0, sticky="nsew")
            self.preview_figure = None
            self.preview_mpl_canvas = None
            self.preview_section_summary = tk.StringVar(value="Sections: rectangle")
            ttk.Label(preview_frame, textvariable=self.preview_section_summary, foreground="#4c4c4c").grid(row=2, column=0, sticky="w", pady=(3, 0))
            self._refresh_preview()

        def _mousewheel(self, event):
            widget = self.winfo_containing(event.x_root, event.y_root)
            widget_path = "" if widget is None else str(widget)
            if widget_path != str(self) and not widget_path.startswith(str(self) + "."):
                return
            if event.num == 4:
                self.canvas.yview_scroll(-1, "units")
            elif event.num == 5:
                self.canvas.yview_scroll(1, "units")
            elif getattr(event, "delta", 0):
                self.canvas.yview_scroll(int(-event.delta / 120), "units")

        def _detector_value_changed(self, *_args) -> None:
            self._refresh_preview()

        def _shape_selected(self, prefix: str) -> None:
            shape = self.detector_vars[f"{prefix}_shape"].get()
            # A concrete shape is necessarily a custom active area.  Selecting
            # the full terminal face switches back to the explicit full-face
            # mode, so the selector always has an immediate visual effect.
            self.detector_vars[f"{prefix}_mode"].set("match_exit_surface" if shape == "full_exit_surface" else "custom")
            face_width, face_height = self.geometry_model.face_dimensions("plus_x" if prefix == "plus_x" else "minus_x")
            if shape in {"rectangle", "ellipse"}:
                if not self.detector_vars[f"{prefix}_width_mm"].get().strip():
                    self.detector_vars[f"{prefix}_width_mm"].set(f"{face_width * 1e3:g}")
                if not self.detector_vars[f"{prefix}_height_mm"].get().strip():
                    self.detector_vars[f"{prefix}_height_mm"].set(f"{face_height * 1e3:g}")
            if shape == "ellipse":
                if not self.detector_vars[f"{prefix}_semi_axis_a_mm"].get().strip():
                    self.detector_vars[f"{prefix}_semi_axis_a_mm"].set(f"{face_width * 500:g}")
                if not self.detector_vars[f"{prefix}_semi_axis_b_mm"].get().strip():
                    self.detector_vars[f"{prefix}_semi_axis_b_mm"].set(f"{face_height * 500:g}")
            elif shape == "circle" and not self.detector_vars[f"{prefix}_radius_mm"].get().strip():
                self.detector_vars[f"{prefix}_radius_mm"].set(f"{min(face_width, face_height) * 500:g}")
            self._refresh_preview()

        def _preview_geometry(self) -> GeometryModel:
            """Build a forgiving preview model while an entry is being edited."""

            try:
                return self.get_geometry()
            except (TypeError, ValueError):
                geometry = geometry_from_rows(self.rows, self.geometry_model)
                values = {key: variable.get() for key, variable in self.detector_vars.items()}
                areas = {}
                for prefix, face in (("plus_x", "plus_x"), ("minus_x", "minus_x")):
                    width, height = geometry.face_dimensions(face)
                    shape = values.get(f"{prefix}_shape", "full_exit_surface")
                    values.setdefault(f"{prefix}_width_mm", "")
                    values.setdefault(f"{prefix}_height_mm", "")
                    values.setdefault(f"{prefix}_radius_mm", "")
                    values.setdefault(f"{prefix}_semi_axis_a_mm", "")
                    values.setdefault(f"{prefix}_semi_axis_b_mm", "")
                    if shape in {"rectangle", "ellipse"} and not values[f"{prefix}_width_mm"].strip():
                        values[f"{prefix}_width_mm"] = str(width * 1e3)
                    if shape in {"rectangle", "ellipse"} and not values[f"{prefix}_height_mm"].strip():
                        values[f"{prefix}_height_mm"] = str(height * 1e3)
                    if shape == "circle" and not values[f"{prefix}_radius_mm"].strip():
                        values[f"{prefix}_radius_mm"] = str(min(width, height) * 500.0)
                    if shape == "ellipse":
                        values[f"{prefix}_semi_axis_a_mm"] = values[f"{prefix}_semi_axis_a_mm"].strip() or str(width * 500.0)
                        values[f"{prefix}_semi_axis_b_mm"] = values[f"{prefix}_semi_axis_b_mm"].strip() or str(height * 500.0)
                    try:
                        areas[face] = detector_area_from_values(values, prefix, width * 1e3, height * 1e3)
                    except (TypeError, ValueError):
                        areas[face] = geometry.active_area(face)
                return replace(geometry, detector_plus_x=areas["plus_x"], detector_minus_x=areas["minus_x"])

        def _refresh_preview(self, *_args) -> None:
            if not hasattr(self, "preview_host") or not self.preview_host.winfo_exists():
                return
            geometry = self._preview_geometry()
            try:
                import matplotlib.pyplot as plt
                from matplotlib.backends.backend_tkagg import FigureCanvasTkAgg
                if self.preview_mpl_canvas is not None:
                    self.preview_mpl_canvas.get_tk_widget().destroy()
                if self.preview_figure is not None:
                    plt.close(self.preview_figure)
                self.preview_figure = _geometry_preview_figure(geometry)
                self.preview_mpl_canvas = FigureCanvasTkAgg(self.preview_figure, master=self.preview_host)
                self.preview_mpl_canvas.draw()
                self.preview_mpl_canvas.get_tk_widget().pack(fill="both", expand=True)
            except Exception:
                # Keep edits usable if a local Tk/Matplotlib backend cannot be
                # embedded. The Preview PNG and dedicated preview window still
                # use the same figure builder and report their own errors.
                self.preview_figure = None
                self.preview_mpl_canvas = None
            if hasattr(self, "preview_section_summary"):
                self.preview_section_summary.set("Geometry pieces: " + " | ".join(f"{index + 1}: {section.shape} / {section.material.preset}" for index, section in enumerate(geometry.sections)))

        def open_3d_viewer(self) -> None:
            """Open the current geometry in a dedicated interactive 3D viewer."""

            if getattr(self, "_preview_3d_window", None) is not None:
                try:
                    if self._preview_3d_window.winfo_exists():
                        self._preview_3d_window.lift()
                        self._preview_3d_window.focus_force()
                        return
                except tk.TclError:
                    self._preview_3d_window = None
            try:
                geometry = self._preview_geometry()
                import matplotlib.pyplot as plt
                from matplotlib.backends.backend_tkagg import FigureCanvasTkAgg, NavigationToolbar2Tk
                figure = _geometry_preview_figure(geometry)
            except Exception as exc:
                if messagebox:
                    messagebox.showerror("3D geometry preview", str(exc), parent=self)
                return
            window = tk.Toplevel(self)
            window.title("Quartz geometry preview: 3D + orthographic views")
            window.geometry("1180x820")
            window.minsize(900, 650)
            window.columnconfigure(0, weight=1)
            window.rowconfigure(0, weight=1)
            canvas = FigureCanvasTkAgg(figure, master=window)
            canvas.draw()
            canvas.get_tk_widget().grid(row=0, column=0, sticky="nsew")
            toolbar = NavigationToolbar2Tk(canvas, window, pack_toolbar=False)
            toolbar.update()
            toolbar.grid(row=1, column=0, sticky="ew")
            window._preview_figure = figure
            window._preview_canvas = canvas
            window._preview_toolbar = toolbar
            self._preview_3d_window = window

            def close_preview() -> None:
                plt.close(figure)
                self._preview_3d_window = None
                window.destroy()

            window.protocol("WM_DELETE_WINDOW", close_preview)

        # Compatibility aliases for callers of the earlier editor API.
        def preview_3d(self) -> None:
            return self.open_3d_viewer()

        def render_rows(self) -> None:
            for child in self.row_frame.winfo_children():
                child.destroy()
            self.vars = []
            headings = ("#", "Name", "Shape / Type", "Length mm", "Width start", "Width end", "Thick start", "Thick end", "Material", "Custom n", "Effective n")
            for col, heading in enumerate(headings):
                ttk.Label(self.row_frame, text=heading).grid(row=0, column=col, padx=2, sticky="w")
            for row_index, row in enumerate(self.rows, start=1):
                values: Dict[str, Any] = {}
                for key in ("name", "shape", "length_mm", "width_start_mm", "width_end_mm", "thickness_start_mm", "thickness_end_mm", "material", "custom_n"):
                    variable = tk.StringVar(value=str(row.get(key, "")))
                    values[key] = variable
                ttk.Label(self.row_frame, text=str(row_index)).grid(row=row_index, column=0, padx=2)
                ttk.Entry(self.row_frame, textvariable=values["name"], width=14).grid(row=row_index, column=1, padx=2, pady=1)
                shape_widget = ttk.Combobox(self.row_frame, textvariable=values["shape"], values=SECTION_SHAPES, state="readonly", width=12)
                shape_widget.grid(row=row_index, column=2, padx=2, pady=1)
                shape_widget.bind("<<ComboboxSelected>>", lambda _event, item=values: self._section_shape_changed(item))
                for col, key in enumerate(("length_mm", "width_start_mm", "width_end_mm", "thickness_start_mm", "thickness_end_mm"), start=3):
                    ttk.Entry(self.row_frame, textvariable=values[key], width=9).grid(row=row_index, column=col, padx=2, pady=1)
                ttk.Combobox(self.row_frame, textvariable=values["material"], values=MATERIAL_PRESETS, state="readonly", width=13).grid(row=row_index, column=8, padx=2)
                ttk.Entry(self.row_frame, textvariable=values["custom_n"], width=8).grid(row=row_index, column=9, padx=2)
                effective = tk.StringVar(value=f"{float(row.get('effective_n_450nm', 1.46)):.5f}")
                ttk.Label(self.row_frame, textvariable=effective, width=9).grid(row=row_index, column=10, padx=2)
                values["effective"] = effective
                for key in ("name", "shape", "length_mm", "width_start_mm", "width_end_mm", "thickness_start_mm", "thickness_end_mm", "material", "custom_n"):
                    values[key].trace_add("write", lambda *_args, item=values: self._refresh_effective(item))
                    values[key].trace_add("write", lambda *_args: self._refresh_preview())
                self.vars.append(values)

        def _section_shape_changed(self, variables: Mapping[str, Any]) -> None:
            """Keep dimensions coherent when a row changes type."""

            shape = variables["shape"].get().strip().lower()
            if shape == "rectangle":
                variables["width_end_mm"].set(variables["width_start_mm"].get())
                variables["thickness_end_mm"].set(variables["thickness_start_mm"].get())
            else:
                if not variables["width_end_mm"].get().strip():
                    variables["width_end_mm"].set(variables["width_start_mm"].get())
                if not variables["thickness_end_mm"].get().strip():
                    variables["thickness_end_mm"].set(variables["thickness_start_mm"].get())
            self._refresh_preview()

        def _refresh_effective(self, variables: Mapping[str, Any]) -> None:
            try:
                material = MaterialSpec(variables["material"].get(), None if not variables["custom_n"].get().strip() else float(variables["custom_n"].get()))
                variables["effective"].set(f"{material.effective_index(450e-9):.5f}")
            except (TypeError, ValueError):
                variables["effective"].set("invalid")

        def add_section(self) -> None:
            self.rows.append({"name": f"section_{len(self.rows) + 1}", "shape": "rectangle", "length_mm": 10.0, "width_start_mm": 20.0, "width_end_mm": 20.0, "thickness_start_mm": 10.0, "thickness_end_mm": 10.0, "material": "fused_silica", "custom_n": "", "effective_n_450nm": 1.46})
            self.render_rows()

        def remove_section(self) -> None:
            if len(self.rows) > 1:
                self.rows.pop()
                self.render_rows()

        def get_geometry(self) -> GeometryModel:
            rows = []
            for variables in self.vars:
                row = {key: variable.get() for key, variable in variables.items() if key != "effective"}
                rows.append(row)
            geometry = geometry_from_rows(rows, self.geometry_model)
            values: Dict[str, Any] = {key: variable.get() for key, variable in self.detector_vars.items()}
            plus = detector_area_from_values(values, "plus_x", *[value * 1e3 for value in geometry.face_dimensions("plus_x")])
            minus = detector_area_from_values(values, "minus_x", *[value * 1e3 for value in geometry.face_dimensions("minus_x")])
            return replace(geometry, detector_plus_x=plus, detector_minus_x=minus)

        def set_geometry(self, geometry: GeometryModel) -> None:
            self.geometry_model = geometry
            self.rows = geometry_to_rows(geometry)
            self.render_rows()
            for prefix, area in (("plus_x", geometry.detector_plus_x), ("minus_x", geometry.detector_minus_x)):
                self.detector_vars[f"{prefix}_enabled"].set(area.enabled)
                self.detector_vars[f"{prefix}_mode"].set(area.mode)
                self.detector_vars[f"{prefix}_shape"].set("full_exit_surface" if area.mode == "match_exit_surface" else area.shape)
                self.detector_vars[f"{prefix}_width_mm"].set("" if area.width_m is None else str(area.width_m * 1e3))
                self.detector_vars[f"{prefix}_height_mm"].set("" if area.height_m is None else str(area.height_m * 1e3))
                self.detector_vars[f"{prefix}_radius_mm"].set("" if area.radius_m is None else str(area.radius_m * 1e3))
                self.detector_vars[f"{prefix}_semi_axis_a_mm"].set("" if area.semi_axis_a_m is None else str(area.semi_axis_a_m * 1e3))
                self.detector_vars[f"{prefix}_semi_axis_b_mm"].set("" if area.semi_axis_b_m is None else str(area.semi_axis_b_m * 1e3))
                self.detector_vars[f"{prefix}_center_y_mm"].set(str(area.center_y_m * 1e3))
                self.detector_vars[f"{prefix}_center_z_mm"].set(str(area.center_z_m * 1e3))

        def load_geometry(self) -> None:
            if filedialog is None:
                return
            path = filedialog.askopenfilename(filetypes=(("Geometry JSON", "*.json"), ("All files", "*.*")))
            if path:
                try:
                    self.set_geometry(GeometryModel.from_json(path))
                    self._refresh_preview()
                except Exception as exc:
                    if messagebox:
                        messagebox.showerror("Load Geometry", str(exc))

        def save_geometry(self) -> None:
            if filedialog is None:
                return
            path = filedialog.asksaveasfilename(defaultextension=".json", filetypes=(("Geometry JSON", "*.json"),))
            if path:
                try:
                    self.get_geometry().to_json(path)
                except Exception as exc:
                    if messagebox:
                        messagebox.showerror("Save Geometry", str(exc))

        def load_json(self) -> None:
            return self.load_geometry()

        def save_json(self) -> None:
            return self.save_geometry()

        def preview_png(self) -> None:
            if filedialog is None:
                return
            path = filedialog.asksaveasfilename(defaultextension=".png", filetypes=(("PNG image", "*.png"),))
            if path:
                try:
                    plot_geometry_preview(self.get_geometry(), path)
                except Exception as exc:
                    if messagebox:
                        messagebox.showerror("Geometry preview", str(exc))
else:
    GeometryEditorPanel = None


__all__ = [
    "GeometryEditorPanel", "geometry_from_rows", "geometry_to_rows", "plot_geometry_preview",
    "MATERIAL_PRESETS", "SECTION_SHAPES", "canonical_assembly_model_from_sections", "canonical_complete_assembly_scene", "canonical_sensor_components", "canonical_mylar_overlay_components", "save_complete_detector_assembly", "load_complete_detector_assembly", "reorder_canonical_sections",
]
