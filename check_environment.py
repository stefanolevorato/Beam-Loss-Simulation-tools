#!/usr/bin/env python3
"""Headless v78 checks for the staged simulator source and Gamma transport."""

from __future__ import annotations

import ast
import copy
import json
import py_compile
import tempfile
import zipfile
import unittest
from pathlib import Path

ROOT = Path(__file__).resolve().parent
PYTHON_FILES = (
    "check_environment.py", "quartz_geometry.py", "quartz_geometry_editor.py",
    "quartz_cherenkov_sim.py", "quartz_cherenkov_gui.py", "gamma_transport.py",
    "test_gamma_transport.py",
)


def _assert(condition: bool, message: str) -> None:
    if not condition:
        raise AssertionError(message)


def _sensor_items(scene):
    return [item for item in scene.get("components", ()) if item.get("kind") == "sensor"]


def windows_launcher_structure_checks() -> None:
    """Check the Windows launch flow statically; this sandbox is not Windows."""

    start_bat = (ROOT / "Start_GUI.bat").read_text(encoding="utf-8")
    check_bat = (ROOT / "Check_Dependencies.bat").read_text(encoding="utf-8")
    launcher = (ROOT / "setup_and_start.ps1").read_text(encoding="utf-8")
    dependency_check = (ROOT / "check_dependencies.ps1").read_text(encoding="utf-8")
    _assert("setup_and_start.ps1" in start_bat and "%~dp0" in start_bat, "Start_GUI.bat does not invoke the colocated PowerShell launcher")
    _assert("check_dependencies.ps1" in check_bat and "%~dp0" in check_bat, "Check_Dependencies.bat does not invoke the colocated checker")
    _assert("Scripts\\python.exe" in launcher and "py -3 -m venv" in launcher, "Windows launcher does not create/use the project venv")
    _assert("Test-VenvPython" in launcher and "Move-Item" in launcher and ".venv.broken-" in launcher, "broken venv recovery does not preserve the previous folder")
    _assert("[switch]$InstallDependencies" in launcher and "$InstallDependencies" in launcher, "dependency installation is not opt-in")
    check_index = launcher.index("& $python $checker")
    gui_index = launcher.index("& $python $gui")
    _assert(check_index < gui_index and "Environment checks failed" in launcher[check_index:gui_index], "launcher does not require the environment suite before opening the GUI")
    _assert(".venv\\Scripts\\python.exe" in dependency_check and "check_environment.py" in dependency_check, "dependency checker uses the wrong Python or omits the suite")
    _assert("pip check" not in dependency_check, "dependency checker can fail on unrelated global package conflicts")
    print("Windows launcher structure: correct venv, preserved broken-venv recovery, opt-in installs, required checker-before-GUI PASS (Windows execution unavailable)")


def package_zip_integrity_checks() -> None:
    """Verify the complete available v78 project payload archive."""

    archive = ROOT / "quartz_cherenkov_simulator_v78.zip"
    _assert(archive.is_file(), "v78 package ZIP is missing")
    required = {
        "check_environment.py", "quartz_geometry.py", "quartz_geometry_editor.py",
        "quartz_cherenkov_sim.py", "quartz_cherenkov_gui.py", "gamma_transport.py",
        "test_gamma_transport.py", "gamma_range_model_diagnostic.md",
        "gamma_source_profile_collimated_v1.json", "gamma_reference_lines_v1.csv",
        "gamma_results.json", "gamma_results_example.json", "Assemblycompleto.json",
        "01Assemblycompleto.json", "01DefaultCherenkovAssembly.json", "DefaultScintillatorAssembly.json",
        "README_v69.md", "README_v70.md", "README_v71.md", "README_v72.md",
        "README_v73.md", "README_v74.md", "README_v76.md", "README_v77.md", "README_v78.md",
        "README_windows_setup.txt", "README_gamma_transport_stage.md",
        "README_two_assemblies_v62.md", "README_two_assemblies_v63.md",
        "implementation_summary_v63.txt", "implementation_summary_v74.txt",
        "implementation_summary_v76.txt", "implementation_summary_v77.txt", "implementation_summary_v78.txt",
        "v72_gamma_event_list_geometry_summary.txt", "v73_gamma_results_scroll_summary.txt",
        "requirements.txt", "requirementsroot.txt",
    }
    root = "quartz_cherenkov_simulator_v78/"
    with zipfile.ZipFile(archive, "r") as package:
        bad_member = package.testzip()
        _assert(bad_member is None, f"v78 package ZIP CRC failed for {bad_member}")
        members = set(package.namelist())
    missing = sorted(root + filename for filename in required if root + filename not in members)
    _assert(not missing, f"v78 package ZIP is missing required files: {', '.join(missing)}")
    print(f"v78 available-project ZIP: integrity PASS; required files present ({len(required)})")


def source_checks() -> None:
    for filename in PYTHON_FILES:
        py_compile.compile(str(ROOT / filename), doraise=True)
    release_files = (
        "gamma_source_profile_collimated_v1.json", "gamma_reference_lines_v1.csv", "gamma_results.json",
        "Assemblycompleto.json", "README_v76.md", "README_v77.md", "README_v78.md", "gamma_range_model_diagnostic.md",
        "implementation_summary_v76.txt", "implementation_summary_v77.txt", "implementation_summary_v78.txt", "test_gamma_transport.py", "requirements.txt",
        "requirementsroot.txt", "README_gamma_transport_stage.md", "README_two_assemblies_v62.md",
        "README_two_assemblies_v63.md", "implementation_summary_v63.txt",
    )
    for filename in release_files:
        _assert((ROOT / filename).is_file(), f"staged release source file is missing: {filename}")
    import quartz_geometry as geometry

    _assert(geometry.TERMINAL_SURFACE_CHOICES == ("Open", "Mirror", "Detector/Photosensor coupling"), "terminal choices are not exactly the v25 three choices")
    gui_source = (ROOT / "quartz_cherenkov_gui.py").read_text(encoding="utf-8")
    _assert("validate_canonical_run_assemblies" in gui_source, "canonical run validation entry point missing")
    _assert("build-20260929-v78-gamma-charged-boundaries" in gui_source, "v78 build identity is missing")
    _assert('APP_VERSION = "0.46.0"' in gui_source, "v78 application version is missing")
    gamma_source = (ROOT / "gamma_transport.py").read_text(encoding="utf-8")
    _assert("class ChargedTrackSegment" in gamma_source and "interaction_site" in gamma_source, "charged material-segment accounting is missing")
    _assert("remaining_areal_range" in gamma_source and "_next_volume(configured_volumes" in gamma_source, "charged transport does not query successive configured volumes")
    _assert('text="Open interactive global 3D Gamma view"' in gui_source, "Gamma Results is missing its separate-window viewer button")
    _assert("def _open_gamma_global_3d_viewer(self)" in gui_source and "tk.Toplevel(self)" in gui_source, "Gamma global viewer does not open an independent top-level window")
    _assert("build_gamma_global_view_scene(" in gui_source and "def _draw_gamma_global_3d_view(self)" in gui_source, "Gamma global viewer scene/render path is missing")
    _assert("_gamma_global_3d_event_selected" in gui_source and "_gamma_global_3d_drag_motion" in gui_source and "_gamma_global_3d_zoom_event" in gui_source, "Gamma global viewer event/camera controls are missing")
    _assert("terminal_surfaces_from_payload" in (ROOT / "quartz_geometry.py").read_text(encoding="utf-8"), "canonical and legacy terminal normalization is missing")
    _assert("def _has_instance_attr(instance, name: str) -> bool:" in gui_source and "object.__getattribute__(instance, \"__dict__\")" in gui_source, "safe Tk instance-state check is missing")
    _assert("hasattr(self, " not in gui_source, "unsafe Tk self hasattr lookup remains")
    _assert("getattr(self, " not in gui_source and "def _get_instance_attr(instance, name: str," in gui_source, "unsafe Tk self getattr lookup remains")
    _assert('tkfont.Font(root=canvas, family="TkDefaultFont", size=9, weight="bold")' in gui_source, "real global viewer lost its bold Tk label font")
    _assert("Readout radius (mm)" in gui_source and "Semi-axis A (mm)" in gui_source and "Semi-axis B (mm)" in gui_source, "shape-specific readout controls are missing")
    _assert("_sensor_shape_changed" in gui_source and "sensor_readout_dimensions_from_values" in gui_source, "shape callback/migration helper is missing")
    _assert("Scintillation decay time (ns)" in gui_source, "editable scintillation decay-time field is not visible")
    _assert("scintillator_decay_time_s_from_ns" in gui_source, "GUI decay-time conversion helper is missing")
    _assert('ROOT_NOTEBOOK_TAB_LABELS = (\n    "Run",\n    "Results",' in gui_source, "Results is not the peer tab immediately after Run")
    _assert('cget("background")' not in gui_source and "cget('background')" not in gui_source, "ttk background still queried with unsupported cget")
    _assert(gui_source.count("_ttk_background(info_frame)") == 2, "both 3D popup detail panes must use the ttk-safe background helper")
    import quartz_cherenkov_gui as gui

    class _NoBackgroundCget:
        def cget(self, option):
            raise AssertionError(f"unsafe cget called for {option!r}")

    class _FakeStyle:
        def __init__(self, widget):
            self.widget = widget

        def lookup(self, style_name, option):
            _assert(style_name == "TLabelframe" and option == "background", "unexpected ttk style lookup")
            return ""

    class _FakeTtk:
        Style = _FakeStyle

    original_ttk = gui.ttk
    try:
        gui.ttk = _FakeTtk
        _assert(gui._ttk_background(_NoBackgroundCget()) == "#f0f0f0", "ttk background helper fallback is empty or invalid")
    finally:
        gui.ttk = original_ttk
    _assert("beam-monitor-global-assembly-profile-v1" in gui_source, "complete Global Assembly profile schema is missing")
    _assert("Save complete assembly profile" in gui_source and "Load complete assembly profile" in gui_source, "complete profile controls are missing")
    _assert("Maximum photon bounces (selected detector(s))" in gui_source, "Run bounce control is missing")
    _assert("Inclination from +Z / track tilt (deg)" in gui_source, "Run incidence-angle control is missing")
    _assert("positive tilts toward +X" in gui_source, "Run incidence-angle convention is not visible")
    _assert("build-20260924-v50-global-assembly-geometry-clarity" not in gui_source and "build-20260924-v48-event-global-frame-alignment" not in gui_source and "build-20260924-v47" not in gui_source and "v43" not in gui_source, "stale GUI build label remains")
    _assert("Building assemblies" in gui_source and 'kind == "status"' in gui_source, "run phase/status queue path missing")
    for phase in ("Closing ROOT output", "Aggregating results", "Writing plots", "Saving summary", "Updating GUI"):
        _assert(phase in gui_source, f"post-ROOT phase missing: {phase}")
    _assert("_append_persistent_run_log" in gui_source and "run.log" in gui_source, "persistent run log path missing")
    _assert("compact_event_view_payload" in gui_source and "_event_view_scene" in gui_source, "compact/in-memory event scene split missing")
    _assert('run_directory / "gamma_transport.png"' in gui_source, "Gamma event visualization output is missing")
    _assert("gamma_transport._write_plot" in gui_source, "Gamma event visualization renderer is missing")
    _assert("terminal_sent" in gui_source and "Worker terminated before reporting a terminal state" in gui_source, "finally-safe worker terminal route missing")
    _assert("_refresh_loaded_assembly_views" in gui_source and "_suppress_assembly_refresh" in gui_source, "atomic post-load preview refresh contract missing")
    run_node = next(node for node in ast.walk(ast.parse(gui_source)) if isinstance(node, ast.FunctionDef) and node.name == "_run")
    run_text = ast.get_source_segment(gui_source, run_node) or ""
    _assert("validate_detector_active_areas" not in run_text, "Run still calls the obsolete generic portal validator")
    sim_source = (ROOT / "quartz_cherenkov_sim.py").read_text(encoding="utf-8")
    _assert("status_callback" in sim_source, "simulator status callback path missing")
    _assert("if not getattr(geometry, \"sensor_components\", ())" in sim_source, "canonical transport did not bypass legacy portal validation")
    _assert("TERMINAL_SURFACE_CHOICES" in gui_source and "Optical/light-guide coupling" not in gui_source, "legacy fourth terminal choice remains user-facing")
    _assert("Direct vs guide coupling" not in gui_source and "Through complete light guide" not in gui_source, "editable Sensor-tab route remains")
    _assert("Coupled component (derived status)" in gui_source and "Optical path (derived, read-only)" in gui_source, "derived Sensor-tab fields missing")
    _assert("Apply ideal sensor settings" in gui_source, "ideal sensor button missing")
    _assert('"detector_valid_sensor_hit"' in sim_source and "outcome_reason" in sim_source, "sensor-window outcome accounting fields missing")
    _assert("dimensions_mm" in (ROOT / "quartz_geometry.py").read_text(encoding="utf-8"), "canonical sensor dimensions are not handed to transport")
    _assert("responsive_three_view_panel_bounds" in gui_source and 'canvas.bind("<Configure>"' in gui_source, "responsive three-view contract missing")
    _assert("collision_aware_callout_layout" in gui_source and "detector_callout_label" in gui_source, "shared collision-aware 3D label layout missing")
    _assert("interactive_event_view_layout" in gui_source and "status_max_lines" in gui_source, "responsive event header layout contract is missing")
    tree = ast.parse(gui_source)
    functions = {node.name: node for node in ast.walk(tree) if isinstance(node, ast.FunctionDef)}
    for function_name in ("_draw_global_assembly_3d_view", "_draw_global_assembly_3d_view_impl", "_draw_interactive_event_view"):
        text = ast.get_source_segment(gui_source, functions[function_name]) or ""
        if function_name != "_draw_global_assembly_3d_view":
            _assert("collision_aware_callout_layout" in text, f"{function_name} does not use shared label layout")
    _assert("_draw_global_assembly_3d_view_impl" in (ast.get_source_segment(gui_source, functions["_draw_global_assembly_3d_view"]) or ""), "safe global draw wrapper is not connected")
    _assert("global_assembly_camera_basis" in gui_source and "global_assembly_orientation_triad" in gui_source, "shared camera/triad projection contract is missing")
    _assert("global XYZ • camera-coupled" in gui_source, "popup orientation gizmo is not identified as camera-coupled")
    _assert('global_preview_canvas.bind("<Configure>"' in gui_source, "embedded preview resize binding is missing")
    _assert("_schedule_global_assembly_3d_draw" in gui_source and "after_idle" in gui_source, "mapped global viewer redraw scheduler is missing")
    _assert("global_viewer_render_error" in gui_source and "Global assembly render error" in gui_source, "global viewer errors are not visible")
    _assert("global_assembly_orthographic_layout" in gui_source and "XZ • longitudinal / beam" in gui_source and "YZ • terminal / transverse" in gui_source, "global orthographic pane renderer is missing")
    _assert("Show compact callouts" in gui_source, "3D compact-label toggle is missing")
    for label in ("Save Cherenkov 3D Geometry JSON", "Load Cherenkov 3D Geometry JSON", "Save Scintillator 3D Geometry JSON", "Load Scintillator 3D Geometry JSON"):
        _assert(label in gui_source, f"missing complete 3D command: {label}")
    _assert("complete_detector_assembly_geometry_to_dict" in gui_source and "complete_detector_assembly_geometry_from_dict" in gui_source, "complete assembly serializer is not wired")
    _assert("_save_assembly_model(\"scintillator\")" in gui_source and "_load_assembly_model(\"scintillator\")" in gui_source, "scintillator geometry JSON callbacks missing")
    for field in ("detector_transport_diagnostic", "Cherenkov Transport Diagnostic", "Scintillator Transport Diagnostic", "terminal_surfaces", "light_guide", "mylar", "sensor_readout_state", "optical_transport_path", "warnings_errors"):
        _assert(field in gui_source, f"detector-specific diagnostic field missing: {field}")
    tree = ast.parse(gui_source)
    functions = {node.name: node for node in ast.walk(tree) if isinstance(node, ast.FunctionDef)}
    for function_name in ("_build_cherenkov_sensor_tab", "_build_scintillator_sensor_tab_v13"):
        text = ast.get_source_segment(gui_source, functions[function_name]) or ""
        if function_name == "_build_cherenkov_sensor_tab":
            _assert("ttk.Combobox" not in text, f"editable route combobox remains in {function_name}")
        _assert("Direct vs guide coupling" not in text and "Through complete light guide" not in text, f"editable route remains in {function_name}")
        _assert("sensor_1_component" not in text and "sensor_2_component" not in text, f"editable component remains in {function_name}")


def detector_selection_checks() -> None:
    """Verify the central detector selection, isolation, and mode-labelled viewers."""

    import quartz_cherenkov_gui as gui
    import quartz_cherenkov_sim as simulator
    from quartz_geometry import AssemblyGeometryModel, AssemblySection

    expected_labels = ("Cherenkov only", "Scintillator only", "All")
    _assert(tuple(gui.DETECTOR_RUN_MODE_LABELS.values()) == expected_labels, "Run detector choices changed")
    source = (ROOT / "quartz_cherenkov_gui.py").read_text(encoding="utf-8")
    _assert("Detectors to simulate" in source and "key == \"detector_run_mode\"" in source, "central Run detector selector missing")
    _assert(
        source.count("values=tuple(DETECTOR_RUN_MODE_LABELS.values())") == 1
        or source.count("values=RUN_MODE_COMBO_VALUES") == 1,
        "detector selector is not a single Run combobox",
    )
    _assert("Enable scintillator" not in source, "editable scintillator enable control remains in source")
    _assert(gui.default_values()["detector_run_mode"] == "Cherenkov only", "default is not Cherenkov-only")

    values = gui.default_values()
    values.update({
        "events": "1", "photons_per_event": "4", "scintillator_photon_count": "5",
        "scintillator_deposited_energy_MeV": "1", "display_photons": "4",
    })
    # An explicit Run choice wins over the old compatibility flag.
    values.update({"detector_run_mode": "Cherenkov only", "scintillator_enabled": "1"})
    _assert(gui.validate_settings(values).detector_run_mode == "cherenkov", "legacy checkbox overrode explicit Run mode")
    values["detector_run_mode"] = "Scintillator only"
    settings = gui.validate_settings(values)
    _assert(settings.detector_run_mode == "scintillator", "visible Scintillator-only choice did not normalize")

    particle = gui._particle_for_settings(settings)
    optical = gui._optical_for_settings(settings)
    geometry = simulator.BarGeometry()
    local = gui.scintillator_local_geometry_from_values(values)
    global_geometry = gui.scintillator_global_geometry_from_values(values)
    assembly = simulator.compose_scintillator_assembly(local, global_geometry, geometry)
    config = gui.scintillator_config_from_values(values)
    canonical = gui.global_assembly_scene(values, geometry)
    for mode in ("cherenkov", "scintillator", "combined"):
        statuses, progress = [], []
        result = simulator.simulate_events(
            particle, 1, geometry=geometry, optical=optical, seed=2718,
            progress_callback=lambda *item: progress.append(item),
            status_callback=statuses.append,
            scintillator_assembly=assembly if mode != "cherenkov" else None,
            scintillator_config=config if mode != "cherenkov" else None,
            detector_mode=mode,
        )[0]
        _assert(progress and progress[-1] == (0, 1, 1.0), f"{mode}: progress did not terminate at 100%")
        _assert(any(mode == "scintillator" and "scintillator" in item.lower() or mode != "scintillator" and "Cherenkov" in item for item in statuses), f"{mode}: detector status missing")
        if mode == "cherenkov":
            _assert(not result.scintillator and len(result.photon_records) > 0, "Cherenkov-only produced scintillator records")
        elif mode == "scintillator":
            _assert(not result.photon_records and result.scintillator is not None and result.scintillator.photon_records, "Scintillator-only mixed in Cherenkov records")
        else:
            _assert(result.photon_records and result.scintillator is not None and result.scintillator.photon_records, "All did not produce both detector families")
        payload = gui.canonical_event_view_payload(canonical, mode)
        visible = [mode] if mode != "combined" else ["cherenkov", "scintillator"]
        _assert(payload["visible_assembly_types"] == visible, f"{mode}: payload visibility mismatch")
        _assert([item["assembly_type"] for item in payload["assemblies"]] == visible, f"{mode}: non-participating assembly persisted in event payload")
        event_scene = gui.build_interactive_event_scene(result, payload, detector_mode=mode, display_photons=4)
        _assert(event_scene["visible_assembly_types"] == visible and set(event_scene["per_detector"]) == set(visible), f"{mode}: viewer ownership mismatch")
        _assert(gui.event_table_rows([result])[0]["detector_mode"] == mode, f"{mode}: event table mode missing")

    scint_model = AssemblyGeometryModel(
        "Scintillator", [AssemblySection(name="Active", material="scintillator", axis="X")],
        terminal_plus_surface="Detector/Photosensor coupling",
    )
    single = gui.validate_canonical_run_assemblies(values, {"scintillator": scint_model}, "scintillator")
    _assert(single["scintillator"]["participates"] and not single["cherenkov"]["participates"], "single-detector validation still requires the unselected assembly")


def gamma_viewer_persistence_checks() -> None:
    """Exercise Gamma viewer wiring and independent run/startup persistence."""

    import quartz_cherenkov_gui as gui
    import quartz_cherenkov_sim as simulator
    import quartz_geometry as geometry
    import numpy as np
    from matplotlib.figure import Figure

    source = (ROOT / "quartz_cherenkov_gui.py").read_text(encoding="utf-8")
    _assert("Gamma Results" in gui.ROOT_NOTEBOOK_TAB_LABELS, "dedicated Gamma Results tab is missing")
    _assert("self.gamma_event_figure, self.gamma_event_canvas" in source and "draw_gamma_event_view" in source, "Gamma event viewer is not wired to the tab")
    _assert("gamma_tracks" in source and "gamma_selected_event" in source, "Gamma event selection or track data flow is missing")
    _assert("Per-event results" in source and "gamma_results_event_table" in source, "selectable Gamma per-event table is missing")
    _assert("gamma-global-assembly-geometry-v1" in source, "persisted Gamma geometry snapshot schema is missing")
    _assert("Per-event results" in source and "gamma_results_event_table" in source, "selectable Gamma per-event table is missing")
    _assert("geometry_snapshot" in source and "gamma-global-assembly-geometry-v1" in source, "persisted Gamma Global Assembly geometry snapshot is missing")
    _assert(gui.detector_run_mode("All") == "combined", "All no longer maps to combined optical transport")
    _assert(gui.detector_run_mode("Gamma transport") == gui.GAMMA_RUN_MODE, "Gamma remains unavailable as an explicit run mode")

    event = {
        "event_id": 7,
        "source_energy_keV": 511.0,
        "source_component": "line",
        "interaction_counts": {"compton": 1},
        "gamma_tracks": [{
            "start_position_m": (0, 0, -0.001), "end_position_m": (0, 0, 0.002),
            "volume": "quartz", "interaction": "compton", "energy_keV": 511,
        }],
        "charged_tracks": [{
            "position_m": (0, 0, 0.002), "direction": (1, 0, 0),
            "path_mm": 4, "kind": "electron", "volume": "quartz",
        }],
        "deposits": [{
            "position_m": (0, 0, 0.002), "interaction": "compton",
            "energy_keV": 250, "volume": "quartz",
        }],
        "light": [{"mechanism": "cherenkov", "expected_photons": 12.5}],
    }
    scene = gui.gamma_event_view_data(event)
    _assert(scene["gamma_segments"][0]["end_mm"] == (0.0, 0.0, 2.0), "Gamma track did not map into the global-mm event view")
    _assert(scene["charged_segments"][0]["end_mm"] == (4.0, 0.0, 2.0), "charged-secondary track endpoint was not reconstructed")
    figure = Figure(figsize=(4, 3))
    axis = gui.draw_gamma_event_view(figure, event)
    _assert(len(axis.lines) == 2 and len(axis.collections) == 1, "Gamma event viewer omitted photon/secondary paths or interaction deposits")

    empty_event = {
        "event_id": 0, "source_energy_keV": 835.0, "source_component": "continuum",
        "interaction_counts": {}, "deposits": [], "charged_tracks": [], "gamma_tracks": [], "light": [],
    }
    mixed_events = [empty_event, event]
    rows = gui.gamma_event_table_rows(mixed_events)
    _assert([row["status"] for row in rows] == ["EMPTY", "INTERACTION"], "mixed empty/populated Gamma events were not classified correctly")
    _assert(rows[0]["deposit_record_count"] == 0 and rows[0]["charged_track_count"] == 0, "empty Gamma event counts are wrong")
    _assert(rows[1]["interaction_process_summary"] == "compton ×1" and rows[1]["local_light_expected_photons"] == 12.5, "populated Gamma event summary omitted processes or local light")

    class Var:
        def __init__(self, value=""):
            self.value = value
        def get(self):
            return self.value
        def set(self, value):
            self.value = value

    with tempfile.TemporaryDirectory(prefix="v71_gamma_restore_") as directory:
        root = Path(directory)
        gamma_values = gui.default_values()
        gamma_values.update({
            "output_directory": str(root), "detector_run_mode": "Gamma transport",
            "events": "6", "seed": "1741", "gamma_max_interactions": "30",
            "gamma_max_light_samples": "100",
            "global_offset_x_mm": "5", "global_offset_y_mm": "-3", "global_offset_z_mm": "11",
            "detector2_offset_x_mm": "-4", "detector2_offset_y_mm": "6", "detector2_offset_z_mm": "48",
            "assembly_1_rotation_axis": "Z", "assembly_1_rotation_angle_deg": "15",
            "assembly_2_rotation_axis": "Y", "assembly_2_rotation_angle_deg": "-10",
        })
        gamma_settings = gui.validate_settings(gamma_values)
        gamma_models = {
            "cherenkov": geometry.AssemblyGeometryModel(
                "Cherenkov", [geometry.AssemblySection(
                    name="Quartz active", material="fused silica", axis="X",
                    length_x_mm=32.0, width_y_mm=14.0, thickness_z_mm=8.0,
                )], [geometry.AssemblySection(
                    name="Quartz guide", material="light guide", axis="X",
                    length_x_mm=9.0, width_y_mm=18.0, thickness_z_mm=12.0, center_x_mm=20.5,
                )], True, "Detector/Photosensor coupling", "Mirror",
            ),
            "scintillator": geometry.AssemblyGeometryModel(
                "Scintillator", [geometry.AssemblySection(
                    name="Scintillator active", material="scintillator", axis="X",
                    length_x_mm=24.0, width_y_mm=18.0, thickness_z_mm=10.0,
                )], [geometry.AssemblySection(
                    name="Scintillator guide", material="light guide", axis="X",
                    length_x_mm=8.0, width_y_mm=20.0, thickness_z_mm=12.0, center_x_mm=16.0,
                )], True, "Detector/Photosensor coupling", "Open",
            ),
        }
        gamma_run = root / "quartz_run_gamma"
        gamma_result = gui.run_simulation_to_directory(
            gamma_settings, gamma_run, geometry=simulator.BarGeometry(),
            assembly_models=gamma_models,
        )
        _assert(gamma_result["summary"]["detector_mode"] == "gamma", "successful Gamma run was not integrated into run results")
        artifact = Path(gamma_result["gamma_results_path"])
        saved_gamma = json.loads(artifact.read_text(encoding="utf-8"))
        _assert(len(saved_gamma["events"]) == 6 and all("gamma_tracks" in item for item in saved_gamma["events"]), "completed Gamma event data was not persisted for the viewer")
        geometry_snapshot = saved_gamma.get("geometry_snapshot")
        geometry_data = gui.gamma_geometry_view_data(geometry_snapshot)
        _assert(geometry_data["assembly_types"] == ["cherenkov", "scintillator"], "Gamma result did not preserve both run assemblies")
        _assert(geometry_snapshot["assembly_count"] == 2 and set(geometry_snapshot["profiles"]) == {"cherenkov", "scintillator"}, "Gamma geometry snapshot did not include both complete profiles")
        scene_values = gui._canonical_scene_values_from_settings(gamma_settings)
        global_scene = gui.global_assembly_scene(scene_values, simulator.BarGeometry(), gamma_models)
        saved_transforms = {item["assembly_type"]: item["transform"] for item in geometry_snapshot["scene"]["assemblies"]}
        _assert(saved_transforms == {key: global_scene["assembly_by_type"][key]["transform"] for key in ("cherenkov", "scintillator")}, "saved Gamma transforms differ from Global Assembly")
        gamma_volumes = gui._gamma_snapshot_volumes(gamma_models, scene_values)
        _assert({item.assembly for item in gamma_volumes} == {"cherenkov", "scintillator"} and len(gamma_volumes) == 4, "Gamma transport did not use the two active/guide assembly profiles")
        for assembly_type in ("cherenkov", "scintillator"):
            expected_points = [
                point
                for component in global_scene["assembly_by_type"][assembly_type]["components"]
                if component.get("kind") in {"geometry", "guide"}
                for face in component.get("faces_mm", ())
                for point in face
            ]
            actual_volumes = [item for item in gamma_volumes if item.assembly == assembly_type]
            actual_low = np.min(np.asarray([item.solid.bounds_min_m for item in actual_volumes], dtype=float), axis=0) * 1000.0
            actual_high = np.max(np.asarray([item.solid.bounds_max_m for item in actual_volumes], dtype=float), axis=0) * 1000.0
            expected = np.asarray(expected_points, dtype=float)
            _assert(np.allclose(actual_low, expected.min(axis=0), atol=1e-8) and np.allclose(actual_high, expected.max(axis=0), atol=1e-8), f"{assembly_type} Gamma transport bounds differ from saved Global Assembly faces")

        class Text:
            def __init__(self):
                self.value = ""
            def configure(self, **_kwargs):
                pass
            def delete(self, *_args):
                self.value = ""
            def insert(self, _index, value):
                self.value = str(value)

        class Tree:
            def __init__(self):
                self.items = {}
                self.selected = []
            def get_children(self):
                return tuple(self.items)
            def delete(self, *items):
                for item in items or tuple(self.items):
                    self.items.pop(item, None)
                self.selected = [item for item in self.selected if item in self.items]
            def insert(self, _parent, _index, iid, values):
                self.items[str(iid)] = tuple(values)
            def selection(self):
                return tuple(self.selected)
            def selection_set(self, item):
                self.selected = [str(item)]
            def see(self, _item):
                pass

        class Selector:
            def configure(self, **kwargs):
                self.values = kwargs.get("values", ())

        class Canvas:
            def draw_idle(self):
                pass

        probe = type("GammaTableProbe", (), {})()
        probe.gamma_selected_event = Var("0")
        probe.latest_gamma_payload = {
            "model_status": saved_gamma["model_status"], "summary": saved_gamma["summary"],
            "events": mixed_events, "geometry_snapshot": geometry_snapshot,
        }
        probe.gamma_results_status = Var()
        probe.gamma_results_summary_text = Var()
        probe.gamma_results_aggregate_text = Text()
        probe.gamma_results_event_text = Text()
        probe.gamma_results_event_selector = Selector()
        probe.gamma_results_event_table = Tree()
        probe.gamma_event_figure = Figure(figsize=(4, 3))
        probe.gamma_event_canvas = Canvas()
        probe._set_readonly_text = gui.QuartzCherenkovGUI._set_readonly_text
        probe._refresh_gamma_results_view = lambda: gui.QuartzCherenkovGUI._refresh_gamma_results_view(probe)
        gui.QuartzCherenkovGUI._refresh_gamma_results_view(probe)
        _assert(probe.gamma_results_event_table.get_children() == ("0", "7"), "Gamma event table did not list every event")
        _assert(probe.gamma_results_event_table.items["0"][1] == "EMPTY" and probe.gamma_results_event_table.items["7"][1] == "INTERACTION", "Gamma event table statuses are not visible per row")
        zero_axis = probe.gamma_event_figure.axes[0]
        _assert(len(zero_axis.collections) == len(geometry_data["components"]), "zero-interaction event viewer omitted saved detector geometry")
        _assert("EMPTY" in zero_axis.get_title(), "zero-interaction event is not labeled in the Gamma viewer")
        probe.gamma_results_event_table.selection_set("7")
        gui.QuartzCherenkovGUI._gamma_results_table_selected(probe)
        _assert(probe.gamma_selected_event.get() == "7" and "Event 7 — INTERACTION" in probe.gamma_results_event_text.value, "Gamma row selection did not synchronize the detail panel")
        _assert("compton ×1" in probe.gamma_results_event_text.value and "12.5 expected photons" in probe.gamma_results_event_text.value, "selected Gamma event details omitted process/light values")
        _assert(probe.gamma_event_figure.axes[0].get_title().startswith("Gamma event 7:"), "Gamma row selection did not synchronize the event viewer")

        legacy_probe = type("LegacyGammaProbe", (), {})()
        legacy_probe.gamma_selected_event = Var("0")
        legacy_probe.latest_gamma_payload = {"model_status": "legacy", "summary": {}, "events": [empty_event]}
        legacy_probe.gamma_results_status = Var()
        legacy_probe.gamma_results_summary_text = Var()
        legacy_probe.gamma_results_aggregate_text = Text()
        legacy_probe.gamma_results_event_text = Text()
        legacy_probe.gamma_results_event_selector = Selector()
        legacy_probe.gamma_results_event_table = Tree()
        legacy_probe.gamma_event_figure = Figure(figsize=(4, 3))
        legacy_probe.gamma_event_canvas = Canvas()
        legacy_probe._set_readonly_text = gui.QuartzCherenkovGUI._set_readonly_text
        gui.QuartzCherenkovGUI._refresh_gamma_results_view(legacy_probe)
        _assert("snapshot unavailable" in legacy_probe.gamma_results_status.get(), "legacy Gamma result did not explain missing geometry snapshot")
        _assert("no v72 run geometry snapshot" in legacy_probe.gamma_results_event_text.value, "legacy Gamma event detail did not explain missing geometry")
        _assert(not legacy_probe.gamma_event_figure.axes[0].collections, "legacy Gamma result guessed or substituted detector geometry")
        gamma_loaded, gamma_loaded_run, gamma_loaded_artifact = gui.load_latest_gamma_result(root)
        _assert(gamma_loaded is not None and gamma_loaded_run == gamma_run and gamma_loaded_artifact == artifact, "successful Gamma result was not indexed independently")
        _assert(gui.retain_latest_completed_gamma_payload(gamma_loaded, {"summary": {"detector_mode": "combined"}}) is gamma_loaded, "All optical run erased the in-memory completed Gamma result")
        _assert(gui.retain_latest_completed_gamma_payload(gamma_loaded, {"stopped": True, "_gamma_result": {"events": []}}) is gamma_loaded, "stopped Gamma run replaced the last completed result")

        optical_values = gui.default_values()
        optical_values.update({
            "output_directory": str(root), "detector_run_mode": "All",
            "events": "1", "photons_per_event": "2", "scintillator_photon_count": "2",
            "display_photons": "2", "seed": "1742",
        })
        optical_settings = gui.validate_settings(optical_values)
        optical_run = root / "quartz_run_all"
        optical_result = gui.run_simulation_to_directory(
            optical_settings, optical_run, geometry=simulator.BarGeometry(),
        )
        _assert(optical_result["summary"]["detector_mode"] == "combined", "All no longer executes both optical detector families")
        _assert(optical_result.get("_gamma_result") is None, "All optical run unexpectedly reports an integrated Gamma simulation")
        gamma_after_optical, restored_gamma_run, restored_artifact = gui.load_latest_gamma_result(root)
        _assert(gamma_after_optical is not None and restored_gamma_run == gamma_run and restored_artifact == artifact, "All optical run overwrote the saved Gamma result")

        # Exercise the same startup restore method used by a real Tk session,
        # while replacing only widget operations with light headless doubles.
        probe = type("StartupProbe", (), {})()
        probe.values = {"output_directory": Var(str(root))}
        probe.gamma_selected_event = Var("0")
        probe.last_run_directory = None
        probe.last_run_text = Var()
        probe.latest_run_text = Var()
        probe.latest_summary = None
        probe.latest_gamma_payload = None
        probe.latest_readout_payload = None
        probe.latest_event_global_scene = None
        probe.total_summary_text = Var()
        probe._append_log = lambda _message: None
        probe._populate_event_table = lambda **_kwargs: None
        probe._refresh_results_tab = lambda: None
        probe._refresh_gamma_results_view = lambda: None
        gui.QuartzCherenkovGUI._load_latest_summary(probe)
        _assert(probe.last_run_directory == optical_run, "startup did not keep the newest All optical run as the latest general run")
        _assert(probe.latest_gamma_payload is not None and probe.latest_gamma_payload["events"] == saved_gamma["events"], "startup failed to restore Gamma results independently of the latest All run")
        print("v73 Gamma event table/selection, mixed event statuses, zero-event complete geometry, legacy snapshot messaging, and independent result restore PASS")


def gamma_results_scroll_checks() -> None:
    """Check the Gamma tab scroll container and wheel routing without a display."""

    from types import SimpleNamespace
    import quartz_cherenkov_gui as gui

    source = (ROOT / "quartz_cherenkov_gui.py").read_text(encoding="utf-8")
    _assert(
        "self._build_gamma_results_tab(self._make_scrollable_tab(" in source,
        "Gamma Results is not hosted in the vertical scroll container",
    )
    _assert(
        'ttk.Scrollbar(viewport, orient="vertical", command=canvas.yview)' in source,
        "scrollable tab is missing its vertical scrollbar",
    )
    _assert(
        'ttk.Scrollbar(event_table_frame, orient="horizontal"' in source,
        "Gamma event list lost its horizontal scrollbar",
    )
    _assert(
        "self._scrollable_mousewheel_passthrough.extend((" in source,
        "Gamma event table and text panes do not retain their own wheel scrolling",
    )

    class FakeWidget:
        def __init__(self, path):
            self.path = path

        def __str__(self):
            return self.path

    class FakeCanvas:
        def __init__(self):
            self.scrolls = []

        def yview_scroll(self, units, kind):
            self.scrolls.append((units, kind))

    inner = FakeWidget(".app.gamma_viewport.content")
    event_table = FakeWidget(".app.gamma_viewport.content.event_table")
    canvas = FakeCanvas()
    pointer_widget = {"value": FakeWidget(".app.gamma_viewport.content.viewer")}
    probe = SimpleNamespace(
        winfo_containing=lambda _x, _y: pointer_widget["value"],
        _scrollable_mousewheel_passthrough=[event_table],
        _scrollable_regions=[(canvas, inner)],
    )

    result = gui.QuartzCherenkovGUI._scrollable_mousewheel(
        probe, SimpleNamespace(x_root=10, y_root=20, delta=120, num=None),
    )
    _assert(result == "break" and canvas.scrolls == [(-1, "units")], "Gamma mouse wheel did not scroll the containing results tab")

    pointer_widget["value"] = FakeWidget(".app.gamma_viewport.content.event_table")
    result = gui.QuartzCherenkovGUI._scrollable_mousewheel(
        probe, SimpleNamespace(x_root=10, y_root=20, delta=-120, num=None),
    )
    _assert(result is None and len(canvas.scrolls) == 1, "Gamma tab wheel handler intercepted the event table's own scrolling")

    pointer_widget["value"] = FakeWidget(".app.other_tab.control")
    result = gui.QuartzCherenkovGUI._scrollable_mousewheel(
        probe, SimpleNamespace(x_root=10, y_root=20, delta=0, num=5),
    )
    _assert(result is None and len(canvas.scrolls) == 1, "Gamma wheel handler affected a different tab")
    print("v73 Gamma Results vertical scrollbar and scoped mouse-wheel routing PASS")


def scintillator_transport_debug_checks() -> None:
    """Reproduce the v44 Scintillator-only path and ROOT boundary headlessly."""

    import numpy as np
    import quartz_cherenkov_gui as gui
    import quartz_cherenkov_sim as simulator

    values = gui.default_values()
    values.update({
        "events": "1", "photons_per_event": "8", "display_photons": "8",
        "detector_run_mode": "Scintillator only",
        "scintillator_guide_enabled": "1",
        "scintillator_photon_count": "160",
        "scintillator_deposited_energy_MeV": "1.0",
        "scintillator_sampling_mode": "Fixed Monte Carlo test sample",
        "scintillator_max_bounces": "80",
        "scintillator_readout_index": "1.49",
        "scintillator_readout_acceptance_deg": "90",
        "scintillator_readout_qe": "1.0",
        "scintillator_readout_collection_efficiency": "1.0",
    })
    geometry = simulator.BarGeometry()
    local = gui.scintillator_local_geometry_from_values(values)
    assembly = local.to_assembly(mylar_enabled=True, mylar_reflectivity=0.95)
    config = gui.scintillator_config_from_values(values)
    particle = simulator.Particle(
        mass_kg=0.51099895 * simulator.MEV / simulator.C_LIGHT ** 2,
        kinetic_energy_J=2.0 * simulator.MEV,
    )

    results = {}
    for mode in ("cherenkov", "scintillator", "combined"):
        result = simulator.simulate(
            particle, geometry=geometry, optical=simulator.OpticalConfig(photon_count=8),
            seed=2024, detector_mode=mode,
            scintillator_assembly=assembly if mode != "cherenkov" else None,
            scintillator_config=config if mode != "cherenkov" else None,
        )
        results[mode] = result
        if mode == "cherenkov":
            _assert(result.scintillator is None and result.photon_records, "Cherenkov-only ownership changed")
            continue
        _assert(result.scintillator is not None, f"{mode}: scintillator result was not produced")
        if mode == "scintillator":
            _assert(not result.photon_records, "Scintillator-only mixed in Cherenkov records")
        else:
            _assert(result.photon_records, "combined: Cherenkov records were lost")
        records = result.scintillator.photon_records
        _assert(records and all(record.photon_kind == "scintillator" for record in records), f"{mode}: photon kind is not detector-owned")
        _assert(all(record.path_coordinate_frame == "scintillator_local_m" for record in records), f"{mode}: scintillator path frame is not explicit")
        active = assembly.scintillator
        guide = assembly.light_guide
        _assert(guide is not None, "focused reproduction did not configure a light guide")
        for record in records:
            source = np.asarray(record.emission_position_m, dtype=float)
            _assert(active.contains(source, tolerance=1e-10), "scintillation source is outside the active solid")
            points = [np.asarray(point, dtype=float) for point in record.path_points_m]
            _assert(any(active.contains(point, tolerance=1e-9) for point in points), "path lacks an active-volume point")
            if record.reached_guide:
                _assert(any(guide.contains(point, tolerance=1e-9) for point in points), "guide entry was recorded without a guide-frame path point")
                _assert(any(abs(point[0] - active.bounds_max_m[0]) <= 1e-9 for point in points), "active/guide boundary endpoint is missing")
                _assert(any(abs(point[0] - guide.bounds_min_m[0]) <= 1e-9 for point in points), "guide entry endpoint is missing")
            if record.reached_readout:
                _assert(abs(points[-1][0] - guide.bounds_max_m[0]) <= 1e-9, "readout path does not end at the guide output boundary")
        counts = simulator.detector_diagnostic_counts(records)
        _assert(sum(record.reached_guide for record in records) > 0, f"{mode}: no deterministic guide entries")
        _assert(sum(record.reached_readout for record in records) > 0, f"{mode}: no deterministic guide readout crossings")
        _assert(counts["detector_attempts"] > 0, f"{mode}: no detector/readout attempts were recorded")

        canonical = gui.global_assembly_scene(values, geometry)
        payload = gui.canonical_event_view_payload(canonical, mode)
        event_scene = gui.build_interactive_event_scene(result, payload, detector_mode=mode, display_photons=160)
        expected_visible = ["scintillator"] if mode == "scintillator" else ["cherenkov", "scintillator"]
        _assert(event_scene["visible_assembly_types"] == expected_visible, f"{mode}: event scene visibility changed")
        if mode == "scintillator":
            legend = "\n".join(event_scene["legend"]).lower()
            _assert("cherenkov" not in legend, "Scintillator-only legend names Cherenkov")
            _assert(all(line["detector_type"] == "scintillator" for line in event_scene["photon_lines"]), "Scintillator-only viewer mixed photon ownership")
        if mode == "combined":
            _assert(any(line["detector_type"] == "cherenkov" for line in event_scene["photon_lines"]), "All viewer lost Cherenkov paths")
            _assert(any(line["detector_type"] == "scintillator" for line in event_scene["photon_lines"]), "All viewer lost scintillator paths")

        # Every local point is transformed once by the typed event-view path
        # transform.  Inverting that transform must recover the transport point.
        if mode in {"scintillator", "combined"}:
            transform = gui.RigidAssemblyTransform(**payload["path_transforms"]["scintillator"])
            rotation = transform.rotation_matrix()
            origin = np.asarray(transform.rotation_origin_mm, dtype=float)
            translation = np.asarray((transform.translation_x_mm, transform.translation_y_mm, transform.translation_z_mm), dtype=float)
            for record in result.scintillator.photon_records:
                global_points = gui._event_path_points_global_mm(record, "scintillator", payload)
                for local_point, global_point in zip(record.path_points_m, global_points):
                    recovered = rotation.T @ (np.asarray(global_point, dtype=float) - origin - translation) + origin
                    _assert(np.allclose(recovered, np.asarray(local_point) * 1e3, atol=1e-9), "local scintillator path was not transformed exactly once")

    # The original false dialog came from asking the Cherenkov geometry for a
    # particle entry in a schema-only ROOT row.  The error is now output-only.
    bad_solid = simulator.SolidGeometry(
        shape="box", material="fused_silica", center_x_m=0.5,
        length_x_m=0.1, width_y_m=0.02, thickness_z_m=0.01,
    )
    bad_result = simulator.SimulationResult(
        geometry=simulator.CompoundGeometry((bad_solid,)), particle=particle,
        optical=simulator.OpticalConfig(), detector_mode="cherenkov",
    )
    try:
        simulator._root_event_entry_point(bad_result)
    except RuntimeError as exc:
        message = str(exc).lower()
        _assert("simulation results remain valid" in message and "optional root schema" in message, "ROOT geometry mismatch was not classified as output-only")
        _assert("intersection" not in message, "ROOT geometry mismatch retained the false physics/intersection wording")
    else:
        raise AssertionError("ROOT geometry mismatch unexpectedly produced an entry point")

    # Exercise the full GUI finalization route: ROOT may be unavailable, but
    # the completed state, summary, run log, and PNG products must remain.
    with tempfile.TemporaryDirectory(prefix="v44_headless_run_") as directory:
        values["output_directory"] = directory
        settings = gui.validate_settings(values)
        run_directory = Path(directory) / "quartz_run_debug"
        payload = gui.run_simulation_to_directory(settings, run_directory)
        _assert(payload["summary"]["detector_mode"] == "scintillator", "headless GUI run changed detector mode")
        _assert(payload["summary"]["generated_photons"] == 0, "Scintillator-only GUI run generated Cherenkov photons")
        _assert((run_directory / "summary.json").exists() and (run_directory / "run.log").exists(), "GUI finalization did not preserve JSON/log")
        _assert((run_directory / "event_display_event_0.png").exists(), "GUI finalization did not preserve event PNG")
        if payload["root_error"] is not None:
            root_message = str(payload["root_error"]).lower()
            _assert("optional output" in root_message and "simulation results remain valid" in root_message, "missing ROOT dependency was reported as a physics failure")
            _assert("root output unavailable" in (run_directory / "run.log").read_text(encoding="utf-8").lower(), "run.log did not identify the ROOT-only failure")
        else:
            _assert((run_directory / "quartz_cherenkov.root").exists(), "ROOT succeeded but the output file is missing")
        disk = json.loads((run_directory / "summary.json").read_text(encoding="utf-8"))
        _assert(disk["event_view_scene"]["visible_assembly_types"] == ["scintillator"], "saved event scene retained a non-participating assembly")


def scintillator_decay_time_checks() -> None:
    """Verify the editable decay mean without changing generation or spectrum."""

    import numpy as np
    import quartz_cherenkov_gui as gui
    import quartz_cherenkov_sim as simulator

    values = gui.default_values()
    _assert(values["scintillator_decay_time_ns"] == "2.0", "scintillator decay-time default is not 2.0 ns")
    values.update({
        "detector_run_mode": "Scintillator only",
        "scintillator_decay_time_ns": "0.8",
        "scintillator_photon_count": "512",
        "scintillator_deposited_energy_MeV": "1.0",
        "scintillator_guide_enabled": "0",
    })
    config = gui.scintillator_config_from_values(values)
    _assert(np.isclose(config.decay_time_s, 0.8e-9, rtol=0.0, atol=1e-24), "0.8 ns was not converted to 0.8e-9 s")
    settings = gui.validate_settings(values)
    _assert(np.isclose(settings.scintillator_decay_time_ns, 0.8, rtol=0.0, atol=1e-15), "Run snapshot did not retain the editable decay value")

    for invalid in ("-0.1", "nan", "inf"):
        invalid_values = dict(values)
        invalid_values["scintillator_decay_time_ns"] = invalid
        try:
            gui.validate_settings(invalid_values)
        except ValueError:
            pass
        else:
            raise AssertionError(f"invalid scintillator decay value {invalid!r} was allowed past pre-worker validation")

    # The standalone Scintillator profile already carries this field.  Keep a
    # round trip here so the v60 JSON envelope remains backward compatible.
    profile = gui.scintillator_profile_from_values(values)
    _assert(profile["scintillator_decay_time_ns"] == "0.8", "Scintillator profile omitted the decay field")
    restored_values = gui.default_values()
    gui.apply_scintillator_profile(restored_values, profile)
    _assert(restored_values["scintillator_decay_time_ns"] == "0.8", "Scintillator profile did not restore the decay field")
    _assert(np.isclose(gui.scintillator_config_from_values(restored_values).decay_time_s, 0.8e-9, rtol=0.0, atol=1e-24), "restored profile changed decay conversion")

    particle = simulator.Particle(
        mass_kg=0.51099895 * simulator.MEV / simulator.C_LIGHT ** 2,
        kinetic_energy_J=2.0 * simulator.MEV,
    )
    active = simulator.SolidGeometry(
        shape="box", material="scintillator", center_x_m=0.0,
        length_x_m=0.020, width_y_m=0.020, thickness_z_m=0.010,
        refractive_index_value=1.58,
    )
    assembly = simulator.ScintillatorAssembly(active, terminal_plus_surface="detector")
    common = dict(enabled=True, sampling_mode="fixed", photon_count=512, emission_wavelength_min_m=350e-9, emission_wavelength_max_m=550e-9, max_bounces=40)
    long_result = simulator.simulate_scintillator(particle, assembly, simulator.ScintillatorConfig(**common, decay_time_s=2.0e-9), seed=24680)
    short_result = simulator.simulate_scintillator(particle, assembly, simulator.ScintillatorConfig(**common, decay_time_s=0.8e-9), seed=24680)
    zero_result = simulator.simulate_scintillator(particle, assembly, simulator.ScintillatorConfig(**common, decay_time_s=0.0), seed=24680)
    _assert(len(long_result.photon_records) == len(short_result.photon_records) == len(zero_result.photon_records) == 512, "decay control changed fixed generated photon count")
    long_wavelengths = np.asarray([record.wavelength_m for record in long_result.photon_records])
    short_wavelengths = np.asarray([record.wavelength_m for record in short_result.photon_records])
    _assert(np.array_equal(long_wavelengths, short_wavelengths), "decay control changed the emission wavelength sample")
    speed = simulator.particle_beta(particle.mass_kg, particle.kinetic_energy_J) * simulator.C_LIGHT

    def delays(result):
        values = []
        for record in result.photon_records:
            source_z = float(record.emission_position_m[2])
            particle_time = (source_z - active.bounds_min_m[2]) / speed
            values.append(record.emission_time_s - particle_time)
        return np.asarray(values)

    long_delays, short_delays, zero_delays = delays(long_result), delays(short_result), delays(zero_result)
    _assert(np.all(long_delays >= 0.0) and np.all(short_delays >= 0.0), "exponential emission delays are negative")
    _assert(np.allclose(short_delays, long_delays * 0.4, rtol=0.0, atol=1e-18), "same-seed decay samples did not scale with the configured mean")
    _assert(float(np.mean(short_delays)) < float(np.mean(long_delays)), "shorter decay mean did not reduce the sampled mean")
    _assert(float(np.median(short_delays)) < float(np.median(long_delays)), "shorter decay mean did not reduce the sampled median")
    _assert(np.allclose(zero_delays, 0.0, rtol=0.0, atol=1e-18), "zero decay time did not mean instantaneous emission")

    # Cherenkov-only deliberately does not consume the Scintillator control;
    # its legacy compatibility field remains isolated from this run mode.
    cherenkov_values = gui.default_values()
    cherenkov_values.update({"detector_run_mode": "Cherenkov only", "scintillator_decay_time_ns": "-1"})
    _assert(gui.validate_settings(cherenkov_values).detector_run_mode == "cherenkov", "Cherenkov-only validation became coupled to scintillator decay")
    print("v61 scintillator decay control: visible default/conversion, profile round-trip, invalid pre-worker rejection, deterministic exponential scaling, zero semantics, wavelength/count invariance, and Cherenkov isolation PASS")


def scintillator_v45_invariant_checks() -> None:
    """Exercise the v45 source/readout/scene invariant with exact counts."""

    from dataclasses import replace
    import numpy as np
    import quartz_cherenkov_gui as gui
    import quartz_cherenkov_sim as simulator

    values = gui.default_values()
    values.update({
        "detector_run_mode": "Scintillator only",
        "scintillator_photon_count": "100",
        "scintillator_sampling_mode": "Fixed Monte Carlo test sample",
        "scintillator_deposited_energy_MeV": "2.0",
        "scintillator_guide_enabled": "1",
        "scintillator_readout_index": "1.49",
        "scintillator_readout_acceptance_deg": "90",
        "scintillator_readout_qe": "1.0",
        "scintillator_readout_collection_efficiency": "1.0",
    })
    geometry = simulator.BarGeometry()
    particle = gui._particle_for_settings(gui.validate_settings(values))
    local = gui.scintillator_local_geometry_from_values(values)
    config = gui.scintillator_config_from_values(values)
    assembly = local.to_assembly(mylar_enabled=True, mylar_reflectivity=0.95)

    # The screenshot-sized run: ten events, 100 fixed rows/event, seed 12345.
    results = simulator.simulate_events(
        particle, 10, geometry=geometry, optical=simulator.OpticalConfig(photon_count=100),
        seed=12345, scintillator_assembly=assembly, scintillator_config=config,
        detector_mode="scintillator",
    )
    _assert(len(results) == 10, "v45: event count is not 10")
    total_records = sum(len(item.scintillator.photon_records) for item in results)
    _assert(total_records == 1000, f"v45: expected 1000 accepted scintillator records, got {total_records}")
    aggregate = simulator.aggregate_summary(results, detector_mode="scintillator")
    _assert(aggregate["scintillator_source_rows"] == 1000, "v45: source-row count is not 1000")
    _assert(aggregate["scintillator_accepted_records"] == 1000, "v45: accepted-record count does not reconcile")
    _assert(aggregate["scintillator_excluded_records"] == 0, "v45: unexpected excluded records")
    _assert(aggregate["scintillator_detected_photons"] == sum(
        record.detected for item in results for record in item.scintillator.photon_records
    ), "v45: aggregate/table detected count disagrees with records")
    _assert(aggregate["scintillator_readout_alignment_mismatches"] == 0, "v45: readout endpoint mismatches were hidden")

    canonical = gui.global_assembly_scene(values, geometry)
    payload = gui.canonical_event_view_payload(canonical, "scintillator")
    event_scene = gui.build_interactive_event_scene(results[0], payload, detector_mode="scintillator", display_photons=100)
    _assert(event_scene["visible_assembly_types"] == ["scintillator"], "v45: single-detector scene leaked another assembly")
    active = local.scintillator
    path_transform = gui.RigidAssemblyTransform(**payload["path_transforms"]["scintillator"])
    rotation = path_transform.rotation_matrix()
    origin = np.asarray(path_transform.rotation_origin_mm, dtype=float)
    translation = np.asarray((path_transform.translation_x_mm, path_transform.translation_y_mm, path_transform.translation_z_mm), dtype=float)
    def inverse(point_mm):
        return rotation.T @ (np.asarray(point_mm, dtype=float) - origin - translation) + origin

    track_mid = inverse(np.mean(np.asarray(event_scene["track_mm"], dtype=float), axis=0)) * 1e-3
    _assert(active.contains(track_mid, tolerance=1e-9), "v45: displayed Scintillator track misses active volume")
    _assert(np.allclose(event_scene["impact_markers"][0]["point_mm"], [path_transform.apply_point_mm((0.0, 0.0, 0.0))[i] for i in range(3)], atol=1e-9), "v45: impact marker used the quartz frame")
    _assert(event_scene["per_detector"]["scintillator"]["record_count"] == 100, "v45: viewer record count changed")
    _assert(event_scene["per_detector"]["scintillator"]["selected_record_count"] == 100, "v45: viewer selection count changed")

    sensor_components = [item for item in event_scene["geometry_components"] if item["kind"] == "sensor"]
    _assert(len(sensor_components) == 1, "v45: Scintillator viewer has no unique sensor")
    detected_lines = [line for line in event_scene["photon_lines"] if line["outcome"].startswith("detected_")]
    detected_records = [record for record in results[0].scintillator.photon_records if record.detected]
    _assert(len(detected_lines) == len(detected_records), "v45: displayed detected paths do not reconcile with records")
    detected_lines_by_id = {int(line["photon_id"]): line for line in detected_lines}
    readout_faces = payload["assemblies"][0]["components"]
    sensor = next(item for item in readout_faces if item.get("kind") == "sensor")
    window = np.asarray(sensor["readout_faces_mm"][0], dtype=float)
    window_normal = np.cross(window[1] - window[0], window[2] - window[0])
    window_normal /= np.linalg.norm(window_normal)
    # Rebase the canonical sensor face through the same event geometry transform.
    rebased_sensor = gui._rebase_scene_component(sensor, payload)
    rebased_window = np.asarray(rebased_sensor["readout_faces_mm"][0], dtype=float)
    for record in detected_records:
        line = detected_lines_by_id[int(record.photon_id)]
        endpoint = np.asarray(line["points_mm"][-1], dtype=float)
        local_endpoint = inverse(endpoint)
        distance = simulator.scintillator_readout_endpoint_distance_m(record, assembly)
        _assert(distance <= results[0].scintillator.readout_endpoint_tolerance_m, "v45: detected endpoint missed physical readout")
        _assert(np.allclose(local_endpoint, np.asarray(record.path_points_m[-1]) * 1e3, atol=1e-9), "v45: photon path transformed more than once")
        plane_distance = abs(float(np.dot(endpoint - rebased_window[0], window_normal)))
        _assert(plane_distance <= 1e-7, "v45: displayed detected endpoint is not on displayed sensor window")

    # Active-only, rotated/nonzero placement, and the supported -X readout.
    for guide_enabled in (False, True):
        values["scintillator_guide_enabled"] = "1" if guide_enabled else "0"
        values["global_offset_x_mm"] = "13"
        values["global_offset_y_mm"] = "-4"
        values["global_offset_z_mm"] = "7"
        values["detector2_offset_x_mm"] = "21"
        values["detector2_offset_y_mm"] = "3"
        values["detector2_offset_z_mm"] = "-2"
        values["assembly_2_rotation_axis"] = "Y"
        values["assembly_2_rotation_angle_deg"] = "17"
        local_variant = gui.scintillator_local_geometry_from_values(values)
        variant = local_variant.to_assembly(mylar_enabled=False)
        one = simulator.simulate(
            particle, geometry=geometry, optical=simulator.OpticalConfig(photon_count=0), seed=77,
            scintillator_assembly=variant, scintillator_config=config, detector_mode="scintillator",
        )
        _assert(all(local_variant.scintillator.contains(record.emission_position_m, tolerance=1e-10) for record in one.scintillator.photon_records), "v45: source row escaped active volume")

    minus_values = dict(values)
    minus_values.update({"scintillator_terminal_plus_surface": "Open", "scintillator_terminal_minus_surface": "Detector/Photosensor coupling", "scintillator_guide_enabled": "0"})
    minus_local = gui.scintillator_local_geometry_from_values(minus_values)
    minus_assembly = replace(minus_local.to_assembly(mylar_enabled=False), terminal_plus_surface="Open", terminal_minus_surface="Detector/Photosensor coupling")
    minus_config = gui.scintillator_config_from_values(minus_values)
    minus_result = simulator.simulate(
        particle, geometry=geometry, optical=simulator.OpticalConfig(photon_count=0), seed=123,
        scintillator_assembly=minus_assembly, scintillator_config=minus_config, detector_mode="scintillator",
    )
    _assert(any(record.detected and record.detector_face == "minus_x" for record in minus_result.scintillator.photon_records), "v45: -X terminal readout is not supported")
    _assert(minus_result.scintillator.summary()["scintillator_readout_alignment_mismatches"] == 0, "v45: -X endpoint invariant failed")

    # Ownership matrix: Scintillator-only has no Cherenkov rows; All has both.
    for mode in ("cherenkov", "combined"):
        mode_result = simulator.simulate(
            particle, geometry=geometry, optical=simulator.OpticalConfig(photon_count=4), seed=19,
            scintillator_assembly=assembly, scintillator_config=config, detector_mode=mode,
        )
        if mode == "cherenkov":
            _assert(mode_result.scintillator is None, "v45: Cherenkov-only created scintillator records")
        else:
            _assert(mode_result.scintillator is not None and mode_result.photon_records, "v45: All lost one detector owner")

    all_records = [
        record for item in results for record in item.scintillator.photon_records
    ]
    detected_records_all = [record for record in all_records if record.detected]
    source_failures = sum(
        not active.contains(record.emission_position_m, tolerance=1e-10)
        for record in all_records
    )
    endpoint_distances_um = [
        simulator.scintillator_readout_endpoint_distance_m(record, assembly) * 1e6
        for record in detected_records_all
    ]
    max_endpoint_distance_um = max(endpoint_distances_um, default=0.0)
    max_bounce_records = sum(record.outcome == "max_bounces" for record in all_records)
    print(
        "v45 metrics: events=10 source_rows={source} accepted={accepted} excluded={excluded} "
        "detected={detected} table_detected={table} endpoint_mismatches={mismatches} "
        "max_endpoint_distance_mm={max_distance_mm:.12g} max_endpoint_distance_um={max_distance_um:.12g} "
        "tolerance_um={tolerance:.12g} source_containment_count={source_count} "
        "source_containment_failures={source_failures} active_to_guide={active_to_guide} "
        "guide_to_readout={guide_to_readout} interface_transmissions={transmissions} "
        "interface_reflections={reflections} tir_reflections={tir} mylar_reflections={mylar_reflections} "
        "mylar_losses={mylar_losses} absorbed={absorbed} escaped={escaped} max_bounce_records={max_bounces} "
        "event0_stored={stored} event0_selected={selected} event0_displayed_paths={paths} "
        "excluded_reasons=none mismatch_reasons=none "
        "transport_frame=scintillator_local_m transport_units=m scene_units=mm".format(
            source=aggregate["scintillator_source_rows"],
            accepted=aggregate["scintillator_accepted_records"],
            excluded=aggregate["scintillator_excluded_records"],
            detected=aggregate["scintillator_detected_photons"],
            table=aggregate["scintillator_detected_photons"],
            mismatches=aggregate["scintillator_readout_alignment_mismatches"],
            max_distance_mm=max_endpoint_distance_um / 1e3,
            max_distance_um=max_endpoint_distance_um,
            tolerance=aggregate["scintillator_readout_endpoint_tolerance_um"],
            source_count=len(all_records) - source_failures,
            source_failures=source_failures,
            active_to_guide=aggregate["scintillator_active_to_guide_crossings"],
            guide_to_readout=aggregate["scintillator_guide_to_readout_crossings"],
            transmissions=aggregate["scintillator_interface_transmissions"],
            reflections=aggregate["scintillator_interface_reflections"],
            tir=aggregate["scintillator_tir_reflections"],
            mylar_reflections=aggregate["scintillator_mylar_reflections"],
            mylar_losses=aggregate["scintillator_mylar_losses"],
            absorbed=aggregate["scintillator_absorbed"],
            escaped=aggregate["scintillator_escaped"],
            max_bounces=max_bounce_records,
            stored=event_scene["per_detector"]["scintillator"]["record_count"],
            selected=event_scene["per_detector"]["scintillator"]["selected_record_count"],
            paths=len(event_scene["photon_lines"]),
        )
    )


def scintillator_v46_complete_guide_transport_checks() -> None:
    """Verify canonical multi-section transport, records, and viewer counts."""

    import numpy as np
    import quartz_cherenkov_gui as gui
    import quartz_cherenkov_sim as simulator
    from quartz_geometry import AssemblyGeometryModel, AssemblySection

    values = gui.default_values()
    values.update({
        "detector_run_mode": "Scintillator only",
        "scintillator_guide_enabled": "1",
        "scintillator_photon_count": "120",
        "scintillator_sampling_mode": "Fixed Monte Carlo test sample",
        "scintillator_deposited_energy_MeV": "1.0",
        "scintillator_readout_index": "1.0",
        "scintillator_readout_acceptance_deg": "90",
        "scintillator_readout_qe": "1.0",
        "scintillator_readout_collection_efficiency": "1.0",
        "sensor_2_width_mm": "16.0", "sensor_2_height_mm": "8.0",
        "detector2_offset_x_mm": "13.0", "detector2_offset_y_mm": "-4.0", "detector2_offset_z_mm": "7.0",
        "assembly_2_rotation_axis": "Y", "assembly_2_rotation_angle_deg": "17.0",
    })
    active = AssemblySection(
        name="Active scintillator", material="scintillator", axis="X",
        center_x_mm=0.0, length_x_mm=20.0, width_y_mm=20.0, thickness_z_mm=10.0,
    )
    taper = AssemblySection(
        name="Guide taper", material="light guide", axis="X",
        shape="tapered_rectangular_prism", center_x_mm=30.0, length_x_mm=40.0,
        width_y_start_mm=20.0, width_y_end_mm=16.0,
        thickness_z_start_mm=10.0, thickness_z_end_mm=8.0,
    )
    terminal = AssemblySection(
        name="Guide terminal", material="fused silica", axis="X",
        center_x_mm=70.0, length_x_mm=40.0, width_y_mm=16.0, thickness_z_mm=8.0,
    )
    model = AssemblyGeometryModel(
        "Scintillator", [active], [taper, terminal], True,
        "Detector/Photosensor coupling", "Open", True, 0.95, 0.0, 0.0,
    )
    model.validate()
    particle = simulator.Particle(
        mass_kg=0.51099895 * simulator.MEV / simulator.C_LIGHT ** 2,
        kinetic_energy_J=2.0 * simulator.MEV,
    )
    assembly = gui.scintillator_assembly_from_model(model, values=values, particle=particle)
    config = gui.scintillator_config_from_values(values)
    result = simulator.simulate(
        particle, geometry=simulator.BarGeometry(), optical=simulator.OpticalConfig(photon_count=0),
        seed=4601, detector_mode="scintillator", scintillator_assembly=assembly,
        scintillator_config=config,
    )
    _assert(not result.photon_records, "v46: Scintillator-only produced Cherenkov photon records")
    scint = result.scintillator
    _assert(scint is not None, "v46: Scintillator result is missing")
    records = list(scint.photon_records)
    summary = scint.summary()
    _assert(scint.source_row_count == 120 and len(records) == 120, "v46: source/accepted record count mismatch")
    _assert(scint.readout_mismatch_count == 0, "v46: readout endpoint exclusions occurred")
    _assert(all(active.to_simulation_geometry().contains(record.emission_position_m, tolerance=1e-10) for record in records), "v46: source escaped active section")
    guide_sections = assembly.guide_sections
    guide_report = summary["scintillator_guide_sections"]
    material_indices_450nm = [float(section.refractive_index(450e-9)) for section in assembly.ordered_sections]
    _assert(len(guide_report) == 2, "v46: configured guide sections were not retained")
    for position, section in enumerate(guide_sections):
        report = guide_report[position]
        _assert(report["entries"] == sum(bool(record.guide_section_entry_counts[position]) for record in records), f"v46: {report['name']} entry count is not record-derived")
        _assert(report["path_segments_inside"] == sum(record.guide_section_path_segment_counts[position] for record in records), f"v46: {report['name']} path count is not record-derived")
        _assert(report["entries"] > 0 and report["path_segments_inside"] > 0, f"v46: photons did not physically traverse {report['name']}")
        for record in records:
            if position < len(record.guide_section_path_segment_counts) and record.guide_section_path_segment_counts[position] > 0:
                points = [np.asarray(point, dtype=float) for point in record.path_points_m]
                _assert(any(section.contains((left + right) / 2.0, tolerance=0.0) for left, right in zip(points, points[1:])), f"v46: a {report['name']} path count lacks an in-section segment")
    _assert(summary["scintillator_active_to_guide_crossings"] > 0, "v46: no active-to-guide transmitted interface")
    _assert(summary["scintillator_guide_to_readout_crossings"] == sum(record.guide_to_readout_crossings for record in records), "v46: guide-to-sensor count is not record-derived")
    _assert(summary["scintillator_guide_to_readout_crossings"] == sum(record.detector_transmitted for record in records), "v46: transmitted guide output count disagrees with sensor boundary records")
    detected = [record for record in records if record.detected]
    _assert(len(detected) == sum(record.detected for record in records), "v46: detected records do not reconcile")
    max_endpoint_m = max((simulator.scintillator_readout_endpoint_distance_m(record, assembly) for record in detected), default=0.0)
    _assert(max_endpoint_m <= scint.readout_endpoint_tolerance_m, "v46: detected endpoint missed the physical terminal sensor window")

    canonical = gui.global_assembly_scene(values, simulator.BarGeometry(), {"scintillator": model})
    payload = gui.canonical_event_view_payload(canonical, "scintillator")
    scene = gui.build_interactive_event_scene(result, payload, detector_mode="scintillator", display_photons=120)
    table = gui.event_table_rows([result])[0]
    _assert(scene["per_detector"]["scintillator"]["record_count"] == len(records), "v46: viewer record count disagrees with records")
    _assert(sum(line["outcome"].startswith("detected_") for line in scene["per_detector"]["scintillator"]["photon_lines"]) == len(detected), "v46: viewer detected count disagrees with records")
    _assert(table["scintillator_detected_photons"] == len(detected), "v46: event table detected count disagrees with records")
    _assert(all(line["detector_type"] == "scintillator" for line in scene["photon_lines"]), "v46: Scintillator-only viewer contains Cherenkov paths")

    transform = gui.RigidAssemblyTransform(**payload["path_transforms"]["scintillator"])
    rotation = transform.rotation_matrix()
    origin = np.asarray(transform.rotation_origin_mm, dtype=float)
    translation = np.asarray((transform.translation_x_mm, transform.translation_y_mm, transform.translation_z_mm), dtype=float)
    for record in records[:8]:
        global_points = gui._event_path_points_global_mm(record, "scintillator", payload)
        for local_point, global_point in zip(record.path_points_m, global_points):
            recovered = rotation.T @ (np.asarray(global_point, dtype=float) - origin - translation) + origin
            _assert(np.allclose(recovered, np.asarray(local_point) * 1e3, atol=1e-9), "v46: nonzero global transform was applied more than once")
    _assert(transform.rotation_angle_deg == 17.0 and any(abs(value) > 0.0 for value in (transform.translation_x_mm, transform.translation_y_mm, transform.translation_z_mm)), "v46: nonzero offset/rotation case was not exercised")

    print(
        "v46 verification: mode=Scintillator-only transport_frame=scintillator_local_m "
        "transport_units=m scene_units=mm source_photons={source} accepted_records={accepted} "
        "excluded={excluded} excluded_reasons=none guide_sections={guides} "
        "active_to_guide_transmitted={active_to_guide} interface_transmissions={transmissions} "
        "interface_reflections={reflections} tir_reflections={tir} absorbed={absorbed} escaped={escaped} "
        "guide_to_sensor={guide_to_sensor} detected={detected} max_endpoint_distance_um={endpoint:.12g} "
        "endpoint_tolerance_um={tolerance:.12g} source_containment={containment} "
        "no_cherenkov_records=True table_detected={table} viewer_records={viewer} "
        "global_offset_rotation=verified indices_at_450nm={indices} config=Active scintillator[box 20x20x10 mm scintillator] -> "
        "Guide taper[tapered_rectangular_prism 40 mm, 20x10 to 16x8 mm light_guide] -> "
        "Guide terminal[box 40x16x8 mm fused_silica] -> +X sensor".format(
            source=scint.source_row_count, accepted=len(records), excluded=scint.readout_mismatch_count,
            guides=";".join(f"{item['name']}:entries={item['entries']},path_segments={item['path_segments_inside']}" for item in guide_report),
            active_to_guide=summary["scintillator_active_to_guide_crossings"],
            transmissions=summary["scintillator_interface_transmissions"],
            reflections=summary["scintillator_interface_reflections"], tir=summary["scintillator_tir_reflections"],
            absorbed=summary["scintillator_absorbed"], escaped=summary["scintillator_escaped"],
            guide_to_sensor=summary["scintillator_guide_to_readout_crossings"], detected=len(detected),
            endpoint=max_endpoint_m * 1e6, tolerance=scint.readout_endpoint_tolerance_m * 1e6,
            containment=sum(active.to_simulation_geometry().contains(record.emission_position_m, tolerance=1e-10) for record in records),
            table=table["scintillator_detected_photons"], viewer=scene["per_detector"]["scintillator"]["record_count"],
            indices=",".join(f"{value:.9g}" for value in material_indices_450nm),
        )
    )


def run_status_layout_contract_checks() -> None:
    import quartz_cherenkov_gui as gui

    long_text = "Building global geometry and event view · Scintillator only · elapsed 123.4 s · waiting for readout reconciliation"
    layout = gui.run_status_panel_layout_contract(text=long_text)
    _assert(layout["fixed_width_cells"] and not layout["overlap_with_simulation_parameters"], "v45: Run status cells are not fixed")
    _assert(layout["status_cell_width_px"] == layout["progress_cell_width_px"], "v45: status/progress cells have different widths")
    _assert(len(layout["wrapped_lines"]) > 1, "v45: long Run status text was not wrapped")
    source = (ROOT / "quartz_cherenkov_gui.py").read_text(encoding="utf-8")
    _assert("action_frame.grid_propagate(False)" in source and "status_cell.grid_propagate(False)" in source, "v52: fixed Tk Run panel layout guard is missing")
    _assert(layout["progress_row"] == 2 and layout["status_row_height_px"] == 34 and not layout["dynamic_wrap_height"], "v52: status/progress rows are not fixed")


def complete_global_profile_checks() -> None:
    """Round-trip both canonical models and prove invalid loads are atomic."""

    import quartz_cherenkov_gui as gui
    from quartz_geometry import AssemblyGeometryModel, AssemblySection

    values = gui.default_values()
    values.update({
        "sensor_1_width_mm": "18", "sensor_1_height_mm": "8", "sensor_1_depth_mm": "3", "sensor_1_qe": "0.73",
        "sensor_2_width_mm": "16", "sensor_2_height_mm": "8", "sensor_2_depth_mm": "2", "sensor_2_qe": "0.82",
        "global_offset_x_mm": "1.25", "global_offset_y_mm": "2.5", "global_offset_z_mm": "3.75",
        "detector2_offset_x_mm": "-4", "detector2_offset_y_mm": "5", "detector2_offset_z_mm": "6",
        "assembly_1_rotation_axis": "Y", "assembly_1_rotation_angle_deg": "17",
        "assembly_2_rotation_axis": "Z", "assembly_2_rotation_angle_deg": "-9",
        "scintillator_gap_mm": "12", "detector_separation_mm": "33",
    })
    cherenkov = AssemblyGeometryModel(
        "Cherenkov", [AssemblySection("Q1", length_x_mm=50, center_x_mm=-25), AssemblySection("Q2", length_x_mm=50, center_x_mm=25)],
        [AssemblySection("G1", material="light guide", length_x_mm=20, center_x_mm=60), AssemblySection("G2", material="fused silica", length_x_mm=20, center_x_mm=80)],
        True, "detector", "mirror", True, .91, 3.25, -1.5,
    )
    scintillator = AssemblyGeometryModel(
        "Scintillator", [AssemblySection("S1", material="scintillator", length_x_mm=20, center_x_mm=-10), AssemblySection("S2", material="scintillator", length_x_mm=20, center_x_mm=10)],
        [AssemblySection("SG1", material="light guide", length_x_mm=20, center_x_mm=30), AssemblySection("SG2", material="fused silica", length_x_mm=20, center_x_mm=50)],
        True, "detector", "open", True, .88, -2.5, 2.25,
    )
    models = {"cherenkov": cherenkov, "scintillator": scintillator}
    profile = gui.complete_global_assembly_profile_from_state(values, models)
    _assert(profile["format"] == "beam-monitor-global-assembly-profile-v1" and profile["schema_version"] == 1, "v52 complete profile marker/version missing")
    _assert(profile["assembly_order"] == ["cherenkov", "scintillator"] and profile["units"]["length"] == "mm", "v52 profile identity/units missing")
    restored_values = gui.default_values()
    restored_models = {
        "cherenkov": AssemblyGeometryModel("Cherenkov", [AssemblySection(name="default")]),
        "scintillator": AssemblyGeometryModel("Scintillator", [AssemblySection(name="default", material="scintillator", length_x_mm=20, width_y_mm=20, thickness_z_mm=10, axis="X")]),
    }
    gui.apply_complete_global_assembly_profile(restored_values, restored_models, json.loads(json.dumps(profile)))
    _assert(restored_models["cherenkov"].to_dict() == cherenkov.to_dict(), "v52 Cherenkov complete profile round-trip changed geometry")
    _assert(restored_models["scintillator"].to_dict() == scintillator.to_dict(), "v52 Scintillator complete profile round-trip changed geometry")
    for key in ("sensor_1_qe", "sensor_2_qe", "global_offset_x_mm", "detector2_offset_x_mm", "assembly_1_rotation_axis", "assembly_2_rotation_angle_deg"):
        _assert(restored_values[key] == values[key], f"v52 profile did not recover exact value {key}")
    _assert(restored_models["scintillator"].readout_status()["valid"], "v52 loaded scintillator readout placement is not valid")
    baseline = (restored_models["cherenkov"].to_dict(), restored_models["scintillator"].to_dict(), dict(restored_values))
    invalid_profiles = []
    bad_units = copy.deepcopy(profile); bad_units["units"]["length"] = "m"; invalid_profiles.append(bad_units)
    bad_geometry = copy.deepcopy(profile); bad_geometry["assemblies"]["cherenkov"]["geometry"]["sections"] = []; invalid_profiles.append(bad_geometry)
    bad_assemblies = copy.deepcopy(profile); del bad_assemblies["assemblies"]["scintillator"]; invalid_profiles.append(bad_assemblies)
    for invalid in invalid_profiles:
        try:
            gui.apply_complete_global_assembly_profile(restored_values, restored_models, invalid)
        except (TypeError, ValueError):
            pass
        else:
            raise AssertionError("v52 invalid complete profile was accepted")
        _assert((restored_models["cherenkov"].to_dict(), restored_models["scintillator"].to_dict(), dict(restored_values)) == baseline, "v52 invalid profile partially mutated state")
    with tempfile.TemporaryDirectory() as directory:
        corrupt = Path(directory) / "corrupt.json"
        corrupt.write_text("{bad", encoding="utf-8")
        try:
            gui.load_complete_global_assembly_profile(corrupt, restored_values, restored_models)
        except (OSError, json.JSONDecodeError, TypeError, ValueError):
            pass
        else:
            raise AssertionError("v52 corrupt profile was accepted")
    _assert((restored_models["cherenkov"].to_dict(), restored_models["scintillator"].to_dict(), dict(restored_values)) == baseline, "v52 corrupt profile changed state")
    # Exercise the actual file round-trip into an independent fresh session.
    # This protects the v66 contract that both complete assemblies and their
    # placements are restored when the simulator starts with defaults.
    with tempfile.TemporaryDirectory() as directory:
        saved_profile = Path(directory) / "complete-global-assembly.json"
        gui.save_complete_global_assembly_profile(saved_profile, values, models)
        fresh_values = gui.default_values()
        fresh_models = {
            "cherenkov": AssemblyGeometryModel("Cherenkov", [AssemblySection(name="fresh quartz")]),
            "scintillator": AssemblyGeometryModel("Scintillator", [AssemblySection(name="fresh scintillator", material="scintillator", length_x_mm=20, width_y_mm=20, thickness_z_mm=10, axis="X")]),
        }
        gui.load_complete_global_assembly_profile(saved_profile, fresh_values, fresh_models)
        _assert(fresh_models["cherenkov"].to_dict() == cherenkov.to_dict(), "v66 fresh-session restore changed the complete Cherenkov assembly")
        _assert(fresh_models["scintillator"].to_dict() == scintillator.to_dict(), "v66 fresh-session restore changed the complete Scintillator assembly")
        for key in ("global_offset_x_mm", "global_offset_y_mm", "global_offset_z_mm", "detector2_offset_x_mm", "detector2_offset_y_mm", "detector2_offset_z_mm", "assembly_1_rotation_axis", "assembly_1_rotation_angle_deg", "assembly_2_rotation_axis", "assembly_2_rotation_angle_deg"):
            _assert(fresh_values[key] == values[key], f"v66 fresh-session restore lost position/rotation {key}")
    print("v66 complete profile: fresh-session file round-trip=both assemblies/ordered guides/readouts/positions/transforms exact; invalid loads atomic PASS")


def run_bounce_and_layout_checks() -> None:
    import quartz_cherenkov_gui as gui

    _assert(int(gui.default_values()["maximum_photon_bounces"]) >= 1000, "v52 bounce default is below 1000")
    for mode in ("Cherenkov only", "Scintillator only", "All"):
        values = gui.default_values()
        values.update({"detector_run_mode": mode, "maximum_photon_bounces": "123", "events": "1", "photons_per_event": "1", "scintillator_photon_count": "1", "scintillator_deposited_energy_MeV": "1"})
        settings = gui.validate_settings(values)
        _assert(settings.maximum_photon_bounces == 123 and settings.scintillator_max_bounces == 123, f"v52 {mode} typed run override missing")
        _assert(gui._optical_for_settings(settings).max_bounces == 123, f"v52 {mode} Cherenkov override missing")
        _assert(gui.scintillator_config_from_values(values, max_bounces_override=settings.maximum_photon_bounces).max_bounces == 123, f"v52 {mode} Scintillator override missing")
    for invalid in ("0", "-1", "1.5", "not-a-number"):
        values = gui.default_values(); values["maximum_photon_bounces"] = invalid
        try:
            gui.validate_settings(values)
        except ValueError:
            pass
        else:
            raise AssertionError(f"v52 invalid bounce value accepted: {invalid}")
    for phase in ("Building global geometry and event view | Scintillator only | elapsed 1234 s", "Completed (All)", "Error"):
        layout = gui.run_status_panel_layout_contract(text=phase)
        _assert(layout["progress_row"] == 2 and layout["row_heights_px"] == (32, 34, 22, 34, 30, 28), f"v52 layout changed for {phase}")
    print("v52 Run controls/layout: effective_bounces=C:123,S:123,All:123 invalid=0,-1,1.5,text rejected fixed_status_rows=verified PASS")


def responsive_layout_checks() -> None:
    import quartz_cherenkov_gui as gui

    for width, height in ((240, 360), (900, 420), (2400, 720)):
        bounds = gui.responsive_three_view_panel_bounds(width, height)
        _assert(len(bounds) == 3, f"{width}: expected three view rectangles")
        _assert(all(rect[2] >= rect[0] and rect[3] >= rect[1] for rect in bounds), f"{width}: invalid rectangle")
        widths = [rect[2] - rect[0] for rect in bounds]
        _assert(max(widths) - min(widths) < 1e-9, f"{width}: columns are not equal")
        gaps = [bounds[index + 1][0] - bounds[index][2] for index in range(2)]
        _assert(all(gap >= 0.0 for gap in gaps), f"{width}: overlapping columns")
        _assert(abs(bounds[0][0] - 2.0) < 1e-9 and abs(bounds[-1][2] - (width - 2.0)) < 1e-9, f"{width}: width allocation changed")
        _assert(abs(sum(widths) + sum(gaps) - (width - 4.0)) < 1e-9, f"{width}: columns do not cover allocated width")
        _assert(all(abs(rect[1] - 2.0) < 1e-9 and abs(rect[3] - (height - 2.0)) < 1e-9 for rect in bounds), f"{width}: height allocation changed")


def collision_aware_label_layout_checks() -> None:
    """Verify deterministic, collision-free callouts for global/event cameras."""

    import math
    import quartz_cherenkov_gui as gui

    values = gui.default_values()
    values.update({"cherenkov_guide_enabled": "1", "scintillator_guide_enabled": "1"})
    scene = gui.global_assembly_scene(values, gui.simulator.BarGeometry())
    _assert(sum(item.get("kind") == "sensor" for item in scene["render_components"]) == 2, "representative scene lacks two sensors")
    _assert(any(item.get("kind") == "mylar_overlay" for item in scene["render_components"]), "representative scene lacks Mylar")

    def layout_for(yaw, pitch, width, height):
        fit = gui.global_assembly_fit(scene, width, height)
        center, scale = fit["center_mm"], fit["scale_px_per_mm"]
        cy, sy, cp, sp = math.cos(yaw), math.sin(yaw), math.cos(pitch), math.sin(pitch)

        def project(point):
            x, y, z = (float(point[index]) - center[index] for index in range(3))
            rx, ry = cy * x - sy * y, sy * x + cy * y
            return width / 2.0 + rx * scale, height / 2.0 - (cp * ry - sp * z) * scale

        items = []
        for component in scene["render_components"]:
            points = [project(point) for face in component["faces_mm"] for point in face]
            target = (sum(point[0] for point in points) / len(points), sum(point[1] for point in points) / len(points))
            items.append({
                "component_id": component.get("component_id", component.get("name")),
                "assembly_type": component.get("assembly_type"), "kind": component.get("kind"),
                "label": gui.detector_callout_label(component), "detail": component.get("label"),
                "target": target,
            })
        return items, gui.collision_aware_callout_layout(
            items, width, height,
            occupied_rectangles=(gui.scene_legend_rect(width, height), (0, height - 34, width, height)),
        )

    observed_positions = set()
    for width, height in ((1120, 650), (760, 480), (500, 350), (320, 240)):
        for yaw, pitch in ((0.65, 0.55), (0.0, 0.0), (math.pi / 2.0, 1.5), (2.2, -1.0)):
            items, layout = layout_for(yaw, pitch, width, height)
            rectangles = [item["rect"] for item in layout["callouts"]]
            for index, first in enumerate(rectangles):
                for second in rectangles[index + 1:]:
                    _assert(not gui._rectangles_overlap(first, second, gap=0.0), f"overlapping callout rectangles at {width}×{height}")
            positions = [tuple(item["position"]) for item in layout["callouts"]]
            _assert(len(positions) == len(set(positions)), "two callouts share an exact screen anchor")
            observed_positions.update(positions)
            legend_text = "\n".join(gui.compact_scene_legend_lines(items, layout["omitted"]))
            for omitted in layout["omitted"]:
                _assert(str(omitted.get("detail")) in legend_text, "omitted callout detail was not retained in legend")
    _assert(observed_positions, "label layout produced no positions")

    # A deliberately dense test proves the all-or-legend fallback and the
    # labels-off toggle without depending on Tk's font metrics.
    dense = [{"assembly_type": "cherenkov", "kind": "sensor", "label": "C sensor", "detail": f"sensor detail {i}", "target": (100, 100)} for i in range(30)]
    dense_layout = gui.collision_aware_callout_layout(dense, 180, 120, occupied_rectangles=(gui.scene_legend_rect(180, 120),))
    _assert(dense_layout["omitted"], "dense layout did not omit labels when space was exhausted")
    _assert(all(item["detail"] in "\n".join(gui.compact_scene_legend_lines(dense, dense_layout["omitted"])) for item in dense_layout["omitted"]), "dense omitted detail was lost")
    hidden = gui.collision_aware_callout_layout(dense[:3], 640, 360, show_labels=False)
    _assert(not hidden["callouts"] and len(hidden["omitted"]) == 3, "compact label toggle did not hide callouts")


def _model(geometry, *, guide=False, plus="Detector/Photosensor coupling", minus="Open"):
    active = geometry.AssemblySection(name="Active", material="fused silica", axis="X")
    guides = [geometry.AssemblySection(name="Guide 1", material="light guide", axis="X", length_x_mm=20.0, center_x_mm=60.0)] if guide else []
    return geometry.AssemblyGeometryModel("Cherenkov", [active], guides, guide, plus, minus)


def model_checks() -> None:
    import quartz_geometry as geometry
    import quartz_geometry_editor as editor
    import quartz_cherenkov_gui as gui
    import quartz_cherenkov_sim as simulator

    _assert(gui.APP_VERSION == "0.46.0", "v78 application version marker missing")
    _assert(editor.canonical_assembly_model_from_sections([geometry.AssemblySection(name="active")]).has_light_guide is False, "editor default unexpectedly has a guide")

    direct = _model(geometry)
    direct_scene = direct.complete_scene()
    _assert(direct.readout_status()["optical_path"] == ("Active quartz radiator", "Photosensor"), "direct optical path is wrong")
    _assert(len(_sensor_items(direct_scene)) == 1 and _sensor_items(direct_scene)[0]["coupling_position_mm"][0] == 50.0, "direct sensor placement is wrong")

    guided = _model(geometry, guide=True)
    guided_scene = guided.complete_scene()
    guided_sensor = _sensor_items(guided_scene)
    _assert(guided.readout_status()["optical_path"] == ("Active quartz radiator", "Complete light guide", "Photosensor"), "derived guide path is wrong")
    _assert(len(guided_sensor) == 1 and guided_sensor[0]["coupling_position_mm"][0] == 70.0, "guide output sensor placement is wrong")
    _assert(guided_scene["sensor_routing"]["readout_terminal"] == "plus_x", "selected readout terminal missing")

    no_readout = _model(geometry, guide=True, plus="Open", minus="Mirror")
    no_readout_scene = no_readout.complete_scene()
    _assert(no_readout.readout_status()["code"] == "guide_present_no_readout", "guide-without-readout status missing")
    _assert(not _sensor_items(no_readout_scene), "guide without readout invented a sensor")
    _assert({item["state"] for item in no_readout.to_simulation_geometry().terminal_interfaces} == {"open", "mirror"}, "mirror/open transport state changed")

    mirror = _model(geometry, plus="Mirror")
    open_model = _model(geometry, plus="Open")
    _assert(not _sensor_items(mirror.complete_scene()) and not _sensor_items(open_model.complete_scene()), "mirror/open invented a sensor")
    _assert(any(item.get("kind") == "mirror" for item in mirror.complete_scene()["components"]), "mirror marker missing")

    legacy = geometry.AssemblyGeometryModel.from_dict({
        "name": "Legacy", "sections": [geometry.AssemblySection(name="Active", material="fused silica", axis="X").to_dict()],
        "guide_sections": [geometry.AssemblySection(name="Guide 1", material="light guide", axis="X", length_x_mm=20.0, center_x_mm=60.0).to_dict()],
        "guide_enabled": True, "terminal_plus_surface": "Optical/light-guide coupling", "terminal_minus_surface": "Open",
    })
    _assert(legacy.terminal_plus_surface == "Detector/Photosensor coupling" and legacy.readout_status()["mode"] == "light_guide", "old terminal JSON did not migrate")
    migrated_sensor = simulator.DetectorSensor.from_dict({"component": "guide_" + "1A", "coupling_mode": "light_guide", "width_mm": 20, "height_mm": 10})
    _assert(migrated_sensor.component == "scintillator_1", "old sensor component did not migrate to active component")

    transport = guided.to_simulation_geometry()
    plus = next(item for item in transport.terminal_interfaces if item["face"] == "plus_x")
    _assert(plus["state"] == guided_scene["terminal_states"]["plus_x"] == "detector", "scene/transport terminal disagreement")
    _assert(set(guided_scene["views"]) == {"xz", "xy", "yz"}, "responsive three-view scene contract changed")
    _assert(gui.default_values()["terminal_plus_surface"] == "Detector/Photosensor coupling", "default terminal choice changed")
    simulator._self_test()


def sensor_detection_accounting_checks() -> None:
    """Exercise real terminal transport with deterministic finite-window rays."""

    import math
    import numpy as np
    import quartz_cherenkov_sim as simulator
    import quartz_geometry as geometry

    wavelength_m = 450e-9

    def make_model(face: str, guide: bool = False):
        active = geometry.AssemblySection(
            name="Active quartz", material="fused silica", axis="X",
            length_x_mm=100.0, width_y_mm=20.0, thickness_z_mm=10.0,
        )
        guides = [geometry.AssemblySection(
            name="Guide output", material="light guide", axis="X",
            length_x_mm=20.0, width_y_mm=20.0, thickness_z_mm=10.0,
            center_x_mm=60.0,
        )] if guide else []
        plus = "Detector/Photosensor coupling" if face == "plus_x" else "Open"
        minus = "Detector/Photosensor coupling" if face == "minus_x" else "Open"
        model = geometry.AssemblyGeometryModel(
            "Cherenkov", [active], guides, guide, plus, minus,
        )
        transport = model.to_simulation_geometry(
            sensor_shape="rectangular", sensor_width_mm=4.0,
            sensor_height_mm=2.0, sensor_depth_mm=2.0,
        )
        model.validate_complete_assembly(sensor_shape="rectangular", sensor_width_mm=4.0, sensor_height_mm=2.0, sensor_depth_mm=2.0)
        return model, transport

    def optical(transport, face: str, *, acceptance=math.pi, qe=1.0, collection=1.0):
        if isinstance(transport, simulator.CompoundGeometry):
            section_index = len(transport.sections) - 1 if face == "plus_x" else 0
            index = float(transport.refractive_index_at(
                np.array([transport.sections[section_index].center_x_m, 0.0, 0.0]),
                wavelength_m, section_index,
            ))
        else:
            index = float(transport.refractive_index(wavelength_m))
        return simulator.OpticalConfig(
            positive_end_mode="detector" if face == "plus_x" else "open",
            opposite_end_mode="detector" if face == "minus_x" else "open",
            detector_index=index, detector_acceptance_rad=acceptance,
            detector_quantum_efficiency=qe, detector_collection_efficiency=collection,
            opposite_detector_index=index, opposite_detector_acceptance_rad=acceptance,
            opposite_detector_quantum_efficiency=qe, opposite_detector_collection_efficiency=collection,
            # These legacy values intentionally differ from the canonical
            # 4 x 2 mm sensor to prove the adapter uses the physical window.
            detector_active_width_m=0.020, detector_active_height_m=0.010,
            opposite_detector_active_width_m=0.020, opposite_detector_active_height_m=0.010,
            max_bounces=2,
        )

    def shoot(transport, face, position, direction, config, photon_id):
        return simulator._transport_photon(
            photon_id, wavelength_m, np.asarray(position, dtype=float),
            np.asarray(direction, dtype=float), 0.0, 0.99, transport,
            config, np.random.default_rng(1000 + photon_id),
            coordinate_frame="cherenkov_local_m",
        )

    model_plus, plus = make_model("plus_x")
    ideal_plus = optical(plus, "plus_x")
    records = [
        shoot(plus, "plus_x", (0.0, 0.0, 0.0), (1.0, 0.0, 0.0), ideal_plus, 0),
        shoot(plus, "plus_x", (0.0, 0.001, 0.0005), (1.0, 0.0, 0.0), ideal_plus, 1),
        # Exactly on the finite window boundary exercises the numerical margin.
        shoot(plus, "plus_x", (0.0, 0.0020000000005, 0.0), (1.0, 0.0, 0.0), ideal_plus, 2),
        # This reaches the physical portal but misses the configured 4 mm window.
        shoot(plus, "plus_x", (0.0, 0.003, 0.0), (1.0, 0.0, 0.0), ideal_plus, 3),
    ]
    _assert([record.outcome for record in records[:3]] == ["detected_plus_x"] * 3, "ideal +X finite-window rays were not detected")
    _assert(records[3].outcome_reason == "sensor_window_miss" and records[3].detector_transmitted, "finite-window miss was not distinguished from a face crossing")
    _assert(all(record.detector_valid_sensor_hit for record in records[:3]), "valid sensor hit flag disagrees with ideal detections")

    # The opposite terminal uses the same finite-window authority and its own
    # normal/face mapping.
    model_minus, minus = make_model("minus_x")
    ideal_minus = optical(minus, "minus_x")
    minus_record = shoot(minus, "minus_x", (0.0, -0.001, 0.0005), (-1.0, 0.0, 0.0), ideal_minus, 4)
    _assert(minus_record.outcome == "detected_minus_x" and minus_record.detector_face == "minus_x", "ideal -X detector mapping failed")

    # Post-window filters remain physically active under otherwise ideal setup.
    angle_record = shoot(plus, "plus_x", (0.049, 0.0, 0.0), (math.cos(math.radians(10.0)), math.sin(math.radians(10.0)), 0.0), optical(plus, "plus_x", acceptance=math.radians(5.0)), 5)
    collection_record = shoot(plus, "plus_x", (0.0, 0.0, 0.0), (1.0, 0.0, 0.0), optical(plus, "plus_x", collection=0.0), 6)
    qe_record = shoot(plus, "plus_x", (0.0, 0.0, 0.0), (1.0, 0.0, 0.0), optical(plus, "plus_x", qe=0.0), 7)
    _assert(angle_record.outcome_reason == "acceptance_reject" and angle_record.detector_valid_sensor_hit, "angle rejection bypassed valid-window accounting")
    _assert(collection_record.outcome_reason == "collection_reject" and collection_record.detector_valid_sensor_hit, "collection rejection accounting failed")
    _assert(qe_record.outcome_reason == "qe_reject" and qe_record.detector_valid_sensor_hit, "QE rejection accounting failed")

    # Guide output uses the same active window at the actual complete-guide
    # output plane, not the active-section terminal.
    guided_model, guided = make_model("plus_x", guide=True)
    guided_record = shoot(guided, "plus_x", (0.0, 0.001, 0.0), (1.0, 0.0, 0.0), optical(guided, "plus_x"), 8)
    _assert(guided_record.outcome == "detected_plus_x" and guided_record.reached_readout, "guide-output detector did not detect at the physical output")
    _assert(abs(guided_record.final_position_m[0] - 0.070) < 1e-12, "guide detector endpoint is not the complete guide output")

    # Preview and rigid placement use the same canonical window plane. The
    # optical record remains local metres; the scene transform is global mm.
    scene = model_plus.complete_scene(sensor_width_mm=4.0, sensor_height_mm=2.0)
    sensor = next(item for item in scene["components"] if item.get("kind") == "sensor")
    _assert(sensor["window_faces_mm"] and sensor["dimensions_mm"]["width"] == 4.0, "sensor preview window is not the configured finite window")
    transform = geometry.RigidAssemblyTransform(
        translation_x_mm=17.0, translation_y_mm=-9.0, translation_z_mm=23.0,
        rotation_axis="Z", rotation_angle_deg=23.0,
    )
    local_endpoint_mm = np.asarray(records[0].final_position_m) * 1e3
    global_endpoint_mm = transform.apply_point_mm(local_endpoint_mm)
    global_window_mm = transform.apply_point_mm(np.asarray(sensor["coupling_position_mm"]))
    _assert(abs(global_endpoint_mm[0] - global_window_mm[0]) < 1e-9, "global translation/rotation moved endpoint off sensor plane")

    all_records = records + [minus_record, angle_record, collection_record, qe_record, guided_record]
    _assert(all(record.outcome and record.outcome_reason for record in all_records), "a generated record lacks one final outcome reason")
    reason_counts = simulator.detector_outcome_reason_counts(all_records)
    accepted = [record for record in all_records if record.detected]
    excluded = [record for record in all_records if not record.detected]
    counts = simulator.detector_diagnostic_counts(all_records)
    _assert(len(accepted) == sum(record.outcome_reason == "detected" for record in all_records), "detected count/reason count disagrees")
    _assert(counts["detector_valid_sensor_hits"] == sum(record.detector_valid_sensor_hit for record in all_records), "valid-window counter disagrees with records")
    aggregate = simulator.aggregate_summary([
        simulator.SimulationResult(geometry=plus, particle=simulator.Particle(mass_kg=1.0, kinetic_energy_J=1.0), optical=ideal_plus,
                                   photon_records=all_records, generated_photons=len(all_records), detector_mode="cherenkov")
    ])
    table_counts = {key: aggregate.get(key, 0) for key in ("detector_attempts", "detector_face_crossings", "valid_sensor_hits", "detected_photons", "detector_losses")}
    _assert(table_counts["detected_photons"] == len(accepted), "detection table aggregate does not reconcile with records")
    expected_planes = {id(minus_record): -0.05, id(guided_record): 0.07}
    max_endpoint_distance_mm = max(
        abs(float(record.path_points_m[-1][0]) - expected_planes.get(id(record), 0.05)) * 1e3
        for record in accepted if record.path_points_m
    )
    print(
        "v51 detection metrics: source_records={source} accepted={accepted} excluded={excluded} "
        "excluded_rows={rows} coordinate_frame=cherenkov_local_m transport_units=m scene_units=mm "
        "active_window_hit_count={hits} detected_count={detected} other_outcome_counts={reasons} "
        "max_path_endpoint_to_window_distance_mm={distance:.12g} detection_table_aggregate={table} "
        "configured_acceptance_deg=180 configured_qe=1 configured_collection=1 direct_faces=plus_x,minus_x guide_output=verified "
        "preview=verified global_translation_rotation=verified".format(
            source=len(all_records), accepted=len(accepted), excluded=len(excluded),
            rows=[(record.photon_id, record.outcome_reason) for record in excluded],
            hits=counts["detector_valid_sensor_hits"], detected=len(accepted),
            reasons=reason_counts, distance=max_endpoint_distance_mm, table=table_counts,
        )
    )


def complete_v32_json_and_diagnostic_checks() -> None:
    """Round-trip both complete detector assemblies and prove isolation."""

    import copy
    import json
    import quartz_geometry as geometry
    import quartz_cherenkov_gui as gui

    def make(kind: str):
        material = "fused silica" if kind == "cherenkov" else "scintillator"
        active = [
            geometry.AssemblySection(name="Active taper", material=material, shape="tapered_rectangular_prism", axis="X", length_x_mm=20, width_y_start_mm=10, width_y_end_mm=12, thickness_z_start_mm=6, thickness_z_end_mm=8, center_x_mm=-10),
            geometry.AssemblySection(name="Active block", material=material, shape="box", axis="X", length_x_mm=20, width_y_mm=12, thickness_z_mm=8, center_x_mm=10),
        ]
        guides = [
            geometry.AssemblySection(name="Guide taper", material="light guide", shape="tapered_rectangular_prism", axis="X", length_x_mm=10, width_y_start_mm=12, width_y_end_mm=14, thickness_z_start_mm=8, thickness_z_end_mm=9, center_x_mm=25),
            geometry.AssemblySection(name="Guide output", material="light guide", shape="box", axis="X", length_x_mm=10, width_y_mm=14, thickness_z_mm=9, center_x_mm=35),
        ]
        return geometry.AssemblyGeometryModel(kind.title(), active, guides, True, "Detector/Photosensor coupling", "Mirror", True, .91, 2.25, -1.75)

    values = gui.default_values()
    values.update({"sensor_1_shape": "elliptical", "sensor_1_width_mm": "18", "sensor_1_height_mm": "9", "sensor_1_depth_mm": "3", "sensor_1_qe": ".72", "sensor_2_shape": "circular", "sensor_2_width_mm": "16", "sensor_2_height_mm": "16", "sensor_2_depth_mm": "2.5"})
    reports = {}
    for kind in ("cherenkov", "scintillator"):
        model = make(kind)
        sensor = gui.complete_sensor_readout_from_values(values, kind, model)
        payload = geometry.complete_detector_assembly_geometry_to_dict(model, kind, sensor)
        _assert(payload["format"] == "complete-detector-assembly-geometry" and payload["schema_version"] == 32 and payload["detector_type"] == kind, f"{kind}: complete marker/version missing")
        _assert(len(payload["active_sections"]) == 2 and len(payload["light_guide"]["sections"]) == 2, f"{kind}: ordered complete sections missing")
        _assert(set(payload["terminal_surfaces"]) == {"+X", "-X"} and payload["mylar"]["enabled"] and payload["impact_xy_mm"] == {"x": 2.25, "y": -1.75}, f"{kind}: terminal/Mylar/impact fields missing")
        _assert("position_mm" in payload["photosensor_readout"] and "detector_specific_response" in payload["photosensor_readout"], f"{kind}: complete readout fields missing")
        restored, restored_sensor = geometry.complete_detector_assembly_geometry_from_dict(json.loads(json.dumps(payload)), kind)
        _assert(restored.to_dict() == model.to_dict(), f"{kind}: canonical model round-trip changed")
        _assert(restored_sensor["dimensions_mm"] == sensor["dimensions_mm"] and restored_sensor["shape"] == sensor["shape"], f"{kind}: sensor round-trip changed")
        restored_scene = restored.complete_scene(sensor_label="Photosensor" if kind == "cherenkov" else "Scintillator readout", sensor_shape=restored_sensor["shape"], sensor_width_mm=restored_sensor["dimensions_mm"]["width"], sensor_height_mm=restored_sensor["dimensions_mm"]["height"])
        _assert(next(item for item in restored_scene["components"] if item.get("kind") == "sensor")["coupling_position_mm"][0] == 40.0, f"{kind}: canonical guide-output placement was not recomputed")
        for label, candidate in (
            ("wrong type", dict(payload, detector_type="scintillator" if kind == "cherenkov" else "cherenkov")),
            ("missing sections", {key: value for key, value in payload.items() if key != "sections"}),
            ("bad dimensions", dict(payload, sections=[dict(payload["sections"][0], width_y_start_mm=-1), payload["sections"][1]])),
            ("bad taper", dict(payload, sections=[dict(payload["sections"][0], width_y_end_mm="bad"), payload["sections"][1]])),
            ("bad terminal", dict(payload, terminal_plus_surface="not a terminal")),
            ("bad sensor", dict(payload, photosensor_readout=dict(payload["photosensor_readout"], quantum_efficiency=2))),
        ):
            before = model.to_dict()
            try:
                geometry.complete_detector_assembly_geometry_from_dict(candidate, kind)
            except (TypeError, ValueError):
                pass
            else:
                raise AssertionError(f"{kind}: invalid {label} JSON was accepted")
            _assert(model.to_dict() == before, f"{kind}: invalid {label} mutated model")
        reports[kind] = gui.detector_transport_diagnostic(values, model, kind)
    _assert(reports["cherenkov"]["physical_generation_model"]["model"] != reports["scintillator"]["physical_generation_model"]["model"], "physics diagnostics were conflated")
    _assert("Cherenkov generation/transport" in reports["cherenkov"]["physical_generation_model"]["model"], "Cherenkov generation model missing")
    _assert("deposited-energy/scintillation" in reports["scintillator"]["physical_generation_model"]["model"], "scintillator generation model missing")
    _assert(set(reports["cherenkov"]["terminal_surfaces"]) == {"+X", "-X"} and all(key in reports["cherenkov"] for key in ("light_guide", "mylar", "sensor_readout_state", "optical_transport_path", "warnings_errors")), "Cherenkov diagnostic fields incomplete")
    _assert("Scintillator Transport Diagnostic" in gui.detector_transport_diagnostic_text(reports["scintillator"]), "scintillator diagnostic text is not detector-specific")
    scintillator_before = copy.deepcopy(reports["scintillator"])
    changed = make("cherenkov")
    changed.terminal_plus_surface = "Open"
    changed_report = gui.detector_transport_diagnostic(values, changed, "cherenkov")
    scintillator_after = gui.detector_transport_diagnostic(values, make("scintillator"), "scintillator")
    _assert(changed_report["optical_transport_path"] != reports["cherenkov"]["optical_transport_path"] and scintillator_after == scintillator_before, "changing Cherenkov changed the other detector diagnostic")


def photosensor_visibility_checks() -> None:
    """Exercise both assembly types, both X faces, guides, and all render payloads."""

    import quartz_geometry as geometry
    import quartz_geometry_editor as editor
    import quartz_cherenkov_gui as gui

    def make(name, face, guide=False):
        material = "fused silica" if name == "Cherenkov" else "scintillator"
        active = geometry.AssemblySection(name="Active", material=material, axis="X")
        if face == "plus_x":
            plus, minus, guide_center = "Detector/Photosensor coupling", "Open", 60.0
        else:
            plus, minus, guide_center = "Open", "Detector/Photosensor coupling", -60.0
        guides = [geometry.AssemblySection(name="Guide", material="light guide", axis="X", length_x_mm=20.0, center_x_mm=guide_center)] if guide else []
        return geometry.AssemblyGeometryModel(name, [active], guides, guide, plus, minus)

    for name in ("Cherenkov", "Scintillator"):
        for face in ("plus_x", "minus_x"):
            for guide in (False, True):
                model = make(name, face, guide)
                scene = model.complete_scene()
                sensors = _sensor_items(scene)
                _assert(len(sensors) == 1, f"{name} {face} guide={guide}: canonical sensor missing")
                sensor = sensors[0]
                _assert(all(float(value) > 0.0 for value in sensor["dimensions_mm"].values()), "sensor dimensions are not positive")
                expected_x = 70.0 if face == "plus_x" and guide else -70.0 if face == "minus_x" and guide else 50.0 if face == "plus_x" else -50.0
                _assert(abs(sensor["coupling_position_mm"][0] - expected_x) < 1e-9, f"{name} {face} guide={guide}: wrong coupling coordinate")
                _assert(sensor["solid_faces_mm"] and sensor["window_faces_mm"], "sensor body/window faces are missing")
                for view in ("xz", "xy", "yz"):
                    items = [item for item in scene["views"][view]["components"] if item.get("kind") == "sensor"]
                    _assert(items and any(item.get("label") for item in items), f"{view}: sensor primitive/label missing")
                transport = model.to_simulation_geometry()
                _assert(len(transport.sensor_components) == 1, "transport lost canonical sensor component")
                _assert(transport.sensor_components[0]["coupling_position_mm"] == sensor["coupling_position_mm"], "transport transform diverged from scene")

    both = make("Cherenkov", "plus_x")
    both.terminal_minus_surface = "Detector/Photosensor coupling"
    both_scene = both.complete_scene()
    _assert(len(_sensor_items(both_scene)) == 2, "both detector terminals did not render two explicit sensors")
    _assert(both_scene["sensor_routing"]["ambiguous"] and both_scene["sensor_routing"]["sensor_count"] == 2, "multiple-detector status is not explicit")
    _assert(len(both.to_simulation_geometry().sensor_components) == 2, "transport silently collapsed two detector sensors")

    # Exercise the actual renderer-neutral global payload and the Tk canvas
    # helper with a display-free recording canvas.
    models = {"cherenkov": make("Cherenkov", "plus_x"), "scintillator": make("Scintillator", "plus_x")}
    global_scene = gui.global_assembly_scene(gui.default_values(), gui.simulator.BarGeometry(), models)
    global_sensors = [item for item in global_scene["render_components"] if item.get("kind") == "sensor"]
    _assert(len(global_sensors) == 2 and all(len(item["faces_mm"]) >= 6 for item in global_sensors), "global payload lost sensor bodies")
    _assert(all(item.get("readout_faces_mm") for item in global_sensors), "global payload lost sensor windows")

    class RecordingCanvas:
        def __init__(self): self.operations = []
        def winfo_width(self): return 900
        def winfo_height(self): return 420
        def create_rectangle(self, *args, **kwargs): self.operations.append(("rectangle", args, kwargs))
        def create_polygon(self, *args, **kwargs): self.operations.append(("polygon", args, kwargs))
        def create_line(self, *args, **kwargs): self.operations.append(("line", args, kwargs))
        def create_text(self, *args, **kwargs): self.operations.append(("text", args, kwargs))
        def create_oval(self, *args, **kwargs): self.operations.append(("oval", args, kwargs))
    recording = RecordingCanvas()
    gui.QuartzCherenkovGUI._draw_assembly_scene_canvas(recording, make("Cherenkov", "plus_x").complete_scene(), "cherenkov")
    labels = [str(op[2].get("text", "")) for op in recording.operations if op[0] == "text"]
    _assert(any("Photosensor" in label for label in labels), "2D renderer did not draw a sensor label")
    _assert(any(op[0] == "polygon" for op in recording.operations), "2D renderer did not draw sensor/geometry primitives")

    changed = make("Cherenkov", "plus_x")
    open_scene = changed.complete_scene()
    changed.terminal_plus_surface = "Open"
    open_again = changed.complete_scene()
    _assert(len(_sensor_items(open_scene)) == 1 and not _sensor_items(open_again), "terminal state did not change canonical scene")
    _assert(any(item.get("kind") == "sensor" for item in open_scene["views"]["xz"]["components"]), "viewer payload omitted sensor")

    # The Tk callback is tested without creating a display: changing the
    # terminal state must invoke both local and global refresh hooks.
    class Var:
        def __init__(self, value): self.value = value
        def get(self): return self.value
        def set(self, value): self.value = str(value)
    # The callback is deliberately exercised without a display or a Tk
    # interpreter. Use a plain probe object: an uninitialized tk.Tk instance
    # forwards missing attributes through Tkinter's __getattr__, which
    # recursively looks for self.tk on Windows.
    class Probe:
        pass
    probe = Probe()
    probe.assembly_models = {"cherenkov": make("Cherenkov", "plus_x")}
    probe._assembly_model = lambda assembly: probe.assembly_models[assembly]
    values = gui.default_values()
    probe.values = {key: Var(value) for key, value in values.items()}
    calls = []
    probe._assembly_draw_preview = lambda assembly: calls.append(("local", assembly))
    probe._refresh_global_assembly_preview = lambda *args: calls.append(("global",))
    probe._terminal_surface_changed = lambda assembly: gui.QuartzCherenkovGUI._terminal_surface_changed(probe, assembly)
    probe._terminal_surface_changed("cherenkov")
    probe.values["terminal_plus_surface"].set("Open")
    probe._terminal_surface_changed("cherenkov")
    _assert(calls.count(("local", "cherenkov")) == 2 and calls.count(("global",)) == 2, "terminal callback did not refresh local/global viewers")

    gui_source = (ROOT / "quartz_cherenkov_gui.py").read_text(encoding="utf-8")
    _assert("solid_faces_mm" in gui_source and "window_faces_mm" in gui_source, "3D renderer does not consume sensor body/window solids")
    _assert("_refresh_scintillator_viewer" in gui_source and "_refresh_cherenkov_viewer" in gui_source, "live viewer refresh hooks missing")
    _assert("local_viewer_scene" in gui_source and "geometry_viewer_scene" in gui_source, "live viewer scene snapshots missing")
    _assert("_canonical_local_view_data" in gui_source and "complete_local_assembly_scene" in gui_source, "interactive viewers are not based on canonical complete scenes")
    _assert("create_polygon" in gui_source and "item.get(\"label\"" in gui_source, "2D renderer does not draw labeled sensor primitives")
    _assert("canonical_sensor_components" in (ROOT / "quartz_geometry_editor.py").read_text(encoding="utf-8"), "editor lacks canonical sensor helper")


def scintillator_geometry_json_checks() -> None:
    """Exercise strict round trips and transactional failures headlessly."""

    import json
    import tempfile
    from pathlib import Path
    import quartz_geometry as geometry
    import quartz_cherenkov_gui as gui

    active = [
        geometry.AssemblySection(
            name="Scintillator taper", shape="tapered_rectangular_prism", material="scintillator",
            length_x_mm=20.0, width_y_start_mm=12.0, width_y_end_mm=16.0,
            thickness_z_start_mm=6.0, thickness_z_end_mm=8.0, axis="X", center_x_mm=-10.0,
        ),
        geometry.AssemblySection(
            name="Scintillator block", shape="box", material="scintillator",
            length_x_mm=30.0, width_y_mm=16.0, thickness_z_mm=8.0, axis="X", center_x_mm=15.0,
        ),
    ]
    guides = [
        geometry.AssemblySection(
            name="Guide taper", shape="tapered_rectangular_prism", material="light guide",
            length_x_mm=10.0, width_y_start_mm=16.0, width_y_end_mm=18.0,
            thickness_z_start_mm=8.0, thickness_z_end_mm=9.0, axis="X", center_x_mm=35.0,
        ),
        geometry.AssemblySection(
            name="Guide output", shape="box", material="light guide",
            length_x_mm=10.0, width_y_mm=18.0, thickness_z_mm=9.0, axis="X", center_x_mm=45.0,
        ),
    ]
    original = geometry.AssemblyGeometryModel(
        "Scintillator", active, guides, True, "Detector/Photosensor coupling", "Mirror",
        True, 0.91, 2.25, -1.75,
    )
    payload = geometry.assembly_geometry_to_dict(original, "scintillator")
    _assert(payload["format"] == geometry.ASSEMBLY_GEOMETRY_FORMAT, "geometry format marker missing")
    _assert(payload["schema_version"] == 27 and payload["assembly_type"] == "scintillator", "geometry version/type marker missing")
    _assert("impact_z_mm" not in payload and payload["metadata"]["impact_z"] == "derived; no independent impact Z is serialized", "independent impact Z leaked into JSON")
    restored = geometry.assembly_geometry_from_dict(json.loads(json.dumps(payload)), "scintillator")
    _assert(restored.to_dict() == original.to_dict(), "complete scintillator model did not round-trip")
    scene = restored.complete_scene(sensor_label="Scintillator readout")
    _assert(scene["impact"] == {"x_mm": 2.25, "y_mm": -1.75, "z": "derived"}, "impact metadata was not refreshed")
    _assert(scene["sensor_routing"]["readout_terminal"] == "plus_x", "guide terminal routing changed after load")
    sensor = next(item for item in scene["components"] if item.get("kind") == "sensor")
    _assert(sensor["coupling_position_mm"][0] == 50.0, "sensor was not placed at the ordered guide output")
    transport = restored.to_simulation_geometry()
    _assert(transport.sensor_components[0]["coupling_position_mm"] == sensor["coupling_position_mm"], "transport sensor placement diverged from scene")

    failures = []
    try:
        json.loads("{not json")
    except json.JSONDecodeError:
        failures.append("malformed")
    invalid_payloads = [
        ("wrong type", {**payload, "assembly_type": "cherenkov"}),
        ("missing sections", {key: value for key, value in payload.items() if key != "sections"}),
        ("bad dimension", {**payload, "sections": [dict(payload["sections"][0], width_y_start_mm=-1.0), payload["sections"][1]]}),
        ("bad taper", {**payload, "sections": [dict(payload["sections"][0], width_y_end_mm="not-a-number"), payload["sections"][1]]}),
        ("bad terminal", {**payload, "terminal_plus_surface": "Unknown terminal"}),
    ]
    for label, candidate in invalid_payloads:
        before = original.to_dict()
        try:
            geometry.assembly_geometry_from_dict(candidate, "scintillator")
        except (TypeError, ValueError):
            failures.append(label)
        _assert(original.to_dict() == before, f"invalid {label} changed the live model")
    _assert(set(failures) == {"malformed", "wrong type", "missing sections", "bad dimension", "bad taper", "bad terminal"}, "not all invalid JSON cases were rejected")

    cherenkov = geometry.AssemblyGeometryModel("Cherenkov", [geometry.AssemblySection(name="Quartz", material="fused silica", axis="X")])
    ch_payload = geometry.assembly_geometry_to_dict(cherenkov, "cherenkov")
    _assert(geometry.assembly_geometry_from_dict(ch_payload, "cherenkov").to_dict() == cherenkov.to_dict(), "Cherenkov save/load was broken")
    try:
        geometry.assembly_geometry_from_dict(ch_payload, "scintillator")
    except ValueError:
        pass
    else:
        raise AssertionError("Cherenkov JSON was accepted by the scintillator loader")

    class Var:
        def __init__(self, value): self.value = str(value)
        def get(self): return self.value
        def set(self, value): self.value = str(value)
    class Dialog:
        path = ""
        @classmethod
        def askopenfilename(cls, **_kwargs): return cls.path
    class Box:
        @staticmethod
        def showerror(*_args, **_kwargs): pass
    # Keep this as an uninitialized Tk object: the Windows failure was caused
    # by a missing instance attribute falling through Tkinter.Misc.__getattr__
    # before a Tcl interpreter exists.
    probe = gui.QuartzCherenkovGUI.__new__(gui.QuartzCherenkovGUI)
    probe_state = object.__getattribute__(probe, "__dict__")
    _assert("tk" not in probe_state and "_pending_assembly_refresh" not in probe_state, "regression probe unexpectedly has initialized Tk state")
    probe.assembly_models = {"scintillator": geometry.AssemblyGeometryModel("Scintillator", [geometry.AssemblySection(name="Current", material="scintillator", axis="X")])}
    probe.scintillator_geometry_sections = probe.assembly_models["scintillator"].sections
    probe.values = {key: Var(value) for key, value in gui.default_values().items()}
    probe.status = Var("Ready")
    probe._assembly_refresh_tree = lambda *_args: None
    probe._assembly_build_properties = lambda *_args: None
    probe._assembly_model_changed = lambda *_args: None
    probe._terminal_surface_changed = lambda *_args: None
    probe._sensor_status_changed = lambda *_args: None
    probe._refresh_global_assembly_preview = lambda *_args: None
    probe._refresh_loaded_assembly_views = lambda *_args: None
    old_dialog, old_box = gui.filedialog, gui.messagebox
    try:
        gui.filedialog, gui.messagebox = Dialog, Box
        with tempfile.TemporaryDirectory() as directory:
            path = Path(directory) / "geometry.json"
            path.write_text(json.dumps(payload), encoding="utf-8")
            Dialog.path = str(path)
            probe._load_assembly_model("scintillator")
            _assert(probe._pending_assembly_refresh == {"scintillator": False}, "load callback did not initialize pending refresh state without Tk")
            _assert(probe.assembly_models["scintillator"].to_dict() == original.to_dict(), "GUI callback did not commit a valid scintillator load")
            committed = probe.assembly_models["scintillator"].to_dict()
            committed_values = {key: value.get() for key, value in probe.values.items()}
            for invalid in ("{bad", json.dumps({**payload, "assembly_type": "cherenkov"}), json.dumps({**payload, "sections": []}), json.dumps({**payload, "terminal_plus_surface": "bad"})):
                path.write_text(invalid, encoding="utf-8")
                probe._load_assembly_model("scintillator")
                _assert(probe.assembly_models["scintillator"].to_dict() == committed, "invalid GUI load changed model")
                _assert({key: value.get() for key, value in probe.values.items()} == committed_values, "invalid GUI load changed controls")
    finally:
        gui.filedialog, gui.messagebox = old_dialog, old_box


def cherenkov_terminal_roundtrip_checks() -> None:
    """Round-trip both terminal mirrors and check the local 3D labels."""

    import copy
    import json
    import tempfile
    from pathlib import Path
    import quartz_geometry as geometry
    import quartz_cherenkov_gui as gui

    class Var:
        def __init__(self, value=""):
            self.value = str(value)
        def get(self):
            return self.value
        def set(self, value):
            self.value = str(value)

    class Canvas:
        def __init__(self):
            self.width, self.height = 900, 420
            self.idle = []
            self.operations = []
        def winfo_exists(self): return True
        def winfo_width(self): return self.width
        def winfo_height(self): return self.height
        def delete(self, *args): self.operations.append(("delete", args, {}))
        def create_rectangle(self, *args, **kwargs): self.operations.append(("rectangle", args, kwargs))
        def create_polygon(self, *args, **kwargs): self.operations.append(("polygon", args, kwargs))
        def create_line(self, *args, **kwargs): self.operations.append(("line", args, kwargs))
        def create_text(self, *args, **kwargs): self.operations.append(("text", args, kwargs))
        def create_oval(self, *args, **kwargs): self.operations.append(("oval", args, kwargs))
        def after_idle(self, callback):
            self.idle.append(callback)
            return len(self.idle)
        def after_cancel(self, _token): pass

    class Dialog:
        path = ""
        @classmethod
        def askopenfilename(cls, **_kwargs): return cls.path

    class Box:
        @staticmethod
        def showerror(*_args, **_kwargs): pass

    active = geometry.AssemblySection(
        name="Quartz radiator", material="fused silica", axis="X",
        length_x_mm=80.0, width_y_mm=24.0, thickness_z_mm=12.0,
    )
    scintillator = geometry.AssemblyGeometryModel(
        "Scintillator", [geometry.AssemblySection(name="Scintillator", material="scintillator", axis="X")],
    )
    probe = gui.QuartzCherenkovGUI.__new__(gui.QuartzCherenkovGUI)
    probe.values = {key: Var(value) for key, value in gui.default_values().items()}
    probe.assembly_models = {
        "cherenkov": geometry.AssemblyGeometryModel("Cherenkov", [copy.deepcopy(active)]),
        "scintillator": scintillator,
    }
    probe.geometry_sections = probe.assembly_models["cherenkov"].sections
    probe.scintillator_geometry_sections = scintillator.sections
    probe.geometry_viewer_canvas = Canvas()
    probe.geometry_viewer_yaw, probe.geometry_viewer_pitch, probe.geometry_viewer_zoom = .75, .45, 1.0
    probe.local_viewer_canvas = Canvas()
    probe.local_viewer_yaw, probe.local_viewer_pitch, probe.local_viewer_zoom = .75, .45, 1.0
    canvas = Canvas()
    probe._assembly_editor_state = {
        "cherenkov": {
            "preview": canvas, "preview_text": Var(), "preview_redraw_after": None,
            "selected": ("main", 0),
        },
    }
    probe.status = Var("Ready")
    probe._assembly_refresh_tree = lambda *_args: None
    probe._assembly_build_properties = lambda *_args: None
    probe._sensor_status_changed = lambda *_args: None
    global_scenes = []
    probe._refresh_global_assembly_preview = lambda: global_scenes.append(
        gui.global_assembly_scene(probe.values, gui.simulator.BarGeometry(), probe.assembly_models)
    )

    old_dialog, old_box = gui.filedialog, gui.messagebox
    try:
        gui.filedialog, gui.messagebox = Dialog, Box
        with tempfile.TemporaryDirectory() as directory:
            path = Path(directory) / "cherenkov-complete.json"
            for mirror_face in ("plus_x", "minus_x"):
                plus = "Mirror" if mirror_face == "plus_x" else "Detector/Photosensor coupling"
                minus = "Mirror" if mirror_face == "minus_x" else "Detector/Photosensor coupling"
                model = geometry.AssemblyGeometryModel(
                    "Cherenkov", [copy.deepcopy(active)], terminal_plus_surface=plus,
                    terminal_minus_surface=minus, mylar_enabled=True,
                    mylar_reflectivity=0.91, impact_x_mm=1.25, impact_y_mm=-2.5,
                )
                sensor = gui.complete_sensor_readout_from_values(probe.values, "cherenkov", model)
                payload = geometry.complete_detector_assembly_geometry_to_dict(model, "cherenkov", sensor)
                _assert(payload["schema_version"] == 32, "Cherenkov save did not use the canonical complete v32 envelope")
                _assert(payload["terminal_surfaces"] == {"+X": plus, "-X": minus}, "saved terminal block lost the configured endpoint")
                _assert(payload["opposite_end_mode"] == geometry.normalize_terminal_surface(minus), "legacy -X alias was not exported")
                path.write_text(json.dumps(payload, indent=2), encoding="utf-8")
                Dialog.path = str(path)
                probe._load_assembly_model("cherenkov")

                loaded = probe.assembly_models["cherenkov"]
                _assert(loaded.terminal_plus_surface == geometry.normalize_terminal_surface(plus), f"{mirror_face}: +X terminal changed on callback load")
                _assert(loaded.terminal_minus_surface == geometry.normalize_terminal_surface(minus), f"{mirror_face}: -X terminal changed on callback load")
                expected_alias = geometry.normalize_terminal_surface(minus)
                _assert(probe.values["opposite_end_mode"].get() == expected_alias, f"{mirror_face}: legacy Run control was not synchronized from the -X terminal")
                _assert(loaded.mylar_enabled and abs(loaded.mylar_reflectivity - 0.91) < 1e-12, f"{mirror_face}: Mylar configuration was not restored")
                _assert(probe.values["terminal_plus_surface"].get() == plus and probe.values["terminal_minus_surface"].get() == minus, f"{mirror_face}: editor controls did not restore endpoint labels")

                local_scene = probe._assembly_editor_state["cherenkov"]["scene"]
                face_key = "plus_x" if mirror_face == "plus_x" else "minus_x"
                _assert(local_scene["terminal_states"][face_key] == "mirror", f"{mirror_face}: editor preview terminal state is not Mirror")
                _assert(any(item.get("kind") == "mirror" and item.get("face") == face_key for item in local_scene["physical_components"]), f"{mirror_face}: visible editor mirror component is missing")
                _assert(local_scene["mylar"]["enabled"] and abs(local_scene["mylar"]["reflectivity"] - 0.91) < 1e-12, f"{mirror_face}: editor Mylar configuration did not refresh")
                _assert(face_key in local_scene["mylar"]["excluded_terminal_faces"], f"{mirror_face}: Mylar still covers the mirror face")

                global_scene = global_scenes[-1]
                _assert(global_scene["terminal_states"]["cherenkov"][face_key] == "mirror", f"{mirror_face}: Global Assembly state is not Mirror")
                global_components = global_scene["render_components"]
                _assert(any(item.get("kind") == "mirror" and face_key in item.get("component_id", "") for item in global_components), f"{mirror_face}: visible Global Assembly mirror is missing")

                expected_label = f"Mirror at {'+X' if mirror_face == 'plus_x' else '-X'}"
                view_data = probe._canonical_local_view_data("cherenkov")
                _assert(expected_label in [label for label, _center in view_data["mirror_labels"]], f"{mirror_face}: canonical local viewer mirror label is missing")
                probe._draw_interactive_geometry_view()
                viewer_text = [str(kwargs.get("text", "")) for name, _args, kwargs in probe.geometry_viewer_canvas.operations if name == "text"]
                _assert(expected_label in viewer_text, f"{mirror_face}: interactive Cherenkov 3D viewer did not render a readable mirror label")
                _assert(any(name == "rectangle" and kwargs.get("fill") == "#ffffff" for name, _args, kwargs in probe.geometry_viewer_canvas.operations), f"{mirror_face}: Cherenkov mirror legend lacks a high-contrast backing")

                scintillator_model = probe.assembly_models["scintillator"]
                scintillator_model.terminal_plus_surface = plus
                scintillator_model.terminal_minus_surface = minus
                probe.values["scintillator_terminal_plus_surface"].set(plus)
                probe.values["scintillator_terminal_minus_surface"].set(minus)
                probe._draw_local_3d_view()
                scintillator_text = [str(kwargs.get("text", "")) for name, _args, kwargs in probe.local_viewer_canvas.operations if name == "text"]
                _assert(expected_label in scintillator_text, f"{mirror_face}: interactive Scintillator 3D viewer did not render a readable mirror label")
                _assert(any(name == "rectangle" and kwargs.get("fill") == "#ffffff" for name, _args, kwargs in probe.local_viewer_canvas.operations), f"{mirror_face}: Scintillator mirror legend lacks a high-contrast backing")

                committed_model = copy.deepcopy(loaded.to_dict())
                committed_values = {key: value.get() for key, value in probe.values.items()}
                committed_local_scene = copy.deepcopy(local_scene)
                committed_global_scene = copy.deepcopy(global_scene)
                committed_refresh_count = len(global_scenes)
                committed_canvas_operation_count = len(canvas.operations)
                invalid = copy.deepcopy(payload)
                invalid["terminal_surfaces"]["-X"] = "Open" if minus == "Mirror" else "Mirror"
                path.write_text(json.dumps(invalid), encoding="utf-8")
                probe._load_assembly_model("cherenkov")
                _assert(loaded.to_dict() == committed_model, f"{mirror_face}: invalid/conflicting JSON changed the current geometry")
                _assert({key: value.get() for key, value in probe.values.items()} == committed_values, f"{mirror_face}: invalid/conflicting JSON changed the GUI variables")
                _assert(probe._assembly_editor_state["cherenkov"]["scene"] == committed_local_scene, f"{mirror_face}: invalid/conflicting JSON changed the editor preview")
                _assert(len(global_scenes) == committed_refresh_count and global_scenes[-1] == committed_global_scene, f"{mirror_face}: invalid/conflicting JSON refreshed or changed the global view")
                _assert(len(canvas.operations) == committed_canvas_operation_count, f"{mirror_face}: invalid/conflicting JSON redrew the editor")

            # Nested vXX surfaces restore the -X Mirror even if a stale legacy
            # alias says Open.  Legacy-only v27 inputs still migrate the alias.
            nested_only = copy.deepcopy(payload)
            nested_only.pop("terminal_plus_surface")
            nested_only.pop("terminal_minus_surface")
            nested_only["opposite_end_mode"] = "open"
            restored, _sensor = geometry.complete_detector_assembly_geometry_from_dict(nested_only, "cherenkov")
            _assert(restored.terminal_minus_surface == "mirror", "structured terminal_surfaces did not outrank a stale legacy alias")
            legacy_v27 = geometry.assembly_geometry_to_dict(model, "cherenkov")
            legacy_v27.pop("terminal_plus_surface")
            legacy_v27.pop("terminal_minus_surface")
            legacy_v27.pop("terminal_surfaces")
            legacy_v27["opposite_end_mode"] = "mirror"
            migrated = geometry.assembly_geometry_from_dict(legacy_v27, "cherenkov")
            _assert(migrated.terminal_minus_surface == "mirror", "legacy v27 opposite_end_mode did not migrate into -X")
    finally:
        gui.filedialog, gui.messagebox = old_dialog, old_box


def live_viewer_refresh_checks() -> None:
    """Prove both interactive viewer redraw paths are live without Tk/display."""

    import json
    import tempfile
    from pathlib import Path
    import quartz_geometry as geometry
    import quartz_cherenkov_gui as gui

    class Var:
        def __init__(self, value): self.value = str(value)
        def get(self): return self.value
        def set(self, value): self.value = str(value)

    class RecordingCanvas:
        def __init__(self): self.operations = []
        def winfo_exists(self): return True
        def winfo_width(self): return 900
        def winfo_height(self): return 650
        def delete(self, *args): self.operations.append(("delete", args, {}))
        def create_polygon(self, *args, **kwargs): self.operations.append(("polygon", args, kwargs))
        def create_line(self, *args, **kwargs): self.operations.append(("line", args, kwargs))
        def create_text(self, *args, **kwargs): self.operations.append(("text", args, kwargs))
        def create_oval(self, *args, **kwargs): self.operations.append(("oval", args, kwargs))
        def create_rectangle(self, *args, **kwargs): self.operations.append(("rectangle", args, kwargs))

    class LiveWindow:
        def winfo_exists(self): return True

    def make_probe():
        active = geometry.AssemblySection(name="Scintillator active", material="scintillator", axis="X", length_x_mm=20.0, width_y_mm=20.0, thickness_z_mm=10.0)
        guide = geometry.AssemblySection(name="Guide", material="light guide", axis="X", length_x_mm=20.0, width_y_mm=16.0, thickness_z_mm=8.0, center_x_mm=20.0)
        probe = gui.QuartzCherenkovGUI.__new__(gui.QuartzCherenkovGUI)
        probe_state = object.__getattribute__(probe, "__dict__")
        _assert("tk" not in probe_state and "_assembly_editor_state" not in probe_state, "live-viewer regression probe unexpectedly initialized Tk/editor state")
        probe.assembly_models = {
            "scintillator": geometry.AssemblyGeometryModel("Scintillator", [active], [guide], True, "Detector/Photosensor coupling", "Open"),
            "cherenkov": geometry.AssemblyGeometryModel("Cherenkov", [geometry.AssemblySection(name="Quartz", material="fused silica", axis="X")]),
        }
        probe.values = {key: Var(value) for key, value in gui.default_values().items()}
        probe.status = Var("Ready")
        probe.values["scintillator_guide_enabled"].set("1")
        probe.values["sensor_2_automatic_position"].set("1")
        probe.local_viewer = LiveWindow()
        probe.local_viewer_canvas = RecordingCanvas()
        probe.local_viewer_scene = None
        probe.local_viewer_yaw, probe.local_viewer_pitch, probe.local_viewer_zoom = .75, .45, 1.0
        probe.local_viewer_drag = None
        probe.geometry_viewer = LiveWindow()
        probe.geometry_viewer_canvas = RecordingCanvas()
        probe.geometry_viewer_scene = None
        probe.geometry_viewer_bindings = []
        probe.geometry_viewer_yaw, probe.geometry_viewer_pitch, probe.geometry_viewer_zoom = .75, .45, 1.0
        probe.geometry_viewer_drag = None
        return probe

    probe = make_probe()
    probe._refresh_scintillator_viewer()
    before = json.dumps(probe.local_viewer_scene, sort_keys=True, default=str)
    before_components = [(item.get("name"), item.get("kind"), item.get("bounds_mm"), item.get("coupling_position_mm")) for item in probe.local_viewer_scene["components"]]
    probe.assembly_models["scintillator"].sections[0].width_y_mm = 34.0
    probe.values["scintillator_terminal_plus_surface"].set("Mirror")
    probe.assembly_models["scintillator"].terminal_plus_surface = "Mirror"
    probe._refresh_scintillator_viewer()
    after = json.dumps(probe.local_viewer_scene, sort_keys=True, default=str)
    after_components = [(item.get("name"), item.get("kind"), item.get("bounds_mm"), item.get("coupling_position_mm")) for item in probe.local_viewer_scene["components"]]
    _assert(before != after and before_components != after_components, "already-open Scintillator viewer did not receive changed canonical geometry/terminal scene")
    _assert(any(item.get("kind") == "mirror" for item in probe.local_viewer_scene["components"]), "Scintillator mirror disappeared from live viewer scene")
    _assert(probe.local_viewer_scene["guide_enabled"] and probe.local_viewer_scene["guide_section_count"] == 1, "Scintillator active-only guide gating changed")
    probe.assembly_models["scintillator"].guide_enabled = False
    probe.values["scintillator_guide_enabled"].set("0")
    probe._refresh_scintillator_viewer()
    _assert(not probe.local_viewer_scene["guide_present"] and probe.local_viewer_scene["guide_section_count"] == 0, "disabled guide still appeared in the live viewer")
    probe.assembly_models["scintillator"].guide_enabled = True
    probe.values["scintillator_guide_enabled"].set("1")
    probe.assembly_models["scintillator"].mylar_enabled = True
    probe.assembly_models["scintillator"].impact_x_mm = 7.5
    probe.assembly_models["scintillator"].impact_y_mm = -3.25
    probe._refresh_scintillator_viewer()
    _assert(probe.local_viewer_scene["mylar"]["enabled"] and probe.local_viewer_scene["impact"] == {"x_mm": 7.5, "y_mm": -3.25, "z": "derived"}, "Mylar/impact changes did not reach the live canonical scene")

    # Sensor dimensions and placement are canonical scene inputs too.
    probe.values["sensor_2_width_mm"].set("31")
    probe.values["sensor_2_height_mm"].set("12")
    probe.assembly_models["scintillator"].terminal_plus_surface = "Detector/Photosensor coupling"
    probe._refresh_scintillator_viewer()
    sensor = next(item for item in probe.local_viewer_scene["components"] if item.get("kind") == "sensor")
    _assert(sensor["dimensions_mm"]["width"] == 31.0, "live viewer did not refresh sensor dimensions")

    # The same viewer object is redrawn in place, preserving camera values.
    probe.local_viewer_yaw, probe.local_viewer_pitch, probe.local_viewer_zoom = 1.2, .2, 1.7
    probe._refresh_scintillator_viewer()
    _assert((probe.local_viewer_yaw, probe.local_viewer_pitch, probe.local_viewer_zoom) == (1.2, .2, 1.7), "live refresh reset the camera orientation/zoom")

    # Cherenkov keeps the same live callback architecture and canonical scene.
    probe._refresh_cherenkov_viewer()
    _assert(probe.geometry_viewer_scene["assembly"] == "Cherenkov", "Cherenkov viewer stopped building its canonical scene")
    probe.assembly_models["cherenkov"].sections[0].length_x_mm = 140.0
    probe._refresh_cherenkov_viewer()
    _assert(probe.geometry_viewer_scene["components"][0]["bounds_mm"][3] - probe.geometry_viewer_scene["components"][0]["bounds_mm"][0] == 140.0, "Cherenkov live viewer did not refresh geometry")

    # Valid JSON commits and refreshes the same open object; invalid JSON is
    # rejected before the model/snapshot transaction boundary.
    payload = geometry.assembly_geometry_to_dict(geometry.AssemblyGeometryModel(
        "Scintillator", [geometry.AssemblySection(name="Loaded", material="scintillator", axis="X", length_x_mm=44.0)],
        [], False, "Open", "Open"), "scintillator")
    probe.local_viewer_scene = probe.local_viewer_scene
    probe._assembly_refresh_tree = lambda *_args: None
    probe._assembly_build_properties = lambda *_args: None
    probe._assembly_model_changed = lambda *_args: None
    probe._terminal_surface_changed = lambda *_args: None
    probe._sensor_status_changed = lambda *_args: None
    probe._refresh_global_assembly_preview = lambda *_args: None
    probe.scintillator_geometry_sections = probe.assembly_models["scintillator"].sections
    old_dialog, old_box = gui.filedialog, gui.messagebox
    try:
        class Dialog:
            path = ""
            @classmethod
            def askopenfilename(cls, **_kwargs): return cls.path
        class Box:
            @staticmethod
            def showerror(*_args, **_kwargs): pass
        gui.filedialog, gui.messagebox = Dialog, Box
        with tempfile.TemporaryDirectory() as directory:
            path = Path(directory) / "scintillator.json"
            path.write_text(json.dumps(payload), encoding="utf-8")
            Dialog.path = str(path)
            probe._load_assembly_model("scintillator")
            _assert(probe._pending_assembly_refresh == {"scintillator": False}, "live-viewer load did not initialize pending state safely")
            _assert(probe.local_viewer_scene["components"][0]["bounds_mm"][3] - probe.local_viewer_scene["components"][0]["bounds_mm"][0] == 44.0, "valid Scintillator JSON did not refresh open viewer")
            committed_model = probe.assembly_models["scintillator"].to_dict()
            committed_scene = json.dumps(probe.local_viewer_scene, sort_keys=True, default=str)
            path.write_text("{bad", encoding="utf-8")
            probe._load_assembly_model("scintillator")
            _assert(probe.assembly_models["scintillator"].to_dict() == committed_model, "invalid JSON changed the committed Scintillator model")
            _assert(json.dumps(probe.local_viewer_scene, sort_keys=True, default=str) == committed_scene, "invalid JSON changed the open viewer scene")
    finally:
        gui.filedialog, gui.messagebox = old_dialog, old_box

    probe.local_viewer = LiveWindow()
    probe.local_viewer_canvas = RecordingCanvas()
    probe.local_viewer_scene = {"sentinel": True}
    probe.local_viewer = type("DeadWindow", (), {"winfo_exists": lambda _self: False})()
    probe._refresh_scintillator_viewer()
    _assert(probe.local_viewer is None and probe.local_viewer_canvas is None and probe.local_viewer_scene is None, "destroyed Scintillator viewer references were not cleared")


def canonical_global_scene_checks() -> None:
    """Verify the v30 global scene has two live, rigidly transformed assemblies."""

    import json
    import math
    import tempfile
    import quartz_geometry as geometry
    import quartz_cherenkov_gui as gui

    def model(name, material):
        active = geometry.AssemblySection(
            name=f"{name} active", material=material, shape="tapered_rectangular_prism", axis="X",
            length_x_mm=30.0 if name == "Cherenkov" else 22.0,
            width_y_start_mm=10.0, width_y_end_mm=18.0,
            thickness_z_start_mm=6.0, thickness_z_end_mm=12.0,
            center_x_mm=0.0,
        )
        guide = geometry.AssemblySection(
            name=f"{name} guide", material="light guide", shape="box", axis="X",
            length_x_mm=12.0, width_y_mm=18.0, thickness_z_mm=12.0,
            center_x_mm=21.0 if name == "Cherenkov" else 17.0,
        )
        return geometry.AssemblyGeometryModel(
            name, [active], [guide], True, "Detector/Photosensor coupling", "Mirror",
            True, 0.93, 2.0, -1.0,
        )

    models = {"cherenkov": model("Cherenkov", "fused silica"), "scintillator": model("Scintillator", "scintillator")}
    values = gui.default_values()
    scene = gui.global_assembly_scene(values, gui.simulator.BarGeometry(), models)
    _assert(scene["assembly_count"] == 2, "global scene does not have exactly two top-level assemblies")
    _assert(scene["assembly_types"] == ["cherenkov", "scintillator"], "global assembly ordering/types changed")
    _assert([item["assembly_type"] for item in scene["assemblies"]] == ["cherenkov", "scintillator"], "global assembly list is not one-of-each")
    _assert(len([item for item in scene["assemblies"] if item["assembly_type"] == "scintillator"]) == 1, "scintillator was duplicated")
    _assert(len(scene["component_ids"]) == len(set(scene["component_ids"])), "global component IDs are not unique")
    _assert(not any("detector_" in item or "scintillator_2" in item for item in scene["component_ids"]), "stale duplicate detector path leaked into global component IDs")
    _assert(scene["status_summary"].splitlines()[:3] == ["Cherenkov assembly: 1", "Scintillator assembly: 1", "Total detector assemblies: 2"], "explicit global validation summary is missing")
    _assert(any(item["kind"] == "guide" and item["shape"] == "box" for item in scene["render_components"]), "current guide geometry was not used")
    _assert(any(item["kind"] == "sensor" and item["assembly_type"] == "scintillator" for item in scene["render_components"]), "scintillator sensor solid was not used")
    _assert(any(item["kind"] == "mirror" and item["assembly_type"] == "cherenkov" for item in scene["render_components"]), "terminal mirror was not used")

    # A local edit is visible in a newly built scene and in the same scene
    # contract used by an already-open viewer refresh callback.
    before = scene["assembly_by_type"]["scintillator"]["components"][0]["faces_mm"]
    models["scintillator"].sections[0].length_x_mm = 47.0
    models["scintillator"].sections.append(geometry.AssemblySection(name="Scintillator current", material="scintillator", axis="X", length_x_mm=9.0, center_x_mm=28.0))
    models["scintillator"].guide_sections[0].center_x_mm = 38.5
    changed = gui.global_assembly_scene(values, gui.simulator.BarGeometry(), models)
    after = changed["assembly_by_type"]["scintillator"]["components"][0]["faces_mm"]
    _assert(before != after and len(changed["assembly_by_type"]["scintillator"]["sections"]) == 2, "Scintillator global scene did not refresh from the current editor model")
    models["cherenkov"].sections[0].length_x_mm = 91.0
    models["cherenkov"].guide_sections[0].center_x_mm = 51.5
    changed_again = gui.global_assembly_scene(values, gui.simulator.BarGeometry(), models)
    _assert(changed_again["assembly_by_type"]["cherenkov"]["sections"][0]["length_x_mm"] == 91.0, "Cherenkov global scene did not refresh from the current editor model")

    # One rigid transform is applied to every component of each assembly.
    values["global_offset_x_mm"] = "13"
    values["global_offset_y_mm"] = "-4"
    values["global_offset_z_mm"] = "7"
    values["assembly_1_rotation_axis"] = "Z"
    values["assembly_1_rotation_angle_deg"] = "37"
    transformed = gui.global_assembly_scene(values, gui.simulator.BarGeometry(), models)
    local = gui.solid_surface_faces_mm(models["cherenkov"].sections[0].to_simulation_geometry())
    global_faces = transformed["assembly_by_type"]["cherenkov"]["components"][0]["faces_mm"]
    local_distance = math.dist(local[0][0], local[0][1])
    global_distance = math.dist(global_faces[0][0], global_faces[0][1])
    _assert(abs(local_distance - global_distance) < 1e-9, "rigid transform changed internal component distances")
    _assert(transformed["transformed_impact_markers_mm"][0] != [2.0, -1.0, 0.0], "impact marker was not transformed")

    # Valid assembly JSON commits and causes a global rebuild; invalid JSON
    # does neither, preserving the previous canonical scene.
    class Var:
        def __init__(self, value): self.value = str(value)
        def get(self): return self.value
        def set(self, value): self.value = str(value)
    probe = gui.QuartzCherenkovGUI.__new__(gui.QuartzCherenkovGUI)
    probe.assembly_models = models.copy()
    probe.values = {key: Var(value) for key, value in gui.default_values().items()}
    probe.status = Var("Ready")
    probe._assembly_refresh_tree = lambda *_args: None
    probe._assembly_build_properties = lambda *_args: None
    probe._assembly_model_changed = lambda *_args: None
    probe._terminal_surface_changed = lambda *_args: None
    probe._sensor_status_changed = lambda *_args: None
    refreshes = []
    probe._refresh_live_assembly_viewer = lambda *_args: None
    probe._refresh_global_assembly_preview = lambda *_args: refreshes.append(gui.global_assembly_scene(probe.values, gui.simulator.BarGeometry(), probe.assembly_models))
    old_dialog, old_box = gui.filedialog, gui.messagebox
    try:
        class Dialog:
            path = ""
            @classmethod
            def askopenfilename(cls, **_kwargs): return cls.path
        class Box:
            @staticmethod
            def showerror(*_args, **_kwargs): pass
        gui.filedialog, gui.messagebox = Dialog, Box
        with tempfile.TemporaryDirectory() as directory:
            path = Path(directory) / "scintillator.json"
            payload = geometry.assembly_geometry_to_dict(models["scintillator"], "scintillator")
            path.write_text(json.dumps(payload), encoding="utf-8")
            Dialog.path = str(path)
            probe._load_assembly_model("scintillator")
            _assert(refreshes and refreshes[-1]["assembly_count"] == 2, "valid JSON load did not refresh global scene")
            committed = probe.assembly_models["scintillator"].to_dict()
            count = len(refreshes)
            path.write_text("{bad", encoding="utf-8")
            probe._load_assembly_model("scintillator")
            _assert(probe.assembly_models["scintillator"].to_dict() == committed and len(refreshes) == count, "invalid JSON changed/refreshed global scene")
    finally:
        gui.filedialog, gui.messagebox = old_dialog, old_box


def rigid_assembly_transform_checks() -> None:
    """Exercise the model-layer rigid transform for all required configurations."""

    import math
    import numpy as np
    import quartz_geometry as geometry
    import quartz_cherenkov_gui as gui

    def make_model(name, material, detector_face, guide):
        active = geometry.AssemblySection(
            name=f"{name} active", material=material, shape="tapered_rectangular_prism", axis="X",
            length_x_mm=30.0, width_y_start_mm=10.0, width_y_end_mm=18.0,
            thickness_z_start_mm=6.0, thickness_z_end_mm=12.0,
        )
        guide_sections = []
        if guide:
            guide_sections.append(geometry.AssemblySection(
                name=f"{name} guide", material="light guide", shape="box", axis="X",
                length_x_mm=12.0, width_y_mm=18.0, thickness_z_mm=12.0,
                center_x_mm=21.0 if detector_face == "plus_x" else -21.0,
            ))
        plus = "Detector/Photosensor coupling" if detector_face == "plus_x" else "Mirror"
        minus = "Detector/Photosensor coupling" if detector_face == "minus_x" else "Mirror"
        return geometry.AssemblyGeometryModel(
            name, [active], guide_sections, guide,
            plus, minus, True, 0.91, 2.25, -1.75,
        )

    def values_for_transform():
        values = gui.default_values()
        values.update({
            "global_offset_x_mm": "13.5", "global_offset_y_mm": "-4.25", "global_offset_z_mm": "7.75",
            "detector2_offset_x_mm": "-8.0", "detector2_offset_y_mm": "5.5", "detector2_offset_z_mm": "-2.25",
            "assembly_1_rotation_axis": "Z", "assembly_1_rotation_angle_deg": "37",
            "assembly_2_rotation_axis": "Y", "assembly_2_rotation_angle_deg": "-23",
        })
        return values

    def points(value):
        return np.asarray(value, dtype=float)

    def compare_face_lists(before, after, transform, label):
        _assert(len(before) == len(after), f"{label}: face count changed")
        for face_before, face_after in zip(before, after):
            expected = [transform.apply_point_mm(point) for point in face_before]
            _assert(np.allclose(points(face_after), points(expected), atol=1e-9, rtol=0.0), f"{label}: face was not rigidly transformed")

    def relative_transform(before_transform, after_transform):
        """Return T_after * inverse(T_before) for the canonical baseline."""

        # The baseline used by this check has no rotation.  Keeping the
        # calculation explicit avoids conflating the new physical placement
        # translation with the user-supplied rigid transform.
        before_transform.validate()
        after_transform.validate()
        if before_transform.normalized_axis != "None":
            raise AssertionError("rigid transform baseline unexpectedly rotated")
        zero_after = geometry.RigidAssemblyTransform(
            rotation_axis=after_transform.rotation_axis,
            rotation_angle_deg=after_transform.rotation_angle_deg,
            rotation_origin_x_mm=after_transform.rotation_origin_x_mm,
            rotation_origin_y_mm=after_transform.rotation_origin_y_mm,
            rotation_origin_z_mm=after_transform.rotation_origin_z_mm,
        )
        rotated_before_translation = np.asarray(zero_after.apply_point_mm((
            before_transform.translation_x_mm,
            before_transform.translation_y_mm,
            before_transform.translation_z_mm,
        )))
        return geometry.RigidAssemblyTransform(
            translation_x_mm=float(after_transform.translation_x_mm - rotated_before_translation[0]),
            translation_y_mm=float(after_transform.translation_y_mm - rotated_before_translation[1]),
            translation_z_mm=float(after_transform.translation_z_mm - rotated_before_translation[2]),
            rotation_axis=after_transform.rotation_axis,
            rotation_angle_deg=after_transform.rotation_angle_deg,
            rotation_origin_x_mm=after_transform.rotation_origin_x_mm,
            rotation_origin_y_mm=after_transform.rotation_origin_y_mm,
            rotation_origin_z_mm=after_transform.rotation_origin_z_mm,
        )

    def component_anchor(component):
        all_points = [point for face in component.get("faces_mm", ()) for point in face]
        _assert(all_points, f"{component.get('component_id')}: missing physical faces")
        return np.mean(points(all_points), axis=0)

    for name, material in (("Cherenkov", "fused silica"), ("Scintillator", "scintillator")):
        for detector_face in ("plus_x", "minus_x"):
            for guide in (False, True):
                model = make_model(name, material, detector_face, guide)
                models = {"cherenkov": make_model("Cherenkov", "fused silica", detector_face, guide),
                          "scintillator": make_model("Scintillator", "scintillator", detector_face, guide)}
                values = gui.default_values()
                # Drive the assembly under test through its actual canonical
                # key; the other assembly remains at the identity transform.
                if name == "Cherenkov":
                    values.update(values_for_transform())
                    assembly_key = "cherenkov"
                else:
                    transformed_values = values_for_transform()
                    for key in ("global_offset_x_mm", "global_offset_y_mm", "global_offset_z_mm"):
                        transformed_values[key] = values_for_transform()[key]
                    values.update(transformed_values)
                    assembly_key = "scintillator"
                # Test the full global scene with one-of-each assemblies.
                identity = gui.global_assembly_scene(gui.default_values(), gui.simulator.BarGeometry(), models)
                transformed = gui.global_assembly_scene(values, gui.simulator.BarGeometry(), models)
                before_assembly = identity["assembly_by_type"][assembly_key]
                after_assembly = transformed["assembly_by_type"][assembly_key]
                transform = relative_transform(before_assembly["rigid_transform"], after_assembly["rigid_transform"])
                before = {item["component_id"]: item for item in before_assembly["components"]}
                after = {item["component_id"]: item for item in after_assembly["components"]}
                _assert(set(before) == set(after), f"{name} face={detector_face} guide={guide}: component IDs changed")
                _assert(sum(item.get("kind") == "sensor" for item in after.values()) == 1, f"{name} face={detector_face} guide={guide}: sensor count")
                anchors_before, anchors_after = [], []
                for component_id, before_item in before.items():
                    after_item = after[component_id]
                    compare_face_lists(before_item["faces_mm"], after_item["faces_mm"], transform, component_id)
                    for key in ("interface_faces_mm", "readout_faces_mm", "solid_faces_mm", "window_faces_mm"):
                        if key in before_item:
                            compare_face_lists(before_item.get(key, ()), after_item.get(key, ()), transform, f"{component_id}/{key}")
                    if "position_mm" in before_item:
                        _assert(np.allclose(points(after_item["position_mm"]), points(transform.apply_point_mm(before_item["position_mm"])), atol=1e-9, rtol=0.0), f"{component_id}: body position not transformed")
                    if "coupling_position_mm" in before_item:
                        _assert(np.allclose(points(after_item["coupling_position_mm"]), points(transform.apply_point_mm(before_item["coupling_position_mm"])), atol=1e-9, rtol=0.0), f"{component_id}: coupling position not transformed")
                    if "normal_vector_mm" in before_item:
                        _assert(np.allclose(points(after_item["normal_vector_mm"]), points(transform.apply_vector_mm(before_item["normal_vector_mm"])), atol=1e-9, rtol=0.0), f"{component_id}: sensor normal not rotated")
                    anchors_before.append(component_anchor(before_item))
                    anchors_after.append(component_anchor(after_item))
                # Every pairwise vector and distance is invariant under one
                # rigid transform, including active, guide, mirror and sensor.
                for index in range(len(anchors_before)):
                    for other in range(index + 1, len(anchors_before)):
                        local_vector = anchors_before[other] - anchors_before[index]
                        global_vector = anchors_after[other] - anchors_after[index]
                        _assert(np.allclose(global_vector, transform.apply_vector_mm(local_vector), atol=1e-9, rtol=0.0), f"{name}: pairwise relative vector changed")
                        _assert(abs(np.linalg.norm(local_vector) - np.linalg.norm(global_vector)) < 1e-9, f"{name}: pairwise relative distance changed")
                absolute_transform = after_assembly["rigid_transform"]
                _assert(np.allclose(points(after_assembly["impact_marker_global_mm"]), points(absolute_transform.apply_point_mm(after_assembly["impact_marker_local_mm"])), atol=1e-9, rtol=0.0), f"{name}: impact marker was not transformed")
                # Sensor window must meet the transformed terminal/output and
                # its body centre must remain outside the assembly along the
                # transformed outward normal.
                sensor = next(item for item in after.values() if item.get("kind") == "sensor")
                terminal = next(item for item in after_assembly["terminal_interfaces"] if item["face"] == detector_face)
                _assert(np.allclose(points(sensor["coupling_position_mm"]), points(terminal["global_position_mm"]), atol=1e-9, rtol=0.0), f"{name} face={detector_face} guide={guide}: sensor window is detached")
                published_sensor = transformed["canonical_sensor_transforms"][assembly_key][0]
                _assert(np.allclose(points(published_sensor["position_mm"]), points(sensor["position_mm"]), atol=1e-9, rtol=0.0), f"{name}: published global sensor transform is stale/local")
                outward = points(sensor["position_mm"]) - points(sensor["coupling_position_mm"])
                _assert(float(np.dot(outward, points(sensor["normal_vector_mm"]))) > 0.0, f"{name} face={detector_face}: sensor body is not outward")
                _assert(set(after_assembly["mylar"].keys()) >= {"enabled", "reflectivity"}, f"{name}: Mylar metadata missing")

    # The two top-level assemblies remain exactly one of each, with one sensor
    # per readout even when both are transformed.
    both = gui.global_assembly_scene(values_for_transform(), gui.simulator.BarGeometry(), {
        "cherenkov": make_model("Cherenkov", "fused silica", "plus_x", True),
        "scintillator": make_model("Scintillator", "scintillator", "minus_x", True),
    })
    _assert(both["assembly_types"] == ["cherenkov", "scintillator"] and both["assembly_count"] == 2, "global scene assembly cardinality changed")
    _assert(sum(item.get("kind") == "sensor" for item in both["render_components"]) == 2, "global scene duplicated or omitted sensors")
    _assert(len(both["component_ids"]) == len(set(both["component_ids"])), "global component IDs are not unique after transform")

    # Fake an already-open viewer: an offset/rotation edit must rebuild the
    # same canonical scene that the viewer would draw.
    class FakeVar:
        def __init__(self, value): self.value = str(value)
        def get(self): return self.value
        def set(self, value): self.value = str(value)

    probe = gui.QuartzCherenkovGUI.__new__(gui.QuartzCherenkovGUI)
    probe.values = {key: FakeVar(value) for key, value in gui.default_values().items()}
    probe.assembly_models = {"cherenkov": make_model("Cherenkov", "fused silica", "plus_x", True),
                             "scintillator": make_model("Scintillator", "scintillator", "plus_x", False)}
    probe.global_preview_canvas = None
    probe.global_viewer = object()
    probe.global_viewer_canvas = object()
    captured = []
    probe._draw_global_assembly_3d_view = lambda: captured.append(gui.global_assembly_scene(probe.values, gui.simulator.BarGeometry(), probe.assembly_models))
    probe._refresh_global_assembly_preview()
    initial = captured[-1]
    probe.values["global_offset_x_mm"].set("44")
    probe.values["assembly_1_rotation_axis"].set("X")
    probe.values["assembly_1_rotation_angle_deg"].set("31")
    probe._refresh_global_assembly_preview()
    refreshed = captured[-1]
    _assert(initial["assembly_by_type"]["cherenkov"]["transform"] != refreshed["assembly_by_type"]["cherenkov"]["transform"], "fake open global viewer did not refresh transform")
    initial_sensor = next(item for item in initial["assembly_by_type"]["cherenkov"]["components"] if item.get("kind") == "sensor")
    refreshed_sensor = next(item for item in refreshed["assembly_by_type"]["cherenkov"]["components"] if item.get("kind") == "sensor")
    _assert(initial_sensor["position_mm"] != refreshed_sensor["position_mm"], "fake open global viewer left sensor stationary")


def mylar_overlay_checks() -> None:
    """Verify explicit eligible-face overlays and display-only renderer flow."""

    import numpy as np
    import quartz_geometry as geometry
    import quartz_cherenkov_gui as gui

    def make(name, guide, plus, minus, enabled=True):
        material = "fused silica" if name == "Cherenkov" else "scintillator"
        active = geometry.AssemblySection(name="Active", material=material, axis="X", length_x_mm=40.0)
        guides = [geometry.AssemblySection(name="Guide", material="light guide", axis="X", length_x_mm=20.0, center_x_mm=30.0)] if guide else []
        return geometry.AssemblyGeometryModel(name, [active], guides, guide, plus, minus, enabled, 0.91)

    for guide in (False, True):
        for plus in geometry.TERMINAL_SURFACE_CHOICES:
            for minus in geometry.TERMINAL_SURFACE_CHOICES:
                model = make("Scintillator", guide, plus, minus)
                scene = model.complete_scene()
                overlays = scene["overlay_components"]
                _assert(all(item["kind"] == "mylar_overlay" for item in overlays), "overlay kind is not explicit")
                _assert(all(item["display_only"] and item["display_offset_mm"] > 0 for item in overlays), "display-only epsilon metadata missing")
                _assert(all(item["color"] == geometry.MYLAR_DISPLAY_COLOR and item["hatch"] == geometry.MYLAR_DISPLAY_HATCH for item in overlays), "Mylar style metadata missing")
                for view in ("xz", "xy", "yz"):
                    view_overlays = [item for item in scene["views"][view]["components"] if item.get("kind") == "mylar_overlay"]
                    _assert(all("Mylar" in item.get("label", "") and item.get("hatch") for item in view_overlays), f"{view}: Mylar hatch/label missing")
                active_indices = {int(item["source_face_index"]) for item in overlays if item["section_name"] == "Active"}
                guide_indices = {int(item["source_face_index"]) for item in overlays if item["section_name"] == "Guide"}
                active_faces = geometry.canonical_section_faces_mm(model.sections[0])
                active_low = next(index for index, face in enumerate(active_faces) if max(point[0] for point in face) - min(point[0] for point in face) < 1e-9 and sum(point[0] for point in face) / len(face) < 0)
                active_high = next(index for index, face in enumerate(active_faces) if max(point[0] for point in face) - min(point[0] for point in face) < 1e-9 and sum(point[0] for point in face) / len(face) > 0)
                if guide:
                    guide_faces = geometry.canonical_section_faces_mm(model.guide_sections[0])
                    guide_low = next(index for index, face in enumerate(guide_faces) if max(point[0] for point in face) - min(point[0] for point in face) < 1e-9 and sum(point[0] for point in face) / len(face) < 30)
                    guide_high = next(index for index, face in enumerate(guide_faces) if max(point[0] for point in face) - min(point[0] for point in face) < 1e-9 and sum(point[0] for point in face) / len(face) > 30)
                    _assert(active_high not in active_indices and guide_low not in guide_indices, "active/guide coupling cap was wrapped")
                else:
                    guide_high = None
                if plus != "Open":
                    _assert((guide_high if guide else active_high) not in (guide_indices if guide else active_indices), "+X optical/mirror terminal was wrapped")
                if minus != "Open":
                    _assert(active_low not in active_indices, "-X optical/mirror terminal was wrapped")
                disabled = make("Scintillator", guide, plus, minus, False).complete_scene()
                _assert(not disabled["overlay_components"], "disabled Mylar still has overlay components")
                _assert(not any(item.get("kind") == "mylar_overlay" for view in disabled["views"].values() for item in view.get("components", ())), "disabled Mylar still has 2D overlays")

    enabled_model = make("Scintillator", True, "Detector/Photosensor coupling", "Open", True)
    disabled_model = make("Scintillator", True, "Detector/Photosensor coupling", "Open", False)
    _assert(enabled_model.to_simulation_geometry().sections == disabled_model.to_simulation_geometry().sections, "Mylar overlay changed transport sections")
    _assert(enabled_model.to_simulation_geometry().terminal_interfaces == disabled_model.to_simulation_geometry().terminal_interfaces, "Mylar overlay changed terminal transport metadata")

    values = gui.default_values()
    values.update({"global_offset_x_mm": "17", "global_offset_y_mm": "-4", "global_offset_z_mm": "9", "assembly_1_rotation_axis": "Z", "assembly_1_rotation_angle_deg": "28", "detector2_offset_x_mm": "-12", "detector2_offset_y_mm": "6", "detector2_offset_z_mm": "3", "assembly_2_rotation_axis": "Y", "assembly_2_rotation_angle_deg": "-19"})
    models = {"cherenkov": make("Cherenkov", True, "Detector/Photosensor coupling", "Mirror"), "scintillator": enabled_model}
    global_scene = gui.global_assembly_scene(values, gui.simulator.BarGeometry(), models)
    for assembly in global_scene["assemblies"]:
        transform = assembly["rigid_transform"]
        local_overlays = assembly["local_scene"]["overlay_components"]
        global_overlays = [item for item in assembly["components"] if item.get("kind") == "mylar_overlay"]
        _assert(len(local_overlays) == len(global_overlays), "global overlay count changed")
        for local_item, global_item in zip(local_overlays, global_overlays):
            expected = [[transform.apply_point_mm(point) for point in face] for face in local_item["faces_mm"]]
            _assert(np.allclose(np.asarray(expected), np.asarray(global_item["faces_mm"]), atol=1e-9, rtol=0), "global Mylar overlay did not follow rigid transform")

    class RecordingCanvas:
        def __init__(self): self.operations = []
        def winfo_width(self): return 900
        def winfo_height(self): return 420
        def create_rectangle(self, *args, **kwargs): self.operations.append(("rectangle", args, kwargs))
        def create_polygon(self, *args, **kwargs): self.operations.append(("polygon", args, kwargs))
        def create_line(self, *args, **kwargs): self.operations.append(("line", args, kwargs))
        def create_text(self, *args, **kwargs): self.operations.append(("text", args, kwargs))
        def create_oval(self, *args, **kwargs): self.operations.append(("oval", args, kwargs))
    canvas = RecordingCanvas()
    gui.QuartzCherenkovGUI._draw_assembly_scene_canvas(canvas, models["cherenkov"].complete_scene(), "cherenkov")
    _assert(any(op[2].get("stipple") == "gray25" for op in canvas.operations if op[0] == "polygon"), "2D renderer did not draw Mylar stipple")
    _assert(any("Mylar" in str(op[2].get("text", "")) for op in canvas.operations if op[0] == "text"), "2D renderer did not draw Mylar label")
    _assert(any(op[2].get("dash") for op in canvas.operations if op[0] == "line"), "2D renderer did not draw Mylar stripe edge")

    gui_source = (ROOT / "quartz_cherenkov_gui.py").read_text(encoding="utf-8")
    geometry_source = (ROOT / "quartz_geometry.py").read_text(encoding="utf-8")
    _assert("derive_mylar_overlay_surfaces" in geometry_source and "overlay_components" in geometry_source, "canonical Mylar overlay derivation missing")
    _assert("MYLAR_DISPLAY_HATCH" in gui_source and "mylar_overlay" in gui_source, "renderer Mylar style route missing")


def event_viewer_scene_checks() -> None:
    """Verify detector-separated global event scenes and camera-fit inputs."""

    import quartz_cherenkov_gui as gui
    import quartz_cherenkov_sim as simulator

    values = gui.default_values()
    values.update({
        "detector_run_mode": "combined", "scintillator_enabled": "1",
        "scintillator_sampling_mode": "fixed", "scintillator_photon_count": "3",
        "scintillator_deposited_energy_MeV": "1", "photons_per_event": "3",
        "display_photons": "3", "cherenkov_guide_enabled": "1",
        "scintillator_guide_enabled": "1", "global_offset_x_mm": "0",
        "detector2_offset_x_mm": "0", "detector2_offset_z_mm": "96", "assembly_1_rotation_axis": "Z",
        "assembly_1_rotation_angle_deg": "17", "assembly_2_rotation_axis": "Y",
        "assembly_2_rotation_angle_deg": "-11",
    })
    settings = gui.validate_settings(values)
    geometry = simulator.BarGeometry()
    canonical = gui.global_assembly_scene(values, geometry)
    _assert(canonical["assembly_types"] == ["cherenkov", "scintillator"], "canonical scene does not contain exactly one assembly of each type")
    _assert(canonical["assembly_by_type"]["cherenkov"]["transform"]["rotation_angle_deg"] == 17.0, "Cherenkov global rotation was not retained")
    _assert(canonical["assembly_by_type"]["scintillator"]["transform"]["translation_z_mm"] == 96.0, "Scintillator global offset was not retained")
    payload = gui.canonical_event_view_payload(canonical, settings.detector_run_mode)
    _assert(payload["coordinate_system"] == "global_mm" and payload["visible_assembly_types"] == ["cherenkov", "scintillator"], "combined event scene payload is not typed/global")

    particle = gui._particle_for_settings(settings)
    optical = gui._optical_for_settings(settings)
    local = gui.scintillator_local_geometry_from_values(values)
    global_geometry = gui.scintillator_global_geometry_from_values(values)
    assembly = simulator.compose_scintillator_assembly(local, global_geometry, geometry)
    scint_config = gui.scintillator_config_from_values(values)
    results = {}
    for mode, scintillator in (("cherenkov", None), ("scintillator", True), ("combined", True)):
        result = simulator.simulate(
            particle, geometry, optical, seed=123, scintillator_assembly=assembly if scintillator else None,
            scintillator_config=scint_config if scintillator else None, detector_mode=mode,
        )
        results[mode] = result
        event_scene = gui.build_interactive_event_scene(result, payload, detector_mode=mode, display_photons=3)
        expected = [mode] if mode != "combined" else ["cherenkov", "scintillator"]
        _assert(event_scene["visible_assembly_types"] == expected, f"{mode}: wrong visible assembly selection")
        _assert(set(event_scene["per_detector"]) == set(expected), f"{mode}: detector records are not separated")
        _assert(len(event_scene["track_mm"]) >= 2, f"{mode}: shared particle trajectory missing")
        bounds = event_scene["bounds_mm"]
        for assembly_type in expected:
            for component in event_scene["geometry_components"]:
                if component["assembly_type"] != assembly_type:
                    continue
                for face in component["faces_mm"]:
                    for point in face:
                        _assert(all(bounds[i] <= point[i] <= bounds[i + 3] for i in range(3)), f"{mode}: camera fit clips {assembly_type}")
    combined_scene = gui.build_interactive_event_scene(results["combined"], payload, detector_mode="combined", display_photons=3)
    _assert(combined_scene["per_detector"]["cherenkov"]["record_count"] == len(results["combined"].photon_records), "Cherenkov ownership changed")
    _assert(combined_scene["per_detector"]["scintillator"]["record_count"] == len(results["combined"].scintillator.photon_records), "Scintillator ownership changed")
    _assert(combined_scene["bounds_mm"][3] >= canonical["assembly_high_x_mm"], "combined camera fit clips downstream assembly")
    _assert(all("assembly_type" in line and "color" in line for line in combined_scene["photon_lines"]), "detector-specific photon colors/ownership missing")


def static_selected_event_renderer_checks() -> None:
    """Regression for the v48 static image dropping the S path family."""

    import tempfile
    import numpy as np
    import quartz_cherenkov_gui as gui
    import quartz_cherenkov_sim as simulator

    values = gui.default_values()
    values.update({
        # Reproduce the supplied static-image contract: All mode, 40 shown
        # per detector, and deterministic 100/100 source records.
        "detector_run_mode": "All", "photons_per_event": "100", "scintillator_photon_count": "100",
        "scintillator_deposited_energy_MeV": "1.0", "display_photons": "40",
        "scintillator_guide_enabled": "1",
        # Straight global geometry with the requested canonical Z separation.
        "global_offset_x_mm": "0", "global_offset_y_mm": "0", "global_offset_z_mm": "0",
        "detector2_offset_x_mm": "0", "detector2_offset_y_mm": "0", "detector2_offset_z_mm": "30",
        "assembly_1_rotation_axis": "Z", "assembly_1_rotation_angle_deg": "0",
        "assembly_2_rotation_axis": "Y", "assembly_2_rotation_angle_deg": "0",
    })
    geometry = simulator.BarGeometry()
    local = gui.scintillator_local_geometry_from_values(values)
    assembly = local.to_assembly(mylar_enabled=True, mylar_reflectivity=0.95)
    config = gui.scintillator_config_from_values(values)
    particle = gui._particle_for_settings(gui.validate_settings(values))
    canonical = gui.global_assembly_scene(values, geometry)
    output_root = Path(tempfile.mkdtemp(prefix="v49_static_event_"))

    for mode in ("cherenkov", "scintillator", "combined"):
        result = simulator.simulate(
            particle, geometry=geometry, optical=simulator.OpticalConfig(photon_count=100),
            seed=404, detector_mode=mode,
            scintillator_assembly=assembly if mode != "cherenkov" else None,
            scintillator_config=config if mode != "cherenkov" else None,
        )
        payload = gui.canonical_event_view_payload(canonical, mode)
        scene = gui.build_interactive_event_scene(
            result, payload, detector_mode=mode, display_photons=40,
            selected_categories={"photon_paths"}, track_margin_mm=20.0,
        )
        expected = ["cherenkov"] if mode == "cherenkov" else ["scintillator"] if mode == "scintillator" else ["cherenkov", "scintillator"]
        _assert(scene["visible_assembly_types"] == expected, f"{mode}: static scene visible assemblies changed")
        component_types = [str(item["assembly_type"]) for item in scene["geometry_components"]]
        _assert(set(component_types) == set(expected), f"{mode}: static scene geometry ownership leaked")
        _assert(all(item in expected for item in component_types), f"{mode}: non-participating geometry present")
        _assert(scene["legend"][0] == {
            "cherenkov": "MODE: CHERENKOV ONLY",
            "scintillator": "MODE: SCINTILLATOR ONLY",
            "combined": "MODE: ALL",
        }[mode], f"{mode}: mode-specific legend is wrong")
        _assert(set(scene["per_detector"]) == set(expected), f"{mode}: detector record groups are not isolated")
        _assert(all(line["detector_type"] in expected for line in scene["photon_lines"]), f"{mode}: photon ownership leaked")
        _assert(scene["schema"].endswith("v50"), f"{mode}: static scene schema is stale")
        if mode == "combined":
            _assert(component_types.count("cherenkov") > 0 and component_types.count("scintillator") > 0, "combined: one assembly was dropped")
            _assert(any(item["kind"] == "guide" and item["assembly_type"] == "scintillator" for item in scene["geometry_components"]), "combined: scintillator guide geometry was dropped")
            _assert(any(item["kind"] == "sensor" and item["assembly_type"] == "scintillator" for item in scene["geometry_components"]), "combined: scintillator readout was dropped")
            _assert(any(line["detector_type"] == "cherenkov" for line in scene["photon_lines"]), "combined: Cherenkov paths were dropped")
            _assert(any(line["detector_type"] == "scintillator" for line in scene["photon_lines"]), "combined: scintillator paths were dropped")
            _assert(scene["per_detector"]["cherenkov"]["record_count"] == 100, "combined: Cherenkov source record count is not 100")
            _assert(scene["per_detector"]["scintillator"]["record_count"] == 100, "combined: Scintillator source record count is not 100")
            _assert(scene["static_rendered_count"]["cherenkov"] == 40, "combined: static renderer did not show 40 Cherenkov paths")
            _assert(scene["static_rendered_count"]["scintillator"] == 40, "combined: static renderer did not show 40 Scintillator paths")
            _assert(sum(scene["static_rendered_count"].values()) == 80, "combined: static renderer still applied one combined cap")
            _assert(np.isclose(float(canonical["assembly_by_type"]["scintillator"]["transform"]["translation_z_mm"]), 30.0), "combined: Scintillator Z transform is not 30 mm")
            # The records remain owned and counted even when a selected
            # category excludes one detector family.  This prevents the
            # renderer from misattributing an empty S selection to missing
            # transport records.
            escaped_only = gui.build_interactive_event_scene(
                result, payload, detector_mode=mode, display_photons=40,
                selected_categories={"escaped"}, track_margin_mm=20.0,
            )
            _assert(escaped_only["per_detector"]["scintillator"]["record_count"] == 100, "category filter changed Scintillator source count")
            _assert(escaped_only["per_detector"]["scintillator"]["selected_record_count"] == 0, "escaped-only filter selected Scintillator records")
            _assert(escaped_only["static_rendered_count"]["scintillator"] == 0, "escaped-only filter rendered unselected Scintillator paths")
            _assert(escaped_only["static_rendered_count"]["cherenkov"] > 0, "escaped-only filter lost selected Cherenkov paths")
        bounds = np.asarray(scene["bounds_mm"], dtype=float)
        all_points = [point for line in scene["photon_lines"] for point in line["points_mm"]] + list(scene["track_mm"])
        _assert(all(np.all(np.asarray(point) >= bounds[:3] - 1e-9) and np.all(np.asarray(point) <= bounds[3:] + 1e-9) for point in all_points), f"{mode}: camera bounds exclude a path")
        target = output_root / mode
        target.mkdir()
        output = gui.create_selected_event_plot(
            [result], 0, target, display_photons=40, track_margin_mm=20.0,
            selected_categories={"photon_paths"}, canonical_scene=payload,
        )
        _assert(output.exists() and output.stat().st_size > 1000, f"{mode}: canonical static PNG was not written")

    # A typed local path receives the same transform as the canonical
    # component scene, exactly once.  This is the static/interactive parity
    # invariant for the separated All-mode case.
    transformed_payload = gui.canonical_event_view_payload(canonical, "combined")
    ch_result = simulator.simulate(
        particle, geometry=geometry, optical=simulator.OpticalConfig(photon_count=1),
        seed=405, detector_mode="cherenkov",
    )
    record = ch_result.photon_records[0]
    global_points = gui._event_path_points_global_mm(record, "cherenkov", transformed_payload)
    transform = gui.RigidAssemblyTransform(**transformed_payload["path_transforms"]["cherenkov"])
    expected_point = transform.apply_point_mm(np.asarray(record.path_points_m[0]) * 1e3)
    _assert(np.allclose(global_points[0], expected_point, atol=1e-9), "Cherenkov path was not transformed once into canonical global coordinates")


def event_global_frame_alignment_checks() -> None:
    """Regression for v48's canonical global/event frame and shared beam."""

    import dataclasses
    import math
    import numpy as np
    import quartz_cherenkov_gui as gui
    import quartz_cherenkov_sim as simulator

    screenshot_values = gui.default_values()
    screenshot_values.update({
        "detector_run_mode": "All", "scintillator_enabled": "1",
        # This is the displayed Quartz -> D1 gap diagnostic.  It is not a
        # signed assembly offset and must never be folded into X/Y placement.
        "scintillator_gap_mm": "-80",
        "global_offset_x_mm": "0", "global_offset_y_mm": "0", "global_offset_z_mm": "0",
        "detector2_offset_x_mm": "0", "detector2_offset_y_mm": "0", "detector2_offset_z_mm": "30",
        "scintillator_impact_x_mm": "0", "scintillator_impact_y_mm": "0",
        "entry_x_mm": "0", "entry_y_mm": "0", "display_photons": "4",
    })
    geometry = simulator.BarGeometry()
    canonical = gui.global_assembly_scene(screenshot_values, geometry)
    _assert(canonical["configured_scintillator_gap_mm"] == -80.0, "screenshot D1 gap was not reproduced")
    _assert(canonical["relative_placement"]["cherenkov_signed_offset_mm"] == [0.0, 0.0, 0.0], "Cherenkov offset changed")
    _assert(canonical["relative_placement"]["scintillator_signed_offset_mm"] == [0.0, 0.0, 30.0], "Scintillator Z offset changed axis")
    _assert(canonical["beam_alignment_mismatches"] == [], "default screenshot impacts are not on one shared beam")
    _assert(canonical["beam_alignment"]["scintillator"]["inside_track"], "Scintillator impact is outside the finite shared track")

    payload = gui.canonical_event_view_payload(canonical, "combined")
    _assert(payload["coordinate_system"] == "global_mm", "event payload changed coordinate system")
    _assert(payload["event_view_schema"].endswith("v50"), "event payload schema is not v50")
    for assembly_type in ("cherenkov", "scintillator"):
        canonical_transform = canonical["assembly_by_type"][assembly_type]["transform"]
        _assert(payload["path_transforms"][assembly_type] == canonical_transform, f"{assembly_type} path transform is not canonical")
        _assert(payload["event_geometry_transforms"][assembly_type] == canonical_transform, f"{assembly_type} event geometry transform diverged")

    # Use valid physical run settings for transport; the screenshot's -80 is
    # intentionally a display diagnostic and is not consumed as a placement.
    run_values = dict(screenshot_values)
    run_values["scintillator_gap_mm"] = "10"
    settings = gui.validate_settings(run_values)
    models, _ = gui._canonical_global_models(run_values, geometry)
    particle = gui._particle_for_settings(settings)
    optical = dataclasses.replace(gui._optical_for_settings(settings), photon_count=4)
    assembly = gui.scintillator_assembly_from_model(models["scintillator"], values=run_values, particle=particle)
    _assert(math.isclose(assembly.particle_impact_x_m, 0.0) and math.isclose(assembly.particle_impact_y_m, 0.0), "local scintillator impact was not threaded into transport")
    result = simulator.simulate(
        particle, geometry=geometry, optical=optical, seed=4801,
        detector_mode="combined", scintillator_assembly=assembly,
        scintillator_config=gui.scintillator_config_from_values(run_values),
    )
    _assert(result.photon_records and result.scintillator and result.scintillator.photon_records, "All-mode transport lost a detector-owned photon family")
    event_scene = gui.build_interactive_event_scene(result, payload, detector_mode="combined", display_photons=4)

    # No pre-projection geometry rebase is permitted: each participating
    # vertex in the event scene is the corresponding Global Assembly vertex.
    global_components = {
        str(item["component_id"]): item
        for assembly_item in canonical["assemblies"]
        for item in assembly_item["components"]
    }
    event_components = {str(item["component_id"]): item for item in event_scene["geometry_components"]}
    _assert(set(event_components) == set(global_components), "interactive event scene changed participating component ownership")
    for component_id, global_item in global_components.items():
        _assert(np.allclose(np.asarray(event_components[component_id]["faces_mm"], dtype=float), np.asarray(global_item["faces_mm"], dtype=float), atol=1e-10, rtol=0.0), f"{component_id}: event vertices differ before projection")
    _assert(event_scene["canonical_assembly_transforms"] == {key: canonical["assembly_by_type"][key]["transform"] for key in ("cherenkov", "scintillator")}, "interactive assembly transforms differ from Global Assembly")

    for assembly_type, alignment in canonical["beam_alignment"].items():
        _assert(alignment["aligned"] and alignment["inside_track"], f"{assembly_type}: shared track does not intersect impact point")
        marker = np.asarray(alignment["impact_point_mm"], dtype=float)
        projected = np.asarray(alignment["track_interpolation_point_mm"], dtype=float)
        _assert(float(np.linalg.norm(marker - projected)) <= alignment["tolerance_mm"], f"{assembly_type}: impact interpolation error")
    _assert(event_scene["track_mm"] == payload["particle_path_mm"], "interactive viewer changed the shared black track")

    # Fit data remains valid for the initial camera and each named preset;
    # this is a displayless contract check, not a GUI screenshot assertion.
    _assert(set(gui.global_assembly_view_preset_names()) >= {"Reset view", "Isometric", "View along X", "View along Y", "View along Z", "Fit assembly"}, "camera preset contract changed")
    fit = gui.global_assembly_fit(canonical, 900, 600)
    _assert(fit["extent_mm"] > 0 and all(math.isfinite(float(value)) for value in fit["center_mm"]), "global camera fit is invalid")
    _assert(event_scene["extent_mm"] > 0 and all(math.isfinite(float(value)) for value in event_scene["center_mm"]), "interactive camera fit is invalid")
    _assert(all(len(event_scene["per_detector"][name]["photon_lines"]) >= 0 for name in ("cherenkov", "scintillator")), "detector photon groups missing")

    for mode in ("cherenkov", "scintillator"):
        single_payload = gui.canonical_event_view_payload(canonical, mode)
        single_scene = gui.build_interactive_event_scene(result, single_payload, detector_mode=mode, display_photons=4)
        _assert(single_scene["visible_assembly_types"] == [mode], f"{mode}-only viewer includes another assembly")
        _assert(set(single_scene["per_detector"]) == {mode}, f"{mode}-only photon ownership changed")
        _assert(all(item["assembly_type"] == mode for item in single_scene["geometry_components"]), f"{mode}-only geometry ownership leaked")

    mismatch_values = dict(screenshot_values)
    mismatch_values["scintillator_impact_x_mm"] = "4"
    mismatch_scene = gui.global_assembly_scene(mismatch_values, geometry)
    try:
        gui.canonical_event_view_payload(mismatch_scene, "combined")
    except ValueError as exc:
        text = str(exc)
        _assert("Scintillator impact" in text and "shared-beam contract" in text and "shared Z track" in text, "impact mismatch was not reported precisely")
    else:
        raise AssertionError("separate local impact could not be represented by one combined shared trajectory")


def run_validation_and_preview_regression_checks() -> None:
    """Exercise the two v37 regressions without requiring a Tk display."""

    import copy
    import dataclasses
    import json
    import tempfile
    import quartz_geometry as geometry
    import quartz_cherenkov_gui as gui
    import quartz_cherenkov_sim as simulator

    def model(name, material):
        return geometry.AssemblyGeometryModel(
            name, [geometry.AssemblySection(name="Active", material=material, axis="X", width_y_mm=10.0, thickness_z_mm=5.0)],
            terminal_plus_surface="Detector/Photosensor coupling", terminal_minus_surface="Open",
        )

    values = gui.default_values()
    # Deliberately retain the old generic fields at 20×10 while the separate
    # canonical photosensors are valid 10×5 windows for the actual terminal.
    values.update({
        "detector_active_width_mm": "20", "detector_active_height_mm": "10",
        "sensor_1_width_mm": "10", "sensor_1_height_mm": "5",
        "sensor_2_width_mm": "10", "sensor_2_height_mm": "5",
    })
    models = {"cherenkov": model("Cherenkov", "fused silica"), "scintillator": model("Scintillator", "scintillator")}
    for mode in ("cherenkov", "scintillator", "combined"):
        values["detector_run_mode"] = mode
        values["scintillator_enabled"] = "1" if mode != "cherenkov" else "0"
        settings = gui.validate_settings(values)
        gui.validate_canonical_run_assemblies(values, models, settings.detector_run_mode)
        optical = dataclasses.replace(gui._optical_for_settings(settings), photon_count=0)
        transport = models["cherenkov"].to_simulation_geometry(sensor_shape="rectangular", sensor_width_mm=10, sensor_height_mm=5, sensor_depth_mm=2)
        assembly = simulator.compose_scintillator_assembly(
            gui.scintillator_local_geometry_from_values(values),
            gui.scintillator_global_geometry_from_values(values), transport,
        )
        scintillator_config = gui.scintillator_config_from_values(values)
        result = simulator.simulate(
            gui._particle_for_settings(settings), transport, optical, seed=17,
            scintillator_assembly=assembly if mode != "cherenkov" else None,
            scintillator_config=scintillator_config if mode != "cherenkov" else None,
            detector_mode=mode,
        )
        _assert("10 x 5" not in "".join(getattr(result, "notes", ())), f"{mode}: obsolete portal wording leaked into run notes")
    try:
        models["cherenkov"].validate_complete_assembly(sensor_width_mm=20, sensor_height_mm=10, sensor_depth_mm=2)
    except ValueError as exc:
        text = str(exc)
        _assert("actual plus_x readout terminal" in text and "10 × 5 mm" in text, "actual-terminal mismatch did not report assembly dimensions")
        _assert("rectangular portal (10 x 5 mm)" not in text, "legacy portal error text returned")
    else:
        raise AssertionError("actual terminal mismatch was accepted")

    class Var:
        def __init__(self, value=""): self.value = str(value)
        def get(self): return self.value
        def set(self, value): self.value = str(value)

    class Canvas:
        def __init__(self, width=900, height=420): self.width, self.height, self.operations, self.idle = width, height, [], []
        def winfo_exists(self): return True
        def winfo_width(self): return self.width
        def winfo_height(self): return self.height
        def delete(self, *args): self.operations.append(("delete", args, {}))
        def create_rectangle(self, *args, **kwargs): self.operations.append(("rectangle", args, kwargs))
        def create_polygon(self, *args, **kwargs): self.operations.append(("polygon", args, kwargs))
        def create_line(self, *args, **kwargs): self.operations.append(("line", args, kwargs))
        def create_text(self, *args, **kwargs): self.operations.append(("text", args, kwargs))
        def create_oval(self, *args, **kwargs): self.operations.append(("oval", args, kwargs))
        def after_idle(self, callback): self.idle.append(callback); return len(self.idle)
        def after_cancel(self, _token): pass
        def run_idle(self):
            callbacks, self.idle = self.idle, []
            for callback in callbacks: callback()

    loaded = geometry.AssemblyGeometryModel(
        "Cherenkov", [geometry.AssemblySection(name="Loaded", material="fused silica", axis="X", width_y_mm=30, thickness_z_mm=12)],
        terminal_plus_surface="Detector/Photosensor coupling", terminal_minus_surface="Open",
    )
    values["sensor_1_width_mm"], values["sensor_1_height_mm"] = "30", "12"
    sensor = gui.complete_sensor_readout_from_values(values, "cherenkov", loaded)
    payload = geometry.complete_detector_assembly_geometry_to_dict(loaded, "cherenkov", sensor)
    probe = gui.QuartzCherenkovGUI.__new__(gui.QuartzCherenkovGUI)
    probe.values = {key: Var(value) for key, value in gui.default_values().items()}
    probe.values["sensor_1_width_mm"].set("20"); probe.values["sensor_1_height_mm"].set("10")
    probe.assembly_models = {"cherenkov": model("Cherenkov", "fused silica"), "scintillator": model("Scintillator", "scintillator")}
    probe.geometry_sections = probe.assembly_models["cherenkov"].sections
    probe.geometry_canvases = {}
    probe._assembly_editor_state = {"cherenkov": {"preview": Canvas(), "preview_text": Var(), "preview_redraw_after": None}}
    probe._assembly_refresh_tree = lambda *_args: None
    probe._assembly_build_properties = lambda *_args: None
    probe._refresh_global_assembly_preview = lambda *_args: None
    probe._refresh_transport_diagnostic = lambda *_args: None
    probe.status = Var("Ready")
    probe.geometry_viewer = None; probe.geometry_viewer_canvas = None
    probe.local_viewer = None; probe.local_viewer_canvas = None

    class Dialog:
        path = ""
        @staticmethod
        def askopenfilename(**_kwargs): return Dialog.path
    class Box:
        @staticmethod
        def showerror(*_args, **_kwargs): pass
    old_dialog, old_box = gui.filedialog, gui.messagebox
    try:
        gui.filedialog, gui.messagebox = Dialog, Box
        with tempfile.TemporaryDirectory() as directory:
            path = Path(directory) / "complete.json"
            path.write_text(json.dumps(payload), encoding="utf-8")
            Dialog.path = str(path)
            probe._load_assembly_model("cherenkov")
            canvas = probe._assembly_editor_state["cherenkov"]["preview"]
            canvas.run_idle()
            _assert(probe.assembly_models["cherenkov"].sections[0].name == "Loaded", "valid atomic load did not commit model")
            scene = probe._assembly_editor_state["cherenkov"].get("scene")
            _assert(scene and all(view in scene["views"] for view in ("xz", "xy", "yz")), "atomic load did not rebuild all orthographic payloads")
            titles = [op[2].get("text") for op in canvas.operations if op[0] == "text"]
            _assert(any("Longitudinal XZ" in str(value) for value in titles) and any("Transverse XY" in str(value) for value in titles) and any("Terminal YZ" in str(value) for value in titles), "orthographic canvases were not drawn after load")
            committed_model = copy.deepcopy(probe.assembly_models["cherenkov"].to_dict())
            committed_values = {key: value.get() for key, value in probe.values.items()}
            path.write_text("{bad", encoding="utf-8")
            probe._load_assembly_model("cherenkov")
            _assert(probe.assembly_models["cherenkov"].to_dict() == committed_model, "invalid JSON changed the model")
            _assert({key: value.get() for key, value in probe.values.items()} == committed_values, "invalid JSON changed sensor/model values")
            # A first unmapped canvas is harmless; the queued configure/map
            # redraw uses its eventual size and does not need a tab switch.
            canvas.width, canvas.height = 1, 1
            probe._assembly_draw_preview("cherenkov")
            canvas.width, canvas.height = 900, 420
            probe._schedule_assembly_preview_redraw("cherenkov")
            canvas.run_idle()
            titles = [op[2].get("text") for op in canvas.operations if op[0] == "text"]
            _assert(any("Longitudinal XZ" in str(value) for value in titles), "mapped/resized preview did not recover from 1×1 layout")
    finally:
        gui.filedialog, gui.messagebox = old_dialog, old_box


def run_progress_diagnostics_checks() -> None:
    """Verify honest progress, phase heartbeats, and worker error routing without Tk."""

    import queue
    import threading
    import time
    from dataclasses import asdict, replace
    import quartz_cherenkov_gui as gui
    import quartz_cherenkov_sim as simulator

    values = gui.default_values()
    values.update({
        "events": "2", "photons_per_event": "8", "scintillator_photon_count": "8",
        "scintillator_enabled": "1",
    })
    settings = gui.validate_settings(values)
    particle = gui._particle_for_settings(settings)
    geometry = simulator.BarGeometry()
    optical = replace(gui._optical_for_settings(settings), photon_count=8, max_bounces=8)
    settings_values = asdict(settings)
    assembly = simulator.compose_scintillator_assembly(
        gui.scintillator_local_geometry_from_values(settings_values),
        gui.scintillator_global_geometry_from_values(settings_values), geometry,
    )
    scintillator_config = replace(gui.scintillator_config_from_values(settings_values), enabled=True, photon_count=8)

    for mode in ("cherenkov", "scintillator", "combined"):
        statuses, progress = [], []
        results = simulator.simulate_events(
            particle, 2, geometry=geometry, optical=optical, seed=41,
            progress_callback=lambda *item: progress.append(item),
            status_callback=statuses.append, scintillator_assembly=assembly,
            scintillator_config=scintillator_config, detector_mode=mode,
        )
        _assert(len(results) == 2, f"{mode}: event count changed")
        _assert(progress and progress[-1] == (1, 2, 1.0), f"{mode}: no honest 100% terminal progress")
        _assert(any("Starting" in item and "event 1/2" in item for item in statuses), f"{mode}: missing immediate event phase")
        _assert(any("Transporting" in item for item in statuses), f"{mode}: missing transport status")
        if mode == "combined":
            _assert(any("Cherenkov" in item for item in statuses) and any("scintillator" in item for item in statuses), "combined: detector phase missing")

    cancelled = threading.Event()
    cancelled.set()
    cancel_status = []
    cancelled_results = simulator.simulate_events(
        particle, 1, geometry=geometry, optical=optical, seed=43,
        status_callback=cancel_status.append, stop_event=cancelled,
        scintillator_assembly=assembly, scintillator_config=scintillator_config,
        detector_mode="combined",
    )
    _assert(cancelled_results == [] and any("Cancellation requested" in item for item in cancel_status), "cancellation path changed")

    class Var:
        def __init__(self, value=""): self.value = str(value)
        def set(self, value): self.value = str(value)
        def get(self): return self.value

    class Progress:
        def __init__(self): self.value = None
        def configure(self, **kwargs):
            if "value" in kwargs: self.value = kwargs["value"]
        def stop(self): pass

    class Button:
        def configure(self, **_kwargs): pass

    probe = gui.QuartzCherenkovGUI.__new__(gui.QuartzCherenkovGUI)
    probe.queue = queue.Queue()
    probe.progress_started_at = time.monotonic()
    probe.progress_total_events = 1
    probe.progress_event_index = 0
    probe.progress_event_fraction = 0.0
    probe.progress = Progress()
    probe.progress_text = Var()
    probe.status = Var()
    probe.run_state = "running"
    probe.after = lambda *_args: None
    probe._append_log = lambda _message: None
    probe._update_status("Building assemblies")
    _assert(probe.progress.value is None and "Building assemblies" in probe.status.get(), "phase status fabricated percentage")
    probe.queue.put(("status", "Transporting Cherenkov photons"))
    probe._poll_queue()
    _assert(probe.progress.value is None and "Transporting" in probe.status.get(), "heartbeat did not preserve honest percentage")

    original_runner = gui.run_simulation_to_directory
    try:
        def slow_fake(*_args, progress_callback=None, status_callback=None, **_kwargs):
            status_callback("Starting event 1/1")
            time.sleep(0.02)
            status_callback("Transporting Cherenkov photons")
            progress_callback(0, 1, 0.0)
            time.sleep(0.02)
            progress_callback(0, 1, 1.0)
            return {"summary": {"events": 1}, "run_directory": str(Path(tempfile.gettempdir()) / "fake")}

        gui.run_simulation_to_directory = slow_fake
        probe.queue = queue.Queue()
        probe.stop_requested = threading.Event()
        probe.assembly_models = {}
        fake_run_directory = Path(tempfile.gettempdir()) / "fake"
        probe._worker_run(settings, fake_run_directory, geometry)
        queued = []
        while True:
            try: queued.append(probe.queue.get_nowait())
            except queue.Empty: break
        _assert(any(kind == "status" and "Transporting" in str(payload) for kind, payload in queued), "slow worker heartbeat was not queued")
        _assert(any(kind == "done" for kind, _payload in queued), "slow worker did not finish")

        def raising_fake(*_args, status_callback=None, **_kwargs):
            status_callback("Building assemblies")
            raise RuntimeError("deliberate worker failure")

        gui.run_simulation_to_directory = raising_fake
        probe.queue = queue.Queue()
        probe._worker_run(settings, fake_run_directory, geometry)
        queued = []
        while True:
            try: queued.append(probe.queue.get_nowait())
            except queue.Empty: break
        _assert(any(kind == "error" and "deliberate worker failure" in payload for kind, payload in queued), "worker exception was swallowed")
        probe.run_button = Button()
        probe.stop_button = Button()
        old_messagebox = gui.messagebox
        class MessageBox:
            @staticmethod
            def showerror(*_args, **_kwargs): pass
        gui.messagebox = MessageBox
        try:
            probe._finish_error("RuntimeError: deliberate worker failure")
        finally:
            gui.messagebox = old_messagebox
        _assert(probe.run_state == "error" and probe.status.get() == "Error", "worker error did not exit running state")
    finally:
        gui.run_simulation_to_directory = original_runner


def global_assembly_geometry_clarity_checks() -> None:
    """v50 face, shared-track, measurement, and reserved-callout checks."""

    import math
    import quartz_cherenkov_gui as gui

    values = gui.default_values()
    values.update({
        "detector_run_mode": "All", "detector2_offset_z_mm": "40",
        "entry_x_mm": "0", "entry_y_mm": "0",
        "scintillator_impact_x_mm": "0", "scintillator_impact_y_mm": "0",
    })
    scene = gui.global_assembly_scene(values, gui.simulator.BarGeometry())
    _assert(scene["beam_alignment_mismatches"] == [], "v50 aligned screenshot-like All scene failed")
    _assert(all(item["active_volume_contains"] and item["line_intersects_active_volume"] for item in scene["beam_alignment"].values()), "v50 active-volume segment check failed")
    _assert(scene["relative_global_assembly_offset_mm"] == [0.0, 0.0, 40.0], "v50 relative offset semantics changed")
    _assert(scene["cher_d1_gap_mm"] != scene["d1_d2_separation_mm"], "v50 distance definitions still alias")
    _assert("current assemblies" not in str(scene["measurement_definitions"]), "ambiguous current-assemblies label survived")

    readout_metrics = []
    for assembly in scene["assemblies"]:
        sensor = next(item for item in assembly["components"] if item.get("kind") == "sensor")
        interface = next(item for item in assembly["components"] if item.get("kind") == "readout_interface")
        sensor_plane = {round(float(point[0]), 9) for face in sensor["window_faces_mm"] for point in face}
        interface_plane = {round(float(point[0]), 9) for face in interface["faces_mm"] for point in face}
        _assert(sensor_plane == interface_plane, f"{assembly['assembly_type']}: sensor/readout planes diverged")
        readout_metrics.append(f"{assembly['assembly_type']} plane={sorted(sensor_plane)}")

    for width, height in ((1120, 650), (760, 480), (500, 350), (320, 240)):
        items = []
        for index, component in enumerate(scene["render_components"]):
            items.append({
                "component_id": component["component_id"], "assembly_type": component.get("assembly_type"),
                "kind": component.get("kind"), "label": gui.detector_callout_label(component),
                "detail": component.get("label"), "target": (80.0 + (index % 5) * 75.0, 160.0 + (index // 5) * 32.0),
            })
        layout = gui.collision_aware_callout_layout(items, width, height, occupied_rectangles=gui.scene_reserved_rectangles(width, height))
        rectangles = [item["rect"] for item in layout["callouts"]]
        for index, first in enumerate(rectangles):
            for second in rectangles[index + 1:]:
                _assert(not gui._rectangles_overlap(first, second, gap=0.0), f"v50 callout overlap at {width}x{height}")

    print(
        "v50 metrics: mode=All active_sources=C:{c},S:{s} accepted=C:{ca},S:{sa} "
        "excluded_reasons=none unit/frame=global_mm track_residual_mm=C:{cr:.12g},S:{sr:.12g} "
        "readout_planes={planes} cher_d1_gap_mm={gap:g} d1_d2_separation_mm={sep:g} "
        "relative_offset_mm={offset} labels_rendered={shown} omitted/aggregated={omitted}".format(
            c=len(scene["beam_alignment"]["cherenkov"]["active_sections"]),
            s=len(scene["beam_alignment"]["scintillator"]["active_sections"]),
            ca=scene["beam_alignment"]["cherenkov"]["aligned"],
            sa=scene["beam_alignment"]["scintillator"]["aligned"],
            cr=scene["beam_alignment"]["cherenkov"]["distance_mm"],
            sr=scene["beam_alignment"]["scintillator"]["distance_mm"],
            planes=";".join(readout_metrics), gap=scene["cher_d1_gap_mm"],
            sep=scene["d1_d2_separation_mm"], offset=scene["relative_global_assembly_offset_mm"],
            shown="headless-layout-tested", omitted="finite-size-dependent",
        )
    )


def global_preview_layout_z_beam_checks() -> None:
    """Regression for the embedded main preview's physical Z beam and layout."""

    import math
    import numpy as np
    import quartz_cherenkov_gui as gui

    values = gui.default_values()
    values.update({
        "detector_run_mode": "All", "global_beam_axis": "Z",
        "global_offset_x_mm": "0", "global_offset_y_mm": "0", "global_offset_z_mm": "0",
        "detector2_offset_x_mm": "0", "detector2_offset_y_mm": "0", "detector2_offset_z_mm": "50",
        "assembly_1_rotation_axis": "None", "assembly_1_rotation_angle_deg": "0",
        "assembly_2_rotation_axis": "None", "assembly_2_rotation_angle_deg": "0",
        "entry_x_mm": "0", "entry_y_mm": "0",
        "scintillator_impact_x_mm": "0", "scintillator_impact_y_mm": "0",
    })
    scene = gui.global_assembly_scene(values, gui.simulator.BarGeometry())
    _assert(scene["global_beam_axis_vector"] == [0.0, 0.0, 1.0], "canonical default beam vector is not global +Z")
    _assert(np.allclose(scene["particle_direction_global"], [0.0, 0.0, 1.0], atol=1e-15), "default particle direction is not global +Z")
    _assert(scene["beam_alignment_mismatches"] == [], "default Z beam does not align with both impacts")
    for assembly_type, alignment in scene["beam_alignment"].items():
        _assert(alignment["distance_mm"] <= 1e-7, f"{assembly_type}: impact residual exceeds 1e-7 mm")
        _assert(alignment["active_volume_contains"] and alignment["line_intersects_active_volume"], f"{assembly_type}: Z track misses active volume")
        _assert(0.0 <= alignment["track_fraction"] <= 1.0, f"{assembly_type}: impact is outside finite track")

    transforms = scene["rigid_transforms"]
    _assert(transforms["cherenkov"]["translation_z_mm"] == 0.0, "Cherenkov transform changed")
    _assert(transforms["scintillator"]["translation_z_mm"] == 50.0, "Scintillator Z offset changed")
    _assert(scene["relative_global_assembly_offset_mm"] == [0.0, 0.0, 50.0], "relative offset is not a Z-only 50 mm translation")

    for width, height in ((760, 480), (520, 360), (420, 260)):
        layout = gui.global_assembly_preview_layout(scene, width, height)
        bounds = layout["geometry_bounds_px"]
        title_band = (0.0, 0.0, layout["width"], 42.0)
        axis_band = (0.0, layout["height"] - 42.0, layout["width"], layout["height"])
        _assert(not gui._rectangles_overlap(bounds, title_band, gap=0.0), f"geometry overlaps title band at {width}x{height}")
        _assert(not gui._rectangles_overlap(bounds, axis_band, gap=0.0), f"geometry overlaps axis band at {width}x{height}")
        _assert(layout["plot_bounds_px"][1] >= title_band[3] and layout["plot_bounds_px"][3] <= axis_band[1], "plot bounds are not outside reserved bands")
        _assert(layout["axis_vectors_projected"] == gui.global_preview_axis_vectors(), "triad basis differs from preview projection")
        _assert(layout["camera_yaw"] == gui.GLOBAL_PREVIEW_DEFAULT_YAW and layout["camera_pitch"] == gui.GLOBAL_PREVIEW_DEFAULT_PITCH, "embedded preview camera defaults drifted")
        content = layout["content_bounds_px"]
        plot = layout["plot_bounds_px"]
        _assert(plot[0] <= content[0] <= content[2] <= plot[2] and plot[1] <= content[1] <= content[3] <= plot[3], "fitted content extends outside plot bounds")
        for point in layout["beam_path_px"]:
            _assert(layout["plot_bounds_px"][0] <= point[0] <= layout["plot_bounds_px"][2] and layout["plot_bounds_px"][1] <= point[1] <= layout["plot_bounds_px"][3], "beam endpoint was omitted from preview fit")
        start = gui._global_preview_raw_point(gui._global_assembly_render_point(scene["particle_path_mm"][0], scene))
        end = gui._global_preview_raw_point(gui._global_assembly_render_point(scene["particle_path_mm"][1], scene))
        path_vector = np.asarray(end) - np.asarray(start)
        z_vector = np.asarray(layout["axis_vectors_projected"]["Z"])
        projected_cross = float(path_vector[0] * z_vector[1] - path_vector[1] * z_vector[0])
        _assert(abs(projected_cross) < 1e-9, "rendered particle path is not projected global Z")

    # Execute the actual embedded callback with a display-free Canvas mock.
    class Var:
        def __init__(self, value=""): self.value = value
        def set(self, value): self.value = str(value)
        def get(self): return self.value

    class Canvas:
        def __init__(self): self.operations = []
        def delete(self, *args): self.operations.append(("delete", args, {}))
        def winfo_width(self): return 760
        def winfo_height(self): return 480
        def create_polygon(self, *args, **kwargs): self.operations.append(("create_polygon", args, kwargs))
        def create_line(self, *args, **kwargs): self.operations.append(("create_line", args, kwargs))
        def create_oval(self, *args, **kwargs): self.operations.append(("create_oval", args, kwargs))
        def create_rectangle(self, *args, **kwargs): self.operations.append(("create_rectangle", args, kwargs))
        def create_text(self, *args, **kwargs): self.operations.append(("create_text", args, kwargs))

    probe = gui.QuartzCherenkovGUI.__new__(gui.QuartzCherenkovGUI)
    probe.values = values
    probe.assembly_models = {}
    probe._geometry_for_simulation = lambda: gui.simulator.BarGeometry()
    probe.global_preview_canvas = Canvas()
    probe.global_preview_text = Var()
    probe.global_geometry_status = Var()
    probe.global_viewer = None
    probe.global_viewer_canvas = None
    probe.global_preview_info_text = None
    probe._refresh_global_assembly_preview()
    operations = probe.global_preview_canvas.operations
    green_lines = [item for item in operations if item[0] == "create_line" and item[2].get("fill") == gui.GLOBAL_PARTICLE_PATH_COLOR]
    _assert(green_lines and any(item[2].get("dash") for item in green_lines), "embedded preview did not draw the green shared particle path")
    texts = "\n".join(str(item[2].get("text", "")) for item in operations if item[0] == "create_text")
    _assert("+Z / beam axis" not in texts and "MEASUREMENTS" not in texts, "main canvas still contains the old overlay annotations")
    _assert("particle axis: global +Z" in probe.global_preview_text.get(), "main diagnostic does not expose global +Z")
    print("v57 embedded preview: camera-coupled global XYZ, global +Z beam, C/S residuals <=1e-7 mm, active-volume intersections, transform preservation, reserved-band fit, and real Canvas callback PASS")


def v59_orthographic_and_run_angle_checks() -> None:
    """Verify the new global panes and the single Run inclination authority."""

    import math
    import numpy as np
    import quartz_cherenkov_gui as gui
    import quartz_cherenkov_sim as sim

    values = gui.default_values()
    values.update({
        "detector_run_mode": "All", "global_beam_axis": "Z",
        "global_offset_x_mm": "0", "global_offset_y_mm": "0", "global_offset_z_mm": "0",
        "detector2_offset_x_mm": "0", "detector2_offset_y_mm": "0", "detector2_offset_z_mm": "50",
        "assembly_1_rotation_axis": "None", "assembly_1_rotation_angle_deg": "0",
        "assembly_2_rotation_axis": "None", "assembly_2_rotation_angle_deg": "0",
        "entry_x_mm": "0", "entry_y_mm": "0", "scintillator_impact_x_mm": "0", "scintillator_impact_y_mm": "0",
    })
    settings = gui.validate_settings(values)
    particle = gui._particle_for_settings(settings)
    _assert(np.allclose(particle.direction, [0.0, 0.0, 1.0], atol=1e-15), "Run 0-degree direction is not +Z")
    values["inclination_deg"] = "10"
    settings10 = gui.validate_settings(values)
    particle10 = gui._particle_for_settings(settings10)
    expected10 = [math.sin(math.radians(10.0)), 0.0, math.cos(math.radians(10.0))]
    _assert(np.allclose(particle10.direction, expected10, atol=1e-15), "+10-degree direction does not use [sin(theta),0,cos(theta)]")
    values["inclination_deg"] = "5"
    scene = gui.global_assembly_scene(values, sim.BarGeometry())
    _assert(np.allclose(scene["particle_direction_global"], [math.sin(math.radians(5)), 0.0, math.cos(math.radians(5))], atol=1e-15), "global scene angle is not the Run direction")
    _assert(not scene["beam_alignment_mismatches"], "default +5-degree All line misses an active volume")
    _assert(scene["beam_alignment"]["scintillator"]["distance_mm"] > 4.0, "tilted scene did not preserve the saved scintillator marker")
    _assert(scene["beam_alignment"]["scintillator"]["line_intersects_active_volume"], "tilted line does not intersect scintillator active volume")

    for width, height, expected_columns, expected_rows in ((900, 420, 3, 1), (520, 420, 2, 2), (420, 300, 2, 2)):
        layout = gui.global_assembly_orthographic_layout(scene, width, height)
        _assert(layout["columns"] == expected_columns and layout["rows"] == expected_rows, f"responsive grid wrong at {width}x{height}")
        _assert(tuple(layout["view_order"]) == ("xz", "xy", "yz"), "orthographic view order changed")
        for name, pane in layout["views"].items():
            left, top, right, bottom = pane["rect"]
            _assert(0 <= left < right <= layout["width"] and 0 <= top < bottom <= layout["height"], f"{name} pane is clipped")
            _assert(pane["axis_labels"] in (("X (mm)", "Z / beam (mm)"), ("X (mm)", "Y (mm)"), ("Y (mm)", "Z / beam (mm)")), f"{name} axis labels are not physical")
            _assert(pane["plot_rect"][0] <= pane["origin_px"][0] <= pane["plot_rect"][2], f"{name} origin not in plot")
        _assert(len(layout["path_mm"]) == len(scene["particle_path_mm"]) and layout["path_mm"] == tuple(scene["particle_path_mm"]), "orthographic panes do not share one particle path")

    class Canvas:
        def __init__(self): self.operations = []
        def delete(self, *args): self.operations.append(("delete", args, {}))
        def create_rectangle(self, *args, **kwargs): self.operations.append(("create_rectangle", args, kwargs))
        def create_polygon(self, *args, **kwargs): self.operations.append(("create_polygon", args, kwargs))
        def create_line(self, *args, **kwargs): self.operations.append(("create_line", args, kwargs))
        def create_oval(self, *args, **kwargs): self.operations.append(("create_oval", args, kwargs))
        def create_text(self, *args, **kwargs): self.operations.append(("create_text", args, kwargs))

    probe = gui.QuartzCherenkovGUI.__new__(gui.QuartzCherenkovGUI)
    canvas = Canvas()
    probe._draw_global_orthographic_preview(canvas, scene, 900, 420)
    titles = "\n".join(str(item[2].get("text", "")) for item in canvas.operations if item[0] == "create_text")
    _assert(all(title in titles for title in ("XZ • longitudinal / beam", "XY • transverse", "YZ • terminal / transverse")), "not all three orthographic panes rendered")
    _assert(sum(1 for item in canvas.operations if item[0] == "create_line" and item[2].get("fill") == gui.GLOBAL_PARTICLE_PATH_COLOR) >= 3, "shared particle path was not rendered in all panes")
    _assert("X (mm)" in titles and "Y (mm)" in titles and "Z / beam (mm)" in titles, "physical axis labels missing")
    print("v59 orthographic preview + Run angle: exact transformed XZ/XY/YZ panes, responsive 3/2x2 layout, shared path, origin/axes, 0° and +10° direction, +5° All intersection PASS")


def camera_coupled_axes_checks() -> None:
    """Prove the gizmo is a projected global basis, not static screen arrows."""

    import math
    import numpy as np
    import quartz_cherenkov_gui as gui

    basis_a = gui.global_assembly_camera_basis(0.0, 0.55)
    basis_b = gui.global_assembly_camera_basis(math.pi / 2.0, 0.55)
    for basis in (basis_a, basis_b, gui.global_assembly_camera_basis(math.pi, -0.4)):
        for axis in ("X", "Y", "Z"):
            _assert(abs(np.linalg.norm(basis[axis]) - 1.0) < 1e-12, f"camera basis {axis} is not unit length")
    triad_a = gui.global_assembly_orientation_triad(0.0, 0.55)
    triad_b = gui.global_assembly_orientation_triad(math.pi / 2.0, 0.55)
    _assert(not np.allclose(triad_a["X"], triad_b["X"], atol=1e-12), "X gizmo endpoint did not change with camera yaw")
    _assert(np.allclose(triad_a["Z"], triad_b["Z"], atol=1e-12), "global Z changed under its own yaw rotation")

    origin = gui.global_assembly_camera_project((0, 0, 0), (0, 0, 0), 7.0, 300, 200, 0.37, 0.61)
    endpoint = gui.global_assembly_camera_project((1, 0, 0), (0, 0, 0), 7.0, 300, 200, 0.37, 0.61)
    expected = gui.global_assembly_camera_basis(0.37, 0.61)["X"]
    _assert(np.allclose((endpoint[0] - origin[0]) / 7.0, expected[0], atol=1e-12), "projected X screen-x differs from geometry transform")
    _assert(np.allclose((origin[1] - endpoint[1]) / 7.0, expected[1], atol=1e-12), "projected X screen-y differs from geometry transform")

    values = gui.default_values()
    values.update({
        "detector_run_mode": "All", "global_beam_axis": "Z",
        "global_offset_x_mm": "0", "global_offset_y_mm": "0", "global_offset_z_mm": "0",
        "detector2_offset_x_mm": "0", "detector2_offset_y_mm": "0", "detector2_offset_z_mm": "50",
        "assembly_1_rotation_axis": "None", "assembly_1_rotation_angle_deg": "0",
        "assembly_2_rotation_axis": "None", "assembly_2_rotation_angle_deg": "0",
    })
    scene = gui.global_assembly_scene(values, gui.simulator.BarGeometry())
    _assert(scene["rigid_transforms"]["cherenkov"]["translation_z_mm"] == 0.0, "camera test changed physical C transform")
    _assert(scene["rigid_transforms"]["scintillator"]["translation_z_mm"] == 50.0, "camera test changed physical S transform")
    path = scene["particle_path_mm"]
    projected_path = np.asarray(gui._global_preview_raw_point(path[-1])) - np.asarray(gui._global_preview_raw_point(path[0]))
    z_vector = np.asarray(gui.global_preview_axis_vectors()["Z"])
    _assert(abs(float(projected_path[0] * z_vector[1] - projected_path[1] * z_vector[0])) < 1e-9, "default global +Z path and projected +Z gizmo diverge")
    for yaw, pitch in ((0.0, 0.0), (0.65, 0.55), (math.pi / 2.0, 0.55), (math.pi, -0.4)):
        layout = gui.global_assembly_preview_layout(scene, 760, 480, yaw=yaw, pitch=pitch)
        _assert(layout["axis_vectors_projected"] == gui.global_preview_axis_vectors(yaw, pitch), "layout did not use requested camera state")
        for point in layout["beam_path_px"]:
            _assert(layout["plot_bounds_px"][0] <= point[0] <= layout["plot_bounds_px"][2] and layout["plot_bounds_px"][1] <= point[1] <= layout["plot_bounds_px"][3], "camera fit clipped physical beam endpoint")
    # A rigid assembly rotation changes its rendered vertices, while changing
    # camera yaw leaves the stored global transform and canonical axes alone.
    rotated = dict(values)
    rotated["assembly_2_rotation_axis"] = "X"
    rotated["assembly_2_rotation_angle_deg"] = "22"
    rotated_scene = gui.global_assembly_scene(rotated, gui.simulator.BarGeometry())
    _assert(rotated_scene["rigid_transforms"]["scintillator"]["rotation_axis"] == "X", "physical rigid rotation was not stored")
    _assert(rotated_scene["global_axes"] == ("X", "Y", "Z"), "global axes were replaced by assembly-local axes")
    print("v57 camera projection: unit basis, yaw-dependent X endpoint, rigid-transform/camera separation, +Z path alignment, and multi-angle fit PASS")


def global_viewer_initial_render_checks() -> None:
    """Execute the real global draw callback against a display-free Canvas mock."""

    import quartz_cherenkov_gui as gui

    class RecordingCanvas:
        def __init__(self, width=920, height=620):
            self.width, self.height, self.operations = width, height, []

        def winfo_exists(self):
            return True

        def winfo_width(self):
            return self.width

        def winfo_height(self):
            return self.height

        def delete(self, *args):
            self.operations.append(("delete", args, {}))

        def __getattr__(self, name):
            if name.startswith("create_"):
                def record(*args, **kwargs):
                    self.operations.append((name, args, kwargs))
                    return len(self.operations)
                return record
            raise AttributeError(name)

    class FakeViewer:
        def __init__(self):
            self.pending = []
            self.cancelled = []

        def winfo_exists(self):
            return True

        def after_idle(self, callback):
            self.pending.append(callback)
            return f"idle-{len(self.pending)}"

        def after(self, _delay, callback):
            self.pending.append(callback)
            return f"after-{len(self.pending)}"

        def after_cancel(self, token):
            self.cancelled.append(token)

    class Status:
        def __init__(self): self.messages = []
        def set(self, message): self.messages.append(str(message))

    probe = gui.QuartzCherenkovGUI.__new__(gui.QuartzCherenkovGUI)
    probe.values = gui.default_values()
    probe.assembly_models = {}
    probe._geometry_for_simulation = lambda: gui.simulator.BarGeometry()
    probe.global_viewer_yaw, probe.global_viewer_pitch, probe.global_viewer_zoom = .65, .55, 1.0
    probe.global_viewer_info = None
    probe.global_viewer_info_text = None
    probe.global_viewer_render_error = ""
    probe.status = Status()
    probe._global_viewer_draw_after = None
    canvas, viewer = RecordingCanvas(), FakeViewer()
    probe.global_viewer, probe.global_viewer_canvas = viewer, canvas

    old_font = gui.tkfont
    gui.tkfont = None
    try:
        _assert(gui.tkfont is None, "headless Canvas test must use the tuple-font fallback")
        # Open and Configure both enqueue one identity-bound idle draw.
        probe._schedule_global_assembly_3d_draw(viewer, canvas)
        probe._schedule_global_assembly_3d_draw(viewer, canvas, object())
        _assert(len(viewer.pending) == 1, "open/configure scheduled duplicate global draws")
        viewer.pending.pop(0)()
        operation_names = [operation[0] for operation in canvas.operations]
        first_red_axis = [args for name, args, kwargs in canvas.operations if name == "create_line" and kwargs.get("fill") == "#c62828"][-1]
        probe.global_viewer_yaw = 1.2
        # Exercise the direct implementation call too while Tk font creation
        # is disabled; the mock has no .tk interpreter/root to pass to Font.
        probe._draw_global_assembly_3d_view_impl(canvas, gui.global_assembly_scene(probe.values, gui.simulator.BarGeometry()), 920, 620)
        second_red_axis = [args for name, args, kwargs in canvas.operations if name == "create_line" and kwargs.get("fill") == "#c62828"][-1]
    finally:
        gui.tkfont = old_font
    _assert(operation_names[0] == "delete", "global callback did not clear the canvas")
    _assert(operation_names.count("create_polygon") > 0, "global callback submitted no geometry polygons")
    _assert(operation_names.count("create_line") > 0, "global callback submitted no geometry/axis lines")
    _assert(operation_names.count("create_text") > 0, "global callback submitted no title/legend/callout text")
    _assert(any("GLOBAL ASSEMBLY 3D" in str(kwargs.get("text", "")) for name, _args, kwargs in canvas.operations if name == "create_text"), "global title was not rendered")
    _assert(any("+X" == str(kwargs.get("text", "")) for name, _args, kwargs in canvas.operations if name == "create_text"), "camera-coupled global +X gizmo label was not rendered")
    _assert(any("camera-coupled" in str(kwargs.get("text", "")) for name, _args, kwargs in canvas.operations if name == "create_text"), "camera-coupled gizmo caption was not rendered")
    _assert(any(kwargs.get("font") == ("TkDefaultFont", 9, "bold") for name, _args, kwargs in canvas.operations if name == "create_text"), "headless renderer did not preserve bold tuple-font style")
    _assert(first_red_axis != second_red_axis, "popup XYZ triad remained screen-fixed after camera yaw")
    _assert(probe.global_viewer_render_error == "", "successful global draw recorded an error")

    # A queued callback from an old popup must not draw into the replacement.
    old_canvas, old_viewer = RecordingCanvas(), FakeViewer()
    probe.global_viewer, probe.global_viewer_canvas = old_viewer, old_canvas
    probe._global_viewer_draw_after = None
    probe._schedule_global_assembly_3d_draw(old_viewer, old_canvas)
    stale_callback = old_viewer.pending.pop(0)
    new_canvas, new_viewer = RecordingCanvas(), FakeViewer()
    probe.global_viewer, probe.global_viewer_canvas = new_viewer, new_canvas
    stale_callback()
    _assert(not old_canvas.operations and not new_canvas.operations, "stale popup callback rendered after window replacement")

    # The wrapper converts a real renderer exception into visible Canvas text.
    failing_canvas = RecordingCanvas()
    probe.global_viewer, probe.global_viewer_canvas = viewer, failing_canvas
    original_scene = gui.global_assembly_scene
    try:
        gui.global_assembly_scene = lambda *_args, **_kwargs: (_ for _ in ()).throw(RuntimeError("scene sentinel"))
        probe._draw_global_assembly_3d_view(920, 620, viewer=viewer, canvas=failing_canvas)
    finally:
        gui.global_assembly_scene = original_scene
    error_text = [str(kwargs.get("text", "")) for name, _args, kwargs in failing_canvas.operations if name == "create_text"]
    _assert(any("scene sentinel" in text for text in error_text), "renderer exception was not shown on the Canvas")
    _assert(any("scene sentinel" in message for message in probe.status.messages), "renderer exception was not sent to GUI status")
    print("v55 global viewer: actual callback mock render operations=%d polygons=%d text=%d; stale callback and visible error assertions PASS" % (
        len(operation_names), operation_names.count("create_polygon"), operation_names.count("create_text"),
    ))


def interactive_event_header_layout_checks() -> None:
    """Exercise the real event draw callback and its responsive header bands."""

    import quartz_cherenkov_gui as gui

    statuses = (
        "",
        "No photon paths match the selected categories.",
        "No photon categories selected; solid and particle track remain visible.",
        ("No photon paths match the selected categories because the detector filter excluded the selected records; "
         "full diagnostic is retained in the scrollable event details panel for resize-safe inspection and logging. " * 12),
        "Exception: event payload could not be rendered because the selected event contained an invalid detector record.",
    )
    for width, height in ((420, 300), (560, 420), (900, 650), (1600, 1000)):
        for status in statuses:
            layout = gui.interactive_event_view_layout(
                width, height, "INTERACTIVE EVENT 3D • event 0", status,
            )
            rectangles = {
                key: layout[key] for key in ("title_rect", "status_rect", "legend_rect", "scene_rect")
            }
            names = tuple(rectangles)
            for index, first_name in enumerate(names):
                for second_name in names[index + 1:]:
                    _assert(not gui._rectangles_overlap(rectangles[first_name], rectangles[second_name], gap=0.0), f"event header bands overlap at {width}x{height}: {first_name}/{second_name}")
            _assert(layout["scene_rect"][1] >= layout["legend_rect"][3], "event scene starts inside legend band")
            _assert(layout["status_rect"][1] >= layout["title_rect"][3], "event status starts inside title band")
            _assert(len(layout["status_lines"]) <= 2, "event status exceeded its two-line cap")
            _assert(layout["full_status"] == status, "full event status was not retained")
        long_layout = gui.interactive_event_view_layout(width, height, "INTERACTIVE EVENT 3D • event 0", statuses[3])
        _assert(long_layout["status_clipped"], f"long status was not capped at {width}x{height}")
        _assert("…" in long_layout["status_lines"][-1] or len(long_layout["status_lines"]) == 2, "long status has no bounded display form")

    class RecordingCanvas:
        def __init__(self, width=900, height=650):
            self.width, self.height, self.operations = width, height, []

        def winfo_exists(self):
            return True

        def winfo_width(self):
            return self.width

        def winfo_height(self):
            return self.height

        def delete(self, *args):
            self.operations.append(("delete", args, {}))

        def __getattr__(self, name):
            if name.startswith("create_"):
                def record(*args, **kwargs):
                    self.operations.append((name, args, kwargs))
                    return len(self.operations)
                return record
            raise AttributeError(name)

    face = [(-35.0, -10.0, -8.0), (35.0, -10.0, -8.0), (35.0, 10.0, 8.0), (-35.0, 10.0, 8.0)]
    scene = {
        "center_mm": (0.0, 0.0, 0.0), "extent_mm": 110.0, "mode": "combined",
        "geometry_components": [{"component_id": "C-active", "assembly_type": "cherenkov", "kind": "geometry", "label": "quartz", "faces_mm": [face]}],
        "geometry_faces": [(face, "#78a9d3", "cherenkov", "geometry")],
        "geometry_lines": [(face, "cherenkov", "geometry")],
        "track_mm": [(-48.0, 0.0, 0.0), (48.0, 0.0, 0.0)],
        "impact_markers": [],
        "photon_lines": [{"points_mm": [(-20.0, 0.0, 0.0), (25.0, 6.0, 5.0)], "color": "#1769aa", "detected_marker": True}],
    }
    old_font = gui.tkfont
    gui.tkfont = None
    try:
        probe = gui.QuartzCherenkovGUI.__new__(gui.QuartzCherenkovGUI)
        probe.interactive_event_yaw, probe.interactive_event_pitch = 0.75, 0.45
        probe.interactive_event_zoom = 1.0
        probe.selected_event = type("Selected", (), {"get": lambda _self: "0"})()
        probe._three_d_labels_enabled = lambda: False
        probe.interactive_event_filter_message = statuses[1]
        probe.interactive_event_scene = scene
        probe.interactive_event_canvas = RecordingCanvas()
        probe._draw_interactive_event_view()
        text_values = [str(kwargs.get("text", "")) for name, _args, kwargs in probe.interactive_event_canvas.operations if name == "create_text"]
        _assert(any("INTERACTIVE EVENT 3D • event 0" in value for value in text_values), "headless event callback omitted its title")
        _assert(any(statuses[1] in value for value in text_values), "headless event callback omitted its status")
        _assert(any(name == "create_polygon" for name, _args, _kwargs in probe.interactive_event_canvas.operations), "headless event callback omitted geometry")
        _assert(any(name == "create_line" for name, _args, _kwargs in probe.interactive_event_canvas.operations), "headless event callback omitted scene/path lines")

        empty_scene = dict(scene)
        empty_scene["photon_lines"] = []
        probe.interactive_event_filter_message = statuses[2]
        probe.interactive_event_scene = empty_scene
        probe.interactive_event_canvas.operations.clear()
        probe._draw_interactive_event_view()
        empty_text_values = [str(kwargs.get("text", "")) for name, _args, kwargs in probe.interactive_event_canvas.operations if name == "create_text"]
        _assert(any(statuses[2] in value for value in empty_text_values), "empty-filter status was not rendered separately")
        _assert(any(name == "create_polygon" for name, _args, _kwargs in probe.interactive_event_canvas.operations), "empty-filter event lost its solid geometry")
    finally:
        gui.tkfont = old_font
    print("v58 event viewer: measured title/status bands, two-line status clamp, resize fit, non-empty and empty-filter mock renders PASS")


def results_arrival_analysis_checks() -> None:
    """Exercise detector ownership, selected-event data, and unit-area bins."""

    import quartz_cherenkov_sim as simulator
    from types import SimpleNamespace

    def record(time_ns, detected=False):
        return SimpleNamespace(
            detector_valid_sensor_hit=True,
            readout_arrival_time_s=float(time_ns) * 1e-9,
            arrival_time_s=float(time_ns) * 1e-9,
            detected=bool(detected), outcome="detected_plus_x" if detected else "qe_loss",
        )

    def result(event_id, mode, cherenkov, scintillator=()):
        return SimpleNamespace(
            event_id=event_id, detector_mode=mode, generated_photons=len(cherenkov),
            photon_records=list(cherenkov),
            scintillator=SimpleNamespace(photon_records=list(scintillator)) if mode in {"scintillator", "combined"} else None,
        )

    combined = [result(10, "combined", [record(1, True), record(3)], [record(2), record(4, True)]), result(20, "combined", [record(5)], [record(6)])]
    analysis = simulator.readout_arrival_analysis(combined, event_id=10, bins=4)
    _assert(analysis["selected"]["detectors"]["cherenkov"]["arrivals"] == 2, "selected Cherenkov arrivals are not exact")
    _assert(analysis["selected"]["detectors"]["cherenkov"]["detected"] == 1, "detected count was conflated with arrival count")
    _assert(analysis["aggregate"]["detectors"]["scintillator"]["arrivals"] == 3, "aggregate scintillator owner stream is wrong")
    edges = analysis["aggregate"]["histograms"]["bin_edges_ns"]
    _assert(len(edges) > 1, "aggregate analysis did not create common bins")
    for detector in ("cherenkov", "scintillator"):
        histogram = analysis["aggregate"]["histograms"]["histograms"][detector]
        _assert(len(histogram["counts"]) == len(edges) - 1, "detector histograms do not share bin edges")
        integral = sum(value * width for value, width in zip(histogram["normalized_density_per_ns"], analysis["aggregate"]["histograms"]["bin_widths_ns"]))
        _assert(abs(integral - 1.0) < 1e-12, "normalized aggregate histogram is not unit area")
    c_only = simulator.readout_arrival_analysis([result(0, "cherenkov", [record(1)])], event_id=0)
    _assert(not c_only["selected"]["detectors"]["scintillator"]["included"], "disabled scintillator was shown as zero")
    rows = simulator.readout_arrival_rows(combined)
    saved = simulator.readout_arrival_analysis_from_rows(rows, [{"event_id": 10, "generated_photons": 2, "scintillator_generated_photons": 2}, {"event_id": 20, "generated_photons": 1, "scintillator_generated_photons": 1}], "combined", event_id=10)
    _assert(saved["selected"]["detectors"]["cherenkov"]["arrivals"] == 2, "persisted rows cannot rebuild selected analysis")
    print("v60 Results readout arrival ownership/count/normalization/persistence checks PASS")


def sensor_radius_geometry_checks() -> None:
    """Validate explicit circular/elliptical apertures for both assemblies."""

    import quartz_cherenkov_sim as simulator
    from quartz_geometry import (
        AssemblyGeometryModel, AssemblySection, complete_detector_assembly_geometry_from_dict,
        complete_detector_assembly_geometry_to_dict,
    )
    import quartz_cherenkov_gui as gui

    def make_model(name):
        return AssemblyGeometryModel(
            name=name,
            sections=[AssemblySection(name="active", shape="box", material="fused silica", length_x_mm=100.0, width_y_mm=20.0, thickness_z_mm=10.0)],
            terminal_plus_surface="Detector/Photosensor coupling",
        )

    for assembly in ("cherenkov", "scintillator"):
        index = "1" if assembly == "cherenkov" else "2"
        model = make_model(assembly.title())
        values = gui.default_values()
        values.update({f"sensor_{index}_shape": "circular", f"sensor_{index}_radius_mm": "0.7"})
        circle = gui.complete_sensor_readout_from_values(values, assembly, model)
        _assert(circle["readout_radius_mm"] == 0.7 and circle["dimensions_mm"]["width"] == 1.4, f"{assembly}: circular diameter migration/serialization is wrong")
        transport = model.to_simulation_geometry(sensor_shape="circular", sensor_width_mm=1.4, sensor_height_mm=1.4, sensor_radius_mm=0.7, sensor_depth_mm=2.0)
        optical = simulator.OpticalConfig()
        _assert(simulator.active_area_contains((0.05, 0.0007, 0.0), transport, optical, "plus_x"), f"{assembly}: circular boundary point rejected")
        _assert(not simulator.active_area_contains((0.05, 0.000701, 0.0), transport, optical, "plus_x"), f"{assembly}: circular point outside radius accepted")

        values.update({f"sensor_{index}_shape": "elliptical", f"sensor_{index}_semi_axis_a_mm": "1.0", f"sensor_{index}_semi_axis_b_mm": "0.4"})
        ellipse = gui.complete_sensor_readout_from_values(values, assembly, model)
        _assert(ellipse["readout_semiaxes_mm"] == {"a": 1.0, "b": 0.4}, f"{assembly}: ellipse semiaxes were not serialized")
        transport = model.to_simulation_geometry(sensor_shape="elliptical", sensor_width_mm=2.0, sensor_height_mm=0.8, sensor_semi_axis_a_mm=1.0, sensor_semi_axis_b_mm=0.4, sensor_depth_mm=2.0)
        _assert(simulator.active_area_contains((0.05, 0.001, 0.0), transport, optical, "plus_x"), f"{assembly}: ellipse A boundary rejected")
        _assert(simulator.active_area_contains((0.05, 0.0, 0.0004), transport, optical, "plus_x"), f"{assembly}: ellipse B boundary rejected")
        _assert(not simulator.active_area_contains((0.05, 0.0009, 0.00035), transport, optical, "plus_x"), f"{assembly}: point outside ellipse but inside bounding rectangle accepted")
        payload = complete_detector_assembly_geometry_to_dict(model, assembly, ellipse)
        _, restored = complete_detector_assembly_geometry_from_dict(payload, assembly, name=assembly.title())
        _assert(restored["readout_semiaxes_mm"] == {"a": 1.0, "b": 0.4}, f"{assembly}: JSON round trip changed ellipse")
        legacy = copy.deepcopy(payload)
        legacy["photosensor_readout"].pop("readout_radius_mm", None)
        legacy["photosensor_readout"].pop("readout_semiaxes_mm", None)
        _, migrated = complete_detector_assembly_geometry_from_dict(legacy, assembly, name=assembly.title())
        _assert(migrated["readout_semiaxes_mm"] == {"a": 1.0, "b": 0.4}, f"{assembly}: v61 legacy ellipse migration changed physical aperture")
    print("v62 sensor radius/ellipse geometry: C/S circular boundary, ellipse boundary/outside, solid, transport, JSON and v61 migration checks PASS")


def gamma_transport_targeted_checks() -> None:
    """Run the focused range-diagnostic, closure, and ancestry regressions."""

    import test_gamma_transport

    suite = unittest.defaultTestLoader.loadTestsFromModule(test_gamma_transport)
    result = unittest.TextTestRunner(verbosity=2).run(suite)
    _assert(result.wasSuccessful(), "targeted Gamma transport checks failed")


def gamma_global_viewer_checks() -> None:
    """Check global Gamma scene contents and detached-window controls headlessly."""

    import quartz_cherenkov_gui as gui

    square_a = [[0.0, 0.0, 0.0], [4.0, 0.0, 0.0], [4.0, 4.0, 0.0], [0.0, 4.0, 0.0]]
    square_b = [[0.0, 0.0, 2.0], [4.0, 0.0, 2.0], [4.0, 4.0, 2.0], [0.0, 4.0, 2.0]]
    geometry_snapshot = {
        "schema": "gamma-global-assembly-geometry-v1",
        "scene": {"assemblies": [
            {"assembly_type": assembly, "components": [{
                "component_id": assembly + "-active", "kind": "geometry", "name": assembly,
                "color": "#8db6ce", "faces_mm": [square_a, square_b],
            }]} for assembly in ("cherenkov", "scintillator")
        ]},
    }
    event = {
        "event_id": 7, "source_energy_keV": 125.0, "source_component": "continuum",
        "gamma_tracks": [{
            "start_position_m": [0.0, 0.0, -0.001], "end_position_m": [0.0, 0.0, 0.002],
            "volume": "quartz", "interaction": "compton", "energy_keV": 125.0,
        }],
        "charged_tracks": [{
            "position_m": [0.0, 0.0, 0.001], "direction": [1.0, 0.0, 0.0],
            "path_mm": 3.0, "kind": "electron", "volume": "quartz",
        }],
        "deposits": [{
            "position_m": [0.001, 0.0, 0.001], "interaction": "compton",
            "energy_keV": 14.0, "volume": "quartz",
        }],
    }
    scene = gui.build_gamma_global_view_scene(event, geometry_snapshot, [0.0, 0.0, 1.0])
    _assert(scene["geometry_available"] and len(scene["geometry_components"]) == 2, "detached Gamma scene omitted saved global geometry")
    _assert(scene["source_position_mm"] == [0.0, 0.0, -1.0] and scene["source_ray_end_mm"] is not None, "Gamma source position/direction was not placed in the global frame")
    _assert(scene["gamma_segments"][0]["end_mm"] == (0.0, 0.0, 2.0), "Gamma path was not converted to global millimetres")
    _assert(scene["charged_segments"][0]["end_mm"] == (3.0, 0.0, 1.0), "charged secondary track was not included in the global scene")
    _assert(scene["deposits"][0]["position_mm"] == (1.0, 0.0, 1.0), "energy deposit was not included in the global scene")
    segmented_view = gui.gamma_event_view_data({"charged_tracks": [{
        "kind": "electron", "volume": "old aggregate",
        "segments": [
            {"start_position_m": [0.0, 0.0, 0.0], "end_position_m": [0.001, 0.0, 0.0], "path_mm": 1.0, "volume": "quartz", "deposited_energy_keV": 25.0},
            {"start_position_m": [0.001, 0.0, 0.0], "end_position_m": [0.001, 0.002, 0.0], "path_mm": 2.0, "volume": "scintillator", "deposited_energy_keV": 75.0},
        ],
    }]})["charged_segments"]
    _assert(len(segmented_view) == 2 and [row["volume"] for row in segmented_view] == ["quartz", "scintillator"], "Gamma viewer collapsed material-specific charged segments")
    _assert(segmented_view[1]["end_mm"] == (1.0, 2.0, 0.0) and segmented_view[1]["deposited_energy_keV"] == 75.0, "Gamma viewer lost charged segment path or deposit energy")
    _assert(gui.gamma_global_view_binding_sequences() == ("<Configure>", "<ButtonPress-1>", "<B1-Motion>", "<MouseWheel>", "<Button-4>", "<Button-5>"), "Gamma camera rotate/zoom bindings changed")
    _assert(gui.gamma_global_view_control_names() == ("Event", "Reset view", "Show geometry", "Show source", "Show tracks / deposits"), "Gamma global viewer control contract changed")
    no_geometry = gui.build_gamma_global_view_scene(event, None, [0.0, 0.0, 1.0])
    _assert(not no_geometry["geometry_available"] and len(no_geometry["gamma_segments"]) == 1, "legacy Gamma result fabricated geometry or discarded supported tracks")
    inferred = gui.build_gamma_global_view_scene({"event_id": 0}, geometry_snapshot, [0.0, 0.0, 1.0])
    _assert(inferred["source_position_mm"] == [0.0, 0.0, -1.0] and inferred["source_position_inferred"], "empty-event source fallback was not clearly marked")

    class Var:
        def __init__(self, value=""):
            self.value = str(value)
        def get(self):
            return self.value
        def set(self, value):
            self.value = str(value)

    class Table:
        def __init__(self):
            self.selected = []
            self.seen = []
        def selection_set(self, item):
            self.selected = [str(item)]
        def see(self, item):
            self.seen.append(str(item))

    class Selector:
        def __init__(self):
            self.value = ""
        def set(self, value):
            self.value = str(value)

    probe = type("GammaGlobalSelectionProbe", (), {})()
    probe.gamma_global_3d_event_choice = Var("7")
    probe.gamma_selected_event = Var("0")
    probe.latest_gamma_payload = {"events": [{"event_id": 0}, {"event_id": 7}]}
    probe.gamma_results_event_table = Table()
    probe.gamma_results_event_selector = Selector()
    probe._refresh_gamma_results_view = lambda: setattr(probe, "refresh_count", getattr(probe, "refresh_count", 0) + 1)
    gui.QuartzCherenkovGUI._gamma_global_3d_event_selected(probe)
    _assert(probe.gamma_selected_event.get() == "7" and probe.gamma_results_event_table.selected == ["7"], "detached event selection did not synchronize Gamma Results")
    _assert(probe.gamma_results_event_selector.value == "7" and probe.refresh_count == 1, "detached event selection did not refresh the result view")

    camera_probe = type("GammaCameraProbe", (), {})()
    camera_probe.gamma_global_3d_yaw = -0.96
    camera_probe.gamma_global_3d_pitch = 0.38
    camera_probe.gamma_global_3d_zoom = 1.0
    camera_probe.gamma_global_3d_drag = None
    camera_probe.redraw_count = 0
    camera_probe._draw_gamma_global_3d_view = lambda: setattr(camera_probe, "redraw_count", camera_probe.redraw_count + 1)
    class Pointer:
        x = 10
        y = 10
        delta = 120
        num = 4
    gui.QuartzCherenkovGUI._gamma_global_3d_press(camera_probe, Pointer())
    Pointer.x, Pointer.y = 35, 1000
    gui.QuartzCherenkovGUI._gamma_global_3d_drag_motion(camera_probe, Pointer())
    _assert(camera_probe.gamma_global_3d_yaw != -0.96 and camera_probe.gamma_global_3d_pitch == 1.45, "Gamma mouse drag did not rotate or clamp the camera")
    gui.QuartzCherenkovGUI._gamma_global_3d_zoom_event(camera_probe, Pointer())
    _assert(camera_probe.gamma_global_3d_zoom > 1.0, "Gamma positive wheel input did not zoom in")
    Pointer.delta = -120
    for _ in range(100):
        gui.QuartzCherenkovGUI._gamma_global_3d_zoom_event(camera_probe, Pointer())
    _assert(camera_probe.gamma_global_3d_zoom == 0.2 and camera_probe.redraw_count >= 100, "Gamma zoom limit or redraw was not enforced")
    print("Detached global Gamma viewer: saved geometry, source/direction, Gamma and charged tracks, deposits, event sync, rotate/zoom controls PASS")


def main() -> int:
    source_checks()
    gamma_transport_targeted_checks()
    gamma_global_viewer_checks()
    launcher_files = (
        "Start_GUI.bat", "setup_and_start.ps1", "Check_Dependencies.bat", "check_dependencies.ps1",
    )
    missing_launchers = [name for name in launcher_files if not (ROOT / name).is_file()]
    if missing_launchers:
        print("Windows launcher checks: SKIP; missing original files: " + ", ".join(missing_launchers))
    else:
        windows_launcher_structure_checks()
    if (ROOT / "quartz_cherenkov_simulator_v78.zip").is_file():
        package_zip_integrity_checks()
    else:
        print("Available-project ZIP checks: SKIP; v78 archive is unavailable")
    detector_selection_checks()
    gamma_viewer_persistence_checks()
    gamma_results_scroll_checks()
    scintillator_decay_time_checks()
    scintillator_transport_debug_checks()
    scintillator_v45_invariant_checks()
    scintillator_v46_complete_guide_transport_checks()
    run_status_layout_contract_checks()
    complete_global_profile_checks()
    run_bounce_and_layout_checks()
    responsive_layout_checks()
    collision_aware_label_layout_checks()
    model_checks()
    sensor_detection_accounting_checks()
    complete_v32_json_and_diagnostic_checks()
    photosensor_visibility_checks()
    scintillator_geometry_json_checks()
    cherenkov_terminal_roundtrip_checks()
    live_viewer_refresh_checks()
    canonical_global_scene_checks()
    rigid_assembly_transform_checks()
    mylar_overlay_checks()
    event_viewer_scene_checks()
    static_selected_event_renderer_checks()
    event_global_frame_alignment_checks()
    run_validation_and_preview_regression_checks()
    run_progress_diagnostics_checks()
    global_assembly_geometry_clarity_checks()
    v59_orthographic_and_run_angle_checks()
    camera_coupled_axes_checks()
    global_viewer_initial_render_checks()
    interactive_event_header_layout_checks()
    results_arrival_analysis_checks()
    sensor_radius_geometry_checks()
    print("v78 staged-source checks: detached global Gamma viewer, charged multi-volume path/energy regressions, Gamma energy closure and ancestry; preserved v77 viewer, v76 transport fixes, v74 geometry roundtrips, v73 Gamma Results, persistence, Windows/Tk fixes, and prior regression suite: PASS")
    print(f"workspace: {ROOT}")
    print("GUI visual runtime check: not run (no display server in sandbox)")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
