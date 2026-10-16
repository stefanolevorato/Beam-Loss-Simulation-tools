# Quartz Cherenkov Simulator v74

Release build: `0.42.0` / `build-20260929-v74-cherenkov-terminal-roundtrip`.

This package contains the Cherenkov and scintillator optical simulator, the
Gamma transport and event viewer, canonical geometry editor, sample Gamma
source assets, Windows launch scripts, and the headless environment/regression
checker.

## Start on Windows

Install Python 3.10 or newer with Tk support. Open PowerShell in this folder
and run `Start_GUI.bat`. The launcher creates a project `.venv`, preserves a
broken environment in a dated `.venv.broken-*` folder, runs the environment
checks, and then opens the GUI. Dependency installation is opt-in:

```powershell
powershell -ExecutionPolicy Bypass -File .\setup_and_start.ps1 -InstallDependencies
```

The GUI displays its version and build identifier in the window title and
Source/Diagnostics tab. Use `Check_Dependencies.bat` to run the installed
environment checker without opening the GUI.

## Geometry JSON

Cherenkov and scintillator geometry exports use the canonical v27 geometry
schema or the complete v32 envelope. Both serialize `terminal_plus_surface`,
`terminal_minus_surface`, and the `terminal_surfaces` map. The legacy
`opposite_end_mode` value remains as a compatibility alias for the `-X`
terminal. Modern endpoint fields are authoritative; conflicting modern
representations are rejected before the GUI changes its current geometry.
Loaded endpoint states refresh the editor, local preview, global assembly
view, Mylar exclusion state, and the legacy `-X` Run controls together.

## Verify

Run the complete headless release checker from this folder:

```powershell
python check_environment.py
```

It includes the v74 mirror roundtrip regression and the v73 Gamma, persistence,
Windows launcher structure, geometry, rendering-payload, and transport checks.
Tk window rendering still requires a desktop session.

## Package files

- `quartz_cherenkov_gui.py`, `quartz_cherenkov_sim.py`
- `quartz_geometry.py`, `quartz_geometry_editor.py`
- `gamma_transport.py`, `check_environment.py`
- `gamma_source_profile_collimated_v1.json`, `gamma_reference_lines_v1.csv`,
  and `gamma_results.json`
- `Assemblycompleto.json`
- `Start_GUI.bat`, `setup_and_start.ps1`, `Check_Dependencies.bat`,
  `check_dependencies.ps1`
- `requirements.txt` and `requirements-root.txt`
- The v62/v63 assembly and Gamma stage notes
