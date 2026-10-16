# Quartz Cherenkov and Scintillator Simulator — v71

**Application version:** `0.39.0`  
**Build ID:** `build-20260928-v71-gamma-viewer-restore`

v71 restores the dedicated Gamma event viewer and preserves the last completed
Gamma result independently from optical runs. Gamma results are written to
their run folder and indexed in `latest_gamma_result.json` in the configured
output directory. Startup restores that Gamma result even when a newer optical
run is the latest general run folder. Existing v70 installations without an
index are migrated by scanning completed Gamma run folders.

## Run mode semantics

- **Cherenkov only**, **Scintillator only**, and **All** are optical modes.
  **All** combines the Cherenkov and Scintillator optical simulations.
- **Gamma transport** is a separate explicit run mode. It transports gamma
  photons and records interactions, charged secondaries, energy deposits, and
  expected local light. It does not transport optical photons or calculate
  fiber transport, sensor arrivals, or detections.
- Running an optical mode, including **All**, leaves the last completed Gamma
  result available in **Gamma Results**. A successful new Gamma run replaces
  that saved result. A stopped Gamma run retains the previous completed result.

The Gamma Results tab provides aggregate and per-event details plus an
interactive 3D event viewer. The viewer displays gamma flight segments,
charged-secondary tracks, and interaction deposits in global millimeters.

## Included release

The archive contains the GUI, simulator and geometry modules, Gamma transport
module and source assets, full regression checker, requirements, Windows setup
and launch scripts, READMEs, and `Assemblycompleto.json`.

## Windows launch

Extract all files to one folder and run `Start_GUI.bat`. The launcher uses that
folder's `.venv`, preserves an unusable existing environment as
`.venv.broken-<timestamp>`, runs `check_environment.py`, and opens the GUI only
after the checks pass. Dependency installation is opt-in; see
`README_windows_setup.txt`.

The Windows scripts were structurally checked in the release sandbox but were
not executed on Windows. The GUI requires an interactive desktop session.
