# Quartz Cherenkov and Scintillator Simulator — v72

**Application version:** `0.40.0`  
**Build ID:** `build-20260928-v72-gamma-event-list-geometry`

v72 retains the v71 independent Gamma result index and restores the visible
v70/v71 Windows, Tk, and complete Global Assembly profile behavior. Gamma
Results now includes a selectable row for every event, with source energy and
component, interaction processes, deposit count and total, charged-track
count, expected local light when present, and an `INTERACTION` or `EMPTY`
status. Selecting a row synchronizes the event viewer and detail panel.

Each new Gamma result stores the run's two complete Cherenkov and Scintillator
assembly profiles and the transformed Global Assembly scene captured from the
detached run snapshot. Gamma transport and the viewer use the same active and
guide solids and the same rigid transforms. The viewer draws the complete
saved detector geometry underneath event paths, including for an event with
no interactions. Older Gamma results without a geometry snapshot are labeled
clearly and are shown without substituted or guessed solids.

## Run mode and output

- **Cherenkov only**, **Scintillator only**, and **All** are optical modes.
  **All** combines Cherenkov and Scintillator optical transport.
- **Gamma transport** is a separate explicit run mode. It records gamma
  interactions, charged secondaries, energy deposits, and expected local
  light. It does not transport optical photons or calculate fiber transport,
  sensor arrivals, or detections.
- A successful Gamma result is indexed independently from optical results.
  Later optical runs do not replace it, and startup restores it even when an
  optical run is newer.

## Windows launch

Extract the complete archive into one folder and run `Start_GUI.bat`. The
launcher uses that folder's `.venv`, runs `check_environment.py`, and starts
the GUI only if the checks pass. An unusable existing environment is preserved
in a timestamped `.venv.broken-*` folder before recreation. Dependency
installation remains opt-in; see `README_windows_setup.txt`.

## Verification

`python -m py_compile` and the full `python check_environment.py` suite pass.
The suite covers mixed empty/populated Gamma events, table selection,
zero-interaction geometry rendering, v71 result restore, complete Global
Assembly save/load, and the retained Windows/Tk regression checks. Windows
scripts were structurally checked but not executed on Windows. Interactive GUI
visual QA was unavailable in the release sandbox.
