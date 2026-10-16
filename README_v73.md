# Quartz Cherenkov and Scintillator Simulator — v73

**Application version:** `0.41.0`  
**Build ID:** `build-20260928-v73-gamma-results-scroll`

v73 carries forward the v72 Gamma event list, selected-event viewer and detail
panel, persisted complete assembly geometry, independent Gamma result restore,
and the v71/v70 fixes. Gamma Results now has a dedicated vertical scrollbar
and mouse-wheel scrolling through the run summary, event list, and selected
event viewer. The table keeps its horizontal scrollbar, while the table and
text panes retain their own wheel scrolling.

Each new Gamma result stores both complete Cherenkov and Scintillator assembly
profiles and the transformed Global Assembly scene captured from the detached
run snapshot. Gamma transport and the viewer use the same active and guide
solids and rigid transforms. The viewer draws the complete saved detector
geometry underneath event paths, including for an event with no interactions.
Older Gamma results without a geometry snapshot are labeled clearly and shown
without substituted or guessed solids.

## Run mode and output

- **Cherenkov only**, **Scintillator only**, and **All** are optical modes.
  **All** combines Cherenkov and Scintillator optical transport.
- **Gamma transport** is a separate explicit run mode. It records gamma
  interactions, charged secondaries, energy deposits, and expected local light.
  It does not transport optical photons or calculate fiber transport, sensor
  arrivals, or detections.
- A successful Gamma result is indexed independently from optical results.
  Later optical runs do not replace it, and startup restores it even when an
  optical run is newer.

## Windows launch

Extract the complete archive into one folder and run `Start_GUI.bat`. The
launcher uses that folder's `.venv`, runs `check_environment.py`, and starts the
GUI only if the checks pass. An unusable existing environment is preserved in
a timestamped `.venv.broken-*` folder before recreation. Dependency
installation remains opt-in; see `README_windows_setup.txt`.

## Verification

`python -m py_compile` and the full `python check_environment.py` suite pass.
The suite includes targeted Gamma Results scrollbar and mouse-wheel routing
checks alongside mixed empty/populated Gamma events, event selection,
zero-interaction geometry rendering, v71 result restore, complete Global
Assembly save/load, and retained Windows/Tk regression checks. Windows scripts
were structurally checked but not executed on Windows. The optional ROOT
writer was skipped because `uproot` and `awkward` are unavailable. Interactive
GUI visual QA was unavailable in the release sandbox.
