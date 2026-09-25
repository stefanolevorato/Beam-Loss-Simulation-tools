# Quartz Cherenkov / Scintillator two-assembly simulator — v62

Build ID: `build-20260925-v62-sensor-radius-geometry`

This release starts from the frozen v61 scintillator-decay baseline and adds
explicit editable readout geometry for both the Cherenkov and Scintillator
Sensor tabs. Existing v46–v61 behavior remains in the five packaged Python
modules, including decay time, Results/readout-arrival analysis, Run
inclination, orthographic views, complete profiles, bounce limits, guide
transport, sensor active-window accounting, canonical terminal placement, and
the v62 build identity. No fiber code or benchmark changes are included.

## Readout conventions

- Rectangular: `Readout width (mm)` and `Readout height (mm)` are full extents.
- Circular: `Readout radius (mm)` is the physical radius. Width and height are
  derived full diameter values and are hidden from editing.
- Elliptical: `Semi-axis A (mm)` and `Semi-axis B (mm)` are physical semiaxes
  in the terminal-face-local Y/Z coordinates. Width and height are derived
  full extents and are hidden from editing.
- `Readout depth (mm)` remains independent and is the extrusion depth along
  the terminal normal.

The active-window test uses the actual local face coordinates. Ellipses use
`(dy/a)^2 + (dz/b)^2 <= 1`; circular windows use `dy^2 + dz^2 <= r^2`.
Circular and elliptical readout bodies are rendered as sampled extruded
elliptical prisms, with the coupling face at the selected `+X`/`-X` terminal
or valid guide output. The same canonical component feeds local, global,
orthographic, event, JSON, and transport views.

## Compatibility

v61 and older complete JSON/profile data remains readable. Legacy
`dimensions_mm.width/height` values were full transverse extents, so a legacy
circular `2 mm × 2 mm` window migrates to radius `1 mm`; legacy ellipse values
migrate to semiaxes `width/2` and `height/2`. Saved v62 data retains normalized
full extents and adds `readout_radius_mm` or `readout_semiaxes_mm`.

All dimensions are validated as finite and strictly positive before preview,
profile commit, or Run. Cherenkov and Scintillator settings are independent.

## Verification

`check_environment.py` passes the historical v46–v61 checks plus v62 tests
for both detectors: circular center/boundary/outside acceptance, elliptical
center/boundary/outside-rectangle rejection, physical solid radius/semiaxes,
transport handoff, JSON round-trips, and v61 migration. The simulator
`--self-test` passes. GUI visual QA was not run because the sandbox has no
display server; headless scene payload generation and regression rendering
checks passed.

## Windows startup

1. Extract all package files into the same folder.
2. Double-click `Start_GUI.bat`.

On the first run, the launcher creates `.venv` in this project folder and
installs the packages listed in `requirements.txt`. Later starts reuse the
same environment and install packages only when the required imports are
missing. The project environment is activated only when it is not already the
active `VIRTUAL_ENV`; an unrelated active environment is switched to the
project `.venv` for this launcher process.

The `.venv` folder is intentionally not included in the ZIP. It is generated
on the user's Windows machine and is specific to that Python installation.
If a `.venv` folder exists without `Scripts\\python.exe`, the launcher stops
without deleting it and prints manual repair instructions.
