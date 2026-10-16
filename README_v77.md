# Quartz Cherenkov Simulator v77

Application version: `0.45.0`  
Build ID: `build-20260929-v77-gamma-global-interactive-viewer`

v77 adds **Open Interactive Global 3D** to Gamma Results. The detached viewer
uses the selected Gamma event and its saved global coordinates. When a
`gamma-global-assembly-geometry-v1` snapshot is present, it renders the
canonical Cherenkov and Scintillator assemblies at their saved global offsets.
Gamma segments, charged secondary tracks, interaction/deposit points, the
emission point, and the global source direction are shown when their event
records provide them. Drag to orbit; use the mouse wheel to zoom; choose an
event in the viewer's selector.

New Gamma event records persist `source_position_m` and
`source_direction_global`. The source marker uses the transport's actual
initial point. The direction arrow is scaled for the view and is labeled as an
illustration of direction, not a physical path length. Older results can use
the recorded primary Gamma segment and run summary when available. If a saved
result lacks a usable geometry snapshot or source coordinates, the viewer says
so and does not substitute current/default geometry.

The existing Cherenkov, Scintillator, Global Assembly, and optical event
viewers remain in place. v76 Gamma ancestry and energy-accounting fixes are
preserved.

## Verify

Run the targeted Gamma transport regressions:

```text
python -m unittest -v test_gamma_transport
```

Run syntax, viewer mapping/invocation, and headless simulator regressions:

```text
python check_environment.py
```

The checks cover saved source-point/direction metadata, global-mm event
mapping, both canonical assemblies, Gamma Results button wiring, and the
viewer popup invocation. Syntax compilation, all six focused Gamma transport
tests, and the full staged-source checker passed in this run. The checker
skipped Windows launcher checks because the four original launcher files
were not staged and skipped its optional ROOT writer self-test because
`uproot` and `awkward` are unavailable. Desktop visual testing was unavailable
in the sandbox.

## Package

`quartz_cherenkov_simulator_v77.zip` contains the complete project payload
available in this staging workspace, including the simulator sources,
geometry/source data, checker, regression tests, and release notes. The
original Windows launcher files were not staged, so they are not included.
