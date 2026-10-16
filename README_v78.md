# Quartz Cherenkov Simulator v78

Application version: `0.46.0`  
Build ID: `build-20260929-v78-gamma-charged-boundaries`

## Gamma charged-particle transport

Electron and positron tracks now query the configured transport volumes as they
cross material boundaries. Each in-solid segment records its start/end points,
material, path length, kinetic energy change, and deposited energy. Deposits
are assigned to the segment's actual volume. The Gamma Results global viewer
renders these material-specific track segments and keeps its older aggregate
track format readable.

Transport uses the current legacy electron range expression as a single
initial areal-range budget. Energy loss is apportioned uniformly per g/cm², so
density changes alter physical path length without resetting the particle's
range. Vacuum gaps between configured solids are crossed without loss. If no
configured solid remains ahead, residual kinetic energy is recorded as escaped.
The result summary reports deposited, escaped charged, escaped gamma, and
sub-threshold terminated gamma energy, with a balance residual. Interaction
counts remain counts of Gamma interactions and do not count the material loss
segments as separate interactions. The same energy-accounting fields are saved
for each individual event.

This is a boundary and energy-accounting improvement within the existing
straight-line model. It does not add multiple scattering or a validated
stopping-power curve. The legacy range expression remains unvalidated and has
a documented nonmonotonicity: its output at 12.5 keV exceeds the output at
100 keV. The available NIST ESTAR search result contained no numeric tables,
and the offline sandbox cannot retrieve them; no range formula change or
reference accuracy claim is made. See
[`gamma_range_model_diagnostic.md`](gamma_range_model_diagnostic.md).

## Preserved behavior and baseline

The source package consolidates the available v76 Gamma transport and tests
with the recovered v77 GUI source. The v77 detached global 3D Gamma viewer is
present in both staged source and the v77 archive. The v76 charged-track parent
IDs and photoelectric energy double-counting fix are retained. The existing
Gamma event table, inline viewer, saved geometry snapshots, and viewer
selection synchronization remain covered by the headless checker.

## Verification

Run:

```text
python -m unittest -v test_gamma_transport
python check_environment.py
```

The focused suite covers the known range limitation, adjacent-material
crossing, per-segment deposition/path data, charged escape, energy closure,
Gamma ancestry, and the v76 interaction accounting fixes. The full checker
also compiles the Python sources and runs the existing simulator and viewer
regressions. No desktop session was available for visual GUI QA. The original
Windows launcher files were not staged, so launcher checks are reported as
skipped and those launchers are absent from the package.

## Package

`quartz_cherenkov_simulator_v78.zip` contains the changed transport, GUI,
focused tests, checker, release notes, and the complete available staged
project payload. It does not contain the unavailable original Windows
launchers.
