# Quartz Cherenkov Simulator v75

Application version: `0.43.0`  
Build ID: `build-20260929-v75-local-mirror-label`

v75 makes terminal mirror labels readable in both interactive complete local
3D viewers. Each active mirror has a red legend swatch and a clear `Mirror at +X`
or `Mirror at -X` label. The local views use the mirror state and position
from the canonical complete assembly scene; the plate geometry, camera
controls, and Global Assembly viewer are unchanged.

The release carries forward the v74 Cherenkov terminal roundtrip fix and the
earlier simulator, Gamma transport, persistence, and geometry editor work.

## Verify

The headless release checker includes the v75 label checks and the existing
regression suite:

```text
python check_environment.py
```

It exercises the canvas drawing code with a recording canvas; it does not open
a desktop GUI window. The staged v74 file set is missing `Start_GUI.bat`, which
the existing package inventory requires, so running the checker as a whole
stops at that inventory check. The remaining Python regression functions,
including both v75 local viewer paths, passed when run directly. Windows
launcher checks and desktop GUI visual verification were unavailable because
the launcher files and a display session are not present in this staged set.

## Files

`quartz_cherenkov_simulator_v75.zip` packages all files staged for v74 plus
this README and the v75 implementation summary.
