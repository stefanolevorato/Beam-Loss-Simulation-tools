# Quartz Cherenkov and Scintillator Simulator — v69

**Application version:** `0.37.0`  
**Build ID:** `build-20260928-v69-windows-tk-live-viewer-fix`

This is the complete simulator release, carrying forward v68 and fixing the
remaining Windows Tkinter recursion path during a headless assembly load and
live-viewer refresh.

## Included features

- Cherenkov, Scintillator, and combined optical runs, geometry editing,
  detector response, optical Results, and the interactive event viewer.
- Gamma transport with the supplied collimated source profile and reference
  lines. Gamma runs retain their event visualization and write
  `gamma_transport.png`. Gamma local light values remain estimates; this
  stage does not transport optical photons to a sensor.
- Complete Global Assembly profiles store both assemblies, their sensors,
  guides, readout settings, translations, and rotations. Loading a saved
  profile in a fresh session restores both complete assemblies and positions.
- The required `check_environment.py` regression suite, Windows launchers,
  dependency list, and setup documentation.

`Assemblycompleto.json` is an example complete profile. Use **Load complete
assembly profile** in the Global Assembly tab to load it.

## Windows launch

Extract the full archive into one folder, then run `Start_GUI.bat`. The
launcher uses that folder's `.venv`, runs the complete checker, and starts the
GUI only if the checks pass. If the existing venv's Python cannot start, the
launcher preserves it in a timestamped `.venv.broken-*` folder and creates a
fresh project environment with `py -3`.

Dependency installation is opt-in. If needed, run this in PowerShell from the
extracted folder:

```powershell
.\setup_and_start.ps1 -InstallDependencies
```

`numpy` and `matplotlib` are required. `uproot` and `awkward` enable optional
ROOT file output. Tkinter is supplied by the Windows Python installation. See
`README_windows_setup.txt` for setup details. Windows scripts were reviewed
statically but were not executed on Windows.

## Verification

From a Python environment with dependencies available:

```powershell
python -m py_compile check_environment.py quartz_cherenkov_gui.py quartz_cherenkov_sim.py quartz_geometry.py quartz_geometry_editor.py gamma_transport.py
python check_environment.py
```

The live-viewer regression uses an uninitialized `tk.Tk` object, leaves
`_assembly_editor_state` absent, then loads geometry through the real callback
and refreshes the open viewer. GUI instance reads and presence checks inspect
`__dict__` directly, so a missing Tcl interpreter cannot trigger Tkinter's
`__getattr__` fallback.
