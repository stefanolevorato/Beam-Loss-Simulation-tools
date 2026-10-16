# Quartz Cherenkov and Scintillator Simulator — v70

**Application version:** `0.38.0`  
**Build ID:** `build-20260928-v70-windows-headless-font-fix`

This is the complete simulator release, carrying forward the v69 Windows
Tkinter fixes and making the Global Assembly renderer safe to exercise with a
headless canvas mock. Real Tk canvases retain the bold Tk font; headless mocks
use the equivalent Tk font tuple without creating a Tk font object.

## Included features

- Cherenkov, Scintillator, and combined optical runs, geometry editing,
  detector response, optical Results, and the interactive event viewer.
- Gamma transport with the supplied collimated source profile and reference
  lines. Gamma runs retain their event visualization and write
  `gamma_transport.png`.
- Complete Global Assembly profiles store both assemblies, sensors, guides,
  readout settings, translations, and rotations. Loading a saved profile in a
  fresh session restores both complete assemblies and positions.
- Full regression checker, Windows launchers and dependency scripts,
  requirements, and the supplied `Assemblycompleto.json` profile.

## Windows launch

Extract the full archive into one folder, then run `Start_GUI.bat`. The
launcher uses that folder's `.venv`, runs the complete checker, and starts the
GUI only if the checks pass. A broken venv is preserved in a timestamped
`.venv.broken-*` folder before recreation. Dependency installation remains
opt-in; see `README_windows_setup.txt`.

Windows scripts were structurally checked but were not executed on Windows.
