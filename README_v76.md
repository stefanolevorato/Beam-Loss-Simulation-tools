# Quartz Cherenkov Simulator v76

Application version: `0.44.0`  
Build ID: `build-20260929-v76-gamma-ancestry-energy-audit`

v76 is a first-phase Gamma transport correction based on the staged v75
workspace. Gamma-created charged tracks now point to their producing gamma
track ID. Event energy accounting no longer counts photoelectric absorption a
second time as terminated gamma energy. Focused regressions cover ancestry,
energy closure for photoelectric, Compton, and pair outcomes, and gamma and
charged-particle escape, including Compton gamma termination below the
transport threshold.

The electron/positron range formula is unchanged. Its outputs decrease from
0.419657799 g/cm² at 12.5 keV to 0.242048532 g/cm² at 100 keV. No authoritative
range reference is present in the staged files, and network access was
unavailable. The diagnostic and required reference are documented in
`gamma_range_model_diagnostic.md`; no reference-point tests are claimed.

## Verify

Run focused Gamma regressions:

```text
python -m unittest -v test_gamma_transport
```

Run the staged-source checker:

```text
python check_environment.py
```

The checker runs syntax compilation, Gamma regressions, and the practical
headless simulator regressions. It reports package and Windows-launcher checks
as skipped when required original files are absent. No desktop display was
available for live GUI verification.

In this staging run, syntax compilation and all five focused Gamma tests
passed, and `python check_environment.py` completed successfully. The checker
skipped the absent launcher/package originals above. Its optional ROOT writer
self-test was also skipped because `uproot` and `awkward` are not installed.

## Package status

No v76 project ZIP is included because the staged workspace lacks the original
Windows launcher files `Start_GUI.bat`, `setup_and_start.ps1`,
`Check_Dependencies.bat`, and `check_dependencies.ps1`. The checker's required
`requirements-root.txt` is also absent; the workspace has a differently named
`requirementsroot.txt`, whose equivalence was not assumed. The existing v75
archive is an earlier baseline and does not contain these v76 changes.

Geometry and UI behavior were not changed; the GUI version/build strings were
incremented for this release.
