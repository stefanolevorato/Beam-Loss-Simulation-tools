# Quartz Cherenkov / Scintillator two-assembly simulator — v63

Build ID: `build-20260925-v63-gamma-gui-integration`

This package preserves the v62 optical GUI and adds an explicit `Gamma
transport` Run mode. Gamma mode uses the current detached GUI geometry
snapshot: the canonical Cherenkov/Scintillator assembly sections for shape,
size, and material, plus the two placed detector assemblies from the current
beam-monitor controls for global positions and separation. It does not
substitute the prototype's arbitrary default geometry.

The approved `gamma_source_profile_collimated_v1.json` is included beside the
sources. Its relative continuum is used unchanged. The initial source is a
thin, centered, collimated beam along global `+Z`. The six reference lines
(511, 835, 1173, 1275, 1332, 1461 keV) default to zero/unset and enter the
sampling profile only when explicit relative weights are typed in Run. The
example weights used during development are not defaults.

Select `Gamma transport` in Run, choose the profile if needed, enter Events,
Seed, and any line weights, then press Run Simulation. Results are written to
`gamma_results.json` and shown in the separate `Gamma Results` tab. The tab
reports source energy/component, active line weights, interaction processes,
energy deposits by assembly/volume, charged secondaries, and expected local
Cherenkov/scintillation light.

Important semantic boundary: those light values are local expected production
only. This version does not transport optical photons to a sensor, does not
transport fibers, and does not report gamma-mode arrivals or detections. It is
not a fully radiation-realistic or Geant4-equivalent application. The status
is explicitly marked `first_stage_effective_cross_sections_no_geant4`.

The v62 optical modes, optical Results tab, geometry editors, profiles, and
global assembly view remain available. `Start_GUI.bat` launches the same GUI
entry point on Windows; Windows runtime testing was not available here.

Gamma currently guards non-zero assembly rotations because the inherited v62
two-detector placement contract does not apply those rotations to its placed
solids. Set both gamma assembly rotations to `None`/`0 deg` for a run; a
non-zero value is rejected explicitly rather than silently ignored.
