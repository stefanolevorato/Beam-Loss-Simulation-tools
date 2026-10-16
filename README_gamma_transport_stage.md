# Gamma transport stage

`gamma_transport.py` provides the project's effective Gamma interaction stage
alongside the frozen v62 optical simulator. It consumes
`gamma_source_profile_collimated_v1.json` and the existing v62 geometry
objects. Results remain marked
`first_stage_effective_cross_sections_no_geant4`.

The default source is the approved thin, centered, collimated beam along global
`+Z`. Continuum energies are sampled from the profile's relative digitized
bins. Reference lines are inactive unless the user supplies relative weights,
for example:

```sh
PYTHONPATH=. python gamma_transport.py \
  --source-profile gamma_source_profile_collimated_v1.json \
  --line-weight 511=1 --line-weight 1173=0.5 --line-weight 1332=0.5 \
  --events 1000 --json-output gamma_transport.json \
  --plot gamma_transport.png
```

The model includes effective Compton, photoelectric, and above-threshold pair
production branches. Charged secondaries follow straight paths. Their initial
range comes from the legacy, unvalidated range expression, used as one areal
range budget with uniform energy loss per g/cm². At every material boundary,
transport queries the configured volumes again and records a segment and its
deposit in the material traversed. The residual kinetic energy is recorded as
escaped after the path leaves all configured solids. Vacuum gaps between solids
are crossed without energy loss; only configured solid segments appear in the
charged path records.

The current range expression is nonmonotonic between 12.5 and 100 keV. NIST
ESTAR numeric tables were not present in the available search result, and the
offline sandbox cannot retrieve them. The implementation therefore preserves
the formula, makes no reference accuracy claim, and documents the limitation in
`gamma_range_model_diagnostic.md`.

The `TransportVolume` adapter covers the quartz radiator and both placed v62
detector assemblies, including optional guide solids. Cherenkov light is
estimated locally in fused silica above threshold; scintillation light is
estimated locally in plastic-like volumes using the configured
10,000 photons/MeV first-stage yield. Optical-photon transport, fiber
transport, atomic relaxation, multiple-scattering electron transport,
energy-dependent stopping-power integration, quenching, and electronics
response are not included. These results are not a validated general-purpose
radiation-transport calculation.
