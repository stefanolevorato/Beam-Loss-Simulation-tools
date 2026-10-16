# Gamma charged-particle range model diagnostic

Status: **unvalidated discrepancy; the physics formula is unchanged through v78.**

The staged model applies this expression to both electrons and positrons, with
energy in MeV and output in g/cm²:

```text
E = max(kinetic_energy_keV / 1000, 1e-6)
R = max(1e-5, 0.412 * E^1.265 - 0.0954 * ln(max(E, 0.001)))
```

The project files do not include an authoritative electron or positron range
table, a derivation for this exact expression, or its material and energy
domain. The available NIST ESTAR search result identifies the database landing
page but contains no numeric stopping-power or range rows. The offline sandbox
cannot retrieve the underlying ESTAR tables, so there are no ESTAR reference
points to validate here. The earlier attribution in a code comment was not
enough to justify changing the model.

Since v78, the expression supplies one initial areal-range budget for a
straight track. That budget is apportioned across configured material segments
with constant energy loss per g/cm². This makes boundary bookkeeping explicit;
it does not validate the legacy range or provide an energy-dependent stopping
power curve.

## Observed discrepancy

The expression is nonmonotonic in the project's source energy band. Its
current outputs include:

| Kinetic energy | Current output |
| ---: | ---: |
| 12.5 keV | 0.419657799 g/cm² |
| 100 keV | 0.242048532 g/cm² |
| 300 keV | 0.204696081 g/cm² |
| 1 MeV | 0.412 g/cm² |

In particular, increasing the initial energy from 12.5 to 100 keV decreases
the returned CSDA range. These values are outputs of the legacy expression,
**not reference values**. Focused tests preserve this known limitation and make
no reference-point accuracy claim.

## Reference needed before a physics change

Provide or identify a traceable electron and positron CSDA reference over
12.5 keV–2.5 MeV for the modeled materials (fused silica and the selected
plastic/PMMA surrogate), including source/version, units, and any density or
composition assumptions. The current transport function is shared by both
particle species, so the reference must also establish whether that is valid.
Once available, add table-backed point tests before replacing the expression.
