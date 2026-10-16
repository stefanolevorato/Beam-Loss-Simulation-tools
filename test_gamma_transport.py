"""Targeted first-phase regressions for Gamma transport bookkeeping.

The range diagnostic intentionally records the current legacy formula's
nonmonotonic behavior. It is not a reference-point validation: no authoritative
electron/positron range table is staged with this project.
"""

from __future__ import annotations

import unittest
from unittest.mock import patch

import numpy as np

import gamma_transport as gamma


class FixedSource:
    source_direction = (0.0, 0.0, 1.0)
    source_profile_path = ""
    lines = ()

    def __init__(self, energy_keV: float):
        self.energy_keV = float(energy_keV)

    def sample_energy_keV(self, _rng):
        return self.energy_keV, "test"


class TestSolid:
    bounds_min_m = np.array((-1.0, -1.0, 0.0))

    def __init__(self, charged_path_m: float = 1e9):
        self.charged_path_m = float(charged_path_m)

    def ray_intersection(self, _origin, _direction):
        return self.charged_path_m, "test-face", (0.0, 0.0, 1.0)

    def ray_interval(self, _origin, _direction):
        return 0.0, self.charged_path_m


class SlabSolid:
    """Infinite X slab used to make exact adjacent-volume ray intervals."""

    bounds_min_m = np.array((-1.0, -1.0, -1.0))

    def __init__(self, x_min: float, x_max: float):
        self.x_min = float(x_min)
        self.x_max = float(x_max)

    def ray_interval(self, origin, direction):
        origin_x = float(origin[0])
        direction_x = float(direction[0])
        if abs(direction_x) < 1e-15:
            if self.x_min <= origin_x <= self.x_max:
                return -1e9, 1e9
            return 1e9, 1e9
        first = (self.x_min - origin_x) / direction_x
        second = (self.x_max - origin_x) / direction_x
        return min(first, second), max(first, second)


def _volume(charged_path_m: float = 1e9) -> gamma.TransportVolume:
    return gamma.TransportVolume(
        "test volume", "test assembly", "test component", TestSolid(charged_path_m), "custom"
    )


def _event(energy_keV: float, process, max_interactions: int = 100, charged_path_m: float = 1e9):
    volume = _volume(charged_path_m)
    source = FixedSource(energy_keV)
    def next_volume(_volumes, position, _direction):
        if float(position[2]) < 0.0:
            return volume, 0.0, 1e8  # Keep the primary gamma ray in the deterministic test volume.
        return volume, 0.0, volume.solid.charged_path_m

    with patch.object(gamma, "_next_volume", side_effect=next_volume), patch.object(
        gamma, "mass_attenuation_coefficients", side_effect=process
    ):
        return gamma.simulate_event(
            source, [volume], 0, np.random.default_rng(24680), max_interactions=max_interactions
        )


def _closure_energy_keV(event: gamma.GammaEvent) -> float:
    # Deposits summarize charged-track deposits. Count charged tracks once.
    charged = sum(track.deposited_energy_keV + track.escaped_energy_keV for track in event.charged_tracks)
    return charged + event.escaped_gamma_energy_keV + event.terminated_gamma_energy_keV


class GammaTransportTargetedTests(unittest.TestCase):
    def test_known_range_model_discrepancy_is_visible_pending_reference(self):
        """Diagnostic only: not a pass claim for physical range accuracy."""
        range_12_5_kev = gamma.electron_range_g_cm2(12.5)
        range_100_kev = gamma.electron_range_g_cm2(100.0)
        self.assertGreater(
            range_12_5_kev,
            range_100_kev,
            "legacy range behavior changed; review against an authoritative reference before updating this diagnostic",
        )
        print(
            "KNOWN RANGE MODEL DISCREPANCY (unvalidated): "
            f"R(12.5 keV)={range_12_5_kev:.9g} g/cm^2 > "
            f"R(100 keV)={range_100_kev:.9g} g/cm^2; no reference-point assertion is made."
        )

    def test_event_energy_closure_for_photoelectric_compton_and_pair(self):
        def photoelectric(_energy, _material):
            return {"compton": 0.0, "photoelectric": 1.0, "pair": 0.0}

        def compton(_energy, _material):
            return {"compton": 1.0, "photoelectric": 0.0, "pair": 0.0}

        def pair_then_photoelectric(energy, _material):
            if energy >= gamma.PAIR_THRESHOLD_KEV:
                return {"compton": 0.0, "photoelectric": 0.0, "pair": 1.0}
            return photoelectric(energy, _material)

        cases = (
            (2000.0, photoelectric, 10, "photoelectric"),
            (1200.0, compton, 4, "compton"),
            (2000.0, pair_then_photoelectric, 10, "pair"),
        )
        for energy_keV, process, limit, expected_process in cases:
            with self.subTest(interaction=expected_process):
                event = _event(energy_keV, process, max_interactions=limit)
                self.assertGreater(event.interaction_counts[expected_process], 0)
                self.assertAlmostEqual(_closure_energy_keV(event), energy_keV, delta=1e-8)
                self.assertAlmostEqual(
                    sum(deposit.energy_keV for deposit in event.deposits),
                    sum(track.deposited_energy_keV for track in event.charged_tracks),
                    delta=1e-8,
                )

    def test_event_energy_closure_includes_charged_and_gamma_escape(self):
        def photoelectric(_energy, _material):
            return {"compton": 0.0, "photoelectric": 1.0, "pair": 0.0}

        event = _event(700.0, photoelectric, charged_path_m=0.0)
        self.assertAlmostEqual(_closure_energy_keV(event), 700.0, delta=1e-8)
        self.assertAlmostEqual(event.charged_tracks[0].escaped_energy_keV, 700.0, delta=1e-8)

        volume = _volume()
        with patch.object(gamma, "_next_volume", return_value=None):
            escaped = gamma.simulate_event(FixedSource(700.0), [volume], 1, np.random.default_rng(1))
        self.assertAlmostEqual(_closure_energy_keV(escaped), 700.0, delta=1e-8)
        self.assertAlmostEqual(escaped.escaped_gamma_energy_keV, 700.0, delta=1e-8)

    def test_subthreshold_compton_gamma_is_counted_once_as_terminated(self):
        def compton(_energy, _material):
            return {"compton": 1.0, "photoelectric": 0.0, "pair": 0.0}

        event = _event(1.1, compton, max_interactions=100)
        self.assertGreater(event.terminated_gamma_energy_keV, 0.0)
        self.assertAlmostEqual(_closure_energy_keV(event), 1.1, delta=1e-8)

    def test_gamma_parent_ids_reference_the_producing_gamma(self):
        calls = 0

        def root_compton_then_child_photoelectric(_energy, _material):
            nonlocal calls
            process = "compton" if calls < 2 else "photoelectric"
            calls += 1
            return {
                "compton": float(process == "compton"),
                "photoelectric": float(process == "photoelectric"),
                "pair": 0.0,
            }

        event = _event(1200.0, root_compton_then_child_photoelectric, max_interactions=5)
        self.assertEqual([segment.track_id for segment in event.gamma_tracks], [0, 1])
        self.assertEqual([segment.parent_track_id for segment in event.gamma_tracks], [None, 0])
        self.assertEqual([track.track_id for track in event.charged_tracks], [1, 2])
        self.assertEqual([track.parent_track_id for track in event.charged_tracks], [0, 1])

    def test_charged_track_crosses_adjoining_volumes_and_segments_deposition(self):
        first = gamma.TransportVolume(
            "first material", "assembly_a", "first", SlabSolid(0.0, 0.0002), "custom"
        )
        second = gamma.TransportVolume(
            "second material", "assembly_b", "second", SlabSolid(0.0002, 0.0008), "fused_silica"
        )
        event = gamma.GammaEvent(0, 100.0, "test")
        with patch.object(gamma, "electron_range_g_cm2", return_value=0.05):
            track, next_id = gamma._charged_track(
                event, first, np.array((0.0001, 0.0, 0.0)), np.array((1.0, 0.0, 0.0)),
                100.0, "electron", 7, 3, 100.0, np.random.default_rng(42), 0,
                (first, second), "photoelectric",
            )

        self.assertEqual(next_id, 4)
        self.assertEqual([segment.volume for segment in track.segments], ["first material", "second material"])
        self.assertEqual([deposit.volume for deposit in event.deposits], ["first material", "second material"])
        self.assertTrue(all(not deposit.interaction_site for deposit in event.deposits))
        self.assertAlmostEqual(track.segments[0].path_mm, 0.1, places=6)
        self.assertAlmostEqual(track.segments[1].path_mm, 0.04 / 2.2 * 10.0, places=6)
        self.assertAlmostEqual(track.path_mm, sum(segment.path_mm for segment in track.segments), places=10)
        self.assertAlmostEqual(sum(segment.deposited_energy_keV for segment in track.segments), track.deposited_energy_keV, places=10)
        self.assertAlmostEqual(sum(deposit.energy_keV for deposit in event.deposits), track.deposited_energy_keV, places=10)
        self.assertAlmostEqual(track.deposited_energy_keV + track.escaped_energy_keV, 100.0, places=10)
        self.assertAlmostEqual(track.escaped_energy_keV, 0.0, places=10)
        self.assertGreater(track.segments[1].start_position_m[0], track.segments[0].start_position_m[0])
        self.assertAlmostEqual(_closure_energy_keV(event), event.source_energy_keV, places=10)
        self.assertAlmostEqual(event.energy_accounting()["balance_residual_keV"], 0.0, places=10)

    def test_charged_track_reports_residual_energy_after_leaving_configured_volumes(self):
        first = gamma.TransportVolume(
            "short solid", "assembly_a", "short", SlabSolid(0.0, 0.0002), "custom"
        )
        event = gamma.GammaEvent(0, 100.0, "test")
        with patch.object(gamma, "electron_range_g_cm2", return_value=0.05):
            track, _ = gamma._charged_track(
                event, first, np.array((0.0001, 0.0, 0.0)), np.array((1.0, 0.0, 0.0)),
                100.0, "positron", 2, 5, 100.0, np.random.default_rng(17), 0,
                (first,), "photoelectric",
            )

        self.assertAlmostEqual(track.deposited_energy_keV, 20.0, places=9)
        self.assertAlmostEqual(track.escaped_energy_keV, 80.0, places=9)
        self.assertAlmostEqual(track.deposited_energy_keV + track.escaped_energy_keV, track.energy_keV, places=10)
        self.assertEqual(len(track.segments), 1)

    def test_charged_segment_summary_separates_energy_loss_from_interaction_counts(self):
        first = gamma.TransportVolume(
            "first material", "assembly_a", "first", SlabSolid(0.0, 0.0002), "custom"
        )
        second = gamma.TransportVolume(
            "second material", "assembly_b", "second", SlabSolid(0.0002, 0.0008), "fused_silica"
        )
        event = gamma.GammaEvent(0, 100.0, "test")
        event.interaction_counts["photoelectric"] = 1
        with patch.object(gamma, "electron_range_g_cm2", return_value=0.05):
            track, _ = gamma._charged_track(
                event, first, np.array((0.0001, 0.0, 0.0)), np.array((1.0, 0.0, 0.0)),
                100.0, "electron", 7, 3, 100.0, np.random.default_rng(42), 0,
                (first, second), "photoelectric",
            )
        event.deposits.append(gamma.EnergyDeposit(
            0, first.name, first.assembly, first.component, first.material, "photoelectric",
            (0.0001, 0.0, 0.0), 0.0, 100.0, 7, (track.track_id,), None, True,
        ))
        result = gamma.GammaTransportResult(FixedSource(100.0), [event], volumes=(first, second))
        summary = result.summary()
        self.assertEqual(summary["interaction_counts"], {"photoelectric": 1})
        self.assertAlmostEqual(summary["by_volume"][first.name]["deposited_energy_keV"], 20.0, places=9)
        self.assertAlmostEqual(summary["by_volume"][second.name]["deposited_energy_keV"], 80.0, places=9)
        self.assertAlmostEqual(summary["energy_accounting"]["source_energy_keV"], 100.0, places=10)
        self.assertAlmostEqual(summary["energy_accounting"]["balance_residual_keV"], 0.0, places=10)
        self.assertAlmostEqual(result.to_dict()["events"][0]["energy_accounting"]["balance_residual_keV"], 0.0, places=10)


if __name__ == "__main__":
    unittest.main(verbosity=2)
