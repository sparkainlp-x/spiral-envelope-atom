"""Tests for spiral-envelope-atoms. Standard library only (unittest).

Every headline number in the README and demo is pinned here, so a
change in the model or the measurement rule fails loudly instead of
silently moving a claim.
"""

import cmath
import math
import os
import random
import subprocess
import sys
import unittest

sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.abspath(__file__))))

from spiral_envelope_atoms import modes as M

R = 1.0
TWO_PI = 2.0 * math.pi

# The seeded candidate grid used by the demo and the README.
SEED_M = [-2, -1, -0.5, 0, 0.5, 1, 2]
SEED_N_FACTORS = [0.5, 1, 1.5, 2, 2.5, 3]

# The mode-number ranges used by the demo and the README.
MODE_M = [-2, -1, 0, 1, 2]
MODE_N = [1, 2, 3, 4, 5, 6]
LAMBDAS = [1.0, 0.5, 0.3, 0.2, 0.15, 0.1]
SWEEP_COUNTS = [30, 30, 20, 15, 10, 0]


def seeded_candidates():
    return [(m, TWO_PI * f / R) for m in SEED_M for f in SEED_N_FACTORS]


class TestEnvelopeAndField(unittest.TestCase):
    def test_envelope_at_zero_is_one(self):
        self.assertAlmostEqual(M.envelope(0.0, 1.0), 1.0, places=12)

    def test_envelope_decays_and_matches_formula(self):
        for r in [0.25, 0.5, 1.0]:
            self.assertAlmostEqual(M.envelope(r, 0.5), math.exp(-r / 0.5), places=12)
        self.assertLess(M.envelope(1.0, 0.5), M.envelope(0.5, 0.5))
        self.assertLess(M.envelope(0.5, 0.5), M.envelope(0.0, 0.5))

    def test_envelope_rejects_bad_inputs(self):
        with self.assertRaises(ValueError):
            M.envelope(-0.1, 1.0)
        with self.assertRaises(ValueError):
            M.envelope(0.5, 0.0)

    def test_field_magnitude_equals_envelope(self):
        # |exp(i * phase)| = 1, so |psi| = A(r) at every grid point.
        radii, thetas, values = M.field_on_grid(2, TWO_PI * 3 / R, 0.5, R)
        for i, r in enumerate(radii):
            for j in range(len(thetas)):
                self.assertAlmostEqual(abs(values[i][j]), M.envelope(r, 0.5), places=12)

    def test_grid_shape_and_half_offset_radii(self):
        radii, thetas, values = M.field_on_grid(1, TWO_PI / R, 1.0, R, n_r=8, n_theta=16)
        self.assertEqual(len(radii), 8)
        self.assertEqual(len(thetas), 16)
        self.assertEqual(len(values), 8)
        self.assertEqual(len(values[0]), 16)
        self.assertGreater(radii[0], 0.0)  # r = 0 avoided
        self.assertLessEqual(radii[-1], R)


class TestClosureGate(unittest.TestCase):
    def test_gate_admits_exact_closing_modes(self):
        for m in [-2, -1, 0, 1, 2]:
            for n in [1, 2, 3, 6]:
                self.assertTrue(M.closes(m, TWO_PI * n / R, R), (m, n))

    def test_gate_refuses_half_integer_m(self):
        for m in [-1.5, -0.5, 0.5, 1.5]:
            self.assertFalse(M.closes(m, TWO_PI * 2 / R, R), m)

    def test_gate_refuses_non_integer_radial_factor(self):
        for f in [0.5, 1.5, 2.5]:
            self.assertFalse(M.closes(1, TWO_PI * f / R, R), f)

    def test_gate_refuses_zero_radial_cycles(self):
        # k = 0 has no pitch to fit the disk.
        self.assertFalse(M.closes(1, 0.0, R))

    def test_gate_tolerance_boundary(self):
        # Within tolerance: admitted. Beyond it: refused.
        self.assertTrue(M.closes(1 + 1e-10, TWO_PI * 2 / R, R))
        self.assertFalse(M.closes(1 + 1e-6, TWO_PI * 2 / R, R))
        self.assertTrue(M.closes(1, TWO_PI * (2 + 1e-10) / R, R))
        self.assertFalse(M.closes(1, TWO_PI * (2 + 1e-6) / R, R))

    def test_azimuthal_closure_ratio(self):
        # Integer m returns the field to itself; m = 0.5 flips its sign.
        self.assertAlmostEqual(abs(M.azimuthal_closure_ratio(2) - 1), 0.0, places=12)
        self.assertAlmostEqual(abs(M.azimuthal_closure_ratio(0.5) + 1), 0.0, places=12)

    def test_radial_phase_total_is_integer_cycles(self):
        for n in [1, 2, 3]:
            k = TWO_PI * n / R
            self.assertAlmostEqual(M.radial_cycles(k, R), n, places=12)

    def test_mode_numbers_and_quantized_radii(self):
        self.assertEqual(M.mode_numbers(2, TWO_PI * 3 / R, R), (2, 3))
        self.assertIsNone(M.mode_numbers(0.5, TWO_PI * 3 / R, R))
        rings = M.quantized_radii(TWO_PI * 3 / R, R)
        self.assertEqual(len(rings), 3)
        self.assertAlmostEqual(rings[0], 1.0 / 3.0, places=12)
        self.assertAlmostEqual(rings[1], 2.0 / 3.0, places=12)
        self.assertAlmostEqual(rings[2], 1.0, places=12)
        self.assertIsNone(M.quantized_radii(TWO_PI * 1.5 / R, R))

    def test_pitch_formula(self):
        self.assertAlmostEqual(M.pitch(TWO_PI * 3 / R), R / 3.0, places=12)
        with self.assertRaises(ValueError):
            M.pitch(0.0)


class TestSeededClosureCount(unittest.TestCase):
    def test_seeded_grid_15_of_42_close(self):
        admitted, refused = M.count_closing(seeded_candidates(), R)
        self.assertEqual(len(seeded_candidates()), 42)
        self.assertEqual(len(admitted), 15)
        self.assertEqual(len(refused), 27)

    def test_admitted_are_exactly_integer_pairs(self):
        admitted, _ = M.count_closing(seeded_candidates(), R)
        for (m, k) in admitted:
            numbers = M.mode_numbers(m, k, R)
            self.assertIsNotNone(numbers)
            self.assertIn(numbers[0], [-2, -1, 0, 1, 2])
            self.assertIn(numbers[1], [1, 2, 3])

    def test_gate_does_not_depend_on_decay_length(self):
        # The gate is pure phase arithmetic; lam never enters closes().
        cands = seeded_candidates()
        first = M.count_closing(cands, R)[0]
        for lam in [0.1, 1.0, 10.0]:
            # Rebuild the same fields at each lam; the gate verdicts
            # are unchanged because closes() takes no lam argument.
            self.assertEqual(M.count_closing(cands, R)[0], first)
            for (m, k) in first:
                M.field(0.5, 0.3, m, k, lam)  # field itself is well defined


class TestResolution(unittest.TestCase):
    def test_mode_enumeration_count(self):
        modes = M.enumerate_modes(MODE_M, MODE_N, R)
        self.assertEqual(len(modes), 30)

    def test_visible_radius_formula_and_cap(self):
        self.assertAlmostEqual(M.visible_radius(1.0, R), 1.0, places=12)  # capped
        self.assertAlmostEqual(M.visible_radius(0.5, R), 1.0, places=12)  # capped
        self.assertAlmostEqual(M.visible_radius(0.2, R), 0.2 * math.log(20), places=12)
        with self.assertRaises(ValueError):
            M.visible_radius(1.0, R, tau=1.0)

    def test_visible_pitches_cases(self):
        # lam = 0.3 -> r_max ~ 0.8987: n=2 shows 1 pitch, n=3 shows 2.
        self.assertEqual(M.visible_pitches(2, 0.3, R), 1)
        self.assertEqual(M.visible_pitches(3, 0.3, R), 2)
        # lam = 0.1 -> r_max ~ 0.2996: even n=6 shows only 1 pitch.
        self.assertEqual(M.visible_pitches(6, 0.1, R), 1)

    def test_lambda_sweep_counts_pinned(self):
        counts = [len(M.resolved_modes(MODE_M, MODE_N, R, lam)) for lam in LAMBDAS]
        self.assertEqual(counts, SWEEP_COUNTS)

    def test_low_n_modes_drop_first(self):
        # At lam = 0.3 the surviving radial numbers are 3..6 only.
        resolved = M.resolved_modes(MODE_M, MODE_N, R, 0.3)
        self.assertEqual(sorted(set(n for (_, n, _) in resolved)), [3, 4, 5, 6])
        # At lam = 0.15 only 5 and 6 survive.
        resolved = M.resolved_modes(MODE_M, MODE_N, R, 0.15)
        self.assertEqual(sorted(set(n for (_, n, _) in resolved)), [5, 6])

    def test_resolved_sets_shrink_monotonically(self):
        previous = None
        for lam in LAMBDAS:
            current = set((m, n) for (m, n, _) in M.resolved_modes(MODE_M, MODE_N, R, lam))
            if previous is not None:
                self.assertTrue(current <= previous, lam)
            previous = current


class TestNegativeControls(unittest.TestCase):
    def test_random_candidates_all_refused(self):
        for seed in [0, 1, 42]:
            cands = M.random_candidates(500, seed)
            admitted, _ = M.count_closing(cands, R)
            self.assertEqual(len(admitted), 0, seed)

    def test_random_phase_field_matches_no_mode(self):
        rng = random.Random(7)
        phases = [[rng.uniform(-math.pi, math.pi) for _ in range(32)] for _ in range(16)]
        modes = M.enumerate_modes(MODE_M, MODE_N, R)
        score = M.best_match_score(phases, modes, R, 1.0)
        # Far below any match threshold (0.5); measured value is ~0.07.
        self.assertLess(score, 0.5)

    def test_demo_runs_and_prints_headlines(self):
        here = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
        out = subprocess.run(
            [sys.executable, os.path.join(here, "demo.py")],
            capture_output=True, text=True, cwd=here, timeout=120,
        )
        self.assertEqual(out.returncode, 0, out.stderr)
        self.assertIn("15/42", out.stdout)
        self.assertIn("30, 30, 20, 15, 10, 0", out.stdout)
        self.assertIn("0/500", out.stdout)


if __name__ == "__main__":
    unittest.main(verbosity=2)
