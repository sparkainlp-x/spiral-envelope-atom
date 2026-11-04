"""Demo: spiral-envelope atoms on a 2D disk (standard library only).

Prints the three measured results:
  1. How many seeded candidate modes the frozen gate admits.
  2. How envelope decay changes the resolved closure set.
  3. The negative controls (random candidates, random-phase field).

Run:  python3 demo.py
Tests: python3 -m unittest discover -s tests -v
"""

import math
import os
import random
import sys

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))

from spiral_envelope_atoms import modes as M

R = 1.0
TWO_PI = 2.0 * math.pi

SEED_M = [-2, -1, -0.5, 0, 0.5, 1, 2]
SEED_N_FACTORS = [0.5, 1, 1.5, 2, 2.5, 3]
MODE_M = [-2, -1, 0, 1, 2]
MODE_N = [1, 2, 3, 4, 5, 6]
LAMBDAS = [1.0, 0.5, 0.3, 0.2, 0.15, 0.1]


def main():
    print("Spiral-envelope atoms (classical wave toy, SYNTHETIC only)")
    print("Disk radius R = 1, field psi = exp(-r/lam) * exp(i*(m*theta + k*r))")
    print()

    # 1. Seeded closure count.
    cands = [(m, TWO_PI * f / R) for m in SEED_M for f in SEED_N_FACTORS]
    admitted, refused = M.count_closing(cands, R)
    print(f"1. Seeded candidates: {len(admitted)}/{len(cands)} close "
          f"({len(refused)} refused by the frozen gate)")
    print("   Admitted mode numbers (m, n): "
          + ", ".join(str(M.mode_numbers(m, k, R)) for (m, k) in admitted))
    print()

    # 2. Envelope decay vs the resolved closure set.
    print("2. Resolved closure set vs decay length lam (tau = 0.05):")
    counts = []
    for lam in LAMBDAS:
        resolved = M.resolved_modes(MODE_M, MODE_N, R, lam)
        n_set = sorted(set(n for (_, n, _) in resolved))
        counts.append(len(resolved))
        print(f"   lam = {lam:<4} -> {len(resolved):>2} modes resolved, "
              f"radial numbers {n_set}, visible radius "
              f"{M.visible_radius(lam, R):.4f}")
    print(f"   Sweep counts: {', '.join(str(c) for c in counts)}")
    print("   Reading: low-n modes drop out first as lam shrinks;")
    print("   high-n modes survive longer because their tight pitch still")
    print("   shows two full turns inside the visible inner disk.")
    print()

    # 3. Negative controls.
    total_admitted = 0
    for seed in [0, 1, 42]:
        rc = M.random_candidates(500, seed)
        a, _ = M.count_closing(rc, R)
        total_admitted += len(a)
    print(f"3. Negative control, random candidates: "
          f"{total_admitted}/1500 admitted (0/500 per seed)")
    rng = random.Random(7)
    phases = [[rng.uniform(-math.pi, math.pi) for _ in range(32)] for _ in range(16)]
    modes = M.enumerate_modes(MODE_M, MODE_N, R)
    score = M.best_match_score(phases, modes, R, 1.0)
    print(f"   Negative control, random-phase field: best mode match "
          f"{score:.4f} (threshold 0.5, no match)")
    print()
    print("Limits: toy rule, not a physical law. No claim about real atoms,")
    print("real electromagnetism, or quantum mechanics. See README.md.")


if __name__ == "__main__":
    main()
