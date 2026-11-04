# spiral-envelope-atoms

A classical wave toy on a 2D disk. Synthetic only.

## The model

A mode is a complex scalar field on a disk of radius `R`:

```
psi(r, theta) = exp(-r / lam) * exp(i * (m * theta + k * r))
```

a decaying radial envelope times a helical spiral phase. The word
"electromagnetic" in the package title is an analogy label only: this
field is **not** a solution of Maxwell's equations, and nothing here
describes a real atom, real electromagnetism, or quantum mechanics.

## The closure condition (a chosen toy rule, not a physical law)

A candidate `(m, k)` is admitted by the frozen gate only if both hold:

1. **Azimuthal closure:** `m` is an integer, so going once around the
   disk (`theta -> theta + 2*pi`) returns the field to itself. The
   spiral closes around the circumference.
2. **Radial fit:** `k * R / (2*pi)` is a positive integer `n`, so the
   spiral pitch `p = 2*pi / k` fits the disk radius an integer number
   of times: `R = n * p`.

An admitted mode is labelled by its mode numbers `(m, n)`. Its
quantized radii are the ring radii `r_j = j * p` for `j = 1..n`, the
radii at which successive spiral turns complete. The gate is frozen:
exact integer arithmetic within a tolerance of `1e-9`. It never learns
and never admits a near miss. Envelope decay does **not** change the
gate — the same candidates close at every decay length.

## Resolution: where envelope decay matters

Whether a closing mode is *resolved* is a separate, stated measurement
rule. With an amplitude floor `tau = 0.05`, the field is only visible
out to `r_max = min(R, lam * ln(1/tau))`. The number of pitches visible
inside that region for mode `n` is `V_n = floor(n * r_max / R)`, and a
mode is resolved when at least two pitches are visible (or all of them,
when the mode has fewer than two).

Consequence: as `lam` shrinks, **low-`n` modes drop out first**,
because their wide pitch no longer shows two full turns inside the
visible inner disk, while high-`n` modes survive longer because their
tight pitch still packs two turns into the visible region. This is a
property of the stated measurement rule, not a physical law.

## Measured results (all pinned by tests)

Seeded candidates: `m` in `{-2, -1, -0.5, 0, 0.5, 1, 2}` crossed with
radial factors `{0.5, 1, 1.5, 2, 2.5, 3}` — **15/42 close**, 27 refused.

Closing modes for `m` in `{-2..2}`, `n` in `{1..6}`: **30 modes**.

Resolved closure set vs decay length (`tau = 0.05`):

| lam | resolved modes | surviving radial numbers |
|-----|---------------:|--------------------------|
| 1.0 | 30 | 1, 2, 3, 4, 5, 6 |
| 0.5 | 30 | 1, 2, 3, 4, 5, 6 |
| 0.3 | 20 | 3, 4, 5, 6 |
| 0.2 | 15 | 4, 5, 6 |
| 0.15 | 10 | 5, 6 |
| 0.1 | 0 | — |

Negative controls: **0/500** random continuous candidates admitted per
seed (seeds 0, 1, 42 — 0/1500 total), and a random-phase field's best
match against any closing mode is **0.07**, far below the 0.5 match
threshold.

## Run it

```sh
python3 demo.py
python3 -m unittest discover -s tests -v
```

Standard library only. No dependencies.

## Honest limits

- Classical wave-envelope analogy toy. Synthetic only.
- No claim about real atoms, real electromagnetism as physics, or
  quantum mechanics. The closure condition is a toy rule chosen for
  this package, not a physical law.
- The resolution rule (two visible pitches, amplitude floor) is a
  measurement convention of this toy, not a detector model.
- The field is evaluated on a polar grid that avoids `r = 0` (the
  phase `m * theta` is undefined there); results are grid-checked,
  not continuum-proven.
- No data from any real device, experiment, or observation is used.

## Licence and authorship

Licence: AGPL-3.0-only (see `LICENSE`).
Author: Jean-François Brisson (Spark AI NLP), sole author. No other
author's code is in this package.
