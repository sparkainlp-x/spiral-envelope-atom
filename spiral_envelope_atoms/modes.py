"""Spiral-envelope atoms: a classical wave toy on a 2D disk.

Model (a chosen toy rule, not a physical law)
--------------------------------------------
A mode is a complex scalar field on a disk of radius ``R``::

    psi(r, theta) = A(r) * exp(i * (m * theta + k * r))

with a decaying radial envelope ``A(r) = exp(-r / lam)`` and a helical
spiral phase ``m * theta + k * r``.  The word "electromagnetic" in the
package title is an analogy label only: this field is NOT a solution of
Maxwell's equations, and nothing here describes a real atom, real
electromagnetism, or quantum mechanics.

Closure condition (the toy rule that selects "atoms")
-----------------------------------------------------
A candidate ``(m, k)`` is admitted by the frozen gate only if BOTH hold:

1. Azimuthal closure: ``m`` is an integer.  Going once around the disk
   (``theta -> theta + 2*pi``) then returns the field to itself, so the
   spiral closes around the circumference.
2. Radial fit: ``k * R / (2*pi)`` is a positive integer ``n``.  The
   spiral pitch ``p = 2*pi / k`` then fits the disk radius an integer
   number of times: ``R = n * p``.

An admitted mode is labelled by its mode numbers ``(m, n)``.  Its
quantized radii are the ring radii ``r_j = j * p`` for ``j = 1..n``:
the radii at which successive spiral turns complete.

The gate is frozen: it is exact integer arithmetic within a tolerance,
it never learns, and it never admits a near miss.  Envelope decay does
NOT change the gate: the same candidates close at every decay length.

Resolution (a separate measurement layer)
------------------------------------------
Whether a closing mode is *resolved* is a different, stated measurement
rule, and this is where envelope decay matters.  With an amplitude floor
``tau``, the field is only visible out to ``r_max = min(R, lam * ln(1/tau))``.
The number of pitches visible inside that region for mode ``n`` is
``V_n = floor(n * r_max / R)``.  A mode is resolved when at least two
pitches are visible, or all of its pitches when it has fewer than two::

    resolved  <=>  V_n >= min(2, n)

Consequence, stated plainly: as ``lam`` shrinks, low-``n`` modes drop out
first, because their wide pitch no longer shows two full turns inside the
visible inner disk, while high-``n`` modes survive longer because their
tight pitch still packs two turns into the visible region.  This is a
property of the stated measurement rule, not a physical law.

Only the Python standard library is used.
"""

import cmath
import math
import random

TWO_PI = 2.0 * math.pi

# Frozen-gate tolerance: a candidate must hit an integer within this
# distance.  1e-9 is far above float64 rounding for the values used
# here (order 1e-15) and far below any genuine near miss (order 1e-2).
GATE_TOL = 1e-9

# Default amplitude floor for the resolution rule.
DEFAULT_TAU = 0.05

# Minimum number of visible pitches required to resolve a mode.
MIN_VISIBLE_PITCHES = 2


def envelope(r, lam):
    """Decaying radial envelope ``exp(-r / lam)``.

    ``r`` must be >= 0 and ``lam`` must be > 0.
    """
    if r < 0:
        raise ValueError("radius r must be >= 0")
    if lam <= 0:
        raise ValueError("decay length lam must be > 0")
    return math.exp(-r / lam)


def pitch(k):
    """Spiral pitch ``p = 2*pi / k`` (radial distance per full turn)."""
    if k == 0:
        raise ValueError("k = 0 has no finite pitch (no radial cycle)")
    return TWO_PI / k


def radial_cycles(k, radius):
    """Number of radial phase cycles across the disk: ``k * R / (2*pi)``."""
    return k * radius / TWO_PI


def _is_integer(value, tol=GATE_TOL):
    """True when ``value`` is within ``tol`` of an integer."""
    return abs(value - round(value)) <= tol


def closes(m, k, radius):
    """Frozen gate: True only for exact closing modes.

    Admits ``(m, k)`` exactly when ``m`` is an integer and
    ``k * R / (2*pi)`` is a positive integer (at least one radial
    cycle).  ``k = 0`` is refused: a mode with no radial cycle has no
    pitch to fit the disk.
    """
    if not _is_integer(m):
        return False
    n = radial_cycles(k, radius)
    if not _is_integer(n):
        return False
    return round(n) >= 1


def mode_numbers(m, k, radius):
    """``(m, n)`` for a closing mode, or ``None`` when the gate refuses."""
    if not closes(m, k, radius):
        return None
    return (int(round(m)), int(round(radial_cycles(k, radius))))


def quantized_radii(k, radius):
    """Ring radii ``r_j = j * p`` for a closing mode, or ``None``.

    Returns the list ``[p, 2p, ..., R]`` (the last ring is the rim).
    """
    numbers = mode_numbers(0, k, radius)
    if numbers is None:
        return None
    n = numbers[1]
    p = pitch(k)
    return [j * p for j in range(1, n + 1)]


def field(r, theta, m, k, lam):
    """Complex field ``psi(r, theta)`` for one mode."""
    return envelope(r, lam) * cmath.exp(1j * (m * theta + k * r))


def field_on_grid(m, k, lam, radius, n_r=32, n_theta=64):
    """Field sampled on a polar grid (half-offset radii, so r = 0 is avoided).

    Returns ``(radii, thetas, values)`` where ``values[i][j]`` is
    ``psi(radii[i], thetas[j])``.
    """
    if n_r < 1 or n_theta < 1:
        raise ValueError("grid sizes must be >= 1")
    radii = [(i + 0.5) * radius / n_r for i in range(n_r)]
    thetas = [TWO_PI * j / n_theta for j in range(n_theta)]
    values = [[field(r, t, m, k, lam) for t in thetas] for r in radii]
    return radii, thetas, values


def azimuthal_closure_ratio(m):
    """``psi(theta + 2*pi) / psi(theta)``: exactly 1 iff ``m`` is an integer.

    This is the numeric form of the azimuthal closure rule.
    """
    return cmath.exp(1j * TWO_PI * m)


def enumerate_modes(m_values, n_values, radius):
    """All closing modes ``(m, n, k)`` for the given mode-number ranges."""
    modes = []
    for m in m_values:
        for n in n_values:
            if n < 1:
                continue
            k = TWO_PI * n / radius
            if closes(m, k, radius):
                modes.append((m, n, k))
    return modes


def visible_radius(lam, radius, tau=DEFAULT_TAU):
    """Largest radius where the envelope stays above ``tau``, capped at ``R``."""
    if not (0.0 < tau < 1.0):
        raise ValueError("tau must lie strictly between 0 and 1")
    return min(radius, lam * math.log(1.0 / tau))


def visible_pitches(n, lam, radius, tau=DEFAULT_TAU):
    """Pitches visible inside the resolved region for mode number ``n``."""
    r_max = visible_radius(lam, radius, tau)
    return int(math.floor(n * r_max / radius))


def is_resolved(n, lam, radius, tau=DEFAULT_TAU):
    """Resolution rule: at least two visible pitches, or all of them if fewer."""
    return visible_pitches(n, lam, radius, tau) >= min(MIN_VISIBLE_PITCHES, n)


def resolved_modes(m_values, n_values, radius, lam, tau=DEFAULT_TAU):
    """Closing modes that are also resolved at this decay length."""
    return [
        (m, n, k)
        for (m, n, k) in enumerate_modes(m_values, n_values, radius)
        if is_resolved(n, lam, radius, tau)
    ]


def count_closing(candidates, radius):
    """How many seeded ``(m, k)`` candidates the gate admits, and which."""
    admitted = [(m, k) for (m, k) in candidates if closes(m, k, radius)]
    refused = [(m, k) for (m, k) in candidates if not closes(m, k, radius)]
    return admitted, refused


def random_candidates(count, seed, m_range=(-4.0, 4.0), n_range=(-1.0, 5.0)):
    """Negative-control candidates: continuous random ``(m, k)`` pairs.

    ``m`` is uniform in ``m_range`` and ``k`` is drawn so that
    ``k * R / (2*pi)`` is uniform in ``n_range`` (radius is fixed to 1
    here; rescale ``k`` by ``1 / radius`` for another disk).  With
    continuous draws, an exact integer hit has probability zero, so the
    gate is expected to refuse every one.
    """
    rng = random.Random(seed)
    out = []
    for _ in range(count):
        m = rng.uniform(*m_range)
        n = rng.uniform(*n_range)
        out.append((m, TWO_PI * n))
    return out


def best_match_score(random_phases, modes, radius, lam, n_r=16, n_theta=32):
    """Largest |correlation| between a random-phase field and any mode.

    ``random_phases`` is a grid of phases (same shape as the mode grid).
    Each mode field is normalised, so the score is 1 for an exact match
    and near 0 for unrelated phases.  Used as a field-level negative
    control: random phases must not match any closing mode.
    """
    radii = [(i + 0.5) * radius / n_r for i in range(n_r)]
    thetas = [TWO_PI * j / n_theta for j in range(n_theta)]
    best = 0.0
    for (m, _n, k) in modes:
        num = 0.0j
        den_mode = 0.0
        den_rand = 0.0
        for i, r in enumerate(radii):
            for j, t in enumerate(thetas):
                mode_val = field(r, t, m, k, lam)
                rand_val = cmath.exp(1j * random_phases[i][j])
                num += rand_val * mode_val.conjugate()
                den_mode += abs(mode_val) ** 2
                den_rand += abs(rand_val) ** 2
        score = abs(num) / math.sqrt(den_mode * den_rand)
        if score > best:
            best = score
    return best
