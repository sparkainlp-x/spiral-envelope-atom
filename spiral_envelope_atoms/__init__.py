"""spiral-envelope-atoms: a classical wave toy on a 2D disk.

Complex field = decaying radial envelope x helical spiral phase.
"Atoms" are the discrete modes admitted by a frozen closure gate:
integer azimuthal number and an integer number of spiral pitches
across the disk radius.  Synthetic only; no claim about real atoms,
real electromagnetism, or quantum mechanics.
"""

from .modes import (
    GATE_TOL,
    DEFAULT_TAU,
    MIN_VISIBLE_PITCHES,
    TWO_PI,
    best_match_score,
    closes,
    count_closing,
    envelope,
    enumerate_modes,
    azimuthal_closure_ratio,
    field,
    field_on_grid,
    is_resolved,
    mode_numbers,
    pitch,
    quantized_radii,
    radial_cycles,
    random_candidates,
    resolved_modes,
    visible_pitches,
    visible_radius,
)

__all__ = [
    "GATE_TOL",
    "DEFAULT_TAU",
    "MIN_VISIBLE_PITCHES",
    "TWO_PI",
    "best_match_score",
    "closes",
    "count_closing",
    "envelope",
    "enumerate_modes",
    "azimuthal_closure_ratio",
    "field",
    "field_on_grid",
    "is_resolved",
    "mode_numbers",
    "pitch",
    "quantized_radii",
    "radial_cycles",
    "random_candidates",
    "resolved_modes",
    "visible_pitches",
    "visible_radius",
]
