"""Paths in the base of a fibration."""

import numpy as np


def pl_path(points, steps=70):
    """Piecewise-linear path through `points`, sampled at roughly `steps` points."""
    pts = [complex(p) for p in points]
    if len(pts) == 1:
        return np.full(steps, pts[0])
    per = max(steps // (len(pts) - 1), 2)
    out = []
    for a, b in zip(pts, pts[1:]):
        out.extend((1 - s) * a + s * b for s in np.linspace(0, 1, per))
    out.extend([pts[-1]] * (steps - len(out)))
    return np.array(out)


def circle(radius, centre=0, steps=200, anticlockwise=True):
    """Circle of the given radius about `centre`."""
    s = np.linspace(0, 2 * np.pi, steps)
    if not anticlockwise:
        s = -s
    return centre + radius * (np.sin(s) + 1j * np.cos(s))
