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


def _bezier(p0, p1, p2, n):
    """Quadratic Bezier from p0 to p2 with control point p1."""
    s = np.linspace(0, 1, n)
    return (1 - s) ** 2 * p0 + 2 * (1 - s) * s * p1 + s ** 2 * p2


def rounded_path(points, steps=70, rounding=0.0):
    """A path through `points`, with corners optionally rounded off.

    `rounding` is either one number applied to every corner, or one number per
    waypoint so that corners can be rounded individually.  Each value in [0, 1]
    is the fraction of the shorter adjacent segment that the corner is cut back
    by, the cut being bridged by a quadratic Bezier; 0 leaves it sharp.  Only
    interior waypoints have corners, so a two-point path is always straight.

    The result is resampled at uniform arc length, so points are spread evenly
    along the path rather than bunching in its short segments.
    """
    pts = [complex(p) for p in points]
    if len(pts) < 2:
        return np.full(max(steps, 1), pts[0] if pts else 0j)

    if np.isscalar(rounding):
        radii = [float(rounding)] * len(pts)
    else:
        radii = [float(r) for r in rounding]
        radii += [0.0] * (len(pts) - len(radii))

    if len(pts) == 2 or not any(radii[1:-1]):
        dense = _polyline(pts, 24)
    else:
        pieces, cursor = [], pts[0]
        for i in range(1, len(pts) - 1):
            a, b, c = pts[i - 1], pts[i], pts[i + 1]
            back = min(abs(b - a), abs(c - b)) * 0.5 * min(max(radii[i], 0.0), 1.0)
            if back <= 0:
                # a corner with no rounding still has to be passed through
                pieces.append(np.array([cursor, b]))
                cursor = b
                continue
            enter = b + (a - b) / abs(a - b) * back
            leave = b + (c - b) / abs(c - b) * back
            pieces.append(np.array([cursor, enter]))
            pieces.append(_bezier(enter, b, leave, 16))
            cursor = leave
        pieces.append(np.array([cursor, pts[-1]]))
        dense = np.concatenate([_polyline(list(p), 8) for p in pieces])

    return _resample(dense, steps)


def _polyline(pts, per):
    out = []
    for a, b in zip(pts, pts[1:]):
        out.append(a + (b - a) * np.linspace(0, 1, per))
    return np.concatenate(out) if out else np.array(pts)


def _resample(dense, steps):
    """Resample a dense polyline at uniform arc length."""
    d = np.abs(np.diff(dense))
    s = np.concatenate([[0.0], np.cumsum(d)])
    if s[-1] == 0:
        return np.full(steps, dense[0])
    want = np.linspace(0, s[-1], steps)
    return np.interp(want, s, dense.real) + 1j * np.interp(want, s, dense.imag)
