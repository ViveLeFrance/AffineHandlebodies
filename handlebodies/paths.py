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


#: samples per smoothed corner.  Odd, so that the node itself -- which the
#: curve passes through at its midpoint -- is one of them.  The path is
#: resampled by linear interpolation along these, so this and not the final
#: `steps` is what decides how smooth a rounded corner comes out.
TURN_SAMPLES = 65


def _cubic(p0, p1, p2, p3, s):
    return ((1 - s) ** 3 * p0 + 3 * (1 - s) ** 2 * s * p1
            + 3 * (1 - s) * s ** 2 * p2 + s ** 3 * p3)


def _turn(a, b, c, radius):
    """A smooth turn through `b`, coming from `a` and leaving towards `c`.

    Two cubics meeting at `b`, the first tangent to the incoming segment and
    the second to the outgoing one, both tangent at `b` to the bisector of the
    two.  So the path still visits every waypoint: rounding a corner rounds the
    turn taken through it, it does not cut the corner off.  `radius` in [0, 1]
    is the fraction of the shorter adjacent segment the turn is spread over; 0
    leaves the corner sharp, and the node alone is returned.

    A reversal has no bisector to turn through and is left sharp.
    """
    a, b, c = complex(a), complex(b), complex(c)
    reach = min(abs(b - a), abs(c - b))
    back = reach * 0.5 * min(max(float(radius), 0.0), 1.0)
    if back <= 0:
        return np.array([b])
    u = (b - a) / abs(b - a)
    w = (c - b) / abs(c - b)
    heading = u + w
    if abs(heading) < 1e-12:
        return np.array([b])
    d = heading / abs(heading)

    half = TURN_SAMPLES // 2
    s = np.linspace(0.0, 1.0, half + 1)
    # Tangent magnitudes equal to each half's chord, which is the usual choice
    # for a Hermite arc and keeps the two halves the same shape.
    into = _cubic(b - back * u, b - (2 * back / 3) * u,
                  b - (back / 3) * d, b, s)
    away = _cubic(b, b + (back / 3) * d,
                  b + (2 * back / 3) * w, b + back * w, s)
    return np.concatenate([into[:-1], away])   # the node itself, once


def _radii(rounding, n):
    if np.isscalar(rounding):
        return [float(rounding)] * n
    radii = [float(r) for r in rounding]
    return radii + [0.0] * (n - len(radii))


def rounded_path(points, steps=70, rounding=0.0, closed=False):
    """A path through `points`, with corners optionally rounded off.

    `rounding` is either one number applied to every corner, or one number per
    waypoint so that corners can be rounded individually.  Each value in [0, 1]
    is the fraction of the shorter adjacent segment the turn is spread over; 0
    leaves it sharp.  The path passes through every waypoint either way -- see
    `_turn` -- so rounding changes how a corner is taken, never where the path
    goes.

    `closed` says the first and last points are the same point and that it is a
    corner of the loop like any other, to be rounded with the first radius.
    The path then begins and ends at that node, in the middle of its turn.

    The result is resampled at uniform arc length, so points are spread evenly
    along the path rather than bunching in its short segments.
    """
    pts = [complex(p) for p in points]
    if len(pts) < 2:
        return np.full(max(steps, 1), pts[0] if pts else 0j)
    radii = _radii(rounding, len(pts))

    turns = [_turn(pts[i - 1], pts[i], pts[i + 1], radii[i])
             for i in range(1, len(pts) - 1)]
    # The seam of a closed path is a corner whose two sides are the last
    # segment and the first; it is opened in the middle of its own turn, which
    # is the node, so the path still starts and ends there exactly.
    seam = (_turn(pts[-2], pts[0], pts[1], radii[0])
            if closed and len(pts) > 2 and pts[0] == pts[-1] else None)
    half = TURN_SAMPLES // 2

    pieces = []
    if seam is None or len(seam) == 1:
        seam, cursor = None, pts[0]
    else:
        pieces.append(seam[half:])
        cursor = seam[-1]
    for turn in turns:
        pieces.append(np.array([cursor, turn[0]]))
        if len(turn) > 1:
            pieces.append(turn)
        cursor = turn[-1]
    if seam is None:
        pieces.append(np.array([cursor, pts[-1]]))
    else:
        pieces.append(np.array([cursor, seam[0]]))
        pieces.append(seam[:half + 1])

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
