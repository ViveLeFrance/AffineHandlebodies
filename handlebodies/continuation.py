"""Numerical path tracking for the critical values of rho along a path in pi.

The exact route computes one eliminant in (W, t) and then finds roots per step;
that elimination is instant for the recorded examples but is out of reach past
about degree 10, and solving each fibre from scratch instead costs ~19s for a
70-step path.  Consecutive steps differ only slightly, though, so the solutions
at one step are excellent starting points for the next.  That is what this
module does: seed exactly at the start of the path, then follow the solutions.

The one thing that makes this delicate is scale.  At degree 10 the solutions
have magnitude ~2800 but sit only ~1 apart, so a corrector tolerance keyed to
the size of the coordinates lets Newton stride from one solution onto another
without anything looking wrong.  Every tolerance here is therefore relative to
the *separation between the tracked points*, never to their magnitude.

Only the affine fibre coordinates are tracked, so a solution escaping to
infinity would be lost rather than followed.  Callers are expected to check the
result against an exact solve; `Bifibration` does.
"""

import numpy as np
import sympy as sp

from . import algebra as A


class TrackingFailure(RuntimeError):
    """Raised when the solutions cannot be followed reliably."""


class System:
    """The critical-point equations of a linear form on a fibre, compiled.

    Two equations in the two fibre coordinates -- the fibre itself and the
    parallel-minor condition -- whose coefficients depend on the base point t.
    """

    def __init__(self, fibre, variables, coeffs, t):
        self.variables = list(variables)
        self.coeffs = [complex(c) for c in coeffs]
        self.t = t
        u, v = self.variables
        minor = A.parallel_minors(fibre, self.variables, coeffs)[0]
        self.equations = (sp.expand(fibre), sp.expand(minor))
        args = (u, v, t)
        self._f = sp.lambdify(args, list(self.equations), 'numpy')
        self._j = sp.lambdify(args, [[sp.diff(e, w) for w in (u, v)]
                                     for e in self.equations], 'numpy')
        self._dt = sp.lambdify(args, [sp.diff(e, t) for e in self.equations],
                               'numpy')

    def f(self, p, t):
        return np.array(self._f(p[0], p[1], t), complex)

    def jacobian(self, p, t):
        return np.array(self._j(p[0], p[1], t), complex)

    def dt(self, p, t):
        return np.array(self._dt(p[0], p[1], t), complex)

    def value(self, p):
        """The value of the linear form at a solution."""
        return self.coeffs[0] * p[0] + self.coeffs[1] * p[1]

    def values(self, points):
        return np.array([self.value(p) for p in points])


def _newton(system, p, t, tol=1e-11, iterations=60):
    """Refine one point; returns (point, converged).

    Convergence is judged on the size of the Newton increment relative to the
    point, never on the residual.  These equations carry a tenth power of a
    coordinate that reaches into the thousands, so their residuals run to
    astronomical values at points that are perfectly converged; an absolute
    residual test simply never succeeds.
    """
    for _ in range(iterations):
        f = system.f(p, t)
        if not np.all(np.isfinite(f)):
            return p, False
        try:
            delta = np.linalg.solve(system.jacobian(p, t), f)
        except np.linalg.LinAlgError:
            return p, False
        p = p - delta
        if np.max(np.abs(delta)) <= tol * max(1.0, float(np.max(np.abs(p)))):
            return p, True
    return p, False


def separation(points):
    """The closest distance between two of the tracked points."""
    if len(points) < 2:
        return np.inf
    d = np.linalg.norm(points[:, None, :] - points[None, :, :], axis=2)
    np.fill_diagonal(d, np.inf)
    return float(d.min())


def seed(fibre, variables, coeffs):
    """The critical points on one fibre, found exactly.

    Eliminating one coordinate gives the other exactly; the companion is picked
    from the roots of the fibre by which of them also satisfies the minor, then
    polished.  The polish matters: the companion arrives with only the accuracy
    that root-finding on a substituted polynomial gives it.
    """
    u, v = variables
    minor = A.parallel_minors(fibre, list(variables), coeffs)[0]
    equations = [sp.expand(fibre), sp.expand(minor)]
    us = A.roots(A.eliminant(equations, [v], extra=u), prec=20)

    t = sp.Symbol('_seed_t')
    system = System(fibre, variables, coeffs, t)

    # Everything below is evaluated numerically: substituting into a degree-10
    # expression symbolically, once per candidate root, dominates the cost.
    coefficients = [sp.lambdify(u, c, 'numpy')
                    for c in sp.Poly(equations[0], v).all_coeffs()]
    minor_at = sp.lambdify((u, v), equations[1], 'numpy')

    points = []
    for uv in us:
        candidates = np.roots([complex(f(uv)) for f in coefficients])
        if not len(candidates):
            raise TrackingFailure('the fibre degenerates at the start of the path')
        resid = [abs(complex(minor_at(uv, c))) for c in candidates]
        p, ok = _newton(system, np.array([uv, candidates[int(np.argmin(resid))]],
                                         complex), 0)
        if not ok:
            raise TrackingFailure('a starting point could not be refined')
        points.append(p)
    return np.array(points, complex)


def _step(system, points, t, h):
    """Advance every point from t to t + h; returns (points, correction).

    The prediction is the tangent to each solution path; the correction is how
    far Newton then had to move, which is what the step size is judged on.
    """
    predicted = np.empty_like(points)
    for i, p in enumerate(points):
        try:
            tangent = -np.linalg.solve(system.jacobian(p, t), system.dt(p, t))
            predicted[i] = p + tangent * h
        except np.linalg.LinAlgError:
            predicted[i] = p

    corrected = np.empty_like(points)
    for i, q in enumerate(predicted):
        r, ok = _newton(system, q, t + h)
        if not ok:
            return None, np.inf
        corrected[i] = r
    return corrected, float(np.abs(corrected - predicted).max())


def _advance(system, points, t0, t1, tolerance, floor, report, max_trials=80):
    """Follow the points across one segment, subdividing where needed.

    Both exits are bounded.  An earlier version accepted a step once the size
    hit the floor and carried on from there, which does not recover: `done`
    then advances by about 1e-11 a time and the segment never finishes.  If the
    step has to collapse that far, the honest answer is that these solutions
    cannot be followed here, and the caller should solve exactly instead.
    """
    span = t1 - t0
    if span == 0:
        return points
    done, h, trials = 0.0, 1.0, 0
    while done < 1.0 - 1e-12:
        trials += 1
        if trials > max_trials:
            raise TrackingFailure(
                f'no progress after {max_trials} attempts near t = '
                f'{t0 + span * done:.6g}')
        h = min(h, 1.0 - done)
        if abs(span) * h < floor:
            raise TrackingFailure(
                f'step collapsed below {floor:g} near t = '
                f'{t0 + span * done:.6g}; the solutions come too close to '
                f'follow here')
        t = t0 + span * done
        before = separation(points)
        moved, correction = _step(system, points, t, span * h)
        # Two guards. The correction has to be small against how far apart the
        # solutions are -- never against how large they are -- and the step
        # must not itself collapse that separation. Without the second, a step
        # can quietly land two tracked points on the same solution while its
        # correction still looks respectable; every later step is then measured
        # against a separation of zero and nothing can be accepted again.
        if (moved is not None and correction <= tolerance * before
                and separation(moved) >= 0.5 * before):
            points = moved
            done += h
            h = min(h * 1.7, 1.0)
        else:
            h /= 2.0
            report['subdivisions'] += 1
    return points


def track(system, start, path, tolerance=0.08, floor=1e-9):
    """Follow `start` along `path`, one frame per point of the path.

    `tolerance` is the fraction of the current inter-point separation that a
    corrector step may consume.  Returns (frames, report), the report carrying
    the residual, the closest the points ever came, and how much subdivision
    and how many unreliable steps it took.
    """
    path = [complex(p) for p in path]
    points = np.array(start, complex)
    frames = [points.copy()]
    report = {'subdivisions': 0, 'unreliable': 0,
              'min_separation': separation(points), 'max_residual': 0.0}

    for a, b in zip(path, path[1:]):
        points = _advance(system, points, a, b, tolerance, floor, report)
        frames.append(points.copy())
        report['min_separation'] = min(report['min_separation'],
                                       separation(points))
        report['max_residual'] = max(report['max_residual'],
                                     _relative_residual(system, points, b))
    return frames, report


def _relative_residual(system, points, t):
    """How far one more Newton step would move the points, relative to them.

    The residuals themselves are useless as a measure here (see `_newton`), so
    convergence is reported as the size of the correction still outstanding.
    """
    worst = 0.0
    for p in points:
        try:
            delta = np.linalg.solve(system.jacobian(p, t), system.f(p, t))
        except np.linalg.LinAlgError:
            return np.inf
        worst = max(worst, float(np.max(np.abs(delta)))
                    / max(1.0, float(np.max(np.abs(p)))))
    return worst
