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


#: a correction this small relative to the points is as accurate as double
#: precision gets on these systems, whatever the step-size rule would prefer
SETTLED = 1e-8

#: how far the separation between solutions may fall in a single step; they
#: genuinely converge as the path runs into a critical value, so this only
#: guards against a step so large that two of them swap places
SHRINK = 0.05


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


def _newton(system, p, t, tol=1e-11, settled=1e-7, iterations=60):
    """Refine one point; returns (point, converged).

    Convergence is judged on the size of the Newton increment relative to the
    point, never on the residual.  These equations carry high powers of
    coordinates of appreciable size, so their residuals run to astronomical
    values at points that are perfectly converged; an absolute residual test
    simply never succeeds.

    It also has to accept a point that has stopped improving.  Those same high
    powers put a floor on how small the increment can get -- on a surface with
    a term like y^8 z^4 and coordinates near 5, the residual carries values
    around 1e8, and roundoff alone leaves the increment bouncing near 1e-10.
    Insisting on `tol` there rejects points that are converged as far as double
    precision allows, and one such point stops the whole path.  So an increment
    that has ceased to shrink counts as converged provided it is already small
    (`settled`) against the point itself -- which is far finer than anything
    the tracking needs to tell two solutions apart.
    """
    previous = np.inf
    for _ in range(iterations):
        f = system.f(p, t)
        if not np.all(np.isfinite(f)):
            return p, False
        try:
            delta = np.linalg.solve(system.jacobian(p, t), f)
        except np.linalg.LinAlgError:
            return p, False
        p = p - delta
        size = float(np.max(np.abs(delta))) / max(1.0, float(np.max(np.abs(p))))
        if size <= tol:
            return p, True
        # Stagnation only counts once the increment is already small.  Newton's
        # early iterations are linear, not quadratic, so a perfectly healthy
        # run can go 1e-1 then 6e-2; treating that as "no longer improving"
        # abandons the point before convergence has even started, and one such
        # point stops the whole path.
        if size <= settled and size >= 0.5 * previous:
            return p, True
        previous = size
    return p, previous <= settled


def separation(points):
    """The closest distance between two of the tracked points."""
    if len(points) < 2:
        return np.inf
    d = np.linalg.norm(points[:, None, :] - points[None, :, :], axis=2)
    np.fill_diagonal(d, np.inf)
    return float(d.min())


def seed(fibre, variables, coeffs, report=None):
    """The critical points on one fibre, found exactly.

    Eliminating one coordinate gives the other exactly; the companion is picked
    from the roots of the fibre by which of them also satisfies the minor, then
    polished.  The polish matters: the companion arrives with only the accuracy
    that root-finding on a substituted polynomial gives it.

    `report`, if given, collects how many points could not be refined.
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
    both_at = sp.lambdify((u, v), equations, 'numpy')
    dboth_at = sp.lambdify((u, v), [sp.diff(e, v) for e in equations], 'numpy')

    points, unrefined = [], 0
    for uv in us:
        candidates = np.roots([complex(f(uv)) for f in coefficients])
        if not len(candidates):
            raise TrackingFailure('the fibre degenerates at the start of the path')
        # Try the companions in order of how well each satisfies the minor, and
        # keep the first that refines.  Picking only the best-looking one and
        # giving up if it fails throws away a solution the eliminant has
        # already found exactly, on the strength of one guess at its partner.
        resid = [abs(complex(minor_at(uv, c))) for c in candidates]
        companion, converged = _refine_companion(
            both_at, dboth_at, uv, candidates[int(np.argmin(resid))])
        if not converged:
            unrefined += 1
        points.append(np.array([uv, companion], complex))
    if report is not None:
        report['unrefined_seeds'] = unrefined
    return np.array(points, complex)


def _refine_companion(both_at, dboth_at, u_value, guess, tol=1e-12,
                      settled=1e-7, iterations=40):
    """Sharpen the second coordinate with the first held fixed.

    The eliminant gives the first coordinate exactly, so it must not be allowed
    to move: refining both at once lets two distinct solutions drift onto the
    same one, and the pair is then lost -- which leaves the tracker with two
    identical points, no separation between them, and nowhere to go.

    Both equations are used, in the least-squares sense.  With one coordinate
    fixed there is one unknown and two conditions, consistent because the fixed
    coordinate came from the eliminant; leaning on the fibre alone fails
    wherever the fibre happens to run tangent to the direction being solved in.
    """
    y = complex(guess)
    previous = np.inf
    for _ in range(iterations):
        f = np.array(both_at(u_value, y), complex)
        d = np.array(dboth_at(u_value, y), complex)
        denom = float(np.vdot(d, d).real)
        if denom == 0 or not np.all(np.isfinite(f)) or not np.all(np.isfinite(d)):
            return y, False
        step = complex(np.vdot(d, f)) / denom
        y -= step
        size = abs(step) / max(1.0, abs(y))
        if size <= tol:
            return y, True
        if size >= 0.5 * previous:       # stopped improving; see _newton
            return y, size <= settled
        previous = size
    return y, previous <= settled


def _step(system, points, t, h):
    """Advance every point from t to t + h.

    Returns (points, correction, stubborn) where `stubborn` lists the points
    Newton could not settle.  Those keep their last iterate rather than sinking
    the step: one solution near a branch point can refuse to converge while the
    other thirty behave perfectly, and abandoning the step for its sake stops
    the path outright.  The correction reported is over the points that did
    settle, since that is what the step size should be judged on.
    """
    predicted = np.empty_like(points)
    for i, p in enumerate(points):
        try:
            tangent = -np.linalg.solve(system.jacobian(p, t), system.dt(p, t))
            predicted[i] = p + tangent * h
        except np.linalg.LinAlgError:
            predicted[i] = p

    corrected = np.empty_like(points)
    stubborn = []
    for i, q in enumerate(predicted):
        r, ok = _newton(system, q, t + h)
        corrected[i] = r
        if not ok:
            stubborn.append(i)
    settled_only = [i for i in range(len(points)) if i not in stubborn]
    if not settled_only:
        return None, np.inf, stubborn
    correction = float(np.abs(corrected[settled_only]
                              - predicted[settled_only]).max())
    return corrected, correction, stubborn


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
        # As two solutions genuinely approach, the allowed correction shrinks
        # with their separation and eventually falls under what roundoff can
        # deliver, so no step could ever be accepted and the path stops right
        # where it gets interesting.  A step whose corrector barely had to move
        # anything is accurate whatever the solutions are doing, so that is
        # accepted too.
        scale = max(1.0, float(np.abs(points).max()))
        allowed = max(tolerance * before, SETTLED * scale)
        moved, correction, stubborn = _step(system, points, t, span * h)
        # Two guards. The correction has to be small against how far apart the
        # solutions are -- never against how large they are -- and the step
        # must not itself collapse that separation. Without the second, a step
        # can quietly land two tracked points on the same solution while its
        # correction still looks respectable; every later step is then measured
        # against a separation of zero and nothing can be accepted again.
        # Only reject for reasons a smaller step can mend.  A point that will
        # not settle will not settle at half the size either, so counting that
        # against the step just halves forever until the floor is hit -- which
        # is what stopped these paths a few steps in, while a fixed step of the
        # same size walked the whole way without difficulty.  Stubborn points
        # are recorded and carried; the correction and the separation are the
        # tests that respond to step size, and they decide.
        if (moved is not None and correction <= allowed
                and separation(moved) >= SHRINK * before):
            points = moved
            done += h
            report['stubborn'] = max(report.get('stubborn', 0), len(stubborn))
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
    report = {'subdivisions': 0, 'unreliable': 0, 'stopped': None,
              'min_separation': separation(points), 'max_residual': 0.0}

    for a, b in zip(path, path[1:]):
        try:
            points = _advance(system, points, a, b, tolerance, floor, report)
        except TrackingFailure as exc:
            # Stop, but keep everything already found.  Each frame is the
            # critical values over one point of the path -- solved for, not
            # guessed -- and they do not stop being true because the next step
            # could not be taken.
            report['stopped'] = str(exc)
            break
        frames.append(points.copy())
        report['min_separation'] = min(report['min_separation'],
                                       separation(points))
        report['max_residual'] = max(report['max_residual'],
                                     _relative_residual(system, points, b))
    report['covered'] = len(frames)
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
