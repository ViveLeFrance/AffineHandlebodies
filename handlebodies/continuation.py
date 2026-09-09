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

#: how far off the minor a root of the fibre may sit and still count as a
#: critical point over the same first coordinate.  Generous: a root arrives
#: with the accuracy np.roots gives it, and a root that is not a critical
#: point misses by order 1, so there is a wide gap to sit in.
ON_MINOR = 1e-6


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

    Eliminating one coordinate gives the values the other takes, exactly; the
    companions are read off the roots of the fibre by which of them also
    satisfy the minor, then polished.  The polish matters: a companion arrives
    with only the accuracy that root-finding on a substituted polynomial gives
    it.

    Every compatible companion is kept, not just the likeliest one.  The
    eliminant records the distinct values the first coordinate takes, not the
    points themselves: a projection can send two critical points to the same
    value, and one companion per root then seeds fewer points than there are.
    Nothing downstream would notice.  The tracking follows whatever it is
    given, and with a point missing it follows the rest perfectly -- residual
    at roundoff, no subdivision, nothing reported; down to a single point it
    cannot even measure a separation, so its own guards go quiet too.  So the
    count is checked here, against the length of the critical scheme, rather
    than left to `verify`, which compares counts only at the far end of the
    path, is off by default, and is skipped outright if the path stopped short.

    `report`, if given, collects how many points could not be refined.
    """
    u, v = variables
    minor = A.parallel_minors(fibre, list(variables), coeffs)[0]
    equations = [sp.expand(fibre), sp.expand(minor)]
    us = A.roots(A.eliminant(equations, [v], extra=u), prec=20)

    # Everything below is evaluated numerically: substituting into a degree-10
    # expression symbolically, once per candidate root, dominates the cost.
    coefficients = [sp.lambdify(u, c, 'numpy')
                    for c in sp.Poly(equations[0], v).all_coeffs()]
    minor_at = sp.lambdify((u, v), equations[1], 'numpy')
    minor_scale = _term_scale(equations[1], (u, v))
    both_at = sp.lambdify((u, v), equations, 'numpy')
    dboth_at = sp.lambdify((u, v), [sp.diff(e, v) for e in equations], 'numpy')

    points, unrefined = [], 0
    for uv in us:
        candidates = np.roots([complex(f(uv)) for f in coefficients])
        if not len(candidates):
            raise TrackingFailure('the fibre degenerates at the start of the path')
        def off(c, uv=uv):
            return _minor_residual(minor_at, minor_scale, uv, c)
        chosen = [c for c in candidates if off(c) <= ON_MINOR]
        if not chosen:
            # The eliminant put a solution over this value, so there is a
            # companion to find; take the likeliest rather than none at all.
            chosen = [min(candidates, key=off)]
        for guess in chosen:
            companion, converged = _refine_companion(
                both_at, dboth_at, uv, guess)
            if not converged:
                unrefined += 1
            point = np.array([uv, companion], complex)
            # Two candidates straddling one companion refine onto it, which is
            # right where the fibre runs tangent to the direction being solved
            # in and the pair is genuinely one point.  Merging is safe because
            # merging too far is caught: the count below would come up short.
            if not any(_same_point(point, q) for q in points):
                points.append(point)
    if report is not None:
        report['unrefined_seeds'] = unrefined

    length = A.quotient_dimension(equations, [u, v])
    if length is not None and len(points) != length:
        raise TrackingFailure(
            f'seeded {len(points)} critical point'
            f'{"" if len(points) == 1 else "s"} on the starting fibre, but its '
            f'critical scheme has length {length}: '
            + _why_short(equations, [u, v]))
    return np.array(points, complex)


def _why_short(equations, variables):
    """Why the seed does not have as many points as the scheme has length.

    Three things can put it there, and they want opposite responses, so it is
    worth the extra elimination to say which.  Both questions are asked of the
    ideal rather than of the numbers, and only ever on the way to raising, so
    nothing on the working path pays for them.
    """
    fibre, minor = equations
    jacobian = sp.Matrix([[sp.diff(e, w) for w in variables]
                          for e in equations]).det()
    degenerate = A.quotient_dimension(equations + [sp.expand(jacobian)],
                                      variables)
    if degenerate == 0:
        return ('the scheme is reduced, so every point is there to be found '
                'and the seeding missed one -- which is a bug, not a choice '
                'of fibration')

    singular = A.quotient_dimension(
        [fibre] + [sp.diff(fibre, w) for w in variables], variables)
    if singular is None or singular > 0:
        return ('the starting fibre is itself singular, and its singular '
                'point is a critical point of every linear rho, so no choice '
                'of rho will help; move the reference value of pi off pi\'s '
                'critical values')
    return ('rho has a degenerate critical point on a smooth fibre, so it is '
            'not Morse there; perturb rho -- almost any other linear form '
            'will do')


def _term_scale(expr, variables):
    """Compile |expr| with every term made positive: what it sums, in size.

    The yardstick a residual has to be read against.  It has to be built from
    the whole expression, in all its variables, and evaluated at the point:
    scaling against the coefficients left after substituting one coordinate
    measures the wrong thing entirely.  Where the substitution is the very
    thing that cancels -- a minor of 3u^2 - 6 at u = sqrt(2), whose terms are
    each worth 6 -- what survives it is the roundoff, about 1.8e-15, and
    dividing that by itself calls a solution a miss by a factor of one.
    """
    p = sp.Poly(expr, *variables)
    absolute = sum(abs(c) * sp.prod([w ** e for w, e in zip(variables, mono)])
                   for mono, c in zip(p.monoms(), p.coeffs()))
    return sp.lambdify(variables, absolute, 'numpy')


def _minor_residual(minor_at, minor_scale, u_value, v_value):
    """How far off the minor a point sits, relative to the terms summed.

    Never absolute: these expressions carry high powers of coordinates of
    appreciable size, so the residual at a perfectly good root runs to
    whatever those powers happen to be worth and says nothing on its own.

    Zero when there is nothing to be off by, the terms themselves having
    vanished.  That is not a degenerate case to guard against but one of the
    two this is here for -- a minor that vanishes identically over this value
    of the first coordinate is what a projection failing to separate the
    critical points looks like, and every root of the fibre then qualifies.
    """
    scale = float(abs(complex(minor_scale(abs(u_value), abs(v_value)))))
    if scale == 0:
        return 0.0
    return abs(complex(minor_at(u_value, v_value))) / scale


def _same_point(a, b, tol=1e-7):
    """Whether two seeds are the same point, relative to their size."""
    scale = max(1.0, float(np.abs(b).max()))
    return bool(np.abs(a - b).max() <= tol * scale)


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
