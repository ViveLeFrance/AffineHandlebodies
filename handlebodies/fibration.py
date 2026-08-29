"""Linear Lefschetz fibrations on affine hypersurfaces, and bifibrations."""

import numpy as np
import sympy as sp

from . import algebra as A
from .algebra import W


def sort_by_angle(points, origin=0, anticlockwise=True):
    """Order points by argument about `origin`, starting from the negative real axis.

    The resulting index is what the paper cites when it names a vanishing
    cycle, so this convention is deliberately identical to the original
    implementation and must not be changed casually.
    """
    pts = [complex(p) for p in points]
    pts.sort(key=lambda p: (-np.pi + np.angle(p - origin)) % (2 * np.pi))
    if not anticlockwise:
        pts.reverse()
    return np.array(pts)


class LinearFibration:
    """A linear map sum(coeffs * variables) restricted to {domain = 0}."""

    def __init__(self, domain, variables, expression, origin=0, solvefor=None):
        self.variables = list(variables)
        self.domain = sp.expand(domain)
        self.coeffs = A.linear_coeffs(expression, self.variables)
        self.origin = A.exact(origin)
        self.solvefor = self._as_symbol(solvefor) if solvefor else self.variables[0]
        self._cache = {}

    def _as_symbol(self, v):
        return sp.Symbol(v) if isinstance(v, str) else v

    @property
    def origin_c(self):
        """The reference fibre as a Python complex, for plotting and sorting."""
        return complex(self.origin)

    @property
    def expression(self):
        return sum(c * v for c, v in zip(self.coeffs, self.variables))

    # -- critical locus ---------------------------------------------------

    def critical_system(self):
        return A.parallel_minors(self.domain, self.variables, self.coeffs) + [self.domain]

    def _eliminant(self):
        """Univariate polynomial in W whose roots are the critical values."""
        if 'eliminant' not in self._cache:
            eqs = self.critical_system() + [sp.expand(self.expression - W)]
            self._cache['eliminant'] = A.eliminant(eqs, self.variables)
        return self._cache['eliminant']

    def critical_values(self, anticlockwise=True):
        """Critical values, ordered and labelled by `sort_by_angle`."""
        return sort_by_angle(A.roots(self._eliminant()), self.origin_c,
                             anticlockwise)

    def n_critical_points(self):
        """Number of critical points counted with multiplicity (exact)."""
        return A.quotient_dimension(self.critical_system(), self.variables)


    # -- fibres -----------------------------------------------------------

    def fibre(self, value, solvefor=None):
        """The fibre over `value`, as a polynomial in the remaining variables."""
        s = self._as_symbol(solvefor) if solvefor else self.solvefor
        value = value if isinstance(value, sp.Expr) else A.exact(value)
        i = self.variables.index(s)
        if self.coeffs[i] == 0:
            raise ValueError(f"cannot eliminate {s}: its coefficient is zero")
        rest = sum(c * v for j, (c, v) in enumerate(zip(self.coeffs, self.variables))
                   if j != i)
        return sp.expand(self.domain.subs(s, sp.cancel((value - rest) / self.coeffs[i])))

    def other_variables(self, solvefor=None):
        s = self._as_symbol(solvefor) if solvefor else self.solvefor
        return [v for v in self.variables if v != s]

    def trace_fibre(self, path):
        """Track the points of the fibre of this fibration along `path`.

        Only meaningful when the domain is a curve, so that the fibre is a
        finite set of points: this is the rho level, where two points of the
        reference fibre colliding over a critical value exhibit the vanishing
        cycle.  Branches are labelled by `sort_by_angle` about the reference
        fibre, matching the labels drawn on the fibre plot.
        """
        rest = self.other_variables()
        if len(rest) != 1:
            raise ValueError(
                f"trace_fibre needs a curve (one variable left after "
                f"eliminating {self.solvefor}); this domain leaves {len(rest)}")
        u, t = rest[0], sp.Symbol('t')
        p = sp.Poly(self.fibre(t), u)
        coeffs = [sp.lambdify(t, c, 'numpy') for c in p.all_coeffs()]

        frames = []
        for s in path:
            c = [complex(f(complex(s))) for f in coeffs]
            r = np.roots(c)
            if len(r) != p.degree():
                raise ValueError(f"fibre degenerates at {s}: leading coefficient vanishes")
            frames.append(r)
        return track(frames, first=sort_by_angle(frames[0], self.origin_c))

    # -- behaviour at infinity --------------------------------------------

    def boundary_points(self, value=None):
        """Points where the closure of the fibre meets the hyperplane at infinity.

        Their number is the number of punctures of the general fibre.
        """
        value = self.origin if value is None else value
        w = sp.Symbol('w')
        V = self.variables
        dom_h = A.homogenize(self.domain, V, w)
        fib_h = A.homogenize(sp.expand(self.expression - A.exact(value)), V, w)

        found = set()
        for k, v in enumerate(V):
            eqs = [e.subs({w: 0, v: 1}) for e in (dom_h, fib_h)]
            rest = [u for u in V if u is not v]
            for sol in sp.solve(eqs, rest, dict=True):
                pt = [sol.get(u, u) for u in V]
                pt[k] = sp.Integer(1)
                if any(getattr(q, 'free_symbols', set()) for q in pt):
                    continue  # positive-dimensional component in this chart
                nz = next((q for q in pt if q != 0), None)
                if nz is None:
                    continue
                point = tuple(sp.simplify(q / nz) for q in pt) + (sp.Integer(0),)
                # Normalising can itself introduce zoo when a chart solution
                # runs off to infinity, so validate the finished point.
                if not all(q.is_number and q.is_finite for q in point):
                    continue
                found.add(point)
        return sorted(found, key=str)

    def transversality(self, value=None):
        """Rank of [grad Xbar ; grad pi ; dw] at each point at infinity.

        Full rank (= number of variables) means the fibre meets the hyperplane
        at infinity transversely, so the fibration extends to a Lefschetz
        pencil and parallel transport is defined.
        """
        value = self.origin if value is None else value
        w = sp.Symbol('w')
        V = self.variables
        coords = V + [w]
        dom_h = A.homogenize(self.domain, V, w)
        fib_h = A.homogenize(sp.expand(self.expression - A.exact(value)), V, w)
        grad_dom = [sp.diff(dom_h, c) for c in coords]
        grad_fib = [sp.diff(fib_h, c) for c in coords]

        out = []
        for pt in self.boundary_points(value):
            sub = dict(zip(coords, pt))
            M = sp.Matrix([[sp.simplify(g.subs(sub)) for g in grad_dom],
                           [sp.simplify(g.subs(sub)) for g in grad_fib],
                           [0] * len(V) + [1]])
            out.append((pt, M.rank()))
        return out



class Bifibration:
    """A pair (pi, rho): rho is a fibration on each fibre of pi.

    The critical values of rho on pi^{-1}(t) are the roots of a single
    polynomial in (W, t), computed once by elimination.  Tracing a path in the
    base of pi then costs one numeric root-find per step instead of a fresh
    symbolic solve, which is what makes matching paths cheap.
    """

    def __init__(self, pi, rho_expression, rho_origin=0, rho_solvefor=None):
        self.pi = pi
        self.rho_vars_all = pi.other_variables()
        self.rho_coeffs = A.linear_coeffs(rho_expression, self.rho_vars_all)
        self.rho_origin = A.exact(rho_origin)
        self.rho_vars = self.rho_vars_all
        self.rho_solvefor = (pi._as_symbol(rho_solvefor) if rho_solvefor
                             else self.rho_vars[0])
        self._cache = {}

    def rho_at(self, t):
        """The fibration rho on the fibre of pi over `t`."""
        expression = sum(c * v for c, v in zip(self.rho_coeffs, self.rho_vars))
        return LinearFibration(self.pi.fibre(t), self.rho_vars, expression,
                               origin=self.rho_origin, solvefor=self.rho_solvefor)

    @property
    def rho_0(self):
        return self.rho_at(self.pi.origin)

    #: how long to spend trying for a closed-form parametric eliminant before
    #: falling back to solving at each step; the recorded examples take <0.1s
    PARAMETRIC_TIMEOUT = 20

    def _parametric(self):
        """Poly in (W, t): roots in W are the critical values of rho on pi^-1(t).

        Raises if no such eliminant can be found in reasonable time; the family
        is not zero-dimensional (t is free), and that elimination gets hard well
        before the zero-dimensional ones do.
        """
        if 'param' not in self._cache:
            t = sp.Symbol('t')
            fibre_t = self.pi.fibre(t)
            eqs = (A.parallel_minors(fibre_t, self.rho_vars, self.rho_coeffs)
                   + [fibre_t,
                      sp.expand(sum(c * v for c, v in zip(self.rho_coeffs, self.rho_vars)) - W)])
            # Eliminate the fibre coordinates, keeping W and the base point t.
            # This family is not zero-dimensional (t is free), so it needs a
            # genuine elimination rather than the shape-lemma route.
            if A.use_singular():
                elim = A.backend.eliminate(eqs, list(self.rho_vars) + [W, t],
                                           self.rho_vars,
                                           timeout=self.PARAMETRIC_TIMEOUT)
                elim = [p for p in elim if p.free_symbols <= {W, t}]
            else:
                basis = sp.groebner(eqs, *self.rho_vars, W, t, order='lex')
                elim = [p for p in basis.exprs if p.free_symbols <= {W, t}]
            if not elim:
                raise ValueError("no eliminant in (W, t); check that rho is generic")
            self._cache['param'] = (sp.Poly(elim[0], W), t)
        return self._cache['param']

    def critical_values_over(self, point):
        """Critical values of rho on the single fibre of pi over `point`."""
        fibre = sp.expand(self.pi.fibre(A.exact_point(point)))
        expression = sum(c * v for c, v in zip(self.rho_coeffs, self.rho_vars))
        rho = LinearFibration(fibre, self.rho_vars, expression,
                              solvefor=self.rho_solvefor)
        return A.roots(rho._eliminant(), prec=16)

    def _frames(self, path):
        """The critical values of rho over each point of `path`.

        The parametric eliminant turns the whole path into one elimination plus
        cheap numeric root-finding, and is what makes this instant for the
        recorded examples.  Where that elimination is out of reach the fibres
        are solved one at a time instead: far slower, but it terminates.
        """
        if self._cache.get('per_step'):
            return [self.critical_values_over(s) for s in path]
        try:
            p, t = self._parametric()
        except Exception:
            self._cache['per_step'] = True
            return [self.critical_values_over(s) for s in path]
        at = A.numeric_root_fn(p, t)
        return [at(s) for s in path]

    def matching_path(self, path):
        """Critical values of rho traced along `path` in the base of pi.

        Returns an array of shape (len(path), n) whose columns are continuous
        branches: entry [k, j] is branch j at path[k].  Two branches meeting at
        the end of a vanishing path is the matching path of rho.
        """
        frames = self._frames(path)
        # Branch j starts at the j-th critical value of rho_0, so branch
        # indices agree with the labels on the critical-value plot.
        return track(frames, first=self.rho_0.critical_values())


def track(frames, first=None):
    """Match roots between consecutive steps into continuous branches.

    Greedy nearest-neighbour is enough here because the steps are small
    relative to the spacing of the roots away from a collision, and at a
    collision the two colliding branches are the ones we want to identify.

    `first` fixes the labelling of the branches: the j-th branch is the one
    starting nearest first[j].  Pass the critical values (or fibre points) in
    their `sort_by_angle` order so that branch indices agree with the labels
    used everywhere else, and hence with the paper.
    """
    frames = [np.asarray(f) for f in frames]
    n = len(frames[0])
    out = np.empty((len(frames), n), complex)
    if first is None:
        out[0] = np.sort_complex(frames[0])
    else:
        pool = list(frames[0])
        for j in range(n):
            i = int(np.argmin(np.abs(np.array(pool) - first[j])))
            out[0, j] = pool.pop(i)
    for k in range(1, len(frames)):
        remaining = list(frames[k])
        for j in range(n):
            i = int(np.argmin(np.abs(np.array(remaining) - out[k - 1, j])))
            out[k, j] = remaining.pop(i)
    return out


def is_closed(path, rtol=1e-9):
    """Whether a path returns to its starting point."""
    path = np.asarray(path)
    scale = max(float(np.abs(path).max()), 1.0)
    return bool(abs(path[-1] - path[0]) < rtol * scale)


def monodromy(branches):
    """The permutation a closed path induces on the tracked branches.

    Entry j is the label of the branch that branch j has moved onto.  A result
    that is not a bijection means the tracking lost a branch, usually because
    the path passes too close to a critical value for the step size.
    """
    first, last = branches[0], branches[-1]
    perm = [int(np.argmin(np.abs(first - last[j]))) for j in range(len(last))]
    return perm if sorted(perm) == list(range(len(perm))) else None


def collision(branches):
    """The pair of branches closest together at the end of a traced path."""
    final = branches[-1]
    d = np.abs(final[:, None] - final[None, :])
    np.fill_diagonal(d, np.inf)
    i, j = np.unravel_index(d.argmin(), d.shape)
    return (int(min(i, j)), int(max(i, j)), float(d[i, j]))
