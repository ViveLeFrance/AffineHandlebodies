"""Exact polynomial algebra underlying the Lefschetz computations.

Everything here works over QQ (or QQ(i)) via Groebner bases, so the results are
exact up to the final root-finding step.  Coefficients supplied by the user are
rationalised on the way in: ``0.3`` means the rational 3/10, not the binary
float nearest to it.
"""

import itertools
import os
import re

import numpy as np
import sympy as sp
from sympy.polys.orderings import ProductOrder, grevlex

from . import singular as backend

W = sp.Symbol('W')  # stands for the value of a fibration


class PositiveDimensional(ValueError):
    """Raised when a system that has to be finite turns out not to be."""


def use_singular():
    """Whether to route Groebner work through the Singular binary.

    Set HANDLEBODIES_NO_SINGULAR=1 to force the pure sympy path.
    """
    if os.environ.get('HANDLEBODIES_NO_SINGULAR'):
        return False
    return backend.available()


_NAME = re.compile(r'[A-Za-z_]\w*')


def _symbols(text):
    """Map every identifier in `text` to a plain Symbol.

    Without this, names such as `beta`, `gamma` or `zeta` would resolve to
    sympy's special functions rather than to the parameters they denote here.
    """
    names = {n: sp.Symbol(n) for n in _NAME.findall(text) if n != 'I'}
    names['I'] = sp.I
    return names


def sympify(expr):
    """sympify a string, treating all identifiers as parameters."""
    if not isinstance(expr, str):
        return expr
    return sp.sympify(expr, locals=_symbols(expr))


def exact(c):
    """Coerce a user-supplied coefficient to an exact sympy number.

    Accepts ints, rationals, floats, and strings such as ``"3/2"`` or
    ``"-0.2*I"``.  Floats are read as decimals, so 0.3 -> 3/10.
    """
    if isinstance(c, complex):
        c = sp.Rational(str(c.real)) + sp.Rational(str(c.imag)) * sp.I
    elif isinstance(c, float):
        c = sp.Rational(str(c))
    else:
        c = sympify(c)
    return sp.nsimplify(c, rational=True)


def poly(expr, subs=None):
    """Parse an expression (string or sympy) and specialise its parameters."""
    e = sympify(expr)
    if subs:
        e = e.subs({sp.Symbol(k) if isinstance(k, str) else k: exact(v)
                    for k, v in subs.items()})
    return sp.expand(sp.nsimplify(e, rational=True))


def linear_coeffs(expr, variables):
    """Read the coefficients of a linear form off an expression.

    Lets a fibration be written the way it is on paper -- "0.3*x + 0.1*y + z"
    -- rather than as a bare tuple whose order has to be remembered.
    """
    e = sp.expand(sympify(expr))
    p = sp.Poly(e, *variables)
    if p.total_degree() > 1:
        raise ValueError(f"the fibration must be linear; got {e}")
    const = e.subs({v: 0 for v in variables})
    if sp.simplify(const) != 0:
        raise ValueError(f"the fibration must have no constant term; got {const}")
    return [exact(e.coeff(v)) for v in variables]


def parallel_minors(G, variables, coeffs):
    """Equations saying grad(G) is parallel to the constant vector `coeffs`.

    The critical points of a linear fibration with these coefficients, on the
    hypersurface {G = 0}, are cut out by these 2x2 minors together with G. The minors are those appearing in the cross product grad pi x grad G.
    """
    gradG = [sp.diff(G, v) for v in variables]
    minors = [sp.expand(gradG[i] * coeffs[j] - gradG[j] * coeffs[i])
              for i, j in itertools.combinations(range(len(variables)), 2)]
    return [m for m in minors if m != 0]


def linear_eliminant(curve, minor, variables, coeffs, extra=W):
    """Eliminate both fibre coordinates from a plane curve, by substitution.

    A general elimination is wasted here.  The fibration is linear, so the
    equation `c.v = extra` can be solved for one coordinate and substituted,
    and what is left is two polynomials in the other coordinate alone --
    whose resultant is the eliminant.  One Sylvester determinant, polynomial
    in the degrees, where a Groebner basis is doubly exponential in the worst
    case: on a quintic surface this is the difference between a tenth of a
    second and a minute of Singular.

    The substitution is not itself the saving, and it is worth not being
    misled about that: handing the substituted pair to a Groebner engine is no
    better than the original system and measurably worse for Singular, which
    was happy to use the linear generator itself and would rather not be given
    a denser one with the denominator cleared.  What the substitution buys is
    the shape -- two polynomials in one unknown -- in which elimination stops
    being a basis computation and becomes a determinant.  So callers falling
    back should hand over the original system, not the substituted one.

    Returns None when the structure does not apply, so the caller can fall
    back to the general route rather than be wrong.
    """
    u, v = variables
    c_u, c_v = [sp.sympify(c) for c in coeffs]
    if c_u == 0:
        u, v, c_u, c_v = v, u, c_v, c_u
    if c_u == 0:
        return None                      # the form involves neither coordinate
    substituted = []
    for equation in (curve, minor):
        e = sp.together(sp.expand(equation.subs(u, sp.together((extra - c_v * v) / c_u))))
        substituted.append(sp.Poly(sp.numer(e), v))
    if any(e.is_zero for e in substituted):
        return None
    res = sp.expand(sp.resultant(*substituted))
    if res == 0:
        # the two curves share a component, so the critical locus is a curve
        raise PositiveDimensional(
            "the critical locus is positive-dimensional, so there is no finite "
            "set of critical values to compute; this normally means the surface "
            "is singular along a curve contained in it")
    return sp.Poly(res, extra)


def elimination_order(n):
    """A block order eliminating the first `n` of the generators.

    Elimination is what a lex basis is usually reached for, but lex builds a
    full triangular form -- far more than an elimination ideal needs, and
    ruinously more once the system is not small.  A block order asks only for
    what the elimination theorem uses: that a monomial touching an eliminated
    variable outranks every monomial in the kept ones.  It is how Singular's
    own `eliminate` works, and it is the difference on a quartic surface
    between twenty seconds and a lex basis that does not finish.
    """
    return ProductOrder((grevlex, lambda m: m[:n]), (grevlex, lambda m: m[n:]))


def eliminant(equations, variables, extra=W):
    """Eliminate `variables`, leaving a univariate polynomial in `extra`.

    Its roots are the values `extra` takes on the solutions of `equations` --
    for a critical system, the critical values of the fibration.
    """
    order = list(variables) + [extra]

    # Check zero-dimensionality first
    if quotient_dimension(equations, order) is None:
        raise PositiveDimensional(
            "the critical locus is positive-dimensional, so there is no finite "
            "set of critical values to compute; this normally means the surface "
            "is singular along a curve contained in it")

    if use_singular():
        polys = backend.eliminate(equations, order, variables)
    else:
        polys = list(sp.groebner(equations, *order,
                                 order=elimination_order(len(variables))).exprs)

    univariate = [p for p in polys if p.free_symbols <= {extra}]
    if not univariate:
        raise ValueError(
            "elimination produced no univariate polynomial; the fibration is "
            "probably not generic enough for its critical values to be finite")
    return sp.Poly(univariate[0], extra)


def quotient_dimension(equations, variables):
    """dim_C C[variables]/I, i.e. the number of solutions with multiplicity.

    Returns None if the ideal is not zero-dimensional.
    """
    if use_singular():
        return backend.vdim(equations, list(variables))
    basis = sp.groebner(equations, *variables, order='grevlex')
    # NB: Poly.monoms() sorts by lex regardless of the basis order, so the
    # leading monomial has to be requested for grevlex explicitly.
    leads = [tuple(sp.Poly(p, *variables).LM(order='grevlex').exponents)
             for p in basis.polys]

    bounds = []
    for i in range(len(variables)):
        pure = [m[i] for m in leads if all(e == 0 for j, e in enumerate(m) if j != i)]
        if not pure:
            return None  # no pure power of this variable => positive-dimensional
        bounds.append(min(pure))

    count = 0
    for m in itertools.product(*[range(b) for b in bounds]):
        if not any(all(m[i] >= l[i] for i in range(len(m))) for l in leads):
            count += 1
    return count


def roots(p, prec=30):
    """Roots of a univariate sympy Poly, to `prec` significant digits."""
    return np.array([complex(r) for r in p.nroots(n=prec, maxsteps=500)])


def exact_point(value, digits=9):
    """An exact rational (or Gaussian rational) close to a complex number.

    Path points arrive as floats, so there is no exactness to preserve, and the
    denominator this produces goes straight into the coefficients an exact
    solver has to carry. 
    """
    c = complex(value)
    return (sp.Rational(str(round(c.real, digits)))
            + sp.Rational(str(round(c.imag, digits))) * sp.I)


def numeric_root_fn(p, param):
    """Compile a Poly in (W, param) into a fast callable param -> roots.

    Used for the path-tracing loops, where accuracy of 1e-12 is ample and the
    cost of exact root-finding at every step is not.
    """
    coeffs = [sp.lambdify(param, c, 'numpy') for c in sp.Poly(p, W).all_coeffs()]
    def at(t):
        return np.roots([complex(f(complex(t))) for f in coeffs])
    return at


def homogenize(expr, variables, w):
    """Homogenize `expr` with respect to the new variable `w`."""
    p = sp.Poly(expr, *variables)
    d = p.total_degree()
    return sp.expand(sum(
        c * sp.prod([v ** e for v, e in zip(variables, mono)]) * w ** (d - sum(mono))
        for mono, c in zip(p.monoms(), p.coeffs())))


def singular_points(G, variables):
    """Points where {G = 0} fails to be smooth.

    Returns [] exactly when the affine variety is smooth: the singular locus is
    empty iff the ideal (G, grad G) is the unit ideal.
    """
    equations = [G] + [sp.diff(G, v) for v in variables]
    if use_singular():
        if backend.is_unit_ideal(equations, list(variables)):
            return []
    else:
        basis = sp.groebner(equations, *variables, order='grevlex')
        if list(basis.exprs) == [sp.Integer(1)]:
            return []
    return sp.solve(equations, variables, dict=True)
