"""Finite non-Morse rho parameters on smooth affine pi fibres.

For P(s,w,t)=0, solve P=P_s=P_ss=0 with P_w nonzero, then remove
parameters having any singular point of the pi fibre. This does not detect
escape to infinity or distinct Morse points sharing a critical value.
Elimination is exact; only the final roots are numerical. A separate process
bounds both Singular and SymPy work, including root finding.
"""
import functools
import json
import math
import os
import signal
import subprocess
import sys
import threading

import sympy as sp

from . import algebra as A

TIMEOUT = 10
_LOCK = threading.Lock()
s, w, t, guard = sp.symbols('s w t guard')


def _projection(equations, drop, fast):
    equations = [sp.expand(e) for e in equations if e != 0]
    if any(not e.free_symbols and e != 0 for e in equations):
        return sp.Poly(1, t)
    variables = [*drop, t]
    if fast:
        basis = A.backend.eliminate(equations, variables, drop)
    else:
        basis = sp.groebner(equations, *variables,
                           order=A.elimination_order(len(drop))).exprs
    polys = [sp.Poly(e, t) for e in basis if e != 0 and e.free_symbols <= {t}]
    if not polys:
        return None  # dense projection, rather than a finite parameter set
    p = polys[0]
    for q in polys[1:]:
        p = sp.gcd(p, q)
    return p.sqf_part().monic()


def _calculate(expression, fast):
    p = sp.Poly(sp.sympify(expression, locals={'s': s, 'w': w, 't': t}),
                s, w, t, extension=sp.I).as_expr()
    ps, pw = sp.diff(p, s), sp.diff(p, w)
    e = _projection([p, ps, sp.diff(ps, s), guard*pw - 1],
                    [s, w, guard], fast)
    if e is None:
        return dict(status='nonfinite', values=[], message=
                    'The non-Morse locus does not give a finite parameter list.')
    if e.degree() > 0:
        # Only test singular fibres at candidate parameters. Computing the
        # entire discriminant of pi can be much more expensive.
        singular = _projection([p, ps, pw, e.as_expr()], [s, w], fast)
        if singular is None:
            return dict(status='nonfinite', values=[], message=
                        'Smooth fibres do not form a generic family.')
        e = e.exquo(sp.gcd(e, singular)).monic()
    roots = A.roots(e) if e.degree() > 0 else []
    values = sorted([[float(z.real), float(z.imag)] for z in roots])
    if not all(math.isfinite(v) for z in values for v in z):
        raise ValueError('Non-finite numerical roots')
    return dict(status='complete', values=values, polynomial=str(e.as_expr()))


@functools.lru_cache(maxsize=16)
def _cached(expression, fast, timeout):
    if sys.platform == 'emscripten':
        return dict(status='unavailable', values=[], message=
                    'This calculation needs the local Python application.')
    proc = subprocess.Popen([sys.executable, '-m', 'handlebodies.exceptional'],
                            stdin=subprocess.PIPE, stdout=subprocess.PIPE,
                            stderr=subprocess.PIPE, text=True,
                            start_new_session=(os.name == 'posix'))
    try:
        out, err = proc.communicate(json.dumps([expression, fast]), timeout=timeout)
    except subprocess.TimeoutExpired:
        if os.name == 'posix':
            os.killpg(proc.pid, signal.SIGKILL)
        else:
            proc.kill()
        proc.communicate()
        return dict(status='timeout', values=[], message=
                    f'Exceptional values were not determined within {timeout:g} seconds.')
    if proc.returncode:
        raise RuntimeError(err.strip() or 'Exceptional-value worker failed')
    return json.loads(out)


def values(bifibration, timeout=TIMEOUT):
    """Return status, values as [real, imag], and (on success) an eliminant.

    Cached by the polynomial family, independently of reference points and
    paths. An incomplete calculation is explicitly distinguished from no values.
    """
    try:
        # Dummies avoid collisions with user coordinate names.
        base, target = sp.Dummy('base'), sp.Dummy('target')
        rho = bifibration.rho_at(base)
        remaining, = rho.other_variables()
        p = rho.fibre(target).xreplace({remaining: s, target: w, base: t})
        with _LOCK:
            report = _cached(str(p), A.use_singular(), timeout)
        return dict(report, values=[list(z) for z in report['values']])
    except Exception as exc:
        return dict(status='unavailable', values=[], message=str(exc))


if __name__ == '__main__':
    try:
        expression, fast = json.load(sys.stdin)
        result = _calculate(expression, fast)
    except Exception as exc:
        result = dict(status='unavailable', values=[], message=str(exc))
    print(json.dumps(result))
