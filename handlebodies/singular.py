"""Optional fast backend: drive the Singular binary for the Groebner work.

sympy's Groebner implementation is pure Python and stalls well before these
problems get interesting -- a degree-10 surface does not finish in minutes,
where Singular takes milliseconds.  Singular is used when the binary is on the
PATH and sympy is the fallback, so nothing here is a hard dependency: the
Painleve examples run either way.

Everything crossing the boundary is exact.  Coefficients are cleared of
denominators on the way out, because Singular reads ``3*x/10`` as polynomial
division rather than as a rational coefficient.
"""

import re
import shutil
import subprocess

import sympy as sp

#: path to the Singular binary, or None
BINARY = shutil.which('Singular')

class SingularError(RuntimeError):
    pass


def available():
    return BINARY is not None


def clear_denoms(expr):
    """An integer-coefficient multiple of `expr`."""
    e = sp.expand(expr)
    dens = []
    for c in e.as_coefficients_dict().values():
        c = sp.sympify(c)
        dens += [sp.denom(sp.re(c)), sp.denom(sp.im(c))]
    dens = [int(d) for d in dens if d]
    return sp.expand(e * sp.ilcm(*dens)) if dens else e


def to_text(expr):
    """Render a sympy polynomial in Singular's syntax."""
    s = str(clear_denoms(expr)).replace('**', '^')
    # the imaginary unit is called i in the ring (0,i)
    return re.sub(r'\bI\b', 'i', s)


def from_text(text, symbols):
    """Parse a Singular polynomial back into sympy."""
    local = {str(s): s for s in symbols}
    local['i'] = sp.I
    return sp.sympify(text.replace('^', '**'), locals=local)


def is_gaussian(expressions):
    return any(sp.expand(e).has(sp.I) for e in expressions)


def _preamble(variables, order, gaussian):
    names = ','.join(str(v) for v in variables)
    base = '(0,i)' if gaussian else '0'
    lines = [f'ring r = {base},({names}),{order};']
    if gaussian:
        lines.append('minpoly = i^2+1;')
    # without this Singular prints x2 for x^2, which will not parse back
    lines.append('short = 0;')
    return '\n'.join(lines)


def run(script, timeout=120):
    if not available():
        raise SingularError('the Singular binary is not on the PATH')
    try:
        proc = subprocess.run([BINARY, '-q', '--no-warn'],
                              input=script + '\nquit;\n',
                              capture_output=True, text=True, timeout=timeout)
    except subprocess.TimeoutExpired:
        raise SingularError(f'Singular did not finish within {timeout}s')
    out = proc.stdout
    if '? ' in out or proc.returncode:
        detail = '\n'.join(l for l in (out + proc.stderr).splitlines()
                           if l.strip().startswith('?'))
        raise SingularError(detail or proc.stderr[:300] or 'unknown failure')
    return out


def _emit(ideal_name):
    """Print one generator per line, tagged so wrapped output can be rejoined."""
    # NB: Singular rejects identifiers with leading underscores.
    return (f'int jdx;\n'
            f'for (jdx=1; jdx<=size({ideal_name}); jdx++) '
            f'{{ "<<", string({ideal_name}[jdx]); }}')


def _collect(out):
    """Rebuild the tagged generators, which Singular may wrap across lines."""
    polys, current = [], None
    for line in out.splitlines():
        line = line.rstrip('\\')
        if line.startswith('<< '):
            if current is not None:
                polys.append(current)
            current = line[3:].strip()
        elif current is not None and line.strip():
            current += line.strip()
    if current is not None:
        polys.append(current)
    return [p for p in polys if p]


def vdim(equations, variables, timeout=120):
    """Number of solutions with multiplicity, or None if positive-dimensional."""
    gaussian = is_gaussian(equations)
    polys = ',\n  '.join(to_text(e) for e in equations)
    script = (f'{_preamble(variables, "dp", gaussian)}\n'
              f'ideal I = {polys};\n'
              f'ideal G = std(I);\n'
              f'"<<", string(vdim(G));')
    value = int(_collect(run(script, timeout))[0])
    return None if value < 0 else value


def is_unit_ideal(equations, variables, timeout=120):
    """Whether the equations have no common solution at all."""
    gaussian = is_gaussian(equations)
    polys = ',\n  '.join(to_text(e) for e in equations)
    script = (f'{_preamble(variables, "dp", gaussian)}\n'
              f'ideal I = {polys};\n'
              f'"<<", string(dim(std(I)));')
    return int(_collect(run(script, timeout))[0]) == -1


def eliminate(equations, variables, drop, timeout=120):
    """Generators of the ideal got by eliminating `drop` from `equations`.

    This does not need the ideal to be zero-dimensional, which
    matters for the parametric family behind a matching path: there the base
    point of pi is a free variable.
    """
    gaussian = is_gaussian(equations)
    polys = ',\n  '.join(to_text(e) for e in equations)
    product = '*'.join(str(v) for v in drop)
    script = (f'{_preamble(variables, "dp", gaussian)}\n'
              f'ideal I = {polys};\n'
              f'ideal E = eliminate(std(I), {product});\n'
              f'{_emit("E")}')
    return [from_text(p, variables) for p in _collect(run(script, timeout))]
