"""Exact exceptional-parameter tests, including non-generic affine families.

Run with python -m tests.exceptional. Both available elimination backends are
checked against independently calculated polynomial conditions.
"""
import json
import time
from unittest.mock import patch

import sympy as sp

from handlebodies import algebra as A
from handlebodies.compute import answer
from handlebodies.exceptional import _calculate, _cached, t


def main():
    cases = [
        ('cubic cusp', 'w-s**3-t*s', 't'),
        ('quadratic', 'w-s**2-t', '1'),
        ('higher ramification', 'w-s**5-t*s', 't'),
        ('quartic surface', 's**4+w**4+s*w+t*s+2*w+1', 't**4-2*t+1'),
        ('quintic surface', 's**5+w**5+s*w+t*s+2*w+1', 't**5+2*t-1'),
        ('degree seven family', 'w-s**7-t*s**3-s', 't**3+7/4'),
        ('multiple parameters', 'w-s**3-(t**4-1)*s', 't**4-1'),
        ('complex coefficients', 'w-s**3-(t-I)*s', 't-I'),
        ('finite image of a curve', 'w-t*s', 't'),
        ('singular point excluded', 'w**2-s**3-t*s', '1'),
        ('singular elsewhere on fibre', 'w**2-1+4*s**3-3*s**4-t*s', 't-16/9'),
        ('degree drop is not a cusp', 't*s**3+s+w', '1'),
    ]
    for fast in [False, True] if A.backend.available() else [False]:
        start = time.monotonic()
        for name, p, expected in cases:
            report = _calculate(p, fast)
            assert report['status'] == 'complete', (name, report)
            actual = sp.Poly(report['polynomial'], t, extension=sp.I).monic()
            target = sp.Poly(sp.sympify(expected), t, extension=sp.I).monic()
            assert actual == target, (name, actual, target)
            assert len(report['values']) == target.degree(), name
            json.dumps(report, allow_nan=False)
        report = _calculate('w-s**4', fast)
        assert report['status'] == 'nonfinite', report
        print(f'{"Singular" if fast else "SymPy"}: {len(cases)+1} families in {time.monotonic()-start:.3f}s')

    # A non-coordinate projection, including rational linear changes of both
    # coordinates: t=x+z, w=y-z, s=z.
    spec = dict(domain='y-z-z**3-(x+z)*z',
                pi=dict(expression='x+z', solvefor='x'),
                rho=dict(expression='y-z', solvefor='y'))
    report = answer('/api/exceptional', spec)
    assert report['status'] == 'complete' and report['values'] == [[0., 0.]], report
    json.dumps(report, allow_nan=False)
    # Origins do not enter the geometric cache, and callers cannot mutate it.
    report['values'][0][0] = 99
    spec['pi']['origin'] = 2
    spec['rho']['origin'] = 3
    with patch('handlebodies.exceptional.subprocess.Popen', side_effect=AssertionError('cache miss')):
        report = answer('/api/exceptional', spec)
    assert report['values'] == [[0., 0.]], report

    start = time.monotonic()
    report = _cached('w-s**3-t*s', False, 0.001)
    assert report['status'] == 'timeout', report
    assert time.monotonic()-start < 3
    assert _calculate('w-s**3-t*s', False)['status'] == 'complete'
    print('API, geometry cache, mutation isolation and worker timeout: ok')


if __name__ == '__main__':
    main()
