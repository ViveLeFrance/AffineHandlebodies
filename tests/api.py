"""Check that every answer the page can be given is one it can be sent.

    python -m tests.api

The page's endpoints are the only things that cross into JSON, and a result
that cannot be encoded is a request that dies rather than one that reports a
failure.  plotly is the way that happens: it encodes the arrays it recognises
and hands back the rest as numpy, so a fibration with no critical values --
a plane, say -- used to draw an empty array and take the request down with it.

Failures are answers too, so a nonsense surface has to encode as cleanly as a
good one; nothing here asserts that a computation succeeded.
"""

import json
import sys

from handlebodies.compute import answer

PI = {'expression': '2*x + y - z', 'solvefor': 'x', 'origin': '0'}
RHO = {'expression': 'z - 3*y', 'solvefor': 'y', 'origin': '0'}

#: surfaces chosen for the shapes of answer they provoke, not for interest
SURFACES = [
    ('a worked one', 'x*y*z - x - y + 1'),
    ('no critical values', 'x + y + z'),
    ('no critical values, again', 'x*y - 1'),
    ('singular at a point', 'x**2 + y**2 + z**2'),
    ('singular along a curve', 'x*y*z'),
    ('higher degree', 'x*y*z + x**4 + y**4 - 1'),
    ('not a surface at all', 'nonsense +'),
    ('empty', ''),
]


def encodes(label, path, body):
    """Whether the page ends up with JSON, the way both back ends arrange it.

    Raising is a legitimate outcome -- app.py and the browser both turn it into
    a reported failure -- so what is checked is the whole contract: whatever
    comes of the call, it can be sent.
    """
    try:
        result = answer(path, dict(body))
    except Exception as exc:
        result = {'ok': False, 'error': f'{type(exc).__name__}: {exc}'}
    try:
        json.dumps(result)
    except Exception as exc:
        return f'{label}: {path} gave something unencodable -- {exc}'
    return None


def main():
    problems = []
    for label, domain in SURFACES:
        spec = {'domain': domain, 'params': None, 'pi': dict(PI), 'rho': dict(RHO)}
        here = [encodes(label, path, body) for path, body in (
            ('/api/surface', {'domain': domain, 'params': None}),
            ('/api/compute', dict(spec)),
            ('/api/compute', dict(spec, matching=True)),
            ('/api/trace', {'spec': spec, 'waypoints': [[0, 0], [2, 1]], 'steps': 8}))]
        print(f'  {label:28} {"ok" if not any(here) else "FAIL"}')
        problems += here
    problems.append(encodes('a path', '/api/path',
                            {'waypoints': [[0, 0], [1, 1]], 'steps': 6}))
    problems = [p for p in problems if p]
    print()
    for line in problems:
        print('FAIL', line)
    print(f'{len(problems)} problems' if problems
          else 'every endpoint answers in JSON, including when it fails')
    return 1 if problems else 0


if __name__ == '__main__':
    sys.exit(main())
