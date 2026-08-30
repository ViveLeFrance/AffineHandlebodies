"""Check the numerical tracker against the exact route.

The tracker follows solutions instead of solving afresh, so its characteristic
failure is landing two branches on the same solution and reporting a matching
path that was never there.  The recorded examples are small enough for the
exact eliminant, so both routes can be run on all of them and compared: same
colliding pair, same values.  That is the evidence the tracker does not jump.
"""

import sys

import numpy as np

import handlebodies as hb
from examples import ALL
from handlebodies.fibration import collision, is_closed, monodromy

TOLERANCE = 1e-6


#: the exact route matches unordered root sets between consecutive samples, so
#: its answer is only as good as the sampling; the reference has to be refined
#: enough to have stopped changing
REFERENCE_SAMPLES = 1200


def refine(path, samples=REFERENCE_SAMPLES):
    """The same route, sampled more finely.

    Interpolates along the given samples rather than drawing a fresh straight
    line: several of the recorded paths turn corners, and a straight
    replacement would be a different path entirely.
    """
    path = np.asarray(path)
    here = np.linspace(0.0, 1.0, len(path))
    want = np.linspace(0.0, 1.0, max(samples, len(path)))
    return np.interp(want, here, path.real) + 1j * np.interp(want, here, path.imag)


def compare(name, tol=TOLERANCE):
    """Every default matching path of one example, both ways.

    The tracker is compared against the exact route on a *finely sampled*
    version of the same path.  Comparing against the exact route at the
    working sampling would be comparing against something that is itself
    aliased on a few of these paths.
    """
    study = hb.load(name)
    problems, checked = [], 0

    for label, fn in study.pi_paths():
        path = np.asarray(fn(study.cv))
        exact = study.bifibration.matching_path(refine(path))
        try:
            # verify=True so the last frame really is the endpoint; with it
            # off the tracker deliberately stops a hair short, and the ends
            # would differ for that reason rather than from any error.
            tracked = study.bifibration.matching_path(path, numeric=True,
                                                      verify=True)
        except Exception as exc:
            problems.append(f'{name}/{label}: tracking failed: '
                            f'{type(exc).__name__}: {exc}')
            continue

        checked += 1
        if tracked.shape[1] != exact.shape[1]:
            problems.append(f'{name}/{label}: {tracked.shape[1]} branches '
                            f'vs {exact.shape[1]}')
            continue

        # The two runs sample the path differently, so only the ends line up.
        # The endpoint is judged loosely on purpose: it is where two values
        # coincide, and a double root is only ever resolved to about the square
        # root of machine precision by either route.  The colliding pair below
        # is the check that matters there.
        for what, a, b, bound in (('start', tracked[0], exact[0], tol),
                                  ('end', tracked[-1], exact[-1], 1e-3)):
            deviation = float(np.abs(a - b).max())
            if deviation > bound * max(1.0, float(np.abs(b).max())):
                problems.append(f'{name}/{label}: {what} values differ by '
                                f'{deviation:.2e}')

        if is_closed(path):
            if monodromy(tracked) != monodromy(exact):
                problems.append(f'{name}/{label}: monodromy '
                                f'{monodromy(tracked)} vs {monodromy(exact)}')
        else:
            a, b = collision(tracked)[:2], collision(exact)[:2]
            if a != b:
                problems.append(f'{name}/{label}: collides {a} vs {b}')
    return checked, problems


def main():
    total, failures = 0, []
    for name in ALL:
        checked, problems = compare(name)
        total += checked
        failures.extend(problems)
        print(f'  {name:6} {checked:2} matching paths '
              f'{"ok" if not problems else "FAIL"}')
    print()
    for line in failures:
        print('FAIL', line)
    print(f'{len(failures)} problems' if failures
          else f'tracker agrees with the exact route on all {total} paths')
    return 1 if failures else 0


if __name__ == '__main__':
    sys.exit(main())
