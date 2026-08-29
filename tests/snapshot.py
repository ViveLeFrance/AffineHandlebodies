"""Snapshot every example's computed results.

Written while the notebooks were still present and verified against them, so
this file remains a valid regression baseline after the notebooks are removed.
Regenerate deliberately with `python -m tests.snapshot --write`.
"""
import json, pathlib, sys

import numpy as np

import handlebodies as hb
from examples import ALL

BASELINE = pathlib.Path(__file__).parent / 'baseline.json'


def _row(m):
    """One traced path: the colliding pair, or the permutation if it is closed."""
    return {'path': m['name'], 'closed': m['closed'],
            'result': list(m['monodromy']) if m['closed']
                      else list(m['collision'][:2])}


def compute():
    out = {}
    for name in ALL:
        s = hb.load(name)
        summary = s.summary()
        summary.pop('name', None)
        out[name] = {
            'summary': summary,
            'pi_values': [[v.real, v.imag] for v in s.cv.pi],
            'rho_values': [[v.real, v.imag] for v in s.cv.rho],
            # lists, not tuples: these round-trip through JSON unchanged
            'matching': [_row(m) for m in s.matching_paths()],
            'vanishing': [_row(m) for m in s.vanishing_cycles()],
        }
    return out


def compare(got, want, tol=1e-6):
    problems = []
    for name in want:
        if name not in got:
            problems.append(f'{name}: missing'); continue
        g, w = got[name], want[name]
        for k, v in w['summary'].items():
            if g['summary'].get(k) != v:
                problems.append(f'{name}.{k}: {g["summary"].get(k)!r} != {v!r}')
        for k in ('pi_values', 'rho_values'):
            a = np.array([complex(*p) for p in g[k]])
            b = np.array([complex(*p) for p in w[k]])
            if a.shape != b.shape:
                problems.append(f'{name}.{k}: {a.shape} != {b.shape}')
            elif np.abs(a - b).max() > tol:
                problems.append(f'{name}.{k}: max deviation {np.abs(a-b).max():.2e}')
        for k in ('matching', 'vanishing'):
            if g[k] != w[k]:
                for x, y in zip(g[k], w[k]):
                    if x != y:
                        problems.append(f'{name}.{k}[{x["path"]}]: '
                                        f'{x["result"]} != {y["result"]}')
    return problems


if __name__ == '__main__':
    got = compute()
    if '--write' in sys.argv:
        BASELINE.write_text(json.dumps(got, indent=1) + '\n')
        print(f'wrote {BASELINE} ({len(got)} examples)')
    else:
        problems = compare(got, json.loads(BASELINE.read_text()))
        for p in problems:
            print('FAIL', p)
        print(f'{len(problems)} problems' if problems
              else f'all {len(got)} examples match the baseline')
        sys.exit(1 if problems else 0)
