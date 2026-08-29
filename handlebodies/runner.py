"""Make an example file runnable: `python examples/p2.py` shows its results.

The example modules are configuration first -- the page imports them
without side effects -- so run() does nothing unless the module is being
executed directly.
"""

import argparse
import math
import sys


def _preset_from(module):
    """The page's starting values, taken from an example file."""
    import sympy as sp

    from .study import Study
    g = vars(module)
    study = Study(name=g.get('NAME', 'example'), domain=g['DOMAIN'],
                  params=g.get('PARAMS'), pi=dict(g['PI']), rho=dict(g['RHO']))
    params = ', '.join(f'{k}={v}' for k, v in (g.get('PARAMS') or {}).items())
    return {
        'name': study.name,
        'domain': g['DOMAIN'],
        'params': params,
        'params_dict': {k: str(v) for k, v in (g.get('PARAMS') or {}).items()},
        'pi': str(study.pi.expression),
        'pisolve': str(study.pi.solvefor),
        'piorigin': str(study.pi.origin),
        'rho': str(study.bifibration.rho_0.expression),
        'rhosolve': str(study.bifibration.rho_solvefor),
        'rhoorigin': str(study.bifibration.rho_origin),
        'auto': True,
    }


def run(name=None, argv=None):
    """Open the interactive page with this example's choices filled in.

    Put `run(__name__)` at the foot of an example file.  It returns at once
    unless that file is the one being executed.
    """
    import argparse
    import sys

    if name is not None and name != '__main__':
        return None
    module = sys.modules['__main__']

    ap = argparse.ArgumentParser(description='Open this example in the browser.')
    ap.add_argument('-p', '--port', type=int, default=8765)
    ap.add_argument('-n', '--no-browser', action='store_true')
    args = ap.parse_args(argv)

    from .app import serve
    g = vars(module)
    preset = _preset_from(module)
    print(f"{preset['name']}", flush=True)
    serve(preset, port=args.port, open_browser=not args.no_browser,
          paths=(g.get('PI_PATHS', ()), g.get('RHO_PATHS', ())))
