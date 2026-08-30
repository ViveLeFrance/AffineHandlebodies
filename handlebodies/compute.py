"""Turn a specification into figures, for the interactive page.

Stateless: the page sends the whole specification each time and gets back
whatever can be computed from it, so there is no session to keep in step.
The pi half is computed alone when rho has not been filled in yet.
"""

import sympy as sp

from .study import Study, _slug


def _fig(f):
    return f.to_plotly_json()


def _err(stage, exc):
    return {'stage': stage, 'error': f'{type(exc).__name__}: {exc}'}


def surface_and_pi(spec, paths=None):
    """Everything that depends only on the surface and pi.

    `paths` carries the extra vanishing paths an example file declares; the
    page itself has no way to express those, so editing a field in the browser
    falls back to the straight-line defaults.
    """
    from . import web as R
    pi_paths, rho_paths = paths or ((), ())
    study = Study(name=spec.get('name', 'surface'),
                  domain=spec['domain'], params=spec.get('params') or None,
                  pi=spec['pi'], rho=spec.get('rho') or dict(expression='0'),
                  pi_paths=pi_paths, rho_paths=rho_paths)
    pi = study.pi
    cv = pi.critical_values()
    ranks = pi.transversality()
    return study, {
        'domain': str(study.domain),
        'pi': str(pi.expression),
        'origin': str(pi.origin),
        'fibre': str(pi.fibre(pi.origin)),
        'tex': _tex_block(pi, r'\pi', study.variables, pi.fibre(pi.origin)),
        'n_values': len(cv),
        'n_points': pi.n_critical_points(),
        'punctures': len(ranks),
        'origin_xy': [pi.origin_c.real, pi.origin_c.imag],
        'ranks': [r for _, r in ranks],
        'values': [[v.real, v.imag] for v in cv],
        'figure': _fig(R.points(cv, title='Critical values of pi',
                                origin=pi.origin_c,
                                origin_name=f'reference fibre {pi.origin}')),
    }


def rho_section(study, numeric=False, verify=False):
    """Everything that needs rho: its critical values, its fibre, the traces.

    `numeric` switches the matching paths to the path tracker.  It is the whole
    reason the flag reaches this far: on a high-degree surface the exact route
    would be attempted here, on every matching path, before the page had drawn
    anything at all.  The vanishing cycles are left alone -- those are numeric
    root-finding already and stay quick whatever the degree.
    """
    from . import web as R
    rho = study.rho
    cv = rho.critical_values()
    out = {
        'rho': str(rho.expression),
        'origin': str(rho.origin),
        'fibre': str(rho.fibre(rho.origin)),
        'tex': _tex_block(rho, r'\rho', rho.variables, rho.fibre(rho.origin)),
        'n_values': len(cv),
        'n_points': rho.n_critical_points(),
        'values': [[v.real, v.imag] for v in cv],
        'origin_xy': [rho.origin_c.real, rho.origin_c.imag],
        'critical_figure': _fig(R.points(
            cv, title='Critical values of rho', origin=rho.origin_c,
            origin_name=f'reference fibre {rho.origin}')),
        'fibre_figure': _fig(R.points(
            study.reference_fibre_points(),
            title=f'Generic fibre of rho over {rho.origin}',
            name='fibre points')),
    }

    def described(m):
        if m['closed']:
            return f"monodromy {m['monodromy']}"
        i, j, gap = m['collision']
        return f"points {i} and {j} collide (gap {gap:.1e})"

    n_rho = len(cv)
    out['vanishing'] = [{
        'name': m['name'], 'key': _slug(m['name']),
        'target': i if i < n_rho else None,
        'result': described(m), 'closed': m['closed'],
        'figure': _fig(R.trace(
            m['branches'], title=f"rho: {m['name']} -- {described(m)}",
            highlight=None if m['closed'] else m['collision'][:2])),
    } for i, m in enumerate(study.vanishing_cycles())]

    n_pi = len(study.cv.pi)
    out['matching'] = [{
        'name': m['name'], 'key': _slug(m['name']),
        'target': i if i < n_pi else None,
        'result': (m['error'] if m.get('error')
                   else described(m).replace('points', 'branches')),
        'failed': bool(m.get('error')),
        'closed': m['closed'],
        'figure': None if m.get('error') else _fig(R.trace3d(
            m['branches'],
            title=f"pi: {m['name']} -- {described(m).replace('points','branches')}",
            highlight=None if m['closed'] else m['collision'][:2],
            origin=rho.origin_c)),
    } for i, m in enumerate(study.matching_paths(numeric=numeric,
                                                 verify=verify))]
    return out


def compute(spec, paths=None):
    """Compute as much as the specification allows; report failures per stage.

    `spec['numeric']` asks for the matching paths to be tracked rather than
    solved afresh at every step.
    """
    result = {'ok': True}
    try:
        study, pi_out = surface_and_pi(spec, paths)
        result['pi'] = pi_out
    except Exception as exc:
        result['ok'] = False
        result['pi'] = _err('pi', exc)
        return result

    try:
        result['smooth'] = study.is_smooth
    except Exception:
        result['smooth'] = None

    if not spec.get('rho'):
        return result
    try:
        result['rho'] = rho_section(study, bool(spec.get('numeric')),
                                    bool(spec.get('verify')))
    except Exception as exc:
        result['ok'] = False
        result['rho'] = _err('rho', exc)
    return result


def _tex_block(fibration, symbol, variables, fibre_expr):
    """The identifying equations for a fibration, as aligned LaTeX.

    Rendered as one `aligned` environment so the three lines line up on their
    equals signs.
    """
    import sympy as sp

    args = ','.join(sp.latex(v) for v in variables)
    sub = '0' if symbol == r'\pi' else '1'
    pt = 'z_0' if symbol == r'\pi' else 'w_0'
    lines = [
        rf'{symbol}({args}) &= {sp.latex(fibration.expression)}',
        rf'{pt} &= {sp.latex(fibration.origin)}',
        rf'{symbol}^{{-1}}({pt}) &= \left\{{\, {sp.latex(fibre_expr)} = 0 \,\right\}}',
    ]
    return r'\begin{aligned}' + r' \\[2pt] '.join(lines) + r'\end{aligned}'


def _described(m):
    if m['closed']:
        return f"monodromy {m['monodromy']}"
    i, j, gap = m['collision']
    return f"points {i} and {j} collide (gap {gap:.1e})"


def retrace(spec, waypoints, steps=70, rounding=0.0, label='edited path',
            level='rho', numeric=False, verify=False):
    """Trace along a path the user has drawn.

    At the rho level this follows the points of rho's fibre; at the pi level it
    follows rho's critical values over the fibres of pi, which is the matching
    path.  `waypoints` are the corners as [re, im] pairs; the rest of the
    specification is the one the page already computed with.
    """
    from . import web as R
    from .fibration import collision, is_closed, monodromy
    from .paths import rounded_path

    study = Study(name=spec.get('name', 'surface'),
                  domain=spec['domain'], params=spec.get('params') or None,
                  pi=spec['pi'], rho=spec['rho'])
    # rounding is one value per waypoint (or a single number for all of them)
    path = rounded_path([complex(a, b) for a, b in waypoints],
                        steps=int(steps), rounding=rounding)

    if level == 'pi':
        branches = study.bifibration.matching_path(path, numeric=numeric,
                                                   verify=verify)
    else:
        branches = study.rho.trace_fibre(path)

    closed = is_closed(path)
    m = {'name': label, 'closed': closed,
         'monodromy': monodromy(branches) if closed else None,
         'collision': None if closed else collision(branches)}
    text = _described(m)
    highlight = None if closed else m['collision'][:2]
    if level == 'pi':
        text = text.replace('points', 'branches')
        figure = R.trace3d(branches, title=f'{label} -- {text}',
                           highlight=highlight, origin=study.rho.origin_c)
    else:
        figure = R.trace(branches, title=f'{label} -- {text}',
                         highlight=highlight)
    return {
        'ok': True,
        'result': text,
        'closed': closed,
        'steps': len(path),
        'figure': _fig(figure),
        'path': [[float(p.real), float(p.imag)] for p in path],
    }


def sample_path(waypoints, steps=70, rounding=0.0):
    """Just the sampled path -- no tracing, so this is cheap enough to preview."""
    from .paths import rounded_path
    path = rounded_path([complex(a, b) for a, b in waypoints],
                        steps=int(steps), rounding=rounding)
    return {'ok': True, 'path': [[float(p.real), float(p.imag)] for p in path]}
