"""Turn a specification into figures, for the interactive page.

Stateless: the page sends the whole specification each time and gets back
whatever can be computed from it, so there is no session to keep in step.
The pi half is computed alone when rho has not been filled in yet.
"""

import json

from .study import Study, _slug


def _fig(f):
    return f.to_plotly_json()


def _err(stage, exc):
    return {'stage': stage, 'error': f'{type(exc).__name__}: {exc}'}


#: the few most recent studies, so pressing a second button does not redo the
#: first one's work.  A Study caches what it computes, and pi's critical values
#: alone can take the best part of a minute on a high-degree surface.
_RECENT = {}
_RECENT_LIMIT = 4


def study_for(spec, paths=None):
    """The Study for this specification, reused if it has been built before."""
    pi_paths, rho_paths = paths or ((), ())
    key = json.dumps({'domain': spec['domain'],
                      'params': spec.get('params') or None,
                      'pi': spec['pi'], 'rho': spec.get('rho'),
                      'declared': bool(pi_paths or rho_paths)}, sort_keys=True)
    if key not in _RECENT:
        if len(_RECENT) >= _RECENT_LIMIT:
            _RECENT.pop(next(iter(_RECENT)))
        _RECENT[key] = Study(
            name=spec.get('name', 'surface'),
            domain=spec['domain'], params=spec.get('params') or None,
            pi=spec['pi'], rho=spec.get('rho') or dict(expression='0'),
            pi_paths=pi_paths, rho_paths=rho_paths)
    return _RECENT[key]


def surface_and_pi(spec, paths=None):
    """Everything that depends only on the surface and pi.

    `paths` carries the extra vanishing paths an example file declares; the
    page itself has no way to express those, so editing a field in the browser
    falls back to the straight-line defaults.
    """
    from . import web as R
    study = study_for(spec, paths)
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


def rho_section(study):
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
        if m.get('branches') is None:
            return m.get('note', 'shown as computed')
        if m['closed']:
            return _monodromy_text(m['monodromy'])
        i, j, gap = m['collision']
        return f"points {i} and {j} collide (gap {gap:.1e})"

    n_rho = len(cv)
    out['vanishing'] = [{
        'name': m['name'], 'key': _slug(m['name']),
        'target': i if i < n_rho else None,
        'result': described(m), 'closed': m['closed'],
        'monodromy': m['monodromy'],
        'figure': _fig(R.trace(
            m['branches'], title=f"rho: {m['name']} -- {described(m)}",
            highlight=None if m['closed'] else m['collision'][:2])),
    } for i, m in enumerate(study.vanishing_cycles())]

    return out


def matching_section(study, numeric=False, verify=False):
    """The matching paths alone.

    Split out from the rest because it is much the most expensive part: rho's
    own critical values, fibre and vanishing cycles are quick, and it should be
    possible to settle on a rho before paying for these.
    """
    from . import web as R
    n_pi = len(study.cv.pi)

    def described(m):
        if m.get('note'):
            return m['note']
        if m.get('branches') is None:
            return 'shown as computed'
        if m['closed']:
            return _monodromy_text(m['monodromy'])
        i, j, gap = m['collision']
        return f"points {i} and {j} collide (gap {gap:.1e})"

    return [_traced(m, i, n_pi, R.trace3d, study.rho.origin_c, described)
            for i, m in enumerate(
                study.matching_paths(numeric=numeric, verify=verify))]


def _traced(m, i, n_default, render, origin, described):
    """One traced path, as the page needs it.

    Three outcomes, and only the first is a failure: the values could not be
    computed at all; they were computed but not sortable into branches, in
    which case they are still drawn; or everything worked.
    """
    row = {'name': m['name'], 'key': _slug(m['name']),
           'target': i if i < n_default else None,
           'closed': m['closed'],
           'monodromy': m.get('monodromy'),
           'failed': bool(m.get('error')),
           'tracked': m.get('branches') is not None}
    if m.get('error'):
        row['result'] = m['error']
        row['figure'] = None
        return row
    text = described(m).replace('points', 'branches')
    row['result'] = text
    title = f"pi: {m['name']} -- {text}"
    zmax = m.get('partial') or 1.0
    if m.get('branches') is None:
        row['figure'] = _fig(render(m['frames'], title=title, loose=True,
                                    origin=origin, zmax=zmax))
    else:
        row['figure'] = _fig(render(
            m['branches'], title=title,
            highlight=None if m['collision'] is None else m['collision'][:2],
            origin=origin, zmax=zmax, returns=m.get('monodromy')))
    return row


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
        result['rho'] = rho_section(study)
    except Exception as exc:
        result['ok'] = False
        result['rho'] = _err('rho', exc)
        return result

    if not spec.get('matching'):
        return result
    try:
        result['matching'] = matching_section(study, bool(spec.get('numeric')),
                                              bool(spec.get('verify')))
    except Exception as exc:
        result['ok'] = False
        result['matching'] = _err('matching', exc)
    return result


def _tex_block(fibration, symbol, variables, fibre_expr):
    """The identifying equations for a fibration, as aligned LaTeX.

    Rendered as one `aligned` environment so the three lines line up on their
    equals signs.
    """
    import sympy as sp

    args = ','.join(sp.latex(v) for v in variables)
    pt = 'z_0' if symbol == r'\pi' else 'w_0'
    lines = [
        rf'{symbol}({args}) &= {sp.latex(fibration.expression)}',
        rf'{pt} &= {sp.latex(fibration.origin)}',
        rf'{symbol}^{{-1}}({pt}) &= \left\{{\, {sp.latex(fibre_expr)} = 0 \,\right\}}',
    ]
    return r'\begin{aligned}' + r' \\[2pt] '.join(lines) + r'\end{aligned}'


def _monodromy_text(perm):
    """How a closed path reads.

    A permutation that is not a bijection means the branches could not be
    matched up across the ends, which is a failure to report rather than a
    monodromy to quote.
    """
    if perm is None:
        return 'closed path; the branches could not be matched up at its end'
    return f'monodromy {perm}'


def _described(m):
    if m['closed']:
        return _monodromy_text(m['monodromy'])
    i, j, gap = m['collision']
    return f"points {i} and {j} collide (gap {gap:.1e})"


def retrace(spec, waypoints, steps=70, rounding=0.0, label='edited path',
            level='rho', numeric=False, verify=False, closed=False):
    """Trace along a path the user has drawn.

    At the rho level this follows the points of rho's fibre; at the pi level it
    follows rho's critical values over the fibres of pi, which is the matching
    path.  `waypoints` are the corners as [re, im] pairs; the rest of the
    specification is the one the page already computed with.
    """
    from . import web as R
    from .fibration import collision, is_closed, monodromy
    from .paths import rounded_path

    # Through study_for, so that redrawing a path does not pay again for the
    # critical values the page has already shown -- which is most of the cost
    # of a recompute on a high-degree surface.
    study = study_for(spec)
    # rounding is one value per waypoint (or a single number for all of them)
    path = rounded_path([complex(a, b) for a, b in waypoints],
                        steps=int(steps), rounding=rounding, closed=closed)

    if level == 'pi':
        frames, branches = study.bifibration.matching_path(
            path, numeric=numeric, verify=verify)
    else:
        frames = branches = study.rho.trace_fibre(path)

    closed = is_closed(path)
    # Branch indices name the values over the fibre the path starts from. Only
    # a path that starts at the reference fibre gets the indices the plots are
    # labelled with; report which it is rather than leaving the numbers to be
    # read against the wrong picture.
    based = _is_based(path, study, level)
    if branches is None:
        # The values are the result; sorting them into branches is not. Say so
        # and draw them as they came out.
        text = ('the values could not be followed as separate branches; '
                'shown as computed')
        figure = (R.trace3d(frames, title=f'{label} -- {text}', loose=True,
                            origin=study.rho.origin_c) if level == 'pi'
                  else R.trace(frames, title=f'{label} -- {text}', loose=True))
        return {'ok': True, 'result': text, 'closed': closed, 'tracked': False,
                'based': based, 'steps': len(path), 'figure': _fig(figure),
                'path': [[float(p.real), float(p.imag)] for p in path]}

    m = {'name': label, 'closed': closed,
         'monodromy': monodromy(branches) if closed else None,
         'collision': None if closed else collision(branches)}
    text = _described(m)
    highlight = None if closed else m['collision'][:2]
    if level == 'pi':
        text = text.replace('points', 'branches')
        figure = R.trace3d(branches, title=f'{label} -- {text}',
                           highlight=highlight, origin=study.rho.origin_c,
                           returns=m['monodromy'])
    else:
        figure = R.trace(branches, title=f'{label} -- {text}',
                         highlight=highlight)
    return {
        'ok': True,
        'result': text,
        'closed': closed,
        'monodromy': m['monodromy'],
        'based': based,
        'tracked': True,
        'steps': len(path),
        'figure': _fig(figure),
        'path': [[float(p.real), float(p.imag)] for p in path],
    }


def _is_based(path, study, level):
    """Whether a path starts at the reference fibre of the base it lives in.

    The page pins a based path there exactly, so this is an equality test with
    room only for the rounding the path picker's coordinates went through.
    """
    origin = study.pi.origin_c if level == 'pi' else study.rho.origin_c
    return bool(abs(complex(path[0]) - origin)
                <= 1e-9 * max(1.0, abs(origin)))


def sample_path(waypoints, steps=70, rounding=0.0, closed=False):
    """Just the sampled path -- no tracing, so this is cheap enough to preview."""
    from .paths import rounded_path
    path = rounded_path([complex(a, b) for a, b in waypoints],
                        steps=int(steps), rounding=rounding, closed=closed)
    return {'ok': True, 'path': [[float(p.real), float(p.imag)] for p in path]}


def answer(path, body, paths=((), ())):
    """Route one of the page's API calls to the function behind it.

    The local server and the built page both come through here, so the same
    request gets the same answer whether the Python is running on the machine
    or inside the browser.  `body` is the decoded request; what comes back is
    what gets sent as JSON.
    """
    if path == '/api/path':
        return sample_path(body['waypoints'], steps=body.get('steps', 70),
                           rounding=body.get('rounding', 0.0),
                           closed=bool(body.get('closed', False)))
    if path == '/api/trace':
        return retrace(body['spec'], body['waypoints'],
                       steps=body.get('steps', 70),
                       rounding=body.get('rounding', 0.0),
                       label=body.get('label', 'edited path'),
                       level=body.get('level', 'rho'),
                       numeric=bool(body.get('numeric', False)),
                       verify=bool(body.get('verify', False)),
                       closed=bool(body.get('closed', False)))
    if path == '/api/compute':
        # Only the untouched preset gets the file's own paths; once a field is
        # edited the spec no longer matches what those paths were written for.
        declared = paths if body.pop('preset', False) else None
        return compute(body, declared)
    raise LookupError(f'no endpoint {path}')
