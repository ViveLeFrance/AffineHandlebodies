"""A worked example: a surface, a bifibration on it, and the chosen paths."""

import importlib
from dataclasses import dataclass

import numpy as np
import sympy as sp

from . import algebra as A
from .fibration import (Bifibration, LinearFibration, collision, is_closed,
                        monodromy)


def _slug(text):
    return ''.join(c if c.isalnum() else '-' for c in text.lower()).strip('-')


@dataclass
class CritValues:
    """Labelled critical values, as referred to by the paths in a config."""
    pi: np.ndarray
    rho: np.ndarray


class Study:
    """Everything computed for one example.

    Results are computed lazily and cached, so a REPL session can ask for
    whatever it needs without paying for the rest.
    """

    def __init__(self, name, domain, params=None, pi=None, rho=None,
                 pi_paths=(), rho_paths=(), notes='', variables='x y z'):
        self.name = name
        self.variables = [sp.Symbol(s) for s in variables.split()]
        self.domain = A.poly(domain, params)
        self.domain_family = domain
        self.params = params or {}
        self.notes = notes
        self._pi_paths = list(pi_paths)
        self._rho_paths = list(rho_paths)

        # PI and RHO are written with the same key names; Bifibration takes
        # its own, so translate here rather than making the configs asymmetric.
        self.pi = LinearFibration(self.domain, self.variables, **pi)
        self.bifibration = Bifibration(
            self.pi,
            rho_expression=rho['expression'],
            rho_origin=rho.get('origin', 0),
            rho_solvefor=rho.get('solvefor'))
        self.rho = self.bifibration.rho_0
        self._cv = None

    # -- basic invariants -------------------------------------------------

    @property
    def singularities(self):
        return A.singular_points(self.domain, self.variables)

    @property
    def is_smooth(self):
        return self.singularities == []

    @property
    def cv(self):
        """Labelled critical values of both fibrations."""
        if self._cv is None:
            self._cv = CritValues(pi=self.pi.critical_values(),
                                  rho=self.rho.critical_values())
        return self._cv

    @property
    def punctures(self):
        return len(self.pi.boundary_points())

    def summary(self):
        ranks = self.pi.transversality()
        return {
            'name': self.name,
            'domain': str(self.domain),
            'smooth': self.is_smooth,
            'fibre': str(self.pi.fibre(self.pi.origin)),
            'pi_values': len(self.cv.pi),
            'pi_points': self.pi.n_critical_points(),
            'rho_fibre': str(self.rho.fibre(self.rho.origin)),
            'rho_values': len(self.cv.rho),
            'rho_points': self.rho.n_critical_points(),
            'punctures': self.punctures,
            'ranks': [r for _, r in ranks],
            'transverse': all(r == len(self.variables) for _, r in ranks),
        }

    # -- paths ------------------------------------------------------------

    def pi_paths(self):
        """Paths in the base of pi: the straight-line vanishing paths, then any
        extra ones the example declares."""
        from .paths import pl_path
        default = [(f'vanishing {i}',
                    (lambda i: lambda cv: pl_path([self.pi.origin_c, cv.pi[i]],
                                                  steps=70))(i))
                   for i in range(len(self.cv.pi))]
        return default + self._pi_paths

    def rho_paths(self):
        """Paths in the base of rho, likewise."""
        from .paths import pl_path
        default = [(f'vanishing {i}',
                    (lambda i: lambda cv: pl_path([self.rho.origin_c, cv.rho[i]],
                                                  steps=70))(i))
                   for i in range(len(self.cv.rho))]
        return default + self._rho_paths

    def matching_paths(self, only=None, numeric=False, verify=False):
        """Trace the critical values of rho along each path in the base of pi.

        Two branches meeting is the matching path of rho corresponding to the
        vanishing cycle of pi over that path.  `only` is a predicate on the
        path name, so a caller after one picture does not pay for the rest --
        including the elimination that the parametric family needs.
        """
        out = []
        for name, fn in self.pi_paths():
            if only is not None and not only(name):
                continue
            path = np.asarray(fn(self.cv))
            try:
                branches = self.bifibration.matching_path(
                    path, numeric=numeric, verify=verify)
            except Exception as exc:
                # One path that cannot be followed must not take the rest with
                # it.  A straight path that happens to run through another
                # critical value of pi is not a vanishing path at all, and on a
                # surface with several real critical values that is a normal
                # thing to meet; the others are still worth having.
                out.append({'name': name, 'path': path, 'branches': None,
                            'closed': False, 'monodromy': None,
                            'collision': None,
                            'error': f'{type(exc).__name__}: {exc}'})
                continue
            out.append(self._describe(name, path, branches))
        return out

    def vanishing_cycles(self, only=None):
        """Trace the points of rho's reference fibre along each path in rho's base.

        Two points meeting over a critical value of rho is the vanishing cycle
        of rho attached to that critical value.  `only` filters by path name.
        """
        out = []
        for name, fn in self.rho_paths():
            if only is not None and not only(name):
                continue
            path = np.asarray(fn(self.cv))
            branches = self.rho.trace_fibre(path)
            out.append(self._describe(name, path, branches))
        return out

    # -- rendering --------------------------------------------------------

    def reference_fibre_points(self):
        """The points of rho's fibre over its reference value, labelled."""
        from .fibration import sort_by_angle
        import numpy as _np
        u = self.rho.other_variables()[0]
        p = sp.Poly(self.rho.fibre(self.rho.origin), u)
        r = _np.roots([complex(c) for c in p.all_coeffs()])
        return sort_by_angle(r, self.rho.origin_c)

    @staticmethod
    def _describe(name, path, branches):
        """Report the permutation for a closed path, the collision otherwise.

        A closed path returns to where it started, so nothing collides: what it
        exhibits is the monodromy.  Only an open path running into a critical
        value produces a colliding pair.
        """
        closed = is_closed(path)
        return {'name': name, 'path': path, 'branches': branches,
                'closed': closed,
                'monodromy': monodromy(branches) if closed else None,
                'collision': None if closed else collision(branches)}


def load(name, package='examples'):
    """Load an example module and build its Study."""
    mod = importlib.import_module(f'{package}.{name}')
    return Study(
        name=getattr(mod, 'NAME', name),
        domain=mod.DOMAIN,
        params=getattr(mod, 'PARAMS', None),
        pi=mod.PI,
        rho=mod.RHO,
        pi_paths=getattr(mod, 'PI_PATHS', ()),
        rho_paths=getattr(mod, 'RHO_PATHS', ()),
        notes=getattr(mod, 'NOTES', ''),
    )
