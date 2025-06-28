from sage.all import *
import numpy as np
from utility import *

class LefschetzFibration:

    def __init__(self, variables, domain_equation, fibration_equation) -> None:
        self.variables = variables
        self.domain = domain_equation
        self.fibration = fibration_equation

    def __call__(self, argument):
        return self.fibration.subs(argument)

    def get_critical_points(self):
        """Solves for when the gradient of the domain
        is parallel to the differential of the fibration."""
        G = self.domain
        f = self.fibration
        a = var('a', domain=CC) # additional variable to solve for parallelity


        constraints = [G==0]
        constraints.extend([G.diff(variable) == a*f.diff(variable) for variable in self.variables])

        points = solve(constraints, self.variables + [a], solution_dict=True)

        for p in points:
            del p[a]

        return points

       
    def get_critical_values(self):
        """Evaluates fibration at critical points."""
        crit_points = self.get_critical_points()
        return list(set([self.__call__(x) for x in crit_points]))

    def get_fibre(self, point, variable=None):
        """Solves for the fibre over a given point in the specified variable.
        If no variable is specified, it will solve for the first variable in the list.
        Note that we assume the fibration to be linear in all variables."""

        if variable is None:
            variable = self.variables[0]

        f = self.fibration
        G = self.domain

        fib_solution = solve(f == point, variable) # f expected to be linear, so we get a single chart

        return G.subs(fib_solution).simplify()


    def get_hessian(self, point):
        """Computes complex hessian to check nondegeneracy. Generic linear sections will be nondegenerate."""
        pass

    def get_fibre_boundary_components(self, point, variable=None):
        """Determines the number of boundary components of the fibre over a specified
        point by computing the number of points at the hyperplane at infinity in
        the projectivization of the fibre."""

        if variable is None:
            variable = self.variables[0]

        fibre = self.get_fibre(point, variable)

        w = var('w', domain=CC)
        variables = [variable for variable in self.variables]
        variables.remove(variable)
        variables.append(w)
        R = PolynomialRing(CC, names=variables)

        fibre_hom = SR(R(fibre).homogenize(var='w'))

        constraints = [w==0, fibre_hom==0]

        intersection = solve(constraints, variables, solution_dict=True)

        return set_free_variable_to_one_list(intersection)


    
class Bifibration:

    def __init__(self, pi: LefschetzFibration, rho: LefschetzFibration) -> None:
        self.pi = pi
        self.rho = rho

    def get_matching_path(self, rho_eq, origin_fibre=None, target_fibre=None, steps=70, solvefor=None, path=None):
        """For a vanishing path of pi, the critical values of rho restricted to the fibres
        over the vanishin path are such that two of them will converge, yielding a matching path
        of rho. This function computes the trace of these critical values."""

        if (origin_fibre is None or target_fibre is None) and path is None:
            raise ValueError("Please provide a path or origin and target fibres.")

        t = var('t', domain=CC)

        if solvefor is None:
            solvefor = self.pi.variables[0]

        fibre_t = self.pi.get_fibre(t, solvefor)

        rho_vars = self.rho.variables

        rho_eq_t = self.rho.fibration.subs(solvefor==fibre_t)
        
        matching_path = {}
        if not path:
            
            for s in np.linspace(0,1,steps):
                fibre_s = fibre_t.subs(t==(1-s)*origin_fibre + s*target_fibre)
                rho_eq_s = rho_eq_t.subs(t==(1-s)*origin_fibre + s*target_fibre)
                rho_s = LefschetzFibration(rho_vars, fibre_s, rho_eq_s)
                matching_path[s] = rho_s.get_critical_values()
        else:
            n = len(path)
            for index, point in enumerate(path):
                fibre_s = fibre_t.subs(t==point)
                rho_eq_s = rho_eq_t.subs(t==point)
                rho_s = LefschetzFibration(rho_vars, fibre_s, rho_eq_s)
                matching_path[index/(n-1)] = rho_s.get_critical_values()

        return matching_path