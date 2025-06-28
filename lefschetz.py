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

        # fibre = self.get_fibre(point, variable)

        w = var('w', domain=CC)
        variables = [variable for variable in self.variables]
        variables.append(w)
        R = PolynomialRing(CC, names=variables)

        domain_hom = SR(R(self.domain).homogenize(var='w'))
        fibre_hom = SR(R(self.fibration).homogenize(var='w'))

        constraints = [w==0, fibre_hom==0, domain_hom==0]

        intersection = solve(constraints, variables, solution_dict=True)

        return set_free_variable_to_one_list(intersection)


    
class Bifibration:

    def __init__(self, pi: LefschetzFibration, rho: LefschetzFibration) -> None:
        self.pi = pi
        self.rho = rho

    def get_matching_path(self, origin_fibre=None, target_fibre=None, steps=70, solvefor=None, path=None):
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
            
            for s in np.linspace(0, 1, steps):
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
    

# def intersection_at_infinity(f: LefschetzFibration):
#     """Determines the intersection of the projectivization of the homogenized domain of f with 
#     its zero set at infinity."""
#     w = var('w', domain=CC)
#     variables = [variable for variable in f.variables]
#     variables.append(w)
#     R = PolynomialRing(CC, names=variables)

#     G_hom = G_hom = SR(R(f.domain).homogenize(var='w'))
#     f_hom = SR(R(f.fibration).homogenize(var='w'))

#     constraints = [w==0, G_hom==0, f_hom==0]

#     intersection = solve(constraints, variables, solution_dict=True)

#     return set_free_variable_to_one_list(intersection)


# def intersection_summary(f: LefschetzFibration):
#     w = var('w', domain=CC)
#     variables = [variable for variable in f.variables]
#     variables.append(w)
#     R = PolynomialRing(CC, names=variables)

#     G_hom = SR(R(f.domain).homogenize(var='w'))
#     f_hom = SR(R(f.fibration).homogenize(var='w'))

#     hyp_at_inf_constraints = [G_hom==0, w==0]
#     f_vanishing = [f_hom==0]

#     constraints = [expression for expression in hyp_at_inf_constraints]
#     constraints.extend(f_vanishing)

#     print(f'The hyperplane at infinity is given by {G_hom.subs(w==0)} == 0.')
#     print(f'The fibration vanishes at {f_vanishing}.')
#     intersection = set_free_variable_to_one_list(solve(constraints, variables, solution_dict=True))

#     print(f'Their intersection consists of {intersection}.')


def kernels(f: LefschetzFibration, point: dict[Any, Any]):
    w = var('w', domain=CC)
    variables = [variable for variable in f.variables]
    variables.append(w)
    R = PolynomialRing(CC, names=variables)

    G_hom = SR(R(f.domain).homogenize(var='w'))
    f_hom = SR(R(f.fibration).homogenize(var='w'))
    print(f_hom)
    
    J_inf = jacobian([G_hom, w], variables).subs(point).transpose()
    J_f = jacobian([G_hom, f_hom], variables).subs(point).transpose()
    # print('*******Jacs*********')
    # print(J_inf)
    # print('---')
    # print(J_f)
    # print('****************')

    ker_inf = J_inf.kernel()
    ker_f = J_f.kernel()

    # print('------Kernels----------')
    # print(ker_inf)
    # print('******************')
    # print(ker_f)
    # print('----------------')

    intersection = ker_inf.intersection(ker_f)

    return intersection


def transversality_at_infinity(pi: LefschetzFibration, origin_fibre=0, variable=None):
    """Checks whether the fibration is transverse to the hyperplane at infinity."""
    
    intersection = pi.get_fibre_boundary_components(point=origin_fibre, variable=variable)

    if len(intersection) == 0:
        print('The fibration does not vanish at infinity.')
    
    else:
        w = var('w', domain=CC)
        variables = [variable for variable in pi.variables]
        variables.append(w)
        R = PolynomialRing(CC, names=variables)

        domain_hom = SR(R(pi.domain).homogenize(var='w'))
        fibre_hom = SR(R(pi.fibration).homogenize(var='w'))

        gradF = [domain_hom.diff(variable) for variable in variables]
        gradpi = [fibre_hom.diff(variable) for variable in variables]

        for point in intersection:
            row1 = [pdev.subs(point) for pdev in gradF]
            row2 = [pdev.subs(point) for pdev in gradpi]
            row3 = [0,0,0,1]

            M = matrix(CC, [row1, row2, row3])
            print(f'The rank of M at {point} is {M.rank()}.')