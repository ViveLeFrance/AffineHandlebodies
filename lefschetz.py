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


        w = var('w', domain=CC)
        variables = [variable for variable in self.variables]
        variables.append(w)
        R = PolynomialRing(CC, names=variables)

        domain_hom = SR(R(self.domain).homogenize(var='w'))
        fibre_hom = SR(R(self.fibration).homogenize(var='w'))

        constraints = [w==0, fibre_hom==0, domain_hom==0]

        intersection = solve(constraints, variables, solution_dict=True)

        solution = set_free_variable_to_one_list(intersection)

        for sol in solution:
            if all(sol[variable].is_zero() for variable in variables):
                solution.remove(sol)

        return solution

    
    def get_homogenized_domain(self, variable=None):
        """Returns the projectivization of the domain.
        The projectivization variable can be specified, default is w."""
        
        variables = [variable for variable in self.variables]
        
        if variable is None:
            w = var('w', domain=CC)
            variable = w
            variables.append(w)
        else:
            variables.append(variable)

        R = PolynomialRing(CC, names=variables)

        domain_hom = SR(R(self.domain).homogenize(var=str(variable)))

        return domain_hom

    def transversality_at_infinity(self, origin_fibre=0, variable=None):
        """Checks whether the fibration is transverse to the hyperplane at infinity.
        If F is a defining polynomial for the projectivized domain, and pi is the (linear)
        fibration, this function computes the rank of the matrix
                        [grad F | grad pi | 0 0 0 1].
        If the rank is 3, then pi^{-1}(0) intersects the hyperplane at infinity w^{-1}(0)
        transversely. Seidel's result says than the rational function pi/w restricted to 
        {w =/= 0}, i.e., pi, is a Lefschetz fibration if {F = 0} is smooth.
        """
        
        intersection = self.get_fibre_boundary_components(point=origin_fibre, variable=variable)

        if len(intersection) == 0:
            print('The fibration does not vanish at infinity.')
        
        else:
            w = var('w', domain=CC)
            variables = [variable for variable in self.variables]
            variables.append(w)
            R = PolynomialRing(CC, names=variables)

            domain_hom = SR(R(self.domain).homogenize(var='w'))
            fibre_hom = SR(R(self.fibration).homogenize(var='w'))

            gradF = [domain_hom.diff(variable) for variable in variables]
            gradpi = [fibre_hom.diff(variable) for variable in variables]

            n = len(variables)-1

            for point in intersection:
                row1 = [pdev.subs(point) for pdev in gradF]
                row2 = [pdev.subs(point) for pdev in gradpi]
                row3 = [0 for index in range(n)]
                row3.append(1)
                M = matrix(CC, [row1, row2, row3])
                print(f'The rank of M at {point} is {M.rank()}.')

    
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
    

def homogenize(expression, variable='w', ignore=None):
    """Homogenizes a polynomial expression by adding a variable w."""
    variables = list(expression.variables())
    if variable not in variables:
        variables.append(variable)
        for var in ignore:
            variables.remove(var)
    else:
        raise ValueError(f"The variable '{variable}' is already in the expression.")
    S = PolynomialRing(CC, names=ignore)
    R = PolynomialRing(S, names=variables)
    homogenized_expr = SR(R(expression).homogenize(var=variable))
    
    return homogenized_expr

def singular_points(expression):
    """Checks whether a polynomial expression defines a smooth affine variety."""
    variables = expression.variables()
    constraints = [expression.diff(var) == 0 for var in variables]
    solutions = solve(constraints, variables, solution_dict=True)
    singularities = [point for point in solutions if expression.subs(point).is_zero()]
    return singularities


def parameterized_fib_crits(fib: LefschetzFibration, fib_param_path: Dict[str, List[complex]]):
    """Compute the critical values of a Lefschetz fibration with parameters specified by 
    a path of complex numbers. """

    if (len(fib_param_path['a']) != len(fib_param_path['b'])):
        raise ValueError("The paths must have the same number of points.")
    
    n = len(fib_param_path['a'])
    crits = {}

    if len(fib.variables) == 2:
        r1, r2 = var('r1, r2', domain=CC)
        domain = fib.domain.subs(fib.variables[0]==r1, fib.variables[1]==r2)
        a_path = fib_param_path['a']
        b_path = fib_param_path['b']
    
        for i in range(n):
            fib_eq = CC(a_path[i])*r1 + CC(b_path[i])*r2
            fib_i = LefschetzFibration([r1,r2], domain, fib_eq)

            crits[i] = fib_i.get_critical_values()
    
    if len(fib.variables) == 3: 
        r1, r2, r3 = var('r1, r2, r3', domain=CC)
        domain = fib.domain.subs(fib.variables[0]==r1, fib.variables[1]==r2, fib.variables[2]==r3)
        a_path = fib_param_path['a']
        b_path = fib_param_path['b']
        c_path = fib_param_path['c']
    
        for i in range(n):
            fib_eq = CC(a_path[i])*r1 + CC(b_path[i])*r2 + CC(c_path[i])*r3
            fib_i = LefschetzFibration([r1,r2,r3], domain, fib_eq)

            crits[i] = fib_i.get_critical_values()   

    return crits

def trace_preimage(rho: LefschetzFibration, path: List[complex], origin_fibre, title=None, solvefor=None, anticlockwise=True):
    """Used to determine the vanishing cycles of the 'small' Lefschetz fibration rho by solving for its fibres
    over a path. For this to work, fibres must be discrete."""
    if solvefor is None:
        solvefor = rho.variables[0]

    t = var('t', domain=CC)

    fibre_rho_t = rho.get_fibre(t, solvefor)

    fibres = []

    for s in path:
        fibre_rho_s = fibre_rho_t.subs(t=s)
        fibres.append(NumericalRoots(fibre_rho_s))

    # fibres = np.array(fibres)

    plot_points_real = []
    plot_points_imag = []
    for preimage in fibres:
        plot_points_real.append([value.real for value in preimage])
        plot_points_imag.append([value.imag for value in preimage])

    fig, ax = plt.subplots()    

    ax.spines['left'].set_position(('data', 0))
    # Move bottom spine to y=0
    ax.spines['bottom'].set_position(('data', 0))

    # Remove the top and right spines
    ax.spines['right'].set_color('none')
    ax.spines['top'].set_color('none')

    ax.set_title(title)

    init_sols = NumericalRoots(fibre_rho_t.subs(t==path[0]))
    init_sols = sort_by_angle(init_sols, origin_fibre=origin_fibre, anticlockwise=anticlockwise)


    init_real = [value.real for value in init_sols]
    init_imag = [value.imag for value in init_sols]
    
    final_sols = NumericalRoots(fibre_rho_t.subs(t==path[-1]))
    final_real = [value.real for value in final_sols]
    final_imag = [value.imag for value in final_sols]

    regular_sols = NumericalRoots(fibre_rho_t.subs(t==path[len(path)//2]))
    regular_real = [value.real for value in regular_sols]
    regular_imag = [value.imag for value in regular_sols]


    for index, point in enumerate(init_sols):
        ax.text(point.real+0.05, point.imag, str(index), fontsize=12, color='blue')
    ax.plot(init_real, init_imag, 'ro', markersize=5)
    
    # ax.plot(init_real, init_imag, 'ro', markersize=5)
    ax.plot(final_real, final_imag, 'co', markersize=5)
    ax.plot(regular_real, regular_imag, 'o', color='purple', markersize=5)


    ax.grid(True)

    x = np.hstack(plot_points_real)    
    y = np.hstack(plot_points_imag)

    N= x.size     

    cmap = plt.get_cmap('viridis')
    gradient = np.arange(N)
    norm = Normalize(vmin=min(gradient), vmax=max(gradient))
    colors = cmap(norm(gradient))

    ax.scatter(plot_points_real, plot_points_imag, color = colors,  s=2)

    return fig, ax