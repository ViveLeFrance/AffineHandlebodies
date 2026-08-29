from handlebodies import circle, pl_path, run

NAME = "Painlevé VI"
DOMAIN = ("x*y*z + x**2 + y**2 + z**2 - (a1*a4+a2*a3)*x - (a2*a4+a3*a1)*y "
          "- (a3*a4+a1*a2)*z + a1*a2*a3*a4 + a1**2 + a2**2 + a3**2 + a4**2 - 4")
PARAMS = {"a1": 1, "a2": 0, "a3": 0, "a4": 0}
NOTES = "The reference fibre of rho is taken at -15i, away from the critical values."
PI  = dict(expression="x/10 - I*y/2 - z/5", solvefor="z")
RHO = dict(expression="9*x/10 + y/2", origin="-15*I", solvefor="x")
RHO_PATHS = [
    ("trivial circle", lambda cv: circle(40, steps=200)),
    ("vanishing 2", lambda cv: pl_path([-15j, cv.rho[2]], steps=400)),
    ("vanishing 5", lambda cv: pl_path([-15j, -15j-2, cv.rho[5]], steps=400)),
]

run(__name__)
