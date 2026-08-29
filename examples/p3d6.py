from handlebodies import pl_path, run

NAME = "Painlevé III (D6)"
DOMAIN = "x*y*z + x**2 + y**2 + (1+alpha*beta)*x + (alpha+beta)*y + alpha*beta"
PARAMS = {"alpha": -1, "beta": 1}
PI  = dict(expression="-7*x/10 + 3*I*y/10 + z", solvefor="z")
RHO = dict(expression="x - I*y/5", solvefor="x")
RHO_PATHS = [
    ("matching 0-2", lambda cv: pl_path([cv.rho[0], -10j, cv.rho[2]], steps=70)),
]

run(__name__)
