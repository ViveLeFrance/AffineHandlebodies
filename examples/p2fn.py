from handlebodies import pl_path, run

NAME = "Painlevé II (Flaschka-Newell)"
DOMAIN = "x*y*z + x - y + z + s1"
PARAMS = {"s1": 1}
PI  = dict(expression="3*x - y - 2*z", solvefor="y")
RHO = dict(expression="3*x + z/2", solvefor="x")
RHO_PATHS = [
    ("matching 0-2", lambda cv: pl_path([cv.rho[0], -10j, cv.rho[2]], steps=70)),
]

run(__name__)
