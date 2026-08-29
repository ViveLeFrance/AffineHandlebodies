from handlebodies import circle, pl_path, run

NAME = "Painlevé III (D7)"
DOMAIN = "x*y*z + x**2 + y**2 + alpha*x + y"   # alpha in C*
PARAMS = {"alpha": 1}
PI  = dict(expression="-7*x/10 + 3*I*y/10 + z", solvefor="z")
RHO = dict(expression="x - I*y/20", solvefor="x")
NOTES = """Two different vanishing paths to critical value 4 of pi are traced and
compared.  They do not give the same matching path: the right-hand path agrees
with the straight-line one (branches 0-3), while the left-hand path gives
branches 1-3.  Reason for the comparison not recorded in the notebook."""
PI_PATHS = [
    ("to pi-crit 4, left", lambda cv: pl_path([-0.1j, -3-0.5j, cv.pi[4]], steps=100)),
    ("to pi-crit 4, right", lambda cv: pl_path([-0.1j,  2-0.1j, cv.pi[4]], steps=100)),
]
RHO_PATHS = [
    ("trivial circle", lambda cv: circle(2, steps=400)),
]

run(__name__)
