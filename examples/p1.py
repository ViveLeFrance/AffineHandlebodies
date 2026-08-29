from handlebodies import circle, pl_path, run

NAME = "Painlevé I"
DOMAIN = "x*y*z - x - y + 1"
PI  = dict(expression="5*x + y*(-5 + I) + z", solvefor="z")
RHO = dict(expression="x + 2*y", solvefor="x")
RHO_PATHS = [
    ("trivial circle", lambda cv: circle(3, centre=-1, steps=200)),
]

run(__name__)
