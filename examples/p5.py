from handlebodies import pl_path, run

NAME = "Painlevé V"
DOMAIN = ("x*y*z + x**2 + y**2 - (s1 + s2*s3)*x - (s2 + s1*s3)*y - s3*z "
          "+ s3**2 + s1*s2*s3 + 1")       # s1, s2 in C, s3 in C*
PARAMS = {"s1": "I", "s2": "-I", "s3": 1}
PI  = dict(expression="3*x - y - 2*z", solvefor="y")
RHO = dict(expression="x + 3*z/2", solvefor="x")
RHO_PATHS = [
    ("matching 0-2", lambda cv: pl_path([cv.rho[0], -10j, cv.rho[2]], steps=70)),
]

run(__name__)
