from handlebodies import pl_path, run

NAME = "Painlevé II"
NOTES = "The surface is singular at a = -1, so a = -2 is chosen."
DOMAIN = "x*y*z - x - a*y - z + a + 1"
PARAMS = {"a": -2}
PI  = dict(expression="3*x/10 + y/10 + z", solvefor="z")
RHO = dict(expression="3*x + 3*y/4", solvefor="x")

run(__name__)
