from handlebodies import pl_path, run

NAME = "Painlevé IV"
DOMAIN = ("x*y*z + x**2 - (s2**2 + s1*s2)*x - s2**2*y - s2**2*z "
          "+ s2**2 + s1*s2**3")           # s1 in C, s2 in C*
PARAMS = {"s1": "1/2", "s2": "-1/2"}
PI  = dict(expression="3*x - y/2 + 2*z/3", solvefor="y")
RHO = dict(expression="-9*x/2 + z", solvefor="z")

run(__name__)
