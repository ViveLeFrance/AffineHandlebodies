from handlebodies import circle, pl_path, run

NAME = "Painlevé III (D8)"
NOTES = "Surface taken from Mazzocco's paper."
DOMAIN = "x*y*z + x**2 + y**2 - y"
PI  = dict(expression="2*I*x + 2*I*y/5 - z", solvefor="z")
RHO = dict(expression="2*x/5 + y/10", solvefor="x")
RHO_PATHS = [
    # 200 samples does not resolve this one: the reading only settles
    # from about 400 up, and disagrees below that.
    ("monodromy hexagon", lambda cv: pl_path([0.8*cv.rho[k] for k in [0, 1, 2, 3, 4, 5, 0]], steps=400)),
    ("trivial circle", lambda cv: circle(8, steps=400)),
    ("upper circuit", lambda cv: pl_path([
        cv.rho[0] + 1j, cv.rho[4] + 1j,
        cv.rho[5].imag*1j + cv.rho[4].real - 1j,
        cv.rho[0].real + cv.rho[5].imag*1j - 1j, cv.rho[0] + 1j], steps=300)),
    ("lower circuit", lambda cv: pl_path([
        cv.rho[0] - 1j, cv.rho[4] - 1j,
        cv.rho[2].imag*1j + cv.rho[4].real + 1j,
        cv.rho[0].real + cv.rho[2].imag*1j + 1j, cv.rho[0] - 1j], steps=300)),
]

run(__name__)
