"""A blank study: run this to work on a surface of your own.

    python model.py

opens the same page the examples do, with nothing filled in -- enter the
surface, pi and rho in the browser and compute from there.

The entries below are the other half of it.  Fill them in and run the file
again and the page opens with them already in place, computed; that is exactly
the shape of the files in examples/, so a study worth keeping is kept simply by
copying this one under a name of its own.  Coefficients are read exactly:
0.3 means 3/10, and the imaginary unit is I.
"""

from handlebodies import run

# NAME   = "my surface"
# NOTES  = "why these choices"
# DOMAIN = "x*y*z + x**2 + y**2 - y"
# PARAMS = {}                        # e.g. {"a": -2}
# PI     = dict(expression="2*I*x + 2*I*y/5 - z", solvefor="z")
# RHO    = dict(expression="2*x/5 + y/10", solvefor="x")

# Paths other than the straight vanishing ones can also be declared here, but
# the page's pickers are usually the easier way to find them.
# PI_PATHS  = []
# RHO_PATHS = []

run(__name__)
