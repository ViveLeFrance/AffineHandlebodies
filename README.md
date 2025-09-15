# AffineHandlebodies
This collection of notebooks implements the code to perform the calculations used to determine a Weinstein handlebody of a smooth complex affine surface, as described in [CM19](https://www.math.ucdavis.edu/~casals/LefsFronts_Reviewed.pdf). 

It is also intended as a companion repo to the paper (...) to document the various choices that lead to the resulting handlebodies.

## Binder

You can try out the code via Binder in your browser (link below) without installing anything on your machine. Note however that binder takes quite long to install the necessary environment (can be of the order of ~20 minutes). 

[![Binder](https://mybinder.org/badge_logo.svg)](https://mybinder.org/v2/gh/ViveLeFrance/AffineHandlebodies/main?urlpath=lab)

## Local Install
A more seamless experience is provided by running the notebooks locally, which requires [SageMath](https://doc.sagemath.org/html/en/installation/index.html), with version at least 10.2.

The notebooks moreover use `%matplotlib widget` in order to pan and zoom to better inspect the plots produced. This requires the package `ipympl`, which is not included in SageMath by default. If you have trouble finding the Python installation shipped with Sage, you can, for example, install it into your Sage Python environment by launching any notebook with the SageMath kernel, and executing a cell containing the command `pip install ipympl`.

### Using the Notebooks

The file `model.ipynb` works as a template to study your variety of choice. It is suggested you follow the instructions contained there. The files `pX.ipynb` record the explicit choices made for the Painlevé moduli spaces, and the other files are libraries for code used in the notebooks.

