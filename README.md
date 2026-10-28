# AffineHandlebodies
This toolkit implements the calculations described in [CM19](https://www.math.ucdavis.edu/~casals/LefsFronts_Reviewed.pdf), used to determine a Weinstein handlebody of a smooth complex affine surface. 

It is also intended as a companion repository for the article [Weinstein Handlebodies for the Painlevé Betti Spaces](https://arxiv.org/abs/2511.17352) to document the various choices made that lead to the resulting handlebodies.


## Run locally
Clone the repo, install the dependencies, and run `model.py`. This will open a page in your browser where you can input your surface analyze your bifibration. Computations will automatically use [Singular](https://www.singular.uni-kl.de/public/) to accelerate computations if it is installed and added to your PATH environment variable, which is helpful for hypersurfaces of high degree.

You can also run the files contained the `examples` folder, displaying the output which was used in the paper.

## Run online
There is a browser-hosted version on [my Website](https://jdbeimler.github.io/handlebodies/) which you can test without installing anything. It does computations without Singular, and does not indicate what Seidel calls "fake" critical values (given a bifibration $(\pi, \rho)$, this is the subset of regular values of $\pi$ $z \in \mathbb{C}$ such that $\rho$ restricted to $\pi^{-1}(z)$ is not Lefschetz).

