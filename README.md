# BH_local_abundance

This repository contains the analysis for the paper available at [https://arxiv.org/abs/2108.10123](https://arxiv.org/abs/2108.10123).

A component of this repository that may interest a broader audience is the Python file `balls.py`. 
While arguably childishly named, it provides a numpy vectorized implementation for calculating the volume of intersection of three spheres.
It mainly follows the procedure described in [K.D. Gibson and H.A. Scheraga. "Exact calculation of the volume and surface area of fused hard-sphere molecules with unequal atomic radii
"](https://doi.org/10.1080/00268978700102951).
It passed all the MCMC tests I produced, so I'm quite confident it is correct.

