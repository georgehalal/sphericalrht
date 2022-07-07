# -*- coding: utf-8 -*-
"""
Spherical Rolling Hough Transform


A fast, efficient implementation of the Rolling Hough
Transform using spherical harmonic convolutions to
perform the algorithm directly on the sphere.

Classes:
    CubeAndStokes: Define an instance of this class with input
        parameters, then use the build_and_save method to run
        the algorthm.
    StokesQU: Handling Stokes linear polarization maps (only use
        if necessary).

Functions:
    set_logger: Log output in the terminal to a file (implemented
        within the CubeAndStokes class).

Author: George Halal
Email: halalgeorge@gmail.com
Date: 07/06/2022
Version: 1.2.2
"""


from .spherical_rht import StokesQU, CubeAndStokes
from .utils import set_logger


__author__ = "George Halal"
__email__ = "halalgeorge@gmail.com"
__version__ = "1.2.1"
__all__ = ["StokesQU", "CubeAndStokes", "set_logger"]
