.. image:: https://github.com/georgehalal/sphericalrht/docs/images/sphericalrht_logo.gif?raw=true
   :align: center
   :width: 450 px
   :alt: sphericalrht Logo


.. image:: https://readthedocs.org/projects/sphericalrht/badge/?version=latest
   :target: https://sphericalrht.readthedocs.io/en/latest/?badge=latest
   :alt: Documentation Status

.. image:: https://github.com/georgehalal/sphericalrht/actions/workflows/tests.yml/badge.svg?branch=main
   :target: https://github.com/georgehalal/sphericalrht/actions/workflows/tests.yml?branch=main
   :alt: Test Status

.. image:: https://coveralls.io/repos/github/georgehalal/sphericalrht/badge.svg?branch=main
   :target: https://coveralls.io/github/georgehalal/sphericalrht?branch=main
   :alt: Coverage Status

.. image:: https://img.shields.io/badge/license-MIT-blue.svg?style=flat
   :target: https://github.com/georgehalal/sphericalrht/LICENSE
   :alt: License


.. class:: center

``sphericalrht`` is a fast, efficient implementation of the
`Rolling Hough Transform <http://seclark.github.io/RHT/>`_
that runs directly on the sphere, using spherical harmonic convolutions.
The documentation can be found at
`readthedocs.org <https://sphericalrht.readthedocs.io/en/latest/>`_
and the package is distributed over
`PyPI <https://pypi.org/project/sphericalrht/>`_.

.. contents::

============
Installation
============

Run the following in a terminal to install:

.. code-block:: bash

   $ pip install sphericalrht

If installing on a computing cluster, you may want to run the following
instead:

.. code-block:: bash
   
    $ pip install sphericalrht --user

Note that the ``sphericalrht`` package requires at least Python 3.7.


=====
Usage
=====

The code runs in parallel on as many CPUs as available, so feel free to
request many CPUs when submitting a job. The only input parameters that
affect the runtime and memory are ``nside`` and ``norients``. Here's an
example of how to run the algorithm and read in the results:

.. code-block:: python

   from sphericalrht import CubeAndStokes
   
   cube_and_stokes = CubeAndStokes(
       in_map="/path/to/input_map.fits",       # Input .fits map
       nside=1024,                             # NSIDE of output maps
       out_dir="/path/to/output_directory",    # Directory to save results in
       wlen=75,                                # window diameter [arcmins]
       fwhm=30,                                # high-pass filter scale [arcmins]
       thresh=0.7,                             # threshold fraction (0-1)
       norients=100)                           # number of orientation angles
   
   cube_and_stokes.build_and_save()
   
   
   # Load the output maps
   import healpy as hp
   
   I, Q, U = hp.read_map("/path/to/output_maps.fits", field=(0,1,2))
   
   
   # Optionally, load the output of all orientation angles for each pixel
   import h5py
   
   with h5py.File("/path/to/output_cube.h5") as cube_file:
       spherical_rht_out = cube_file["spherical_rht_cube"][:, PIXEL_INDEX]


========================
References & Attribution
========================

The paper introducing this package is in preparation. If you make use 
of this code in your research, please contact halalgeorge@gmail.com 
for discussing proper citations.
