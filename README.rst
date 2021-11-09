.. raw:: html

   <div align="center">
   <img src="https://raw.githubusercontent.com/georgehalal/sphericalrht/main/docs/images/sphericalrht_logo.gif" width="450px">
   </img>
   <br/>
   <a href="https://badge.fury.io/py/sphericalrht">
   <img src="https://badge.fury.io/py/sphericalrht.svg" alt="PyPI version" height="18">
   </a>
   <a href='https://sphericalrht.readthedocs.io/en/latest/?badge=latest'>
   <img src='https://readthedocs.org/projects/sphericalrht/badge/?version=latest' alt="Documentation status" />
   </a>
   <a href="https://github.com/georgehalal/sphericalrht/actions/workflows/tests.yml">
   <img src="https://github.com/georgehalal/sphericalrht/actions/workflows/tests.yml/badge.svg" alt="Test status"/>
   </a>
   <a href='https://coveralls.io/github/georgehalal/sphericalrht'>
   <img src='https://coveralls.io/repos/github/georgehalal/sphericalrht/badge.svg' alt='Coverage Status' />
   </a>
   <a href="https://github.com/georgehalal/sphericalrht/LICENSE">
   <img src="https://img.shields.io/badge/License-MIT-blue.svg?style=flat" alt="License"/>
   </a>
   </div>
   <br/>

.. raw:: html

   <div align="center">
   <code>sphericalrht</code> is a fast, efficient implementation of the
   <a href="http://seclark.github.io/RHT/">Rolling Hough Transform</a>
   that runs directly on the sphere, using spherical harmonic convolutions.
   The documentation can be found at <a href="https://sphericalrht.readthedocs.io/en/latest/">readthedocs</a>
   and the package is distributed over <a href="https://pypi.org/project/sphericalrht/">PyPI</a>.
   </div>
   <br/>



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
