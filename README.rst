.. raw:: html

   <div align="center">
   <img src="https://raw.githubusercontent.com/georgehalal/sphericalrht/main/docs/images/sphericalrht_logo.gif" width="450px">
   </img>
   <br/>
   <a href="https://badge.fury.io/py/sphericalrht">
   <img src="https://badge.fury.io/py/sphericalrht.svg?branch=main&kill_cache=1" alt="PyPI version" height="18">
   </a>
   <a href='https://sphericalrht.readthedocs.io/en/latest/?badge=latest'>
   <img src='https://readthedocs.org/projects/sphericalrht/badge/?version=latest' alt="Documentation status" />
   </a>
   <a href="https://github.com/georgehalal/sphericalrht/actions/workflows/tests.yml">
   <img src="https://github.com/georgehalal/sphericalrht/actions/workflows/tests.yml/badge.svg" alt="Test status"/>
   </a>
   <a href="https://coveralls.io/github/georgehalal/sphericalrht">
   <img src="https://coveralls.io/repos/github/georgehalal/sphericalrht/badge.svg?branch=main&kill_cache=1" alt="Coverage status"/>
   </a>
   <a href="https://github.com/georgehalal/sphericalrht/LICENSE">
   <img src="https://img.shields.io/badge/License-MIT-blue.svg?style=flat" alt="License"/>
   </a>
   <a href="https://badge.fury.io/py/sphericalrht">
   <img alt="PyPI - Python Version" src="https://img.shields.io/pypi/pyversions/sphericalrht">
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
   Development of the Spherical RHT was supported in part by the National Science Foundation.
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

Add the ``--upgrade`` option to upgrade to a newer version.

Note that the ``sphericalrht`` package requires at least Python 3.7.


================
Input Parameters
================
Required
--------

``in_map``: Union[str, Tuple[np.ndarray, str]]
   either a path to the input intensity map or a tuple of the input
   intensity map as an array along with its name as a str,
   which is used for saving log file, alms, spherical RHT
   cube, and Stokes Q/U maps. The rest of the input options
   will be appended to this name when saving. The input map
   ordering is assumed to be RING.

``nside``: int
  output NSIDE for intensity and Stokes Q/U maps.

``out_dir``: str
  directory to save log file, alms, spherical RHT cube, and Stokes Q/U
  maps in COSMO/Healpix convention.

Optional
--------

``wlen``: int
  convolution kernel window diameter [arcmins] (the scale at which to
  measure the orientation).

``fwhm``: float
  scale [arcmins] for the unsharp mask applied to pick out filamentary
  structure.

``thresh``: float
  threshold fraction of the window diameter between 0-1 applied to the
  result of the convolution. Higher thresholds focus on the main
  orientations only, while lower thresholds take more orientations into
  account, weighted by their intensity.

``norients``: int
  angular resolution given by the number of orientations to consider. We
  have found empirically that ``norients`` = 25 is sufficient for most
  applications.

``weighting``: Union[str, np.ndarray]
  either a path to the map or an array of the map pixels. This is used as the
  weight for the output Stokes Q/U maps. The map ordering is assumed to be RING.

``overwrite``: bool
  whether to overwrite outputs of same name if they already exist.

``split_factor``: int
  number of data splits to save on runtime memory usage. Default value
  is based on the requested NSIDE (1 for NSIDE < 4096). If your job
  runs out of memory and you can't request more memory per job, increase
  the split factor.


=====
Usage
=====

The code runs in parallel on as many CPUs as available, so feel free to
request many CPUs when submitting a job. The runtime and memory increase
as O(``nside``\ :sup:`2` * ``norients``). The ``split_factor`` can be used
to decrease the memory complexity and increase the time complexity as
O(``split_factor``). The other input parameters have a negligible effect
on the time and memory complexity.

The polarization convention is COSMO with the polarization being
perpendicular to the filaments in the map in the plane of the sky. To
obtain polarization maps corresponding to the IAU convention, multiply U
by -1. To obtain polarization maps parallel to the orientation of
the filaments (e.g. magnetic field orientation) in the COSMO convention,
multiply both Q and U by -1.

Example 1
---------
Here's one way to run the algorithm with all the input parameters:

.. code-block:: python

    from sphericalrht import CubeAndStokes
    
    cube_and_stokes = CubeAndStokes(
        in_map="/path/to/map_name.fits",
        nside=1024,
        out_dir="/path/to/output_dir",
        wlen=75,
        fwhm=30,
        thresh=0.7,
        norients=25,
        weighting="/path/to/weighting_map.fits"
        overwrite=False,
        split_factor=1)

    cube_and_stokes.build_and_save()


Example 2
---------
If your input map is an array instead of a .fits file, you can enter a
tuple with the array in the first entry and the name as the second entry
as shown in this example that uses only the required input parameters:

.. code-block:: python

    import h5py
    
    with h5py.File("/path/to/map_name.h5", "r") as f:
        intensity = f["I"][:, 0]

    from sphericalrht import CubeAndStokes

    cube_and_stokes = CubeAndStokes(
        in_map=(intensity, "map_name"),
        nside=1024,
        out_dir="/path/to/output_dir") 

    cube_and_stokes.build_and_save()


Reading the results
-------------------
.. code-block:: python

    # Load the output maps
    import healpy as hp
    
    out_name = "map_name_nside1024_wlen75_fwhm30_thresh0.7_norients25"

    I, Q, U = hp.read_map(
        f"/path/to/output_dir/IQU_{out_name}.fits", field=(0,1,2))
    
    
    # If you'd like, you can also load the output of
    # all orientation angles for each pixel
    import h5py
    
    with h5py.File("/path/to/output_dir/{out_name}.h5") as cube_file:
        spherical_rht_out = cube_file["spherical_rht_cube"][:, PIXEL_INDEX]


========================
References & Attribution
========================

The paper introducing this package is in preparation. If you make use 
of this code in your research, please contact halalgeorge@gmail.com 
for discussing proper citations.
