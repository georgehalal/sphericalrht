# -*- coding: utf-8 -*-
"""
Quantitative tests of the classes and functions in the sphericalrht
package.

Qualitative tests have been performed separately on toy models.

Author: George Halal
Email: halalgeorge@gmail.com
Date: 11/11/2021
"""


import os
import shutil

import ducc0
import numpy as np
import unittest

from sphericalrht import StokesQU, CubeAndStokes


NSIDE = 2
NPIX = 12 * NSIDE**2
NORIENTS = 6
TEST_MAP_PATH = os.path.abspath("tests/data/test_map.fits")
TEST_MAP2_PATH = os.path.abspath("tests/data/test_map2.fits")
OUT_DIR = os.path.abspath("tests/data/out")


class TestStokesQU(unittest.TestCase):
    """Test of the StokesQU class."""

    def test_update(self):
        """Test of the StokesQU.update method.
        """
        stokes = StokesQU(NPIX)
        cube = np.ones((NORIENTS, NPIX))
        angs = np.zeros((NORIENTS))
        stokes.update(cube, angs)

        self.assertEqual(np.sum(stokes.stokes_q, dtype=int), NORIENTS * NPIX)
        self.assertEqual(np.sum(stokes.stokes_u, dtype=int), 0)

    def test_normalize_and_weight(self):
        """Test of the StokesQU.normalize_and_weight method.
        """
        stokes = StokesQU(NPIX)
        norm = np.ones((NPIX))
        intensity = np.random.randn(NPIX)
        stokes.normalize_and_weight(norm, intensity)

        self.assertEqual(np.sum(stokes.stokes_q), 0)
        self.assertEqual(np.sum(stokes.stokes_u), 0)


class TestCubeAndStokes(unittest.TestCase):
    """Test of the CubeAndStokes class.
    """

    def test_unsharp_mask(self):
        """Test of the CubeAndStokes.unsharp_mask method.
        """
        c_a_s = CubeAndStokes(
            TEST_MAP_PATH, NSIDE, OUT_DIR, weighting=TEST_MAP2_PATH)
        test_map = np.ones((NPIX))
        umask_map = c_a_s.unsharp_mask(test_map)

        self.assertEqual(np.sum(umask_map), NPIX/2)

    def test_prep_intensity(self):
        """Test of the CubeAndStokes.prep_intensity method.
        """
        c_a_s = CubeAndStokes(TEST_MAP_PATH, 4, OUT_DIR)
        out_map = c_a_s.prep_intensity()
        self.assertEqual(out_map.shape[0], 12 * 4**2)

        c_a_s = CubeAndStokes(TEST_MAP_PATH, 1024, OUT_DIR)
        out_map = c_a_s.prep_intensity()
        self.assertEqual(out_map.shape[0], 12 * 1024**2)

    def test_make_ker(self):
        """Test of the CubeAndStokes.make_ker method.
        """
        c_a_s = CubeAndStokes(TEST_MAP_PATH, NSIDE, OUT_DIR)
        ker = c_a_s.make_ker()

        for i in range(4):
            self.assertEqual(ker.count(i), 1)

    def test_get_ker_alm(self):
        """Test of the CubeAndStokes.get_ker_alm method.
        """
        c_a_s = CubeAndStokes(TEST_MAP_PATH, NSIDE, OUT_DIR)
        alms = c_a_s.get_ker_alm()
        mmax = min(c_a_s.mmax, c_a_s.lmax)
        nalm = (((mmax+1)*(mmax+2))//2
                + (mmax+1)*(c_a_s.lmax-mmax))

        self.assertEqual(alms.shape[1], nalm)
        self.assertEqual(alms[0, 0], 1 / 2. / np.sqrt(np.pi))

    def test_get_ptg(self):
        """Test of the CubeAndStokes.get_ptg method.
        """
        c_a_s = CubeAndStokes(TEST_MAP_PATH, NSIDE, OUT_DIR)
        orients = np.zeros((NORIENTS))
        ptg = c_a_s.get_ptg(NPIX, orients)

        self.assertEqual(ptg[0, 2], 0)

    def test_save_cube_and_stokes(self):
        """Test of the CubeAndStokes.save_cube_and_stokes method.
        """
        c_a_s = CubeAndStokes(
            TEST_MAP_PATH, NSIDE, OUT_DIR, norients=NORIENTS,
            weighting=TEST_MAP_PATH)
        lmax = c_a_s.lmax
        mmax = min(c_a_s.mmax, lmax)

        nalm_ker = (((mmax+1)*(mmax+2))//2 + (mmax+1)*(lmax-mmax))
        nalm_map = ((lmax+1)*(lmax+2)) // 2

        alm_ker = (np.random.uniform(-1., 1., (1, nalm_ker))
                   + 1j*np.random.uniform(-1., 1., (1, nalm_ker)))
        alm_map = (np.random.uniform(-1., 1., (1, nalm_map))
                   + 1j*np.random.uniform(-1., 1., (1, nalm_map)))

        # make alms for m==0 real-valued
        alm_ker[:, 0:lmax+1].imag = 0.
        alm_map[:, 0:lmax+1].imag = 0.

        interpolator = ducc0.totalconvolve.Interpolator(alm_map, alm_ker,
                                                        separate=True,
                                                        lmax=lmax,
                                                        kmax=mmax,
                                                        epsilon=1e-4,
                                                        ofactor=1.5,
                                                        nthreads=0)
        intensity = np.ones((NPIX))
        c_a_s.save_cube_and_stokes(intensity, interpolator)

        cube_name = os.path.join(OUT_DIR, c_a_s.out_name + ".h5")
        map_name = os.path.join(OUT_DIR, "IQU_" + c_a_s.out_name + ".fits")
        self.assertTrue(os.path.exists(cube_name))
        self.assertTrue(os.path.exists(map_name))

    def test_build_and_save(self):
        """Test of the CubeAndStokes.build_and_save method.
        """
        c_a_s = CubeAndStokes(
            TEST_MAP_PATH, NSIDE, OUT_DIR, norients=NORIENTS, split_factor=1)
        c_a_s.build_and_save()

        cube_name = os.path.join(OUT_DIR, c_a_s.out_name + ".h5")
        map_name = os.path.join(
            OUT_DIR, "IQU_" + c_a_s.out_name + ".fits")
        self.assertTrue(os.path.exists(cube_name))
        self.assertTrue(os.path.exists(map_name))

        if os.path.exists(OUT_DIR):
            shutil.rmtree(OUT_DIR, ignore_errors=True)
        kernel_alms_dir = os.path.join(
            os.path.expanduser("~"), ".cache/sphericalrht/kernel_alms")
        if os.path.exists(kernel_alms_dir):
            shutil.rmtree(kernel_alms_dir, ignore_errors=True)
        c_a_s = CubeAndStokes((np.ones((NPIX)), c_a_s.name), float(NSIDE),
                              OUT_DIR, norients=float(NORIENTS), wlen=75.,
                              weighting=np.ones((NPIX)), split_factor=1.)
        c_a_s.build_and_save()

        self.assertTrue(os.path.exists(cube_name))
        self.assertTrue(os.path.exists(map_name))


if __name__ == "__main__":
    unittest.main()
