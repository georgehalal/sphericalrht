# -*- coding: utf-8 -*-
"""
Spherical Rolling Hough Transform


A fast, efficient implementation of the Rolling Hough
Transform using spherical harmonic convolutions to
perform the algorithm directly on the sphere.

Classes:
    StokesQU: Handling Stokes linear polarization maps (only use
        if necessary).
    CubeAndStokes: Define an instance of this class with input
        parameters, then use the build_and_save method to run
        the algorthm.

Author: George Halal
Email: halalgeorge@gmail.com
Date: 11/11/2021
Version: 1.1.0
"""


__author__ = "George Halal"
__email__ = "halalgeorge@gmail.com"
__version__ = "1.1.0"
__all__ = ["StokesQU", "CubeAndStokes"]


import os
import time
from collections import deque
import logging
from typing import Union, Tuple
from dataclasses import dataclass

import psutil
import healpy as hp
import numpy as np
import ducc0
import h5py

from .utils import set_logger


@dataclass(order=True)
class StokesQU:
    """Handle Stokes Q and U linear polarization maps.

    Methods:
        update: sum the contribution from different orientations
        normalize_and_weight: normalize and weight by the intensity
    """

    def __init__(self, npix: int) -> None:
        """Initialize linear Stokes maps.

        Args:
            npix (int): Number of map pixels
        """
        assert type(npix) is int, (
            "Number of pixels should bean integer")

        self.stokes_q = np.zeros((npix))
        self.stokes_u = np.zeros((npix))

        return None

    def update(self, spherical_rht_cube: np.ndarray,
               orient_angs: np.ndarray) -> None:
        """Sum the contribution from different orientations.

        Apply the Q/U formula from Clark & Hensley 2019.

        Args:
            spherical_rht_cube (np.ndarray((Norientations, Npix))):
                convolution result over different orientations
            orient_angs (np.ndarray((Norientations,))): orientation
                angles corresponding to the cube
        """
        assert spherical_rht_cube.shape[0] == orient_angs.shape[0], (
            "One of the inputs has the wrong dimensions")

        self.stokes_q += spherical_rht_cube.T.dot(np.cos(2. * orient_angs))
        self.stokes_u += spherical_rht_cube.T.dot(np.sin(2. * orient_angs))

        return None

    def normalize_and_weight(self, norm: np.ndarray,
                             intensity: np.ndarray) -> None:
        """Normalize and weight by the intensity.

        The maps are normalized such that the integral over angles = 1.
        They are then weighted by the intensity map.

        Args:
            norm (np.ndarray((Npix,))): normalization
            intensity (np.ndarray((Npix,))): input map
        """
        assert norm.shape[0] == intensity.shape[0], (
            "One of the inputs has the wrong dimensions")

        self.stokes_q[norm == 0] = np.nan
        self.stokes_u[norm == 0] = np.nan
        norm[norm == 0] = 1
        self.stokes_q = -intensity * np.divide(self.stokes_q, norm)
        self.stokes_u = intensity * np.divide(self.stokes_u, norm)

        return None


@dataclass
class CubeAndStokes:
    """The main class of the sphericalrht algorithm.

    Methods:
        unsharp_mask: high-pass filter the input map and make it binary.
        prep_intensity: process the input map and optionally calculate
            its alms.
        make_ker: select the pixels defining the convolution kernel.
        get_ker_alm: define the convolution kernel as a stick and
            calculate its alms.
        get_ptg: prepare pointing tensor.
        save_cube_and_stokes: run the convolution and save the resulting
            cube and maps.
        build_and_save: main function to use for running the algorithm
            and saving the resulting cube and maps.
    """

    def __init__(self, in_map: Union[str, Tuple[np.ndarray, str]], nside: int,
                 out_dir: str, wlen: int = 75, fwhm: float = 30.,
                 thresh: float = 0.7, norients: int = 100,
                 overwrite: bool = False) -> None:
        """Define necessary variables based on the input arguments and
        make necessary directories.

        Args:
            in_map (str or tuple(np.ndarray((Npix,)), str)): either
                path to input intensity map or a tuple of the input
                intensity map as an array along with its name as a str,
                which is used for saving log file, alms, spherical RHT
                cube, and Stokes Q/U maps. The rest of the input options
                will be appended to this name when saving.
            nside (int): output NSIDE for intensity and Stokes Q/U maps.
            out_dir (str): directory to save log file, alms,
                spherical RHT cube, and Stokes Q/U maps.
            wlen (int): convolution kernel window diameter [arcmins]
                (the scale at which to measure the orientation).
            fwhm (float): scale [arcmins] for the unsharp mask applied
                to pick out filamentary structure.
            thresh (float): threshold fraction of the window diameter
                between 0-1 applied to the result of the convolution.
                Higher thresholds focus on the main orientations only,
                while lower thresholds take more orientations into
                account, weighted by their intensity.
            norients (int): angular resolution given by the number of
                orientations to consider.
            overwrite (bool): whether to overwrite outputs of same name
                if they already exist.
        """
        assert type(in_map) is str or (
            type(in_map) is tuple and type(in_map[0]) is np.ndarray
            and type(in_map[1]) is str), ("Input map must be a path or a"
                                          " tuple(np.ndarray, name)")
        if type(in_map) is str:
            assert os.path.exists(in_map), (
                "Input map does not exist. Check path and try again.")
            self.in_map = in_map
            self.name = in_map.split("/")[-1].split(".fits")[0]
        else:
            self.in_map = in_map[0]
            self.name = in_map[1]
            assert np.sqrt(self.in_map.shape[0]/12) % 1 == 0, (
               "Input map has the wrong shape or number of pixels.")

        if type(nside) is float:
            assert nside % 1 == 0, "NSIDE must be an integer"
            nside = int(nside)
        else:
            assert type(nside) is int, "NSIDE must be an integer"
        self.nside = nside

        assert type(out_dir) is str, "Output directory must be a str"
        if not os.path.exists(out_dir):
            os.makedirs(out_dir)
        self.out_dir = out_dir

        if type(wlen) is float:
            assert wlen % 1 == 0, "Window diameter must be an integer"
            wlen = int(wlen)
        else:
            assert type(wlen) is int, "Window diameter must be an integer"
        assert wlen > 0, "Window diameter must be positive"
        self.wlen = wlen

        assert type(fwhm) is float or type(fwhm) is int, "FWHM type is invalid"
        assert fwhm > 0, "FWHM must be positive"
        self.fwhm = fwhm
        if fwhm % 1 == 0:
            self.fwhm = int(fwhm)

        assert type(thresh) is float or type(thresh) is int, (
            "Threshold type is invalid")
        assert thresh >= 0 and thresh < 1, "Threshold must be between 0-1"
        self.thresh = thresh

        if type(norients) is float:
            assert norients % 1 == 0, (
                "Number of orientations must be an integer")
            norients = int(norients)
        else:
            assert type(norients) is int, (
                "Number of orientations must be an integer")
        assert norients > 0, "Number of orientations must be positive"
        self.norients = norients

        self.overwrite = overwrite

        self.out_name = (f"{self.name}_nside{nside}_wlen{wlen}_fwhm{fwhm}"
                         f"_thresh{thresh}_norients{norients}")

        set_logger(os.path.join(out_dir, self.out_name + ".log"))
        logging.info(f"* Output directory: {out_dir}")

        self.kernel_alms_dir = os.path.join(
            os.path.expanduser("~"), ".cache/sphericalrht/kernel_alms")
        if not os.path.exists(self.kernel_alms_dir):
            os.makedirs(self.kernel_alms_dir)

        self.ker_nside = min(self.nside * 4, 4096)
        self.lmax = min(int(2.5*self.ker_nside - 1), int(1.25*4096 - 1))
        self.mmax = min(50, self.lmax)

        return None

    def unsharp_mask(self, original_map: np.ndarray) -> np.ndarray:
        """High-pass filter the input map and make it binary.

        Args:
            original_map (np.ndarray((Npix,))): map to
                high-pass filter

        Returns:
            high-pass filtered map of 1s and 0s (np.ndarray)
        """
        smoothed_map = hp.smoothing(
            original_map, fwhm=np.radians(self.fwhm / 60.))
        subtracted_map = original_map - smoothed_map

        return (subtracted_map > 0.).astype(int)

    def prep_intensity(self, return_alm: bool = False) -> Union[
            np.ndarray, Tuple[np.ndarray, np.ndarray]]:
        """Process the intensity map for saving with the Stokes Q/U
        maps and optionally calculate alms for the convolution.

        Args:
            return_alm (bool): whether to calculate and return alms

        Returns:
            intensity_out (np.ndarray((Npix,))): intensity map to
                save
            (optional) umask_intensity_alm (np.ndarray((1, Nalms))):
                processed alms used for convolution
        """
        if type(self.in_map) is str:
            intensity = hp.read_map(self.in_map, field=(0))
        else:
            intensity = self.in_map

        intensity_nside = hp.get_nside(intensity)

        if self.nside == intensity_nside:
            intensity_out = intensity
        else:
            intensity_out = hp.ud_grade(intensity, self.nside)

        if return_alm:
            if self.ker_nside == intensity_nside:
                intensity_in = intensity
            elif self.ker_nside == self.nside:
                intensity_in = intensity_out
            else:
                intensity_in = hp.ud_grade(intensity, self.ker_nside)

            umask_intensity = self.unsharp_mask(intensity_in)
            umask_intensity_alm = hp.map2alm(
                umask_intensity, self.lmax).reshape((1, -1))
            return intensity_out, umask_intensity_alm

        return intensity_out

    def make_ker(self):
        """Select line of pixels defining the kernel at the North Pole.

        Returns:
            line (collections.deque): double-ended queue of pixel
                indices defining the kernel
        """
        niters = int(30 * self.wlen/75 * self.ker_nside/1024)

        line = deque([0, 1])
        for i in range(niters):
            line.append(hp.get_all_neighbours(self.ker_nside, line[-1])[-2])
            line.appendleft(hp.get_all_neighbours(self.ker_nside, line[0])[0])
        line.append(2)
        line.appendleft(3)
        for i in range(niters):
            line.append(hp.get_all_neighbours(self.ker_nside, line[-1])[0])
            line.appendleft(hp.get_all_neighbours(self.ker_nside, line[0])[-2])

        return line

    def get_ker_alm(self) -> np.ndarray:
        """Define the convolution kernel as a stick at the North Pole
        and calculate its alms.

        Returns:
            ker_alm (np.ndarray((1, Nalms))): complex-valued kernel
                alms to use in the convolution
        """
        npix = 12 * self.ker_nside**2

        # create a mask at the North Pole with the size of the kernel
        vec = hp.ang2vec(0, 90, lonlat=True)
        disc = hp.query_disc(self.ker_nside, vec,
                             radius=np.radians(self.wlen / 2. / 60.),
                             inclusive=True)
        window = np.zeros((npix))
        window[disc] = 1

        ker = np.zeros((npix))
        line = self.make_ker()
        ker[line] = 1

        ker *= window
        ker_alm = hp.map2alm(ker)

        # smooth to prevent ringing and apply lmax and mmax cuts
        ker_alm = hp.smoothalm(ker_alm, fwhm=hp.nside2resol(self.nside))
        l, m = hp.Alm.getlm(3*self.ker_nside - 1, np.arange(ker_alm.shape[0]))
        ker_alm = ker_alm[np.logical_and(l < self.lmax+1, m < self.mmax+1)]

        # normalize
        ker_alm /= 2. * np.sqrt(np.pi) * ker_alm[0]

        return ker_alm.reshape((1, -1))

    def get_ptg(self, npix, orients) -> np.ndarray:
        """Prepare pointing tensor of co-latitudes, longitudes, and
        kernel orientations.

        Args:
            npix (int): Number of pixels to calculate the convolution on
            orients (np.ndarray((Norientations,))): Angles by which
                to rotate the kernel

        Returns:
            pointing tensor (np.ndarray((N, 3))): tensor of
                co-latitudes, longitudes, and kernel orientations
        """
        thetas, phis = hp.pix2ang(self.nside, np.arange(npix))
        psis = np.repeat(orients, phis.shape[0])
        phis = np.tile(phis, orients.shape[0])
        thetas = np.tile(thetas, orients.shape[0])

        return np.vstack((thetas, phis, psis)).T

    def save_cube_and_stokes(
            self, split_factor: int, intensity: np.ndarray,
            interpolator: ducc0.totalconvolve.Interpolator) -> None:
        """Run the convolution and save the resulting cube and maps.

        Args:
            split_factor (int): number of convolution splits based on
                the requested NSIDE and number of orientations to save
                on memory usage
            intensity (np.ndarray((Npix,)): intensity map to save
            interpolator (ducc0.totalconvolve.Interpolator): object
                encapsulating the convolution functionality
        """
        npix = 12 * self.nside**2
        stokes = StokesQU(npix)
        norm = np.zeros((npix))

        orients = np.linspace(0.0, np.pi, self.norients)
        # Angle values corresponding to the kernel rotation angles
        orient_angs = np.flip(orients)
        orients = np.array_split(orients, split_factor)
        orient_angs = np.array_split(orient_angs, split_factor)

        # Find the size of each split
        bnd = 0
        split_bounds = [0]
        for i in range(len(orients)):
            bnd += orients[i].shape[0]
            split_bounds.append(bnd)

        # Save convolution result as a data cube in hdf5 format
        out_fn = os.path.join(self.out_dir, self.out_name + ".h5")
        if os.path.exists(out_fn):
            os.remove(out_fn)
        f = h5py.File(out_fn, "w")
        spherical_rht_cube = f.create_dataset(name="spherical_rht_cube",
                                              shape=(self.norients, npix),
                                              dtype="f", compression="gzip")

        # Perform convolutions
        for i, (orients_split, orient_angs_split) in enumerate(
                zip(orients, orient_angs)):
            ptg = self.get_ptg(npix, orients_split)
            spherical_rht_temp = interpolator.interpol(ptg).reshape(
                orients_split.shape[0], npix)
            spherical_rht_cube[split_bounds[i]:split_bounds[i + 1], :] = (
                spherical_rht_temp)
            # Apply threshold
            spherical_rht_temp -= self.thresh
            spherical_rht_temp[spherical_rht_temp < 0] = 0.
            stokes.update(spherical_rht_temp, orient_angs_split)
            norm += spherical_rht_temp.sum(axis=0)
            process = psutil.Process(os.getpid())
            logging.info(f"* Using {process.memory_info().rss / 1e9}GB of"
                         " memory to make spherical RHT cube.")

        f.close()

        stokes.normalize_and_weight(norm, intensity)

        hp.write_map(os.path.join(self.out_dir, "IQU_" + self.out_name
                     + ".fits"), (intensity, stokes.stokes_q, stokes.stokes_u),
                     coord="G", column_names=["I", "Q", "U"], overwrite=True)

        return None

    def build_and_save(self) -> None:
        """Run the algorithm and save the resulting orientation cube and
        Stokes maps.
        """
        out_maps_name = os.path.join(
            self.out_dir, "IQU_" + self.out_name + ".fits")
        out_cube_name = os.path.join(self.out_dir, self.out_name + ".h5")
        if not self.overwrite and (os.path.exists(out_maps_name)
                                   and os.path.exists(out_cube_name)):
            logging.info("* Outputs already exist. "
                         "Change overwrite to True to overwrite them.")
            return None

        start_time = time.time()

        intensity_alm_file = os.path.join(
            self.out_dir, self.name
            + f"_alms_nside{self.ker_nside}_fwhm{self.fwhm}.npy")
        if os.path.exists(intensity_alm_file):
            logging.info("* Input map alms exist.")
            intensity_alm = np.load(intensity_alm_file)
            intensity = self.prep_intensity()
        else:
            intensity, intensity_alm = self.prep_intensity(return_alm=True)
            np.save(intensity_alm_file, intensity_alm)
            logging.info("* Created input map alms.")

        ker_alm_file = os.path.join(
            self.kernel_alms_dir, f"alms_nside{self.ker_nside}_"
            f"wlen{self.wlen}.npy")
        if os.path.exists(ker_alm_file):
            logging.info("* Kernel alms exist.")
            ker_alm = np.load(ker_alm_file)
        else:
            ker_alm = self.get_ker_alm()
            np.save(ker_alm_file, ker_alm)
            logging.info("* Created kernel alms.")

        # Use as many threads as available
        interpolator = ducc0.totalconvolve.Interpolator(intensity_alm, ker_alm,
                                                        separate=True,
                                                        lmax=self.lmax,
                                                        kmax=self.mmax,
                                                        epsilon=1e-4,
                                                        ofactor=1.5,
                                                        nthreads=0)

        logging.info("* Interpolator configured.")
        del intensity_alm, ker_alm

        # Decide how many convolution steps to perform in order not to
        # have everything in memory at once
        split_factor = np.ceil(self.nside**2/512**2 * self.norients/300)
        self.save_cube_and_stokes(split_factor, intensity, interpolator)

        logging.info("* Saved cube and maps in output directory.")
        logging.info(
            f"* Total run time = {(time.time()-start_time) / 60.} mins.")

        return None
