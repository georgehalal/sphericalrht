# -*- coding: utf-8 -*-
"""
This script contains tests of the functions in sphericalrht.utils.

Author: George Halal
Email: halalgeorge@gmail.com
Date: 11/11/2021
"""


import os

import unittest
import logging

import sphericalrht as sr


TEST_LOG_FILE = os.path.abspath("tests/data/out/test.log")


class TestSetLogger(unittest.TestCase):
    """Test of the sphericalrht.set_logger function."""

    def test_set_logger(self):
        """Test of the sphericalrht.set_logger function."""
        sr.set_logger(TEST_LOG_FILE)
        logging.info(f"{sr.__version__}")
        logging.info(f"{sr.__author__}")
        logging.info(f"{sr.__email__}")

        self.assertTrue(os.path.exists(TEST_LOG_FILE))


if __name__ == "__main__":
    unittest.main()
