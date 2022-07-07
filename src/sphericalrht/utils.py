# -*- coding: utf-8 -*-
"""
Utility functions for the Spherical Rolling Hough Transform

set_logger: Log output in the terminal to a file.

Author: George Halal
Email: halalgeorge@gmail.com
Date: 07/06/2022
"""


import logging


__author__ = "George Halal"
__email__ = "halalgeorge@gmail.com"
__version__ = "1.2.2"
__all__ = ["set_logger"]


def set_logger(log_path: str) -> None:
    """Log output in the terminal to a file.

    Args:
        log_path (str): where to save the log file
    """
    logger = logging.getLogger()
    logger.setLevel(logging.INFO)

    file_handler = logging.FileHandler(log_path)
    file_handler.setFormatter(
        logging.Formatter("%(asctime)s:%(levelname)s: %(message)s"))
    logger.addHandler(file_handler)

    stream_handler = logging.StreamHandler()
    stream_handler.setFormatter(logging.Formatter("%(message)s"))
    logger.addHandler(stream_handler)
