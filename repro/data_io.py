'''data_io

Tools for reading data files included in this repository.
'''
import csv
from numpy.typing import ArrayLike
from pathlib import Path


def import_reflectance_measurements(band: float, path: Path) -> ArrayLike:
    pass