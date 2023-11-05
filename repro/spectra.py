"""spectra

Definitions of end member mixtures and their spectra.
"""
from collections import namedtuple
from dataclasses import dataclass, field
from pathlib import Path
from typing import List

import numpy as np
from numpy.typing import ArrayLike
import pandas as pd
from pandas import Series
from scipy.interpolate import interp1d

Range = namedtuple('Range', ['min', 'max'])


@dataclass
class Endmember:
    name: str
    density: float
    grain_size: float

    def __init__(self, reflectance_data: Path):
        self.reflectance = spectrum_from_file(reflectance_data)


class Regolith(Endmember):
    density = 1.8


class MatureHighlands(Regolith):
    grain_size = 32e-6


class Pyroxene(Endmember):
    density = 2.8
    grain_size = 69e-6


class HydratedMorbGlass(Endmember):
    """Hydrated mid-ocean-ridge basalt (MORB) glass

    Laboratory reflectance data from step-wise heating experiments are used
    to simulate varying levels of hydration in lunar regolith.
    """

    spectrum_label = 'MORB D38A'
    density = 2.8
    grain_size = 69e-6


@dataclass
class Mixture:
    water_ppm: float
    abundance: Range = field(default_factory=Range)
    reflectance: Series = field(default_factory=Series)
    end_members: List[Endmember] = field(default_factory=list)


def spectrum_from_file(
    path: Path, name='reflectance', header=None, delimiter=',',
) -> Series:
    spectrum = pd.read_csv(
        path,
        index_col=0,
        header=header,
        delimiter=delimiter,
        names=['wavelength', name],
    )
    return spectrum.squeeze()


def get_normalization_factor(data: Series, normalization_index: float) -> float:
    """get value at the index closest to the desired normalization index"""
    wavelength = data.index.values
    norm_index = np.argmin(abs(wavelength - normalization_index))
    return data.iloc[norm_index]


def interpolate_series(data: Series, new_index: ArrayLike) -> Series:
    x = data.index.values
    y = data.values
    f = interp1d(x, y, fill_value='extrapolate', bounds_error=False)
    new_y = f(new_index)
    return Series(data=new_y, index=new_index)
