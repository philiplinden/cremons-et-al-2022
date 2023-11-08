"""spectra

Definitions of end member mixtures and their spectra.
"""
from collections import namedtuple
from dataclasses import dataclass, field
from functools import partial
from pathlib import Path
from stat import FILE_ATTRIBUTE_OFFLINE
from typing import List

import numpy as np
from numpy.typing import ArrayLike
import pandas as pd
from pandas import DataFrame, Series
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
    path: Path,
    name='reflectance',
    header=None,
    delimiter=',',
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


def interpolate_along_series(data: Series, new_index: ArrayLike) -> Series:
    x = data.index.values
    y = data.values
    f = interp1d(x, y, fill_value='extrapolate', bounds_error=False)
    new_y = f(new_index)
    return Series(data=new_y, index=new_index)


def interpolate_hydration_levels(
    ssa: DataFrame, real_levels: List[int], interp_levels: List[int]
) -> DataFrame:
    """Interpolate across give hydration spectra at new hydration levels.

    For each index in the given data:
        1. Fit a line between the values in all columns.
            REAL_LEVELS denotes the x-coordinate to use for each respective col.
            The SSA values in each column denote the y-coordinate.
        2. Create synthetic values for each interpolated level.
            INTERP_LEVELS denotes the interpolation x-coordinates to compute.
            The resulting values denote the synthetic y-coordinate for a given
                index (row) of the SSA dataframe.

    Args:
        ssa (DataFrame): Real observations. All columns be sampled at the same
            indices.
        real_levels (List[int]): The level corresponding to each column in SSA.
        interp_levels (List[int]): The levels to sample the interpolation.

    Returns:
        DataFrame: DataFrame from interpolating "across" the input DF.
    """
    def ssa_to_espat(x):
        return (1-x)/x

    def espat_to_ssa(x):
        return 1/(x+1)
        
    def linear_fit(y_vals, x_vals):
        """y=mx+b"""
        fit_coefficients = np.polyfit(x_vals, y_vals, deg=1)
        m = fit_coefficients[0]
        b = fit_coefficients[1]
        return pd.Series((m, b), index=['m', 'b'])

    # convert SSA to ESPAT values
    espat = ssa.apply(ssa_to_espat)

    # apply linear for each row index to get coefficients
    fit_func = partial(linear_fit, x_vals=real_levels)
    fit_coeffs = espat.apply(fit_func, axis=1)
    
    # apply the coefficients at each row index
    interp_espat = fit_coeffs.apply(
        lambda w: Series(
            w.m * np.array(interp_levels) + w.b,
            index=interp_levels),
        axis=1
    )

    # # convert back to SSA
    return interp_espat.apply(espat_to_ssa)
