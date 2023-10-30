"""spectra

Definitions of end member mixtures and their spectra.
"""
from collections import namedtuple
from dataclasses import dataclass, field
import numpy as np
from numpy.typing import ArrayLike
from pathlib import Path
from typing import List


Range = namedtuple('Range', ['min', 'max'])
Spectrum = namedtuple('Spectrum', ['wavelength', 'value'])


@dataclass
class Endmember:
    name: str
    density: float
    grain_size: float

    def __init__(self, reflectance_data: Path):
        self.reflectance = reflectance_from_file(reflectance_data)


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
    reflectance: Spectrum = field(default_factory=Spectrum)
    end_members: List[Endmember] = field(default_factory=list)


def reflectance_from_file(path: Path) -> Spectrum:
    with open(path, 'r') as file:
        arr = np.loadtxt(file, delimiter=',')
    wavelength = 1e4 / arr[:, 0]  # in cm^1, convert to microns
    value = arr[:, 1] / 100  # convert percent to decimal

    # since we converted frequency to wavelength, need to resort data
    sort_order = np.argsort(wavelength)
    return Spectrum(wavelength[sort_order], value[sort_order])


def normalize_to_wavelength(
        spect: Spectrum, norm_wavelength: float, tolerance: float = 1e-3
        ) -> Spectrum:
    wavelength = spect.wavelength
    values = spect.value

    norm = np.where(abs(wavelength - norm_wavelength) <= tolerance)[0][0]
    new_values = np.asarray(
        [v / norm for v in values]
    )
    return Spectrum(wavelength, new_values)


def concatenate(spectra: List[Spectrum]) -> Spectrum:
    all_wavelengths = [s.wavelength for s in spectra]
    all_values = [s.value for s in spectra]

    new_wavelengths = np.concatenate(all_wavelengths)
    new_values = np.concatenate(all_values)

    # sort by wavelength
    sort_index = np.argsort(new_wavelengths)

    return Spectrum(new_wavelengths[sort_index], new_values[sort_index])


def interpolate(src: Spectrum, new_wavelengths: ArrayLike) -> Spectrum:
    w, v = np.interp(new_wavelengths, src.wavelength, src.value)
    return Spectrum(w, v)