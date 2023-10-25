"""spectra

This module defines classes that represent spectra for each regolith type.
"""

from dataclasses import dataclass
from typing import Tuple


@dataclass
class Endmember:
    name: str
    abundance_range: Tuple[float, float]
    max_abundance: float
    density: float
    grain_size_range: Tuple[float, float]


class HydratedMorbGlass(Endmember):
    """Hydrated mid-ocean-ridge basalt (MORB) glass"""

    sample_label = ""
    spectrum_label = "MORB D38A"
    density = 2.8
    grain_size = 69e-6


class Regolith(Endmember):
    density = 1.8
    grain_size = 32e-6
