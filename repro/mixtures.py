"""mixtures

Definitions of end member mixtures.
"""
from collections import namedtuple
from dataclasses import dataclass, field

Range = namedtuple('Range', 'min max')


@dataclass
class Endmember:
    name: str
    density: float
    grain_size: float
    water_ppm: float
    abundance: Range = field(default_factory=Range)

    def __init__(self, min_abundance: float, max_abundance: float):
        self.abundance = Range(min_abundance, max_abundance)


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
    
