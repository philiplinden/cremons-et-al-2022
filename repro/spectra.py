
from dataclasses import dataclass


@dataclass
class Endmember:
    name: str
    min_abundance: float
    max_abundance: float


class HydratedGlass(Endmember):
    density: float = 2.8
    grain_size: 69e-6


class Regolith(Endmember):
    density: float = 1.8
    grain_size: float = 32e-6
