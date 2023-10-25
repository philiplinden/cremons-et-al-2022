"""main

This module executes the simulation and runs the plotting
functions. It emulates src/Hydrated_Regolith_Spectra_Generator_and_Retrieval.m
"""
from matplotlib import pyplot as plt
import numpy as np
from pprint import pprint

from .simulator import hapke, get_sample_grid


if __name__ == "__main__":
    WLS = get_sample_grid()
    Refl = np.array([np.random.randn() * x for x in WLS])
    R = hapke(Refl, WLS)
    pprint(WLS)
    pprint(Refl)
    pprint(R)
    plt.figure()
    plt.plot(WLS, Refl, WLS, R)
    plt.show()
