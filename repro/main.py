'''main

This module executes the simulation and runs the plotting
functions. 
'''
import numpy as np
from pathlib import Path

from . import spectra


def main():
    # define wavelength range and spacing for simulations
    WLS = np.linspace(1, 4, 601) 
    
    # set the total number of simulations to run
    num_spectra = 100

    # Total water measured from step-heating experiments
    water_ppm=[1522, 762, 176, 22]

    # Low wavelength portion of MORB spectrum (<1.5 microns)
    MORB_D38A_LowLam = spectra.reflectance_from_file(
        Path('../data/Morb_D38A_Low_wavelength.txt'))

    # Use the 650C spectrum for the lower wavelength portion of the spectrum
    MORB_D38A_MidLam = spectra.reflectance_from_file(
        Path('../data/650^oC.csv'))
    
    # MORB step-wise heating spectra
    heated_MORB_spectra = []
    for temperature in [650, 700, 750, 800]:
        MORB_spectrum = spectra.normalize_to_wavelength(
            spectra.reflectance_from_file(
                Path(f'../data/{temperature}^oC.csv')),
            norm_wavelength=2.6  # microns
        )
        full_spectrum = spectra.concatenate([MORB_D38A_LowLam, MORB_spectrum])
        heated_MORB_spectra.append(full_spectrum)


if __name__ == "__main__":
    main()