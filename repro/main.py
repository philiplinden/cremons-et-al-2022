import numpy as np

import figure1
import figure2


def main():
    (endmember_spectra, heated_MORB_spectra,
     MORB_D38A_LowLam) = figure1.import_spectral_data('data')

    simple_spectra = figure1.simplify_spectra(heated_MORB_spectra,
                                              MORB_D38A_LowLam)

    interpolated_spectra = figure1.interpolate_spectra(simple_spectra)

    figure1.make_figure(endmember_spectra, interpolated_spectra)

    # Temperature of step-heating experiment in degrees Celsius
    temperature = [650, 700, 750, 800]
    # Total water measured from step-heating experiments in ppm
    water_ppm = [1522, 762, 176, 22]

    # set the hydration levels to interpolate across
    hydration_levels = np.linspace(0, 1666, 1667)

    interpolated_hydration_MORB = figure2.estimate_hydrated_spectra(
        interpolated_spectra, water_ppm, hydration_levels)
    figure2.make_figure(interpolated_hydration_MORB, hydration_levels)


if __name__ == '__main__':
    main()
