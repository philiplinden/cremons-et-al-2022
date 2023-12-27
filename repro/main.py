import numpy as np

import base
from base import log
import figure1
import figure2
import spectra


def main():
    log.info('Generating all figures...')
    (endmember_spectra, heated_MORB_spectra,
     MORB_D38A_LowLam) = base.import_spectral_data('data')

    simple_spectra = figure1.simplify_spectra(heated_MORB_spectra,
                                              MORB_D38A_LowLam)

    interpolated_spectra = spectra.interpolate_spectra(simple_spectra)

    figure1.make_figure(endmember_spectra, interpolated_spectra)

    # set the hydration levels to interpolate across
    hydration_levels = np.linspace(0, 1666, 1667)

    interpolated_hydration_MORB = figure2.estimate_hydrated_spectra(
        interpolated_spectra, hydration_levels)
    figure2.make_figure(interpolated_hydration_MORB, hydration_levels)


if __name__ == '__main__':
    main()
