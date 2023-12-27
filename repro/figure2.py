from functools import partial

# external modules
from matplotlib import pyplot as plt
import numpy as np
import pandas as pd

import base
from base import log, EXPERIMENT_WATER_PPM
import spectra  # local module for working with spectral data
import hapke  # local module for hapke model


def estimate_hydrated_spectra(interpolated_spectra,
                              hydration_levels, water_ppm = EXPERIMENT_WATER_PPM):
    log.info(f'Estimating SSA from reflectance using Hapke method')
    hapke_model = partial(hapke.ssa_from_reflectance,
                          asymmetry_factor=0.81,
                          emission_angle=0,
                          incident_angle=0,
                          phase_angle=0,
                          filling_factor=0.41,
                          initial_guess=0.5)
    ssa_MORB = pd.DataFrame()
    for ppm in interpolated_spectra:
        log.info(f'Applying Hapke function to hydration level {ppm}')
        ssa_MORB[ppm] = interpolated_spectra[ppm].apply(hapke_model)

    interpolated_hydration_MORB = spectra.interpolate_hydration_levels(
        ssa_MORB, water_ppm, hydration_levels.tolist())

    return interpolated_hydration_MORB


def make_figure(interpolated_hydration_MORB, hydration_levels, show=False):
    log.info('Building Figure 2...')
    subset_hydration_levels = hydration_levels[::100]

    fig, axes = plt.subplots(figsize=(10, 6), nrows=1, ncols=1)
    ax = interpolated_hydration_MORB[subset_hydration_levels].plot(
        ax=axes,
        title='Figure 2',
        xlabel='Wavelength (μm)',
        ylabel='Estimated Single-Scattering Albedo (SSA)',
        grid=True,
        legend=False,
        xlim=(2.5, 3.5),
        color=plt.cm.viridis(np.linspace(0, 1, len(subset_hydration_levels))),
        linewidth=1,
    )
    figfile = 'repro/Figure2.png'
    fig.savefig(figfile)
    log.info(f'Saved {figfile}')
    if show:
        plt.show()

    return fig

def main():
    from . import figure1

    water_ppm = [1522, 762, 176, 22]

    (endmember_spectra, heated_MORB_spectra,
     MORB_D38A_LowLam) = base.import_spectral_data()

    simple_spectra = figure1.simplify_spectra(heated_MORB_spectra,
                                              MORB_D38A_LowLam)

    interpolated_spectra = spectra.interpolate_spectra(simple_spectra)

    # set the hydration levels to interpolate across
    hydration_levels = np.linspace(0, 1666, 1667)

    estimate_hydrated_spectra(interpolated_spectra, water_ppm,
                              hydration_levels)


if __name__ == '__main__':
    main()
