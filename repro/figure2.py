# external modules
from functools import partial
from matplotlib import pyplot as plt
import numpy as np
import pandas as pd
from pathlib import Path

# local module for working with spectral data
import spectra
# local module for hapke model
import hapke


def estimate_hydrated_spectra(interpolated_spectra, water_ppm,
                              hydration_levels):

    hapke_model = partial(hapke.ssa_from_reflectance,
                          asymmetry_factor=0.81,
                          emission_angle=0,
                          incident_angle=0,
                          phase_angle=0,
                          filling_factor=0.41,
                          initial_guess=0.5)
    ssa_MORB = pd.DataFrame()
    for ppm in interpolated_spectra:
        ssa_MORB[ppm] = interpolated_spectra[ppm].apply(hapke_model)

    interpolated_hydration_MORB = spectra.interpolate_hydration_levels(
        ssa_MORB, water_ppm, hydration_levels.tolist())

    return interpolated_hydration_MORB


def make_figure(interpolated_hydration_MORB, hydration_levels, show=False):

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

    fig.savefig('Figure2.png')

    if show:
        plt.show()

    return fig

def main():
    from . import figure1

    water_ppm = [1522, 762, 176, 22]

    (endmember_spectra, heated_MORB_spectra,
     MORB_D38A_LowLam) = figure1.import_spectral_data()

    simple_spectra = figure1.simplify_spectra(heated_MORB_spectra,
                                              MORB_D38A_LowLam)

    interpolated_spectra = figure1.interpolate_spectra(simple_spectra)

    # set the hydration levels to interpolate across
    hydration_levels = np.linspace(0, 1666, 1667)

    estimate_hydrated_spectra(interpolated_spectra, water_ppm,
                              hydration_levels)


if __name__ == '__main__':
    main()
