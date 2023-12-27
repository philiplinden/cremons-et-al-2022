"""Figure 1

Run this script to generate Figure 1.
"""
# external modules
from matplotlib import pyplot as plt
import pandas as pd

import base
from base import log
import spectra  # local module for working with spectral data


def normalize_to_wavelength(heated_MORB_spectra: pd.DataFrame,
                            ppm: int = 22,
                            micron: float = 2.6):
    scale_factor = spectra.get_normalization_factor(
        heated_MORB_spectra[f'{ppm} ppm'], micron)
    normalized_spectra = pd.DataFrame()
    MORB_spectrum = pd.DataFrame()
    for spectrum in heated_MORB_spectra.columns:
        log.info(f'Normalizing {spectrum} to {ppm} ppm at {micron} µm')
        MORB_spectrum = heated_MORB_spectra[spectrum]

        norm_factor = spectra.get_normalization_factor(MORB_spectrum, micron)
        norm_spectrum = MORB_spectrum * (scale_factor / norm_factor)
        normalized_spectra[spectrum] = norm_spectrum

    return normalized_spectra, MORB_spectrum


def simplify_spectra(heated_MORB_spectra: pd.DataFrame,
                     MORB_D38A_LowLam: pd.Series):
    # Normalize all MORB spectra to the reflectance at 2.6 microns at the
    # lowest water amount and highest temperature (22 ppm, 800 C)
    normalized_spectra, MORB_spectrum = normalize_to_wavelength(
        heated_MORB_spectra, ppm=22, micron=2.6)

    # Replace spectrum below 2.6 micron with 650C reflectance
    # to isolate 3 micron feature changes and stitch in low wavelengths
    log.info('Building simple mid spectrum between 1.6 µm and 2.6 µm')
    mid_spectrum = normalized_spectra['22 ppm'].loc[MORB_spectrum.index < 2.6]
    mid_spectrum = mid_spectrum.loc[mid_spectrum.index > 1.6]

    # normalize lower spectrum to meet up with the rest of the spectrum
    log.info('Normalizing MORB D38A reflectance below 1.6 µm')
    stitch_factor = spectra.get_normalization_factor(mid_spectrum, 1.6) / \
                    spectra.get_normalization_factor(MORB_D38A_LowLam, 1.6)
    lower_spectrum = MORB_D38A_LowLam * stitch_factor

    simple_spectra = pd.DataFrame()
    for spectrum in normalized_spectra.columns:
        log.info(f'Stitching together a spectrum for {spectrum} from parts')
        MORB_spectrum = normalized_spectra[spectrum]
        unified_spectrum = mid_spectrum.combine_first(lower_spectrum)
        upper_spectrum = normalized_spectra[spectrum].loc[
            normalized_spectra.index > 2.6]

        simple_spectra[spectrum] = unified_spectrum.combine_first(
            upper_spectrum)

    return simple_spectra


def make_figure(endmember_spectra: pd.DataFrame,
                interpolated_spectra: pd.DataFrame,
                show: bool = False):
    log.info('Building Figure 1...')
    fig, axes = plt.subplots(figsize=(20, 6), nrows=1, ncols=2)

    figure1a = endmember_spectra.plot(
        ax=axes[0],
        title='Figure 1a',
        xlabel='Wavelength (μm)',
        ylabel='Reflectance',
        grid=True,
        legend=True,
        xlim=(1, 4),
    )

    figure1b = interpolated_spectra.plot(
        ax=axes[1],
        title='Figure 1b',
        xlabel='Wavelength (μm)',
        ylabel='Reflectance',
        grid=True,
        legend=True,
        xlim=(1, 4),
    )

    figfile = 'repro/Figure1.png'
    fig.savefig(figfile)
    log.info(f'Saved {figfile}')
    if show:
        plt.show()

    return fig


def main():
    (endmember_spectra, heated_MORB_spectra,
     MORB_D38A_LowLam) = base.import_spectral_data()

    simple_spectra = simplify_spectra(heated_MORB_spectra, MORB_D38A_LowLam)

    interpolated_spectra = spectra.interpolate_spectra(simple_spectra)

    make_figure(endmember_spectra, interpolated_spectra)


if __name__ == '__main__':
    main()
