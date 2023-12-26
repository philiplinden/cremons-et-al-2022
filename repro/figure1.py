"""Figure 1

Run this script to generate Figure 1.
"""

# external modules
from matplotlib import pyplot as plt
import numpy as np
import pandas as pd
from pathlib import Path

# local module for working with spectral data
import spectra

DEFAULT_DATA_DIR = Path('../data')


def import_spectral_data(data_dir: Path | str = DEFAULT_DATA_DIR):
    if isinstance(data_dir, str):
        data_dir = Path(data_dir)

    # Import reflectance spectra (processed to remove hydration and organics)
    samples = {
        'Mature Mare': data_dir / Path('Mare_70181_Spectra.txt'),
        'Mature Highlands': data_dir / Path('Highlands_62231_Spectra.txt'),
        'Pyroxene':
        data_dir / Path('Apollo15Sample15555ReddishBrownPyroxeneB.txt'),
        'Immature Mare': data_dir / Path('Mare_71061_Spectra.txt'),
        'Immature Highlands': data_dir / Path('Highlands_61221_Spectra.txt'),
    }

    lab_spectra = []
    for name, path in samples.items():
        lab_spectra.append(
            spectra.spectrum_from_file(path, name, header=2, delimiter='\t'))

    endmember_spectra = pd.concat(lab_spectra, axis=1)

    # Temperature of step-heating experiment in degrees Celsius
    temperature = [650, 700, 750, 800]
    # Total water measured from step-heating experiments in ppm
    water_ppm = [1522, 762, 176, 22]

    experiments = zip(temperature, water_ppm)

    # Import observations of MORB step-wise heating reflectance spectra
    imported_reflectances = []
    for deg_c, ppm in experiments:
        MORB_spectrum = spectra.spectrum_from_file(data_dir /
                                                   Path(f"{deg_c}^oC.csv"),
                                                   name=f'{ppm} ppm')
        MORB_spectrum.index = 1e4 / MORB_spectrum.index  # 1/cm to microns
        MORB_spectrum = MORB_spectrum / 100  # percent to decimal

        # since we converted frequency to wavelength, need to resort data
        MORB_spectrum = MORB_spectrum.sort_index()
        imported_reflectances.append(MORB_spectrum)

    heated_MORB_spectra = pd.concat(imported_reflectances, axis=1)

    # Low wavelength portion of MORB spectrum (<1.5 microns)
    MORB_D38A_LowLam = spectra.spectrum_from_file(
        data_dir / Path("Morb_D38A_Low_wavelength.txt"),
        name='MORB D38A low wavelengths',
    )

    return endmember_spectra, heated_MORB_spectra, MORB_D38A_LowLam


def normalize_to_wavelength(heated_MORB_spectra: pd.DataFrame,
                            ppm: int = 22,
                            micron: float = 2.6):
    scale_factor = spectra.get_normalization_factor(
        heated_MORB_spectra[f'{ppm} ppm'], micron)
    normalized_spectra = pd.DataFrame()
    MORB_spectrum = pd.DataFrame()
    for spectrum in heated_MORB_spectra.columns:
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
        heated_MORB_spectra)
    # Replace spectrum below 2.6 micron with 650C reflectance
    # to isolate 3 micron feature changes and stitch in low wavelengths
    mid_spectrum = normalized_spectra['22 ppm'].loc[MORB_spectrum.index < 2.6]
    mid_spectrum = mid_spectrum.loc[mid_spectrum.index > 1.6]

    # normalize lower spectrum to meet up with the rest of the spectrum
    stitch_factor = spectra.get_normalization_factor(mid_spectrum, 1.6) / \
                    spectra.get_normalization_factor(MORB_D38A_LowLam, 1.6)
    lower_spectrum = MORB_D38A_LowLam * stitch_factor

    simple_spectra = pd.DataFrame()
    for spectrum in normalized_spectra.columns:
        MORB_spectrum = normalized_spectra[spectrum]
        unified_spectrum = mid_spectrum.combine_first(lower_spectrum)
        upper_spectrum = normalized_spectra[spectrum].loc[
            normalized_spectra.index > 2.6]

        simple_spectra[spectrum] = unified_spectrum.combine_first(
            upper_spectrum)

    return simple_spectra


def interpolate_spectra(simple_spectra: pd.DataFrame):
    # define wavelength range and spacing for simulations as
    # 1-4 micron at 5nm intervals
    WLS = np.linspace(1, 4, 601)

    # interpolate along our sampling grid of interest
    interpolated_spectra = simple_spectra.apply(
        spectra.interpolate_along_series, args=[WLS], axis=0)

    return interpolated_spectra

def make_figure(endmember_spectra: pd.DataFrame,
                interpolated_spectra: pd.DataFrame,
                show: bool = True):
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

    fig.savefig('Figure1.png')

    if show:
        plt.show()

    return fig


def main():
    (endmember_spectra, heated_MORB_spectra,
     MORB_D38A_LowLam) = import_spectral_data()

    simple_spectra = simplify_spectra(heated_MORB_spectra, MORB_D38A_LowLam)

    interpolated_spectra = interpolate_spectra(simple_spectra)

    make_figure(endmember_spectra, interpolated_spectra)


if __name__ == '__main__':
    main()
