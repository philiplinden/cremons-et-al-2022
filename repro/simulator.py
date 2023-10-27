'''simulator

This module emulates the behavior of the following MATLAB functions:
    src/Hapke_Inverse_Function_Passive.m
    src/Hapke_Lidar_R_Function.m
    src/Hapke_Lidar_SSA_Function.m
    src/Hydrated_Regolith_Spectra_Generator_and_Retrieval.m
'''
from collections import namedtuple
from dataclasses import dataclass, field
from functools import partial
from typing import Any, Callable
import numpy as np
from scipy.optimize import fmin


DEFAULT_WLS = np.linspace(1, 4, 601)

Range = namedtuple('Range', 'min max')


@dataclass
class Endmember:
    name: str
    density: float
    grain_size: float
    water_ppm: float
    abundance: Range = field(default_factory=Range)


class HydratedMorbGlass(Endmember):
    '''Hydrated mid-ocean-ridge basalt (MORB) glass'''

    spectrum_label = "MORB D38A"
    density = 2.8
    grain_size = 69E-6
    abundance = Range(0.0, 0.3)


class Regolith(Endmember):
    density = 1.8
    grain_size = 32e-6


class Spectrum:

    def __init__(self, reflectance_data, species: Endmember, grid=DEFAULT_WLS):
        self.reflectance = reflectance_data
        self.species = species
        self.grid = grid

    def ssa(self, phasing,
            emission_angle=0,
            incident_angle=30,
            phase_angle=30,
            filling_factor=0.41):
        return reflectance_to_ssa(
            self.reflectance,
            self.grid,
            phasing,
            emission_angle,
            incident_angle,
            phase_angle,
            filling_factor,
        )


def ordinary_least_squares(x: Any, y: Callable, yx: Any):
    '''Ordinary Least Squares Function.

    y (Callable): The estimator function.
    x (Any): The argument to y.
    yx (Any): The observation. Must be the same type as the result of y.
    '''
    return sum((y(x) - yx) ** 2)


def angular_width(filling_factor):
    '''Angular width parameter, see Equation 3.

    Args:
        filling_factor (float, optional): Filling factor.
    '''
    return (-3 / 8) * np.log(1 - filling_factor)


def backscattering(h, g):
    '''Backscattering function, see Equation 2.

    This describes the opposition effect.

    Args:
        h (float): angular width parameter.
        g (float): phase angle in radians.
    '''
    return 1 / (1 + (1 / h) * np.tan(g / 2))


def bidirectional_reflectance(SSA, P, mu, mu0, B):
    '''Bidirectional reflectance, see Equation 1.

        R = (ω/4) (μ₀ / (μ + μ₀)) {(1 + B)P + H(ω)H₀(ω) - 1}

    I broke this equation down into smaller terms for code readability.
        let:
            R = r1 * r2 * (r3 + r4 - 1)

        where:
            r1 = (ω/4)
            r2 = (μ₀ / (μ + μ₀))
            r3 = (1 + B)P
            r4 = H(ω) * H₀(ω)

    Equivalent to src/Hapke_Lidar_R_function.m

    Args:
        SSA (_type_): single-scattering albedo, aka ω
        P (float): scattering phase function
        mu (float): cosine of emission angle
        mu0 (float): cosine of incident angle
        B (float): backscattering function

    Returns:
        _type_: _description_
    '''

    def h_func(SSA, mu, mu0):
        '''Ambartsumian-Chandrasekhar H functions.

        Computed using the approximation from equation 8.57 from Hapke (2012).

        Args:
            SSA (float): single-scattering albedo, aka ω
            mu (float): cosine of emission angle
            mu0 (float): coside of incident angle
        '''
        gamma = np.sqrt(1 - SSA)
        r0 = (1 - gamma) / (1 + gamma)

        h1 = np.log((1 + mu0) / mu0)
        h2 = (r0 + ((1 - 2 * r0) * mu0) / 2) * h1

        H = (1 - SSA * mu0 * h2) ** -1

        h3 = np.log((1 + mu) / mu)
        h4 = 1 - 2 * r0 * mu
        h5 = r0 + h3 * (h4 / 2)

        H0 = (1 - SSA * mu * h5) ** -1
        return H, H0

    H, H0 = h_func(SSA, mu, mu0)

    r1 = SSA / 4
    r2 = mu0 / (mu0 + mu)
    r3 = (1 + B) * P
    r4 = H * H0
    return r1 * r2 * (r3 + r4 - 1)


def reflectance_to_ssa(
    Refl,
    WLS,
    P=0.15,
    emission_angle=0,
    incident_angle=30,
    phase_angle=30,
    filling_factor=0.41,
):
    '''Convert reflectance spectrum to single-scattering albedo (SSA)

    Uses scipy.optimize.fmin (equivalent to MATLAB fminsearch) to minimize
    ordinary least squares distance between SSA obtained from the supplied
    reflectance, R, and the SSA from the estimated reflectance at each
    sample point in WLS.

    Hapke, B. (2012). Theory of reflectance and emittance spectroscopy (2nd ed.).
        Cambridge University Press.

    Equivalent to src/Hapke_Lidar_SSA_function.m
    Default parameter values replicate src/Hapke_Inverse_Function_Passive.m

    Args:
        R (array[float]): Bidirectional reflectance, see Equation 1.
        WLS: simulation sample grid?
        p (float, optional): Scattering phase function. Defaults to 0.15 for
            ansiotropic scattering on the modeled mean particle phase function
            for lunar soil (Goguen et al., 2010).
        emission_angle (float, optional): Emission angle in degrees.
            Defaults to 0.
        incident_angle (float, optional): Incident angle in degrees.
            Defaults to 30.
        phase_angle (float, optional): Phase angle in degrees. Defaults to 30.
        filling_factor (float, optional: Particle filling factor.
            Defaults to 0.41.
    '''
    mu = np.cos(np.deg2rad(emission_angle))
    mu0 = np.cos(np.deg2rad(incident_angle))
    g = np.deg2rad(phase_angle)
    h = angular_width(filling_factor=filling_factor)
    B = backscattering(h, g)

    w = []
    for m, x in zip(Refl, WLS):
        w0 = 0.5  # initial guess, ω₀
        # turn bidrectional_reflectance() into the form y=f(x)
        y = partial(bidirectional_reflectance, P=P, mu=mu, mu0=mu0, B=B)
        OLS = partial(
            ordinary_least_squares, y=y, yx=m
        )  # turn least squares into the form y=f(x)
        w.append(
            fmin(
                OLS,
                w0,
                args=x,
                disp=False,
                maxiter=10_000,
                maxfun=10_000,
                ftol=1e-7,
                xtol=1e-7,
            )
        )
    return np.concatenate(w)
