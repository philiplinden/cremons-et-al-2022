"""hapke


This module emulates the behavior of the following MATLAB functions:
    src/Hapke_Inverse_Function_Passive.m
    src/Hapke_Lidar_R_Function.m
    src/Hapke_Lidar_SSA_Function.m
    src/Hydrated_Regolith_Spectra_Generator_and_Retrieval.m
"""
from functools import partial
from typing import Callable, Tuple

import numpy as np
from scipy.optimize import fmin


def ordinary_least_squares(x: float, y: Callable, yx: float) -> float:
    """Ordinary Least Squares Function.

    y (Callable): The estimator function.
    x (float): The argument to y.
    yx (float): The observation. Must be the same type as the result of y.
    """
    return sum((y(x) - yx) ** 2)


def angular_width(filling_factor: float) -> float:
    """Angular width parameter, see Equation 3.

    Args:
        filling_factor (float, optional): Filling factor.
    """
    return (-3 / 8) * np.log(1 - filling_factor)


def backscattering(h: float, g: float) -> float:
    """Backscattering function, see Equation 2.

    This describes the opposition effect.

    Args:
        h (float): angular width parameter.
        g (float): phase angle in radians.
    """
    return 1 / (1 + (1 / h) * np.tan(g / 2))


def reflectance_from_ssa(
    SSA: float, P: float, mu: float, mu0: float, B: float
) -> float:
    """Estimate bidirectional reflectance from single-scattering albedo.

    This function applies the Hapke model to estimate reflectance from
    single-scattering albedo and other characteristics of a mixture and
    spectral observation parameters.

    see Equation 1 in the manuscript:

        R = (ω/4) (μ₀ / (μ + μ₀)) {(1 + B)P + H(ω)H₀(ω) - 1}

    Hapke, B. (2012). Theory of reflectance and emittance spectroscopy
        (2nd ed.). Cambridge University Press.

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
        float: bidrectional reflectance
    """

    def h_func(SSA: float, mu: float, mu0: float) -> Tuple[float, float]:
        """Ambartsumian-Chandrasekhar H functions.

        Computed using the approximation from equation 8.57 from Hapke (2012).

        Args:
            SSA (float): single-scattering albedo, aka ω
            mu (float): cosine of emission angle
            mu0 (float): coside of incident angle

        Returns:
            H, H0 (Tuple[float, float]): Ambartsumian-Chandrasekhar H-functions
        """
        gamma = np.sqrt(max(1 - SSA, 0))
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


def ssa_from_reflectance(
    reflectance: float,
    asymmetry_factor: float = 0.81,
    emission_angle: float = 0,
    incident_angle: float = 30,
    phase_angle: float = 30,
    filling_factor: float = 0.41,
    initial_guess: float = 0.5,
) -> float:
    """Estimate single-scattering albedo (SSA) from reflectance

    TLDR; uses ordinary least squares to fit an SSA that produces the observed
    reflectance when the SSA is fed back through the Hapke model.

    Equivalent to src/Hapke_Lidar_SSA_function.m

    Uses scipy.optimize.fmin (equivalent to MATLAB fminsearch) to minimize
    ordinary least squares distance between SSA obtained from the supplied
    reflectance, R, and the SSA from the estimated reflectance at each
    sample point in WLS.

    Default parameter values replicate src/Hapke_Inverse_Function_Passive.m

    Args:
        reflectance (float): Bidirectional reflectance, R. see Equation 1.
        wavelength (float): Wavelength in microns.
        asymmetry_factor (float, optional): Scattering asymmetry factor, p.
        emission_angle (float, optional): Emission angle in degrees.
            Defaults to 0.
        incident_angle (float, optional): Incident angle in degrees.
            Defaults to 30.
        phase_angle (float, optional): Phase angle in degrees. Defaults to 30.
        filling_factor (float, optional): Particle filling factor.
            Defaults to 0.41.
        initial_guess (float, optional): initial guess of SSA, ω₀
    """
    mu = np.cos(np.deg2rad(emission_angle))
    mu0 = np.cos(np.deg2rad(incident_angle))
    g = np.deg2rad(phase_angle)
    h = angular_width(filling_factor=filling_factor)
    B = backscattering(h, g)
    w0 = initial_guess

    # turn reflectance_from_ssa() into the form y=f(x)
    y = partial(reflectance_from_ssa, P=asymmetry_factor, mu=mu, mu0=mu0, B=B)

    # formulate ordinary least squares between estimated reflectance, y
    # and observed reflectance, yx, and arrange into the form y=f(x)
    OLS = partial(ordinary_least_squares, y=y, yx=reflectance)

    # find SSA that minimizes least squares difference between observed
    # reflectance and estimated (Hapke) reflectance
    opts = dict(
        disp=False,
        maxiter=10_000,
        maxfun=10_000,
        ftol=1e-7,
        xtol=1e-7,
    )
    result = fmin(OLS, w0, **opts)
    return result[0]
