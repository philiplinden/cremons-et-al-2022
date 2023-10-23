from functools import partial
from typing import Any, Callable
import numpy as np
from scipy.optimize import fmin as fminsearch


def get_sample_grid(
    min_wavelength: float = 1, max_wavelength: float = 6, spacing: float = 0.005
):
    """Create a linear spectral simulation grid.

    Defines the discrete-space spectrum to be simulated. All units are in microns.

    Args:
        min_wavelength (float, optional): Lowest simulated wavelength. Defaults to 1 um.
        max_wavelength (float, optional): Highest simulated wavelength. Defaults to 6 um.
        spacing (float, optional): Spectral distance between grid samples. Defaults to 0.005 um.

    Returns:
        array: Simulation sample grid

    """
    start = min_wavelength
    stop = max_wavelength
    step = spacing
    grid = np.linspace(start, stop, int((stop - start) / step + 1))
    return grid


def ordinary_least_squares(x: Any, y: Callable, yx: Any):
    '''Ordinary Least Squares Function.

    y (Callable): The estimator function.
    x (Any): The argument to y.
    yx (Any): The observation. Must be the same type as the result of y.
    '''
    return sum((y(x) - yx) ** 2)


def angular_width(filling_factor=0.41):
    '''Angular width parameter, see Equation 3.

    Args:
        filling_factor (float, optional): Filling factor. Defaults to 0.41 for the
            lunar regolith (Bowell et al., 1989).
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


def R(SSA, P, mu, mu0, B):
    """Bidirectional reflectance, see Equation 1.

        R = (ω/4) (μ₀ / (μ + μ₀)) {(1 + B)P + H(ω)H₀(ω) - 1}

    This equation is HUGE so I broke it down into smaller terms for code readability.    
        let:
            R = r1 * r2 * (r3 + r4 - 1)

        where:
            r1 = (ω/4)
            r2 = (μ₀ / (μ + μ₀))
            r3 = (1 + B)P
            r4 = H(ω) * H₀(ω)

    Args:
        SSA (_type_): single-scattering albedo, aka ω
        P (float): scattering phase function
        mu (float): cosine of emission angle
        mu0 (float): cosine of incident angle
        B (float): backscattering function

    Returns:
        _type_: _description_
    """

    def Hfunc(SSA, mu, mu0):
        '''Ambartsumian-Chandrasekhar H functions.
        
        Computed using the approximation from equation 8.57 from Hapke (2012).
        
        Args:
            SSA (float): single-scattering albedo, aka ω
            mu (float): cosine of emission angle
            mu0 (float): coside of incident angle
        '''
        h1 = 1 - np.sqrt(1 - SSA)
        h2 = 1 + np.sqrt(1 - SSA)
        h3 = h1 / h2
        h4 = np.log((1 + mu0) / mu0)
        h5 = (h3 + ((1 - 2 * h3) * mu0) / 2) * h4

        H = (1 - SSA * mu0 * h5) ** -1

        h6 = np.log((1 + mu) / mu)
        h7 = (1 - 2 * h3 * mu)
        h8 = h3 + (h7 / 2) * h6

        H0 = (1 - SSA * mu * h8) ** -1
        return H, H0 
    
    H, H0 = Hfunc(SSA, mu, mu0)

    r1 = SSA / 4
    r2 = mu0 / (mu0 + mu)
    r3 = (1 + B) * P
    r4 = H * H0
    return r1 * r2 * (r3 + r4 - 1)


def Hapke(Refl, WLS, P=0.15, emission_angle=0, incident_angle=30, phase_angle=30, filling_factor=0.41):
    """Convert reflectance spectrum to single-scattering albedo (SSA)

    Uses scipy.optimize.fmin (equivalent to MATLAB fminsearch) to minimize 
    ordinary least squares distance between SSA obtained from the supplied
    reflectance, R, and the SSA from the estimated reflectance at each
    sample point in WLS.

    Hapke, B. (2012). Theory of reflectance and emittance spectroscopy (2nd ed.).
        Cambridge University Press.

    Args:
        R (array[float]): Bidirectional reflectance, see Equation 1.
        WLS: simulation sample grid?
        p (float, optional): Scattering phase function. Defaults to 0.15 for ansiotropic
            scattering on the modeled mean particle phase function for lunar soil
            (Goguen et al., 2010).
        emission_angle (float, optional): Emission angle in degrees. Defaults to 0.
        incident_angle (float, optional): Incident angle in degrees. Defaults to 30.
        phase_angle (float, optional): Phase angle in degrees. Defaults to 30.
        filling_factor (float, optional: Particle filling factor. Defaults to 0.41.
    """
    mu = np.cos(np.deg2rad(emission_angle))
    mu0 = np.cos(np.deg2rad(incident_angle))
    g = np.deg2rad(phase_angle)
    h = angular_width(filling_factor=filling_factor)
    B = backscattering(h, g)

    w = []
    for m, x in zip(Refl, WLS):
        w0 = 0.5  # initial guess, ω₀
        y = partial(R, P=P, mu=mu, mu0=mu0, B=B) # turn R() into the form y=f(x)
        OLS = partial(ordinary_least_squares, y=y, yx=m) # turn least squares into the form y=f(x)
        w.append(fminsearch(OLS, w0, args=x, disp=False, maxiter=10_000, maxfun=10_000, ftol=1e-7, xtol=1e-7))
    return np.concatenate(w)
