from __future__ import division, print_function
import numpy as np
from scipy import special
from ..core import models
from functools import partial
try:
    import matplotlib.pyplot as plt
except ImportError:
    pass

cx_fits = {40: [0.39004112, -34.24186523],
           4: [0.40712338, -33.82360615],
          }


def pcx_couette_velocity_profile(r, mom_dif_length, R_inner, R_outer, V_inner, V_outer):
    """Calculates the torodial velocity profile for PCX

    Args:
        r (Union[np.ndarray, float]): radii to calculate profile on
        mom_dif_length (float): momentum diffusion length scale
        R_inner (float): Radius of inner boundary
        R_outer (float): Radius of outer boundary
        V_inner (float): Velocity at inner boundary
        V_outer (float): Velocity at outer boundary

    Returns:
        np.ndarray: torodial velocity profile as a function of r
    """
    x = np.asarray(r) / mom_dif_length
    xi = R_inner / mom_dif_length
    xo = R_outer / mom_dif_length

    Iv = partial(special.iv, 1)
    Kv = partial(special.kv, 1)

    denom = Iv(xi) * Kv(xo) - Iv(xo) * Kv(xi)
    A = Kv(xo) * V_inner - Kv(xi) * V_outer
    B = Iv(xi) * V_outer - Iv(xo) * V_inner

    A /= denom
    B /= denom

    return A * special.iv(1, x) + B * special.kv(1, x)


def pcx_velocity_profile(r, mom_dif_length, R_outer, V_outer):
    """Calculates the toroidal velocity profile with no inner boundary for PCX

    Args:
        r (Union[np.ndarray, float]): raddi to calculate profile on
        mom_dif_length (float): momentum diffusion length scale
        R_outer (float): Radii for outer boundary
        V_outer (float): Velocity at outer boundary

    Returns:
        np.ndarray: torodial velocity profile as a function of r
    """
    x = np.asarray(r) / mom_dif_length

    xo = R_outer / mom_dif_length
    Iv = partial(special.iv, 1)
    vel = V_outer * Iv(x) / Iv(xo)
    
    if isinstance(r, np.ndarray):
        if any(rr > R_outer for rr in r):
            idx = np.where(r > R_outer)
            vel[idx] = V_outer * np.exp(-(r[idx] - R_outer) ** 2 / 4.0 ** 2)
    else:
        if r > R_outer:
            return V_outer * np.exp(-(r - R_outer) ** 2 / 4.0 ** 2)

    return vel


def density_profile(r, r_edge, gradient_length_scale):
    """Calculates the electron density profile

    Args:
        r (Union[np.ndarray, float]): radii to calculate profile on
        r_edge (float): edge of the electron density profile
        gradient_length_scale (float): length scale of the gradient at r_edge

    Returns:
       np.ndarray: electron density profile as a function of r
    """
    x = np.asarray(r)
    return 0.5 * (1.0 - np.tanh((x - r_edge) / gradient_length_scale))


def calculate_r_theta_x_from_impact_factor(impact_factor, rmax=150.0, npts=101):
    """Calculates the radius array, theta array, and the distance along a chord at the specified impact factor

    Args:
        impact_factor (float): impact factor of the chord
        rmax (float): max radius to include in chord
        npts (int): number of points to use

    Returns:
        tuple: (np.ndarray, np.ndarray, np.ndarray) r, theta, x
    """
    xmax = np.sqrt(rmax ** 2 - impact_factor ** 2)
    x = np.linspace(-1, 1, npts) * xmax
    r = np.sqrt(x ** 2 + impact_factor ** 2)
    theta = np.arctan2(x, impact_factor)

    return r, theta, x


def calculate_line_profile(wavelength, w0, Ti, vel, theta, mu):
    """Calculates the Gaussian line shape for a given set of parameters

    Args:
        wavelength (np.ndarray): wavelength array
        w0 (float): central wavelength
        Ti (float): temperature of emitting species
        vel (float): toroidal velocity in m/s
        theta (float): angle of torodial velocity to line of sight
        mu (float): relative mass in amu

    Returns:
        np.ndarray: gaussian line shape
    """
    vel_dot = vel * np.cos(theta)
    w_shift = models.doppler_shift(w0, vel_dot)
    sigma = models.doppler_broadening(w0, mu, Ti)

    return models.gaussian(wavelength, w_shift, sigma, norm=False)


def calculate_pcx_chord_emission(impact_factor, Ti, w0, mu, Lnu, Vouter, rmax=40.0, nr=101, nlambda=2000,
                                 Lne=2.5, R_outer=35):
    """Calculates PCX emission with only the outer boundary spinning for a given impact factor

    Args:
        impact_factor (float): impact factor for chord
        Ti (float): ion temperature in eV
        w0 (float): central wavelength
        mu (float): mass in amu
        Lnu (float): momentum diffusion length
        Vouter (float): velocity in m/s for outer boundary
        rmax (float): end of the plasma
        nr (int): number of radial points to integrate chord with
        nlambda (int): number of wavelength points
        Lne (float): density gradient scale length at rmax
        R_outer (float): velocity at outer boundary

    Returns:
        tuple: (np.ndarray, np.ndarray) wavelength and spectrum
    """
    r, theta, x = calculate_r_theta_x_from_impact_factor(impact_factor, rmax=rmax, npts=nr)
    vel = pcx_velocity_profile(r, Lnu, R_outer, Vouter)
    # fig, ax = plt.subplots()
    # ax.plot(r, vel)
    # plt.show()
    vel_adjusted = vel * np.cos(theta)

    # ToDo: Should really iterate over w0 to handle the He II complex
    w_shifted_max = models.doppler_shift(w0, np.max(vel_adjusted))
    sigma = models.doppler_broadening(w_shifted_max, mu, Ti)
    wavelength = np.linspace(-1, 1, nlambda) * 10.0 * sigma + w_shifted_max

    # Now to build a big spectrum matrix
    w_shifts = models.doppler_shift(w0, vel_adjusted)
    full_spectrum = models.gaussian(wavelength[np.newaxis, :], w_shifts[:, np.newaxis], sigma, amp=1.0, norm=False)

    # fig, ax = plt.subplots()
    # ax.plot(vel_adjusted, w_shifts)
    # plt.show()

    dens = density_profile(r, rmax, Lne)
    dens = dens[:, np.newaxis]

    full_spectrum *= dens ** 2

    # fig, ax = plt.subplots()
    # for idx, spec in enumerate(full_spectrum):
    #    ax.plot(wavelength, spec, 'C0')
    #    ax.axvline(w_shifts[idx], color='C1')
    # plt.show()

    # print(full_spectrum.shape)
    spectrum = np.trapz(full_spectrum, x=x, axis=0)
    # print(spectrum.shape)

    # fig, ax = plt.subplots()
    # ax.plot(wavelength, spectrum / spectrum.max(), 'C1')
    # plt.show()

    return wavelength, spectrum


def charge_exchange_rate(Ti, mu=40, noise=False):
    mass = int(mu)

    logTi = np.log(Ti)
    cx = np.polyval(cx_fits[mass], logTi)
    cx = np.exp(cx)
    if noise:
        cx = np.random.normal(loc=cx, scale=0.1*cx, size=1)
    return cx

def Lnu(ne_n0, Ti, mu=40, noise=False):
    sigv_cx = charge_exchange_rate(Ti, mu=mu, noise=noise)
    Lnu = np.sqrt(128 * 1e18 * Ti / (np.sqrt(mu) * ne_n0 * sigv_cx))
    return Lnu

