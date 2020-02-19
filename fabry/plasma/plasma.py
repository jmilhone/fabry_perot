from __future__ import division, print_function

from functools import partial

import numpy as np
from scipy import special

from ..core import models

try:
    import matplotlib.pyplot as plt
except ImportError:
    print("Matplotlib is not installed. Moving on...")
    pass

cx_fits = {40: [0.39004112, -34.24186523],
           4: [0.40712338, -33.82360615],
          }


def pcx_vfd_velocity_profile(r, r_anode, V_max):
    """Calculates the toroidal velocity profile for VFD on PCX

    Args:
        r (Union[np.ndarray, float]): radii to calculate profile on
        r_anode (float): location of anodes, assumed location of peak profile starting
        V_max (float): maximum velocity at r_anode
    Returns:
        np.ndarray: torodial velocity profile as a function of r
    """
    x = r / r_anode
    # numpy array case
    if isinstance(x, np.ndarray):
        velocity = np.where(x < 1, V_max*x, V_max)
        return velocity

    # single r value case
    if x < 1:
        return x*V_max
    else:
        return V_max


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
            vel[idx] = V_outer # * np.exp(-(r[idx] - R_outer) ** 2 / 4.0 ** 2)
    else:
        if r > R_outer:
            return V_outer# * np.exp(-(r - R_outer) ** 2 / 4.0 ** 2)

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


def vfd_density_profile(r, r_anode):
    offset = 0.2
    density = 0.85*(1.0 - (r/r_anode)**1.6)**2 + offset
    density = np.where(r > r_anode, offset, density)
    return density


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
    #print('using vfd velocity profile')
    #vel = pcx_vfd_velocity_profile(r, 32.0, Vouter)
    vel = pcx_velocity_profile(r, Lnu, R_outer, Vouter)


    # fig, ax = plt.subplots()
    # ax.plot(r, vel)
    # plt.show()
    vel_adjusted = vel * np.cos(theta)
    #print(vel.max(), vel_adjusted.max())
    # ToDo: Should really iterate over w0 to handle the He II complex
    w_shifted_max = models.doppler_shift(w0, np.max(vel_adjusted))
    sigma = models.doppler_broadening(w_shifted_max, mu, Ti)
    wavelength = np.linspace(-1, 1, nlambda) * 10.0 * sigma + w_shifted_max

    # Now to build a big spectrum matrix
    w_shifts = models.doppler_shift(w0, vel_adjusted)
    #full_spectrum = models.gaussian(wavelength[np.newaxis, :], w_shifts[:, np.newaxis], sigma, amp=1.0)#, norm=False)
    full_spectrum = models.gaussian(wavelength[np.newaxis, :], w_shifts[:, np.newaxis], sigma, amp=1.0)#, norm=False)

    # fig, ax = plt.subplots()
    # ax.plot(vel_adjusted, w_shifts)
    # plt.show()

    #dens = density_profile(r, rmax, Lne)
    #dens = dens[:, np.newaxis]
    dens = vfd_density_profile(r, 32.0)
    dens = dens[:, np.newaxis]
    #print('removing density profile')
    #dens = 1.0

    full_spectrum *= dens ** 2

    #fig, ax = plt.subplots(figsize=(12,9))
    #for idx, spec in enumerate(full_spectrum):
    #   ax.plot(wavelength, spec)#, 'C0')
    #   # ax.axvline(w_shifts[idx], color='C1')
    #ax.axvline(w0, color='k')
    #ax.set_xlim(488-0.03, 488)
    #plt.show()

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


def radiance_from_profiles(x, w0, Ti, V_adjusted, amp, mu=39.948, nlambda=2048, test_plot=False):
    """Calculates radiance for a plasma chord given profiles and position along the chord

    Args:
        x (np.ndarray): position along the chord
        w0 (float): central wavelength
        V_adjusted (np.ndarray): toroidal velocity array but adjusted with the dot product with the chord vector
        amp (np.ndarray): emissivity amplitude along the chord
        mu (float): relative mass of the ion
        nlambda (int): number of wavelengths to calculate

    Returns:
        Tuple(np.ndarray, np.ndarray): wavelength and radiance as a function of wavelength for the chord
    """
    w_shifted_max = models.doppler_shift(w0, np.max(V_adjusted))
    sigma = models.doppler_broadening(w_shifted_max, mu, Ti)
    wavelength = np.linspace(-1, 1, nlambda) * 10.0 * np.max(sigma) + w_shifted_max

    # Now to build a big spectrum matrix
    w_shifts = models.doppler_shift(w0, V_adjusted)

    spectra = models.gaussian(wavelength[np.newaxis, :], w_shifts[:, np.newaxis],
            sigma[:, np.newaxis], amp=amp[:, np.newaxis])#, norm=False)

    if test_plot:
        leg = False 
        fig, ax = plt.subplots(figsize=(10,10))
        for i in range(x.shape[0]):
            ax.plot(wavelength, spectra[i, :], label=f"{x[i]:.2f}")
        ax.set_xlabel("Wavelength (nm)")
        ax.set_ylabel("Spectral Emissivity (AU)")
        ax.axvline(w0, ls='--', color='k')

        # reset x limits
        xlo, xhigh = ax.get_xlim()
        xdiff = xhigh-xlo
        xmid = 0.5*(xhigh+xlo)
        xdiff /= 2.0
        ax.set_xlim(xmid-0.5*xdiff, xmid+0.5*xdiff)

        if leg:
            ax.legend()
        plt.show()
    radiance = np.trapz(spectra, x=x, axis=0)
    return wavelength, radiance


def linear_Ti_profile(r, Ti_center, Ti_edge, r_edge):
    """Calculates a linear Ti profile at the radii r

    Args:
        r (Union[float, np.ndarray]): radius or radius array to calculate Ti
        Ti_center (float): Ti at r=0
        Ti_edge (float): Ti at r=r_edge
        r_edge (float): Edge of the plasma where Ti(r=r_edge)=Ti_ege
    Returns:
        Ti at r=r with the same shape as r
    """
    return Ti_center + (Ti_edge - Ti_center) * r / r_edge


def vfd_chord(impact_factor, Ti_args, V_args, ne_args, w0=487.98634, rmax=40.0, nr=100, nlambda=2048, mu=39.948, test_plot=False):

    r, theta, x = calculate_r_theta_x_from_impact_factor(impact_factor, rmax=rmax, npts=nr)
    
    #for tup in zip(r, x, theta):
    #    print(tup)
    # calculate Ti profile
    Ti_center = Ti_args[0]
    Ti_edge = Ti_args[1]
    Ti_profile = linear_Ti_profile(r, Ti_center, Ti_edge, rmax)

    # calculate V profile
    Vouter = V_args[0]
    Lnu = V_args[1]
    R_outer = V_args[2]
    vel = pcx_velocity_profile(r, Lnu, R_outer, Vouter)
    V_adjusted = vel*np.cos(theta)

    # calculate amplitude profile
    r_anode = ne_args[0]
    dens = vfd_density_profile(r, r_anode)
    amp = dens**2

    if test_plot:
        fig, ax = plt.subplots(3, sharex=True, figsize=(8, 8))
        ax[0].plot(r, Ti_profile)
        ax[1].plot(r, vel)
        ax[1].plot(r, V_adjusted, '--')
        ax[2].plot(r, dens)

        ax[2].set_xlim(0.0, rmax)
        ax[2].set_xlabel("R (cm)")
        ax[0].set_ylabel(r"$T_i$ (eV)")
        ax[1].set_ylabel(r"$V_{\phi}$ (m/s)")
        ax[2].set_ylabel(r"$n_e$ (AU)")

        fig.subplots_adjust(hspace=0.0)
        fig.align_ylabels()
        plt.show()

    wavelength, radiance = radiance_from_profiles(x, w0, Ti_profile, V_adjusted, amp, mu=mu, 
            nlambda=nlambda, test_plot=test_plot)

    return wavelength, radiance
