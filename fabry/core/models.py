from __future__ import absolute_import, division, print_function
from collections import Iterable
import numpy as np
from scipy.integrate import trapz
from .zeeman import zeeman_lambda
from numba import jit
import os.path as path

@jit(nopython=True)
def trapezoidal_integration(y, x):
    """Performs trapezoidal intergration

    Args:
        y (Union[list, np.ndarray]): y data points
        x (Union[list, np.ndarray]): x data points

    Returns:
        float: area under the curve y(x)
    """
    n = len(x)
    area = 0.0
    for i in xrange(n - 1):
        area += (x[i + 1] - x[i]) * (y[i + 1] + y[i])
    return area / 2.0


@jit
def peak_calculator(L, d, w, order):
    """
    Simple peak calculator for ideal Fabry-Perot.
    
    .. math::
        r_j = L  \sqrt{ \left( \\frac{2d/w}{\\rm{Floor}(2d/w)-j}\\right)^2 - 1}

    Args:
        L (float): L value, units will set return units (pixels or mm)
        d (float): etalon spacing in mm
        w (float): peak wavelength in nm
        order (Union[np.ndarray,int]): j order number. 0 gives order nearest center
            of rings, increasing j as you move out.

    Returns:
        Union[np.ndarray, float]: the location of the peak(s)
            in the same units as L with the same length as order
    """
    m = 2.e6 * d / w
    m0 = np.floor(m)
    return L * np.sqrt(m ** 2 / (m0 - order) ** 2 - 1.0)


@jit
def airy_func(wavelength, cos_th, d, F):
    """
    Computes the Airy function (ideal Fabry-Perot instument function)
    as a function of wavelength and cos(theta) with parameters d and F

    .. math::
        A = \left(1 + Q  \sin^2(\pi \\frac{2d}{w} \cos\\theta)\\right)^{-1}
    
    .. math::
        Q = \left(\\frac{2  \mathcal{F} }{ \pi } \\right)^2

    Args:
        wavelength (np.ndarray): wavelength array in nm
        cos_th (np.ndarray): cos(theta) array
        d (float): etalon spacing in mm
        F (float): etalon finesse

    Returns:
        np.ndarray: evaluated airy function
    """
    Q = (2. * F / np.pi) ** 2
    airy = 1.0 / (1.0 + Q * np.sin(np.pi * 2.e6 * d * cos_th / wavelength) ** 2)
    return airy


@jit
def doppler_calc(w0, mu, temp, v):
    """
    Computes the doppler broadening sigma and the new central wavelength
    from the doppler shift

    .. math::
        \sigma = w_0  \sqrt{\\frac{k_B T }{ mc^2}}
    
    .. math::
        w = w_0 (1 - v/c)

    Args: 
        w0 (float): unshifted wavelength in nm
        mu (float): atomic mass in amu
        temp (float): temperature in eV
        v (float): velocity in m/s

    Returns:
        (float, float)): sigma in nm, shifted wavelength in nm
    """
    sigma = w0 * 3.2765e-5 * np.sqrt(temp / mu)
    w = w0 * (1.0 - 3.336e-9 * v)
    return sigma, w


@jit
def doppler_shift(w0, v):
    return w0 * (1.0 - 3.336e-9 * v)


@jit
def doppler_broadening(w0, mu, temp):
    return w0 * 3.2765e-5 * np.sqrt(temp / mu)


@jit
def gaussian(wavelength, w, sigma, amp=1., norm=True):
    """
    Computes a gaussian for a given central wavelength, sigma and amp
    
    .. math::
        G = \\frac{A}{\sigma \sqrt{2 \pi}} \exp{\left( \\frac{ (w - w_0)^2 }{2 \sigma^2 } \\right) }

    Args:
        wavelength (np.ndarray): wavelength array to calculate spec on
        w (float): central wavelength (same units as wavelength array)
        sigma (float): sigma for gaussian (same units as w)
        amp (float): amplitude of spectrum, default=1.0
        norm (bool): if true, the gaussian will be normalized, default=True
            to integrate to 1 over infinity then the amp factor will be multiplied

    Returns:
        np.ndarray: spectrum evaluated on wavelength array
    """
    if norm:
        norm = 1. / (sigma * np.sqrt(2. * np.pi))
    else:
        norm = 1.
    exp = np.exp(-0.5 * (wavelength - w) ** 2 / sigma ** 2)
    return amp * norm * exp


def lorentzian(wavelength, w, gamma, amp=1.):
    """
    Computes a lorentzian for a given central wavelength, gamma and amp
    
    .. math::
        L = \\frac{A}{2 \pi} \\frac{\gamma }{ (w-w_0)^2 + (\gamma/2)^2}

    Args:
        wavelength (np.array): wavelength array to calculate spec on
        w (float): central wavelength (same units as wavelength array)
        gamma (float): lorentzian gamma parameter
        amp (float, default=1.): amplitude in addition to one that integrates
            spec to 1 over infinity

    Returns:
        spec (np.ndarray): spectrum evaluated on wavelength array
    """
    A = (amp * 0.5 * gamma) / np.pi
    return A / ((wavelength - w) ** 2 + (0.5 * gamma) ** 2)


def offset_forward_model(r, L, d, F, w0, mu, amp, temp, v, nlambda=1024, sm_ang=False, coeff=0.15):
    """Forward model with an attempt to model the 'offset' from nuissance lines

    Args:
        r (np.ndarray): array of r values to compute model on 
        L (float): camera lens focal length, same units as r (pixels or mm)
        d (float): etalon spacing (mm)
        F (float): etalon finesse
        w0 (Union[float, list]): central wavelength(s) in nm
        mu (Union[float, list]): mass(es) in amu
        amp (Union[float, list]): amplitude(s) for the lines
        temp (Union[float, list]): temperature(s) in eV
        v (Union[float, list]): velocities in m/s
        nlambda (int): number of points in wavelength array, default=1024
        sm_ang (bool): use the small angle approx or not, default=True
        coeff (float): coefficient to model the relative amplitude of all the nuissance lines

    Returns:
        np.ndarray: array length of r of forward model

    """
    # print(L, d, F, w0, mu, amp, temp, v)
    # print(nlambda, sm_ang, coeff)
    vals = forward_model(r, L, d, F, w0, mu, amp, temp, v, nlambda=nlambda)
    # vals += max(amp) * coeff / (1.0 + F)
    vals += max(amp) * coeff / (1.0 + (2.0 * F / np.pi) ** 2)

    return vals


@jit
def forward_model(r, L, d, F, w0, mu, amp, temp, v, nlambda=1024):
    """
    Convolves the Doppler spectrum with the ideal Fabry-Perot Airy function.

    Args:
        r (np.ndarray): array of r values to compute model on 
        L (float): camera lens focal length, same units as r (pixels or mm)
        d (float): etalon spacing (mm)
        F (float): etalon finesse
        w0 (Union[float, list]): central wavelength(s) in nm
        mu (Union[float, list]): mass(es) in amu
        amp (Union[float, list]): amplitude(s) for the lines
        temp (Union[float, list]): temperature(s) in eV
        v (Union[float, list]): velocities in m/s
        nlambda (int): number of points in wavelength array, default=1024

    Returns:
        np.ndarray: array length of r of forward model
    """
    # if type(w0) in [list, tuple]:
    #    if not all([type(x) in [list,tuple] for x in [mu, amp, temp, v]]):
    #        raise ValueError('need to have a list for all spec params')
    #    if not all([len(x) == len(w0) for x in [mu, amp, temp, v]]):
    #        raise ValueError('spec params are not all the same length')
    if isinstance(w0, Iterable):
        # if not all(isinstance(x, Iterable) for x in [mu, amp, temp, v]):
        #     raise ValueError('Need to have a iterable for all spec params')
        # if not all(len(x) == len(w0) for x in [mu, amp, temp, v]):
        #     raise ValueError('spec params are not all the same length')

        sigma = []
        w = []
        for i, ww in enumerate(w0):
            width, new_w = doppler_calc(ww, mu[i], temp[i], v[i])
            sigma.append(width)
            w.append(new_w)
        # wavelength = np.linspace(min(w) - 10.*max(sigma), max(w) + 10.*max(sigma), nlambda)[:,np.newaxis]
        wavelength = np.linspace(min(w) - 10. * max(sigma), max(w) + 10. * max(sigma), nlambda)  # .reshape(nlambda, 1)
        spec = 0.0
        for idx, ww in enumerate(w):
            spec += gaussian(wavelength, ww, sigma[idx], amp[idx])

    else:
        # if not all([type(x) not in [list,tuple] for x in [mu, amp, temp, v]]):
        #    raise ValueError('need to have a list or not for all spec params')
        # if any(isinstance(x, Iterable) for x in [mu, amp, temp, v]):
        #     raise ValueError('all spec params must be an instance of Iterable or not an instance, no mixing')

        sigma, w = doppler_calc(w0, mu, temp, v)
        wavelength = np.linspace(w - 10. * sigma, w + 10. * sigma, nlambda)  # [:,np.newaxis]
        # wavelength = np.linspace(w - 10.*sigma, w + 10.*sigma, nlambda).reshape(nlambda, 1)
        spec = gaussian(wavelength, w, sigma, amp)

    # sigma, w = doppler_calc(w0, mu, temp, v)
    # wavelength = np.linspace(w - 10.*sigma, w + 10.*sigma, nlambda)#[:,np.newaxis]
    # spec = gaussian(wavelength, w, sigma, amp)
    # if sm_ang:
    #    cos_th = 1.0 - 0.5 * (r/L)**2
    # else:
    #    cos_th = L / np.sqrt(L**2 + r**2)
    cos_th = L / np.sqrt(L ** 2 + r ** 2)

    # cos_th = cos_th.reshape((1,len(r)))
    # cos_th = cos_th[np.newaxis, :]

    # model = trapz(spec*airy_func(wavelength, cos_th, d, F), wavelength, axis=0)
    model = np.zeros_like(cos_th)

    for idx, cos in enumerate(cos_th):
        # print(trapz(spec*airy_func(wavelength, cos, d, F), wavelength).shape)
        # model[idx] = trapz(spec*airy_func(wavelength, cos, d, F), wavelength)
        model[idx] = trapezoidal_integration(spec * airy_func(wavelength, cos, d, F), wavelength)
    return model


def match_finesse_forward(r, L, d, F, temp, v, errtemp=None, w0=487.98634, mu=39.948):
    sigma, w = doppler_calc(w0, mu, temp, v * 1000.)
    if errtemp is not None:
        errsigma, _ = doppler_calc(w0, mu, errtemp, v * 1000.)
        sigma = np.sqrt(sigma ** 2 + errsigma ** 2)
    wavelength = np.linspace(w - 10. * sigma, w + 10. * sigma, 512)[:, np.newaxis]
    spec = gaussian(wavelength, w, sigma, norm=False)

    cos_th = 1.0 - 0.5 * (r / L) ** 2
    model = trapz(spec * airy_func(wavelength, cos_th, d, F), wavelength, axis=0)
    return model


def lyon_temp_forward(r, L, d, F, current, T, V, E=None):
    w0 = 487.98634
    mu = 39.948

    # my previous calculation ###
    # zeeman_fac = [-1.4, -19./15., -17./15., -1., 1., 17./15., 19./15., 1.4]
    # zeeman_amp = [1.0, 3.0, 6.0, 10.0, 6.0, 3.0, 1.0]
    # Victor's calculation ###
    zeeman_fac = [-1., -17. / 15., -19. / 15., -1.4, 1.4, 19. / 15., 17. / 15., 1.]
    zeeman_amp = [20., 12., 6., 2., 2., 6., 12., 20.]

    B = (0.0133 / 80.) * current

    sblah, w = doppler_calc(w0, mu, T, V * 1.e3)
    if E is not None:
        eblah, _ = doppler_calc(w0, mu, E, V * 1.e3)
        sblah = np.sqrt(sblah ** 2 + eblah ** 2)
    lambdas, amps = zeeman_lambda(w, B, zeeman_fac, amps=zeeman_amp)
    mn = w - 10. * sblah
    mx = w + 10. * sblah
    wavelength = np.linspace(mn, mx, 1024, endpoint=True)[:, np.newaxis]
    spec = 0.
    for l, a in zip(lambdas, amps):
        sigma, _ = doppler_calc(l, mu, T, 0.0)  # already did the velocity shift
        if E is not None:
            esigma, _ = doppler_calc(l, mu, E, 0.0)
            sigma = np.sqrt(sigma ** 2 + esigma ** 2)
        spec += gaussian(wavelength, l, sigma, amp=a, norm=False)

    cos_th = 1.0 - 0.5 * (r / L) ** 2
    model = trapz(spec * airy_func(wavelength, cos_th, d, F), wavelength, axis=0)
    return model


# def lyon_temp_forward_prof(r,L,d,F,current,T,V):

def zeeman_with_lyon_profile(r, L, d, F, current, T, V, E=None):
    w0 = 487.98634
    mu = 39.948

    # Victor's calculation ###
    zeeman_fac = [-1., -17. / 15., -19. / 15., -1.4, 1.4, 19. / 15., 17. / 15., 1.]
    zeeman_amp = [20., 12., 6., 2., 2., 6., 12., 20.]

    current_dir = path.abspath(path.join(__file__, ".."))
    b_data = np.genfromtxt(path.join(current_dir, "lyon_magnetic_field.csv"), delimiter=",")

    z = b_data[:, 0]
    bz = b_data[:, 1]

    # I only want to deal with the array where the plasma is emitting
    zbounds = [0.0, 80.0]  # Victor says that plasma exists here
    i_lower = np.abs(z-zbounds[0]).argmin()
    i_higher = np.abs(z-zbounds[1]).argmin()
    sl = slice(i_lower, i_higher)
    z = z[sl]
    bz = bz[sl]

    # Adjust bz for the current in the coil
    bz *= current / 80.0

    nz = len(z)
    nw = 1024
    spectrum = np.zeros((nz, nw))

    sigma_Ti = doppler_broadening(w0, mu, T)
    w = doppler_shift(w0, V)

    # Extra broadening from defocusing the camera lens
    if E:
        sigma_extra = doppler_broadening(w0, mu, E)
        sigma_Ti = np.sqrt(sigma_Ti**2 + sigma_extra**2)

    # Need to loop over reach z location
    for idx, (zz, bb) in enumerate(zip(z, bz)):
        # calculate the spectrum here
        spectrum[idx, :] += 0.0

        w_zee, a_zee = zeeman_lambda(w, bb, zeeman_fac, amps=zeeman_amp)

    final_spectrum = np.trapz(spectrum, z, axis=0)

    return final_spectrum
def general_model(r, L, d, F, wavelength, emission):
    cos_th = L / np.sqrt(L ** 2 + r ** 2)
    cos_th = cos_th.reshape((1, len(r)))
    airy = airy_func(wavelength, cos_th, d, F)

    emis = emission[:, np.newaxis]
    w = wavelength[:, np.newaxis]
    model = trapz(emis * airy, w, axis=0)

    return model
